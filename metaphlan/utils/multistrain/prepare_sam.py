import bz2
import os
import re
import shlex
import shutil
import subprocess as sp
import tempfile

from ...utils import info, warning


# ==============================================================================================
# Input preprocessing: MetaPhlAn's SAM (or BAM) output -> sorted, filtered BAM of the hit markers.
# One pass over the input (an external bzip2, or samtools view for BAM, runs in parallel with the parsing): the
# header is split off and filtered once the hit markers are known, alignments are filtered (mapping quality,
# indels, percent identity, viral markers), then samtools converts, collates, fixes mates and sorts.
# ==============================================================================================
_HEADER_END = re.compile(rb'\n[^@]')
_CO_INDEX = re.compile(rb'^@CO\tindex:([^\n]*)', re.M)
_CHUNK = 1 << 23


class DatabaseMismatch(ValueError):
    pass


class InputError(ValueError):
    """The input cannot be read (truncated or corrupt compressed file)."""
    pass


class NoAlignments(Exception):
    """The input has no alignment to keep (empty file, header only, or every read filtered out)."""
    pass


def _input_kind(input_path):
    """'sam' or 'bam' (optionally .bz2 compressed)."""
    p = input_path.with_suffix('') if input_path.suffix == '.bz2' else input_path
    if p.suffix in ('.sam', '.bam'):
        return p.suffix[1:]
    raise ValueError(f'{input_path}: only SAM and BAM files, optionally bz2-compressed, are supported')


def _open_stream(input_path, kind):
    """SAM text of the input, and the processes producing it. SAM: an external bzip2 when available (else Python's
    bz2 module); BAM: samtools view -h (after bzip2 for .bam.bz2)."""
    if kind == 'bam':
        if input_path.suffix == '.bz2':
            p0 = sp.Popen(['bzip2', '-dc', str(input_path)], stdout=sp.PIPE)
            p1 = sp.Popen(['samtools', 'view', '-h', '-'], stdin=p0.stdout, stdout=sp.PIPE, bufsize=_CHUNK)
            p0.stdout.close()
            return p1.stdout, [p0, p1]
        p1 = sp.Popen(['samtools', 'view', '-h', str(input_path)], stdout=sp.PIPE, bufsize=_CHUNK)
        return p1.stdout, [p1]
    if input_path.suffix == '.bz2':
        exe = shutil.which('bzip2')
        if exe is not None:
            p = sp.Popen([exe, '-dc', str(input_path)], stdout=sp.PIPE, bufsize=_CHUNK)
            return p.stdout, [p]
        return bz2.open(input_path, 'rb'), []
    return open(input_path, 'rb'), []


def _close_stream(stream, procs):
    stream.close()
    for p in procs:
        p.wait()


def _tmp_dir_for(output_dir):
    t = os.environ.get('MS_TMPDIR') or os.environ.get('TMPDIR')
    if t and os.path.isdir(t) and os.access(t, os.W_OK):
        return t
    return str(output_dir)


def _check_db_name(header_data, db_name):
    """@CO index:<database> lines (MetaPhlAn >= 4.2) must name the database used here."""
    for m in _CO_INDEX.finditer(header_data):
        db_sam = m.group(1).decode().strip()
        if db_sam != db_name:
            raise DatabaseMismatch(f'The sample was mapped against the database {db_sam}, not {db_name}')


def get_read_lens_from_bam(sam_file):
    read_lens = []
    for row in sam_file.fetch():
        read_lens.append(row.qlen)

    return read_lens


def _pipe_together_procs(commands, stdout):
    """Pipes the commands together; returns every process (to check each exit status) and the first one's stdin.
    The parent's copies of the intermediate pipes are closed, so a failing process cannot leave the others
    blocked."""
    procs = []
    stdin = sp.PIPE
    for i, c in enumerate(commands):
        last = i == len(commands) - 1
        p = sp.Popen(shlex.split(c), stdin=stdin, stdout=stdout if last else sp.PIPE)
        if procs:
            procs[-1].stdout.close()
        procs.append(p)
        stdin = p.stdout
    return procs, procs[0].stdin


def _filter_alignments(stream, header_raw, body, db_name, input_path, config):
    """Single pass over the SAM text: the leading header block goes to header_raw (database name checked as soon as
    it is read), filtered alignment lines go to body. Returns the markers with hits and the kept read lengths."""
    min_mapq = config['min_mapping_quality']
    filter_indels = config['filter_indels']
    min_pident = config['min_read_pident']
    fix_pairend = config['fix_pairend']
    remove_viral = config['remove_viral_hits']

    selected_markers = set()
    read_lens = []
    read_name_warning_raised = False
    late_sq_names = []  # SN fields of @SQ lines after the first alignment
    with open(header_raw, 'wb') as fh, open(body, 'wb') as fb:
        in_header = True
        carry = b''
        out = []
        while True:
            chunk = stream.read(_CHUNK)
            if not chunk:
                data = carry
                carry = b''
                if not data:
                    break
                if not data.endswith(b'\n'):
                    data += b'\n'
            else:
                data = carry + chunk
                cut = data.rfind(b'\n') + 1
                carry = data[cut:]
                data = data[:cut]
                if not data:
                    continue
            # data now holds whole lines and ends with a newline
            if in_header:
                if not data.startswith(b'@'):
                    in_header = False
                else:
                    m = _HEADER_END.search(data)
                    if m is None:
                        _check_db_name(data, db_name)
                        fh.write(data)
                        if not chunk:
                            break
                        continue
                    _check_db_name(data[:m.start() + 1], db_name)
                    fh.write(data[:m.start() + 1])
                    data = data[m.start() + 1:]
                    in_header = False
            if data.startswith(b'\n') or b'\n\n' in data:
                raise ValueError(f'Blank line among the alignments of {input_path}')
            # lines starting with @ among the alignments: @SQ lines of markers without hits are dropped, any other
            # one goes to samtools unchanged, in place (samtools rejects real header lines, but accepts a record whose
            # name starts with @)
            has_header_lines = data.startswith(b'@') or b'\n@' in data
            # alignment lines
            lines = data.split(b'\n')
            lines.pop()  # empty string after the final newline
            for line in lines:
                if has_header_lines and line.startswith(b'@'):
                    if line.startswith(b'@SQ'):
                        late_sq_names.append(line.split(b'\t')[1])
                    else:
                        _check_db_name(line, db_name)
                        out.append(line)
                    continue
                line_fields = line.split(b'\t')
                read_name, flag, marker, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, *tags = line_fields
                if remove_viral and marker.startswith(b'VDB'):
                    continue
                selected_markers.add(marker)
                if int(mapq) < min_mapq:
                    continue
                if filter_indels and (cigar[-1] != ord('M') or not cigar[:-1].isdigit()):
                    continue
                read_len = len(seq)
                tags = (x.split(b':') for x in tags)
                tags = {x[0]: x[2] for x in tags}
                n_mismatches = int(tags[b'NM'])
                pident = 100 * (read_len - n_mismatches) / read_len
                if pident < min_pident:
                    continue
                read_lens.append(read_len)
                if fix_pairend:
                    i = read_name.rfind(b'__')
                    if i != -1:
                        read_name = read_name[:i]
                    elif not read_name_warning_raised:
                        warning(f'Read {read_name} does not have __ in its name')
                        read_name_warning_raised = True
                    out.append(b'\t'.join((read_name, *line_fields[1:])))
                else:
                    out.append(line)
                if len(out) >= 100000:
                    fb.write(b'\n'.join(out) + b'\n')
                    out = []
            if not chunk:
                break
        if out:
            fb.write(b'\n'.join(out) + b'\n')
    for sn in late_sq_names:
        assert sn.startswith(b'SN:')  # in bowtie2 output SN is always the first tag
    if any(sn[3:] in selected_markers for sn in late_sq_names):
        raise ValueError(f'@SQ line of a marker with hits after the first alignment line in {input_path}')
    return selected_markers, read_lens


def prepare_sam(input_path, output_bam_path, output_dir, db_name, config):
    """Writes the sorted, filtered BAM of the hit markers to output_bam_path; returns the kept read lengths.
    The uncompressed copy of the kept alignments goes to a temporary folder in $MS_TMPDIR (--tmp), else $TMPDIR,
    else output_dir."""
    kind = _input_kind(input_path)
    if input_path.stat().st_size == 0:
        raise NoAlignments('the input file is empty (0 bytes)')
    tmp_dir = output_dir  # samtools -T prefix
    work = tempfile.mkdtemp(prefix='s2m_', dir=_tmp_dir_for(output_dir))
    try:
        return _prepare_sam(input_path, kind, output_bam_path, output_dir, db_name, config, tmp_dir, work)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def _prepare_sam(input_path, kind, output_bam_path, output_dir, db_name, config, tmp_dir, work):
    header_raw = os.path.join(work, 'header.sam')
    header_filt = os.path.join(work, 'header.filt.sam')
    body = os.path.join(work, 'body.sam')

    info('  Single pass: header split and alignment filtering')
    stream, procs = _open_stream(input_path, kind)
    bzip2_failed = False
    try:
        selected_markers, read_lens = _filter_alignments(stream, header_raw, body, db_name, input_path, config)
    except Exception:
        _close_stream(stream, procs)
        # a bzip2 failure (not a signal, which only means we stopped reading) is retried below; anything else raises
        if kind != 'sam' or not procs or procs[0].returncode <= 0:
            raise
        bzip2_failed = True
    except BaseException:
        _close_stream(stream, procs)
        raise
    else:
        _close_stream(stream, procs)
        if kind == 'sam':
            bzip2_failed = bool(procs) and procs[0].returncode != 0
        else:
            failed = [f'"{" ".join(p.args)}" exited with status {p.returncode}' for p in procs if p.returncode != 0]
            if failed:
                raise RuntimeError(f'Reading {input_path} failed: ' + '; '.join(failed))
    if bzip2_failed:
        # bzip2 rejects some files that Python's bz2 module accepts, e.g. trailing bytes that start with a stream
        # signature
        warning(f'bzip2 exited with status {procs[0].returncode} on {input_path}; reading it with Python\'s bz2 module')
        try:
            with bz2.open(input_path, 'rb') as stream:
                selected_markers, read_lens = _filter_alignments(stream, header_raw, body, db_name, input_path,
                                                                 config)
        except (EOFError, OSError) as e:
            raise InputError(f'{input_path} is truncated or not a valid bz2 file ({e})') from None
    info(f'    Selected {len(selected_markers)} markers')
    if not read_lens:
        if not selected_markers:
            raise NoAlignments('the input has no alignments to MetaPhlAn markers (empty or header-only file)')
        raise NoAlignments(f'none of the alignments to {len(selected_markers)} markers passes the read filters '
                           f'(min_mapping_quality, filter_indels, min_read_pident)')

    # filter the header: keep non-@SQ lines and @SQ lines of hit markers, in input order
    with open(header_raw, 'rb') as fh, open(header_filt, 'wb') as fo:
        carry = b''
        while True:
            chunk = fh.read(1 << 24)
            data = carry + chunk
            if not chunk:
                if not data:
                    break
                lines = data.split(b'\n')
                carry = b''
            else:
                cut = data.rfind(b'\n') + 1
                carry = data[cut:]
                lines = data[:cut].split(b'\n')[:-1]
            keep = []
            for line in lines:
                if not line:
                    continue
                if line.startswith(b'@SQ'):
                    f1 = line.split(b'\t', 2)[1]
                    assert f1.startswith(b'SN:')
                    if f1[3:] in selected_markers:
                        keep.append(line)
                else:
                    keep.append(line)
            if keep:
                fo.write(b'\n'.join(keep) + b'\n')
            if not chunk:
                break
    os.unlink(header_raw)

    info('  Converting, collating and sorting')
    commands_2 = ['samtools view -Sb -']
    if config['fix_pairend']:
        commands_2.append(f'samtools collate -O -T {tmp_dir} -')
        commands_2.append('samtools fixmate -m - -')
    commands_2.append(f'samtools sort -T {tmp_dir} -')
    if config['dedup_stats']:
        if config['fix_pairend']:
            markdup_stats_path = output_dir / 'markdup_stats.txt'
            commands_2.append(f'samtools markdup -f {markdup_stats_path} -T {tmp_dir} - -')
        else:
            warning('Dedup stats works only with fix_pairend unfortunately, we need to figure out a workaround for '
                    'already paired data')

    broken_pipe = False
    with open(output_bam_path, 'wb') as fo:
        procs, stdin2 = _pipe_together_procs(commands_2, stdout=fo)
        try:
            for path in (header_filt, body):
                with open(path, 'rb') as f:
                    shutil.copyfileobj(f, stdin2, _CHUNK)
            stdin2.close()
        except BrokenPipeError:  # samtools stopped reading; reported below
            broken_pipe = True
            try:
                stdin2.close()
            except BrokenPipeError:
                pass
        codes = [p.wait() for p in procs]
    failed = [f'"{c}" exited with status {rc}' for c, rc in zip(commands_2, codes) if rc != 0]
    if failed or broken_pipe:
        raise RuntimeError(f'samtools failed while preparing {output_bam_path}: '
                           + ('; '.join(failed) if failed else 'it stopped reading its input'))

    return read_lens
