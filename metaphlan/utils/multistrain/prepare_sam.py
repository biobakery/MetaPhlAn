import re
from typing import Sequence, IO, Any
import bz2
from collections import Counter
import shlex
import subprocess as sp

from ...utils import info, warning, error


md_pat = re.compile(br'(\d+)')


def pipe_together(commands, stdin=sp.PIPE, stdout=sp.PIPE):
    """

    :param Sequence[str] commands:
    :param int|None|IO stdin: file to pass to the first command as input
    :param int|None|IO stdout: file to which the last command will write
    :return: returns the tuple (last process, the stdin of the first process, stdout of the last process)
    """
    ret_stdin = stdin
    ret_p = None
    for i, c in enumerate(commands):
        c = shlex.split(c)
        if i == len(commands) - 1:  # pass stdout to the last command
            p = sp.Popen(c, stdin=stdin, stdout=stdout)
        else:
            p = sp.Popen(c, stdin=stdin, stdout=sp.PIPE)
        stdin = p.stdout
        if ret_stdin == sp.PIPE:
            ret_stdin = p.stdin
        ret_p = p

    return ret_p, ret_stdin, stdin


def prepare_sam(input_path, output_bam_path, output_dir, db_name, config):
    """

    :param pathlib.Path input_path:
    :param pathlib.Path output_bam_path:
    :param pathlib.Path output_dir:
    :param db_name:
    :param dict config:
    :return:
    """
    tmp_dir = output_dir

    # Define the input file
    if input_path.suffix == '.bz2':
        sample_path_uncompressed = input_path.with_suffix('')
        fi = bz2.open(input_path, 'rb')
    else:
        sample_path_uncompressed = input_path
        fi = open(input_path, 'rb')


    ### First pass
    info('  First pass to pre-filter hits')
    commands = []
    if sample_path_uncompressed.suffix == '.bam':
        commands.append('samtools view -')
    elif sample_path_uncompressed.suffix != '.sam':
        raise Exception('Only sam and bam files with possible bz2 compression are supported')

    _, _, read_stdout = pipe_together(commands, fi)

    marker_hits: Counter[bytes] = Counter()
    for line in read_stdout:
        if line.startswith(b'@'):
            if line.startswith(b'@CO\tindex:'):
                db_sam = line[len(b'@CO\tindex:'):].decode().strip()
                if db_sam != db_name:
                    error(f'The database of the sample {db_sam} does not match {db_name}', exit=True)
            continue
        line_fields = line.rstrip(b'\n').split(b'\t')
        marker = line_fields[2]
        if config['remove_viral_hits'] and marker.startswith(b'VDB'):
            continue
        marker_hits[marker] += 1

    selected_markers = set(marker_hits.keys())
    info(f'    Selected {len(selected_markers)} markers')


    ### Second pass
    info('  Second pass to obtain the prepared BAM file')
    fi.seek(0)  # reset the file to the beginning

    commands_1 = []
    # convert BAM to SAM if BAM
    if sample_path_uncompressed.suffix == '.bam':
        commands_1.append('samtools view -')
    elif sample_path_uncompressed.suffix != '.sam':
        raise Exception('Only sam and bam files with possible bz2 compression are supported')

    if sample_path_uncompressed.suffix == '.bam':
        commands_1.append('samtools view -')

    # here goes the filtering step in between commands_1 and commands_2

    commands_2 = []
    # convert SAM to BAM
    if sample_path_uncompressed.suffix == '.sam':
        commands_2.append('samtools view -Sb -')

    # fix the PE information
    if config['fix_pairend']:
        commands_2.append(f'samtools collate -O -T {tmp_dir} -')
        commands_2.append('samtools fixmate -m - -')

    # sort the BAM
    commands_2.append(f'samtools sort -T {tmp_dir} -')

    # mark duplicate reads and produce stats
    if config['dedup_stats']:
        if config['fix_pairend']:
            markdup_stats_path = output_dir / 'markdup_stats.txt'
            commands_2.append(f'samtools markdup -f {markdup_stats_path} -T {tmp_dir} - -')
        else:
            warning('Dedup stats works only with fix_pairend unfortunately, we need to figure out a workaround for '
                    'already paired data')

    p1, stdin1, stdout1 = pipe_together(commands_1, fi)
    fo = open(output_bam_path, 'wb')
    p2, stdin2, stdout2 = pipe_together(commands_2, stdout=fo)

    # Filtering of lines
    read_lens = []
    read_name_warning_raised = False
    for line in stdout1:
        line_fields = line.rstrip(b'\n').split(b'\t')

        if line.startswith(b'@SQ'):
            assert line_fields[1].startswith(b'SN:')  # in bowtie2 output SN is always the first tag
            marker = line_fields[1][3:]
            if marker in selected_markers:
                stdin2.write(line)
        elif line.startswith(b'@'):  # other header lines like @HD, @PG, @CO
            stdin2.write(line)
        else:  # mapping lines
            read_name, flag, marker, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, *tags = line_fields
            if marker not in selected_markers:
                continue

            if int(mapq) < config['min_mapping_quality']:
                continue

            if config['filter_indels'] and (cigar[-1] != ord('M') or not cigar[:-1].isdigit()):  # cigar must be \d+M
                continue

            read_len = len(seq)
            tags = (x.split(b':') for x in tags)
            tags = {x[0]: x[2] for x in tags}
            n_mismatches = int(tags[b'NM'])
            pident = 100 * (read_len - n_mismatches) / read_len
            if pident < config['min_read_pident']:
                continue

            read_lens.append(read_len)

            if config['fix_pairend']:
                i = read_name.rfind(b'__')
                if i != -1:
                    read_name = read_name[:i]
                elif not read_name_warning_raised:
                    warning(f'Read {read_name} does not have __ in its name')
                    read_name_warning_raised = True

                line_w = b'\t'.join((read_name, *line_fields[1:])) + b'\n'
            else:
                line_w = line

            stdin2.write(line_w)

    stdin2.close()


    p2.wait()
    fi.close()
    fo.close()


    return read_lens


def get_read_lens_from_bam(sam_file):
    read_lens = []
    for row in sam_file.fetch():
        read_lens.append(row.qlen)

    return read_lens
