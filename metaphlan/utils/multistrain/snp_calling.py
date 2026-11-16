from collections import Counter, OrderedDict
from dataclasses import dataclass

import numpy as np

DTYPE_ERROR_RATES = np.float64

DTPYE_BASE_FREQUENCY = np.uint16
import numpy.typing as npt
import pandas as pd
import scipy.stats as sps
import statsmodels.stats.multitest as smsm
from typing import Sequence, Any

from ...utils import info, info_debug, warning

ACTGactg = 'ACTGactg'
ACTG = 'ACTG'
ACTG_to_i = dict(zip(ACTG, range(len(ACTG))))

BAM_FUNMAP = 4
BAM_FSECONDARY = 256
BAM_FQCFAIL = 512
BAM_FDUP = 1024


@dataclass
class PileupResult:
    base_frequencies: dict
    avg_error_rates: dict
    marker_to_length: dict
    df_seq_errors: pd.DataFrame | None
    null_err_rate: float | None


def pu_to_r(pu):
    if pu.alignment.is_paired:
        assert (not pu.alignment.is_read1) or (not pu.alignment.is_read2)
        # from metaphlan we can get pairs but we don't know which is R1 and which R2
        return 'R1' if pu.alignment.is_read1 else 'R2' if pu.alignment.is_read2 else 'UN'
    else:
        return 'UN'



# ==============================================================================================
# Pileup: read-centric and vectorised, one marker at a time.
#
# For each marker, the reads that pass the filters (flags, mapping quality) are collected once and their bases
# are counted per column with numpy: counts, mean base error per column (numpy's pairwise summation, so the
# same values as np.mean over the column), and the read-position error statistics of the read-end screen.
# Every read is used at any depth. (MetaPhlAn's original code used pysam's pileup, whose default max_depth of
# 8000 makes htslib skip reads at deep markers, mostly reverse-strand ones; each sample's log reports how many
# markers that concerns.) The one limit is the uint16 count: a marker whose reads would cover a column more
# than 65,535 times keeps every k-th read (thin_reads), with a warning.
#
# Reads the vectorised path does not model go through a per-column pysam loop for that marker
# (pileup_marker_columns): reads flagged as proper pairs or R1/R2 (htslib changes the qualities of overlapping
# mates; the read-position statistics are kept per R1/R2), CIGARs other than a single M (indels and clips,
# only with filter_indels = false), reads without SEQ or QUAL, and reads extending past the marker end.
# MetaPhlAn's bowtie2 output with the default configuration never needs it.
# ==============================================================================================

# 10 ** (-q / 10) for every base quality
ERR_TABLE = 10 ** (-np.arange(256) / 10)

# ASCII -> index in ACTG order, -1 for anything else (N and IUPAC codes are not counted)
_CODE = np.full(256, -1, dtype=np.int8)
for _i, _b in enumerate(ACTG):
    _CODE[ord(_b)] = _i
    _CODE[ord(_b.lower())] = _i

# pysam's default max_depth, used by MetaPhlAn's original sample2markers code (only reported in the log)
HTSLIB_DEFAULT_MAX_DEPTH = 8000
# the largest count DTPYE_BASE_FREQUENCY holds; also the max_depth given to pysam in the per-column loop
PILEUP_MAX_DEPTH = int(np.iinfo(DTPYE_BASE_FREQUENCY).max)

read_types = ['R1', 'R2', 'UN']


class NeedsColumnLoop(Exception):
    pass


# ----------------------------------------------------------------------------------------------
# numpy's pairwise summation (add.reduce on a contiguous float64 array), vectorised over segments
# ----------------------------------------------------------------------------------------------
def segment_pairwise_sum(v, starts, lengths):
    """Sum v[s:s+l] for every segment exactly like np.add.reduce would (same rounding). numpy reduces a contiguous
    array in blocks of its buffer size (np.getbufsize(), 8192 elements): each block is summed pairwise and the block
    sums are added from left to right."""
    starts = np.asarray(starts, dtype=np.int64)
    lengths = np.asarray(lengths, dtype=np.int64)
    block = np.getbufsize()
    out = _pairwise_sum(v, starts, np.minimum(lengths, block))
    big = lengths > block
    if big.any():
        s, l = starts[big], lengths[big]
        acc = out[big]
        i = block
        while True:
            m = l > i
            if not m.any():
                break
            acc[m] = acc[m] + _pairwise_sum(v, s[m] + i, np.minimum(l[m] - i, block))
            i += block
        out[big] = acc
    return out


def _pairwise_sum(v, starts, lengths):
    """numpy's pairwise summation of one block, vectorised over segments."""
    out = np.zeros(len(starts), dtype=np.float64)
    if len(starts) == 0:
        return out

    small = lengths < 8
    if small.any():
        s, l = starts[small], lengths[small]
        res = np.zeros(len(s), dtype=np.float64)
        for j in range(7):
            m = l > j
            res[m] += v[s[m] + j]
        out[small] = res

    mid = (lengths >= 8) & (lengths <= 128)
    if mid.any():
        s, l = starts[mid], lengths[mid]
        r = [v[s + j].copy() for j in range(8)]
        lim = l - (l % 8)
        i = 8
        while True:
            m = i < lim
            if not m.any():
                break
            for j in range(8):
                r[j][m] += v[s[m] + i + j]
            i += 8
        res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]))
        rem = l % 8
        for k in range(7):
            m = k < rem
            res[m] += v[s[m] + lim[m] + k]
        out[mid] = res

    large = lengths > 128
    if large.any():
        s, l = starts[large], lengths[large]
        n2 = l // 2
        n2 -= n2 % 8
        out[large] = _pairwise_sum(v, s, n2) + _pairwise_sum(v, s + n2, l - n2)

    return out


# ----------------------------------------------------------------------------------------------
# per-marker read collection
# ----------------------------------------------------------------------------------------------
def collect_marker_reads(sam_file, marker, length, flag_filter, min_mapq):
    """Reads that pysam's pileup would use for this contig, in BAM order."""
    starts, lens, revs, seqs, quals = [], [], [], [], []
    for read in sam_file.fetch(marker):
        flag = read.flag
        if flag & flag_filter:
            continue
        if read.mapping_quality < min_mapq:
            continue
        # htslib only alters base qualities of overlapping mates flagged as proper pairs (0x2), and the
        # read-position statistics are split by R1 (0x40) / R2 (0x80). MetaPhlAn's unpaired bowtie2 alignments
        # carry neither flag after fixmate (only 0x1), so all reads count as 'UN' and qualities are untouched.
        # Anything else goes to the per-column loop.
        if flag & (2 | 64 | 128):
            raise NeedsColumnLoop('proper pair or R1/R2 flag')
        ct = read.cigartuples
        if ct is None or len(ct) != 1 or ct[0][0] != 0:
            raise NeedsColumnLoop('cigar')
        q = read.query_qualities
        s = read.query_sequence
        if q is None or s is None or len(s) != ct[0][1]:
            raise NeedsColumnLoop('seq')
        if read.reference_start + len(s) > length:
            raise NeedsColumnLoop('read past the marker end')
        starts.append(read.reference_start)
        lens.append(len(s))
        revs.append(read.is_reverse)
        seqs.append(s)
        quals.append(q)
    if not starts:
        return None
    lens = np.array(lens, dtype=np.int64)
    starts = np.array(starts, dtype=np.int64)
    revs = np.array(revs, dtype=bool)
    seq = np.frombuffer(''.join(seqs).encode('ascii'), dtype=np.uint8)
    qual = np.concatenate([np.frombuffer(q, dtype=np.uint8) for q in quals])
    return starts, lens, revs, seq, qual


def _reads_covering(starts, ends, x):
    """Number of reads [start, end) that cover each column in x."""
    return np.searchsorted(np.sort(starts), x, side='right') - np.searchsorted(np.sort(ends), x, side='right')


def original_drops_reads(starts, lens, max_depth=HTSLIB_DEFAULT_MAX_DEPTH):
    """True if MetaPhlAn's original pileup (pysam, max_depth 8000) skips reads of this marker. htslib skips a
    read that starts at the same column as the previous one when its buffer already holds max_depth reads (its node
    count, compared with '> max_depth', includes one empty node). While the reads starting at x are added, the
    buffer holds the reads covering x - 1 (the others were freed when the columns up to x - 1 were processed) plus
    the reads already added at x. (Plain coverage is not enough: two stacks of 4,000 and 4,100 reads, one starting
    where the other ends, fill the buffer to 8,100.) Only reported in the log."""
    if len(starts) <= max_depth:
        return False
    x, n_x = np.unique(starts, return_counts=True)
    held = _reads_covering(starts, starts + lens, x - 1) + n_x - 1  # buffer before the last read at x is added
    return bool(((n_x >= 2) & (held >= max_depth)).any())


def max_read_coverage(starts, lens):
    if len(starts) == 0:
        return 0
    return int(_reads_covering(starts, starts + lens, np.unique(starts)).max())


def thin_reads(reads, max_depth=PILEUP_MAX_DEPTH):
    """Keeps every k-th read (BAM order), with the smallest k that brings the read coverage of every column to at
    most max_depth, so that the base counts fit DTPYE_BASE_FREQUENCY. Returns (reads, k); k = 1 for any marker
    with at most max_depth reads."""
    starts, lens, revs, seq, qual = reads
    if len(starts) <= max_depth:
        return reads, 1
    k = 1
    while max_read_coverage(starts[::k], lens[::k]) > max_depth:
        k += 1
    if k == 1:
        return reads, 1
    keep = np.zeros(len(starts), dtype=bool)
    keep[::k] = True
    bases = np.repeat(keep, lens)
    return (starts[keep], lens[keep], revs[keep], seq[bases], qual[bases]), k


class _Acc:
    """Read-position accumulators, kept in first-insertion order (the order of the Counters they become).
    Read position p is stored at index p + off (the per-column loop gives negative positions for reverse reads
    with soft clips)."""

    def __init__(self, max_len=1024):
        self.off = 0
        self.all = np.zeros(max_len, dtype=np.int64)
        self.minor = np.zeros(max_len, dtype=np.int64)
        self.err = np.zeros(max_len, dtype=np.float64)
        self.order_all = []
        self.seen_all = set()
        self.order_minor = []
        self.seen_minor = set()

    def ensure(self, lo, hi):
        """Make read positions lo..hi addressable."""
        lo_i, hi_i = lo + self.off, hi + self.off
        if lo_i >= 0 and hi_i < len(self.all):
            return
        pad = max(0, -lo_i)
        k = max(hi_i + pad + 1, len(self.all) + pad, 2 * len(self.all))
        for name in ('all', 'minor', 'err'):
            a = getattr(self, name)
            b = np.zeros(k, dtype=a.dtype)
            b[pad:pad + len(a)] = a
            setattr(self, name, b)
        self.off += pad

    def note_order(self, rp, is_minor):
        # first-appearance order of keys, as Counter insertion order would give
        _, first = np.unique(rp, return_index=True)
        for p in rp[np.sort(first)]:
            p = int(p)
            if p not in self.seen_all:
                self.seen_all.add(p)
                self.order_all.append(p)
        rpm = rp[is_minor]
        if len(rpm):
            _, first = np.unique(rpm, return_index=True)
            for p in rpm[np.sort(first)]:
                p = int(p)
                if p not in self.seen_minor:
                    self.seen_minor.add(p)
                    self.order_minor.append(p)

    def counters(self):
        o = self.off
        all_c = Counter({p: int(self.all[p + o]) for p in self.order_all})
        minor_c = Counter({p: int(self.minor[p + o]) for p in self.order_minor})
        err_c = Counter({p: float(self.err[p + o]) for p in self.order_all})
        return all_c, minor_c, err_c


def pileup_marker_fast(reads, length, config, drop_positions_un, acc):
    """Returns (base_frequencies 4xL uint16, avg_error_rates L float64) or None if nothing counted."""
    starts, lens, revs, seq, qual = reads
    n_reads = len(starts)
    read_idx = np.repeat(np.arange(n_reads), lens)
    offs = np.arange(len(seq)) - np.repeat(np.cumsum(lens) - lens, lens)
    pos = np.repeat(starts, lens) + offs
    code = _CODE[seq]
    keep = (qual >= config['min_base_quality']) & (code >= 0) & (pos >= 0) & (pos < length)
    # read position counted from the 5' end of the sequenced read
    rlen = np.repeat(lens, lens)
    rpos = np.where(np.repeat(revs, lens), rlen - offs - 1, offs)
    if drop_positions_un:
        keep &= ~np.isin(rpos, np.fromiter(drop_positions_un, dtype=np.int64))
    if not keep.any():
        return None

    pos, code, q, rpos, read_idx = pos[keep], code[keep], qual[keep], rpos[keep], read_idx[keep]
    order = np.argsort(pos, kind='stable')  # column-major, BAM read order within a column
    pos, code, q, rpos = pos[order], code[order], q[order], rpos[order]

    cols, col_start, col_n = np.unique(pos, return_index=True, return_counts=True)
    if config['min_base_coverage'] > 1:
        ok = col_n >= config['min_base_coverage']
        if not ok.all():
            sel = np.repeat(ok, col_n)
            pos, code, q, rpos = pos[sel], code[sel], q[sel], rpos[sel]
            cols, col_start, col_n = np.unique(pos, return_index=True, return_counts=True)
            if len(cols) == 0:
                return None

    bf = np.zeros((4, length), dtype=DTPYE_BASE_FREQUENCY)
    cnt = np.bincount(pos * 4 + code, minlength=4 * length).reshape(length, 4).T
    bf[:, :] = cnt
    err = ERR_TABLE[q]
    avg = np.zeros(length, dtype=DTYPE_ERROR_RATES)
    avg[cols] = segment_pairwise_sum(err, col_start, col_n) / col_n

    # read-position error statistics at sufficiently covered columns
    stat_cols = col_n >= config['seq_error_stats_min_cov']
    if stat_cols.any():
        col_id = np.repeat(np.arange(len(cols)), col_n)
        in_stat = stat_cols[col_id]
        # majority base, ties broken by first appearance in the column (Counter.most_common)
        counts_c = cnt[:, cols].T  # (n_cols, 4)
        within = np.arange(len(pos)) - np.repeat(col_start, col_n)
        first = np.full((len(cols), 4), np.iinfo(np.int64).max, dtype=np.int64)
        np.minimum.at(first, (col_id, code), within)
        mx = counts_c.max(axis=1, keepdims=True)
        cand = np.where(counts_c == mx, first, np.iinfo(np.int64).max)
        major = cand.argmin(axis=1)
        is_minor = code != major[col_id]

        rp = rpos[in_stat]
        mn = is_minor[in_stat]
        e = err[in_stat]
        acc.ensure(int(rp.min()), int(rp.max()))
        acc.note_order(rp, mn)
        idx = rp + acc.off
        np.add.at(acc.all, idx, 1)
        np.add.at(acc.minor, idx[mn], 1)
        np.add.at(acc.err, idx, e)  # sequential, in column-major / read order

    return bf, avg



# ----------------------------------------------------------------------------------------------
# per-column pysam loop, for markers with reads the vectorised path does not model
# ----------------------------------------------------------------------------------------------
def pileup_marker_columns(sam_file, marker, length, config, flag_filter, drop_positions, acc_by_r):
    """Counts, mean base errors and read-position statistics of one marker from pysam's pileup columns."""
    bfm, aer = None, None
    for base_pileup in sam_file.pileup(contig=marker, stepper='samtools', min_base_quality=config['min_base_quality'],
                                       min_mapping_quality=config['min_mapping_quality'], ignore_orphans=False,
                                       flag_filter=flag_filter, max_depth=PILEUP_MAX_DEPTH):
        pos = base_pileup.reference_pos  # 0 indexed
        if pos >= length:
            continue

        query_sequences = base_pileup.get_query_sequences()
        query_qualities = base_pileup.get_query_qualities()
        query_positions = base_pileup.get_query_positions()
        rs = [pu_to_r(pu) for pu in base_pileup.pileups]

        # positions from the start of the read (if base is lower-case, the read had been reverse-complemented and the
        #   position has to be counted from the end)
        bqpr = [(b.upper(), q, p if b.isupper() else pu.alignment.qlen - p - 1, r)
                for b, q, p, pu, r in zip(query_sequences, query_qualities, query_positions, base_pileup.pileups, rs)
                if len(b) == 1 and b in ACTGactg]  # a deletion gives '' (and '' in 'ACTGactg' is True)

        if drop_positions is not None:
            bqpr = [(b, q, p, r) for b, q, p, r in bqpr if (r, p) not in drop_positions]

        if not bqpr:  # not covered position
            continue

        bases, qualities, positions, rs = zip(*bqpr)
        base_coverage = len(bases)

        if base_coverage < config['min_base_coverage']:
            continue

        base_frequency = Counter(bases)

        qualities = np.array(qualities)
        err_rates = 10 ** (-qualities / 10)
        avg_error_rate = np.mean(err_rates)

        if bfm is None:
            bfm = np.zeros((4, length), DTPYE_BASE_FREQUENCY)
            aer = np.zeros(length, DTYPE_ERROR_RATES)

        for b in base_frequency.keys():
            bfm[ACTG_to_i[b], pos] = base_frequency[b]
        aer[pos] = avg_error_rate

        if base_coverage >= config['seq_error_stats_min_cov']:
            major_base = base_frequency.most_common(1)[0][0]
            for b, p, e, r in zip(bases, positions, err_rates, rs):
                a = acc_by_r[r]
                a.ensure(p, p)
                if p not in a.seen_all:
                    a.seen_all.add(p)
                    a.order_all.append(p)
                a.all[p + a.off] += 1
                a.err[p + a.off] += e
                if b != major_base:
                    if p not in a.seen_minor:
                        a.seen_minor.add(p)
                        a.order_minor.append(p)
                    a.minor[p + a.off] += 1
    if bfm is None:
        return None
    return bfm, aer


# ----------------------------------------------------------------------------------------------
# pileup of a whole sample
# ----------------------------------------------------------------------------------------------
def run_pileup_inner(sam_file, config, drop_read_positions=None, stats=None):
    """Base counts, mean base errors per column and read-position error statistics of every marker in the BAM.
    drop_read_positions: (read type, read position) pairs whose bases are not counted (read-end screen)."""
    all_markers = sam_file.references
    marker_lengths = sam_file.lengths
    marker_to_length = dict(zip(all_markers, marker_lengths))

    flag_filter = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL
    if config['filter_duplicates']:
        flag_filter |= BAM_FDUP

    drop = set(drop_read_positions) if drop_read_positions is not None else None
    drop_un = {p for r, p in drop if r == 'UN'} if drop else set()
    drop_paired = bool(drop) and any(r != 'UN' for r, _ in drop)

    acc_by_r = {r: _Acc() for r in read_types}
    base_frequencies = {}
    avg_error_rates = {}
    n_column_loop = 0
    n_original_drops = 0
    thinned = []

    for marker, length in zip(all_markers, marker_lengths):
        res = None
        try:
            if drop_paired:
                raise NeedsColumnLoop('drop set with paired read types')
            reads = collect_marker_reads(sam_file, marker, length, flag_filter, config['min_mapping_quality'])
            if reads is not None:
                n_original_drops += original_drops_reads(reads[0], reads[1])
                reads, k = thin_reads(reads)
                if k > 1:
                    thinned.append(f'{marker} (1 in {k})')
                res = pileup_marker_fast(reads, length, config, drop_un, acc_by_r['UN'])
        except NeedsColumnLoop:
            n_column_loop += 1
            res = pileup_marker_columns(sam_file, marker, length, config, flag_filter, drop, acc_by_r)
        if res is not None:
            base_frequencies[marker], avg_error_rates[marker] = res

    if drop_read_positions is None:  # the first of the two passes
        if n_original_drops:
            info(f'{n_original_drops} markers have reads that MetaPhlAn\'s original sample2markers skips (pysam '
                 f'pileup, max_depth {HTSLIB_DEFAULT_MAX_DEPTH}); all reads are used')
        if thinned:
            warning(f'Reads used: every k-th read for markers covered by more than {PILEUP_MAX_DEPTH} reads, the '
                    f'largest count stored: {", ".join(thinned)}')
    if stats is not None:
        stats['n_markers_column_loop'] = stats.get('n_markers_column_loop', 0) + n_column_loop
        stats['n_markers_original_drops'] = stats.get('n_markers_original_drops', 0) + n_original_drops
        stats['n_markers_thinned'] = stats.get('n_markers_thinned', 0) + len(thinned)

    minor_positions, all_positions, per_position_eror_rates_sum = {}, {}, {}
    for r in read_types:
        all_positions[r], minor_positions[r], per_position_eror_rates_sum[r] = acc_by_r[r].counters()

    minor_positions_tot = sum(minor_positions.values(), Counter())
    all_positions_tot = sum(all_positions.values(), Counter())

    if all_positions_tot.total() > 0:
        a = config['null_err_rate_pseudocount']
        null_err_rate = (minor_positions_tot.total() + a) / (all_positions_tot.total() + a)
        info_debug(f'Null err. rate: {null_err_rate}')

        if null_err_rate > config['max_null_err_rate']:
            if config['max_null_err_rate'] > 0:
                warning(f'Null error rate is higher than expected ({null_err_rate}), will not use it')
            null_err_rate = None

        dfs = []
        for r in read_types:
            s_minor_positions = pd.Series(minor_positions[r], name='minor_bases')
            s_all_positions = pd.Series(all_positions[r], name='all_bases')
            s_per_position_eror_rates_sum = pd.Series(per_position_eror_rates_sum[r], name='sum_error_rates')
            df_seq_errors_r = pd.concat([s_minor_positions, s_all_positions, s_per_position_eror_rates_sum], axis=1) \
                .rename_axis('read_position').sort_index()
            df_seq_errors_r['r'] = r
            if len(df_seq_errors_r) > 0:
                dfs.append(df_seq_errors_r)

        if len(dfs) > 0:
            df_seq_errors = pd.concat(dfs).reset_index().set_index(['r', 'read_position'])
            a = config['null_err_rate_excess_pseudocount']
            df_seq_errors['null_err_rate_excess'] = (df_seq_errors['minor_bases'] + a) / (df_seq_errors['sum_error_rates'] + a)
        else:
            warning(f'No positions to get seq. errors statistics')
            df_seq_errors = None
    else:
        info('No positions to establish null error rate, will use default')
        df_seq_errors = None
        null_err_rate = None

    return PileupResult(base_frequencies, avg_error_rates, marker_to_length, df_seq_errors, null_err_rate)


def run_pileup(sam_file, config):
    """
    Pileup before and after the read-end screen (read positions with excess errors are left out of the second)
    :param sam_file: pysam opened file
    :param dict config:
    :return: (PileupResult before, PileupResult after)
    """
    pr_before = run_pileup_inner(sam_file, config)

    pr_after = pr_before
    if config['discard_read_ends'] and pr_before.df_seq_errors is not None:
        df_se = pr_before.df_seq_errors.query(f'all_bases>{config["discard_read_ends_support"]}')
        if len(df_se) > 0:
            c = np.log(df_se['null_err_rate_excess'].dropna())
            q25 = c.quantile(0.25)
            q75 = c.quantile(0.75)
            iqr = q75 - q25
            mx = q75 + 1.5 * iqr
            m = c > mx
            drop_read_positions = c.index[m]

            if len(drop_read_positions) > 0:
                warning(f'Will not use the following positions of the reads, '
                        f'there might be adapter contamination '
                        f'or other issue with sequencing: {list(drop_read_positions)}')
                pr_after = run_pileup_inner(sam_file, config, set(drop_read_positions))
        else:
            info_debug('Not dropping any read positions, not enough support')

    return pr_before, pr_after


def mask_iqr(x, q25, q75):
    if q25 is not None and q75 is not None:
        iqr = q75 - q25
        return (x > q75 + 1.5 * iqr) | (x < q25 - 1.5 * iqr)
    else:
        return np.zeros(x.shape, dtype=bool)


def fit_nb(q50, avg_read_len_i, p_t):
    if q50 is None:
        return None, None

    q50s = np.arange(avg_read_len_i + 1) / avg_read_len_i * q50  # taper the median linearly at the ends of the markers

    # convert median to mean assuming poisson dist.
    a = 50
    b = 50 / 3 - 50 * q50s
    c = -1
    ks = (-b + np.sqrt(b ** 2 - 4 * a * c)) / 2 / a

    # parametrize z (overdispersion, i.e. variance:mean ratio) with coverage
    mean = ks[-1]
    cv = 0.2031 # coefficient of variation
    z = cv**2 * mean + 1.4587 # for mean = 1, z = 1.5

    # get Negative binomial parameters, z is the variance:mean ratio
    rs = ks / (z - 1)
    ps = rs / (rs + ks)
    d_nb = sps.nbinom(rs, ps)

    min_cs, max_cs = d_nb.isf(1 - p_t), d_nb.isf(p_t) - 1

    return min_cs, max_cs


def mask_nb(bc, min_cs, max_cs):
    """

    :param npt.NDArray[int] bc:
    :param npt.NDArray[int] min_cs:
    :param npt.NDArray[int] max_cs:
    :return:
    """
    mask = np.ones(len(bc), dtype=bool)
    if min_cs is not None:
        crop = min(len(min_cs) - 1, len(bc))
        mask[:crop] = (bc[:crop] >= min_cs[:crop]) & (bc[:crop] <= max_cs[:crop])
        mask[-crop:] = (bc[-crop:] >= min_cs[:crop:][::-1]) & (bc[-crop:] <= max_cs[:crop][::-1])
        mask[crop:-crop] = (bc[crop:-crop] >= min_cs[-1]) & (bc[crop:-crop] <= max_cs[-1])

    return mask



def filter_loci_snp_call(bfs, err_rates, avg_read_len, config):
    """

    :param dict[str, npt.NDArray[int]] bfs:
    :param dict[str, npt.NDArray[float]] err_rates:
    :param float avg_read_len:
    :param dict config:
    :return:
    """
    result_row = {}
    result_row['n_markers_present'] = len(bfs)

    base_coverages = {m: bf.sum(axis=0) for m, bf in bfs.items()}
    result_row['n_positions_total'] = sum((np.count_nonzero(bc) for bc in base_coverages.values()))
    
    avg_read_len_i = int(avg_read_len)
    assert avg_read_len_i > 0

    base_coverages_fit = {m: bc[avg_read_len_i:-avg_read_len_i] for m, bc in base_coverages.items()
                          if len(bc) > 2 * avg_read_len_i}

    q25, q50, q75 = None, None, None
    if len(base_coverages_fit) > 0:
        base_coverage_fit = np.concatenate(list(base_coverages_fit.values()))
        base_coverage_fit = base_coverage_fit[base_coverage_fit > 0]
        if len(base_coverage_fit) > 0:
            q25 = np.quantile(base_coverage_fit, .25)
            q50 = np.median(base_coverage_fit)
            q75 = np.quantile(base_coverage_fit, .75)


    min_cs, max_cs = fit_nb(q50, avg_read_len_i, config['outlier_coverage_p'])
    position_mask = {m: mask_nb(bc, min_cs, max_cs) for m, bc in base_coverages.items()}
    outlier_coverage = sum(np.count_nonzero(~x) for x in position_mask.values()) /\
                       sum(len(x) for x in position_mask.values())
    marker_to_outlier_coverage = {m: (~x).mean() for m, x in position_mask.items()}
    markers_outlier_covered = set((m for m, oc in marker_to_outlier_coverage.items()
                                   if oc > config['outlier_coverage_marker']))
    # mask all positions of outlier covered markers
    for m in markers_outlier_covered:
        position_mask[m][:] = False


    result_row['coverage_25th'] = q25
    result_row['coverage_50th'] = q50
    result_row['coverage_75th'] = q75
    result_row['perc_positions_outlier_covered'] = outlier_coverage
    result_row['n_markers_outlier_covered'] = len(markers_outlier_covered)


    if config['trim_marker_ends'] > 0:
        t = config['trim_marker_ends']
        for m, pm in position_mask.items():
            pm[:t] = False
            pm[-t:] = False

    marker_covered_positions = [np.count_nonzero(bc[position_mask[m]]) for m, bc in base_coverages.items()]

    result_row['n_positions_retained'] = sum(marker_covered_positions)
    result_row['n_markers_retained'] = result_row['n_markers_present'] - result_row['n_markers_outlier_covered']

    allelism = {m: (bf > 0).sum(axis=0) for m, bf in bfs.items()}
    allelism_filtered = {m: a * position_mask[m] for m, a in allelism.items()}
    result_row['n_positions_biallelic'] = sum((np.count_nonzero(a == 2) for a in allelism_filtered.values()))
    result_row['n_positions_triallelic'] = sum((np.count_nonzero(a == 3) for a in allelism_filtered.values()))
    result_row['n_positions_quadallelic'] = sum((np.count_nonzero(a == 4) for a in allelism_filtered.values()))


    # Polyallelism test with FDR
    m_to_pvals = OrderedDict()
    m_to_poly_mask = {}
    for m in base_coverages.keys():
        polyallelic_mask = allelism_filtered[m] > 1
        if polyallelic_mask.sum() > 0:
            bf_poly = bfs[m][:, polyallelic_mask]
            bc_poly = base_coverages[m][polyallelic_mask]
            er_poly = err_rates[m][polyallelic_mask]
            poly_max_f = bf_poly.max(axis=0)
            p_values = sps.binom.cdf(poly_max_f, bc_poly, 1 - er_poly)
            m_to_pvals[m] = p_values
        m_to_poly_mask[m] = polyallelic_mask

    m_to_qvals = {}
    if len(m_to_pvals) > 0:
        all_pvalues = np.concatenate(list(m_to_pvals.values()))
        _, all_qvalues, _, _ = smsm.multipletests(all_pvalues, method="fdr_bh")

        # split concatenated array to per-marker arrays
        spacer = 0
        for m, ps in m_to_pvals.items():
            m_to_qvals[m] = all_qvalues[spacer: spacer+len(ps)]
            spacer += len(ps)
        assert spacer == len(all_qvalues)

    polyallelic_significant_masks = {}
    for m in base_coverages.keys():
        polyallelic_mask = m_to_poly_mask[m]
        mask_sig_polyallelic_full = np.zeros(len(base_coverages[m]), dtype=bool)
        if m in m_to_qvals:
            q_vals = m_to_qvals[m]
            mask_sig_polyallelic = q_vals < config['polymorphism_q_alpha']
            mask_sig_polyallelic_full[polyallelic_mask] = mask_sig_polyallelic
        polyallelic_significant_masks[m] = mask_sig_polyallelic_full


    # Polyallelism test without FDR (just p-values)
    # polyallelic_significant_masks = {}
    # for m in base_coverages.keys():
    #     polyallelic_mask = allelism_filtered[m] > 1
    #     mask_sig_polyallelic_full = np.zeros(len(base_coverages[m]), dtype=bool)
    #     if polyallelic_mask.sum() > 0:
    #         bf_poly = bfs[m][:, polyallelic_mask]
    #         bc_poly = base_coverages[m][polyallelic_mask]
    #         er_poly = err_rates[m][polyallelic_mask]
    #         poly_max_f = bf_poly.max(axis=0)
    #         p_values = sps.binom.cdf(poly_max_f, bc_poly, 1 - er_poly)
    #         mask_sig_polyallelic = p_values < config['polymorphism_p_alpha']
    #         mask_sig_polyallelic_full[polyallelic_mask] = mask_sig_polyallelic
    #     polyallelic_significant_masks[m] = mask_sig_polyallelic_full

    result_row['n_positions_significantly_polyallelic'] = sum((psm.sum() for psm in polyallelic_significant_masks.values()))
    result_row['n_markers_significantly_polyallelic'] = sum((psm.sum() > 0 for psm in polyallelic_significant_masks.values()))

    return result_row, position_mask, polyallelic_significant_masks
