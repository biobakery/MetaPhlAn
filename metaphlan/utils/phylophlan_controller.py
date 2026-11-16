__author__ = 'Aitor Blanco Miguez (aitor.blancomiguez@unitn.it),' \
             'Michal Puncochar (michal.puncochar@unitn.it)'
__version__ = '4.2.4'
__date__    = '21 Oct 2025'

import abc
import bz2
import configparser
import glob
import os
import shlex
import shutil
import subprocess as sb
from shutil import move

try:
    from .util_fun import info, error, warning
    from .external_exec import generate_phylophlan_config_file, execute_treeshrink
except ImportError:
    from util_fun import info, error, warning
    from external_exec import generate_phylophlan_config_file, execute_treeshrink


class PhylophlanController(abc.ABC):
    """PhyloPhlAnController interface class"""

    @abc.abstractmethod
    def compute_phylogeny(self, *args, **kwargs):
        """Executes PhyloPhlAn to compute phylogeny"""
        pass


    @staticmethod
    def get_phylophlan_configuration():
        """Gets PhyloPhlAn configuration

        Returns:
            dict: The dictionary with the PhyloPhlAn configuration
        """
        configuration = {
            'map': 'blastn',
            'aligner': 'mafft',
            'tree1': 'raxml'
        }
        return configuration


class Phylophlan3Controller(PhylophlanController):
    """PhyloPhlAn3Controller class"""

    # PhyloPhlAn options that StrainPhlAn sets and --phylophlan_params must not change (input/output locations,
    # StrainPhlAn mode, configuration file: use --phylophlan_configuration)
    PROTECTED_OPTIONS = {'-i', '--input', '-c', '--clean', '-o', '--output', '--output_folder', '--data_folder',
                         '--input_folder', '-f', '--config_file', '--strainphlan', '-t', '--db_type',
                         '--genome_extension'}
    # PhyloPhlAn options with a value that StrainPhlAn sets and --phylophlan_params may replace (the last value wins)
    FIXED_VALUES = {'--diversity': 'low', '--min_num_entries': '1', '--min_num_markers': '1',
                    '--fragmentary_threshold': '1.0', '--not_variant_threshold': '1.0', '--gap_perc_threshold': '0.8'}


    def compute_phylogeny(self, samples_markers_dir, num_samples, tmp_dir):
        """Executes PhyloPhlAn to compute phylogeny

        Args:
            samples_markers_dir (str): Directory containing the samples as FASTA
            num_samples (int): Minimum number of samples containing a marker
            tmp_dir (str): Temporal directory
        """
        n_files = self.n_to_gaps(samples_markers_dir)
        if n_files:
            info('\tN replaced by gaps in the marker sequences of {} samples'.format(n_files))
        default_configuration = not self.phylophlan_configuration
        if default_configuration:
            info("\tGenerating PhyloPhlAn configuration file...")
            self.phylophlan_configuration = generate_phylophlan_config_file(tmp_dir, self.get_phylophlan_configuration())
            info("\tDone.")
        info("\tExecuting PhyloPhlAn...")
        r = self.execute_phylophlan(samples_markers_dir, tmp_dir)
        tree = self.final_tree()
        if r.returncode != 0 or tree is None:
            if not default_configuration:
                error('PhyloPhlAn failed (exit status {}) with the configuration {}:\n{}'.format(
                    r.returncode, self.phylophlan_configuration, r.stderr.decode(errors='replace')), exit=True)
            tree = self.fallback_tree(tmp_dir)
        elif default_configuration:
            self.tree_note = 'RAxML on the concatenated alignment'
        else:
            self.tree_note = 'PhyloPhlAn with {}: {}'.format(self.phylophlan_configuration, os.path.basename(tree))
        if self.mutation_rates and os.path.exists(os.path.join(self.output_dir, "mutation_rates.tsv")):
            move(os.path.join(self.output_dir, "mutation_rates.tsv"), os.path.join(
                self.output_dir, "{}.mutation".format(self.clade)))
            move(os.path.join(self.output_dir, "mutation_rates"), os.path.join(
                self.output_dir, "{}_mutation_rates".format(self.clade)))
        info("\tDone.")
        if self.treeshrink:
            info("\tExecuting TreeShrink...")
            execute_treeshrink(tree, self.output_dir, tmp=tmp_dir, centroid=True if num_samples >= 100 else False)
            info("\tDone.")
            info('\tThis StrainPhlAn analysis ran TreeShrink, do not forget to cite:')
            info('\tMai, Uyen, and Siavash Mirarab. 2018. “TreeShrink: Fast and Accurate Detection of Outlier Long'
                 ' Branches in Collections of Phylogenetic Trees.” BMC Genomics 19 (S5): 272.'
                 ' https://doi.org/10.1186/s12864-018-4620-2.')


    def execute_phylophlan(self, samples_markers_dir, tmp_dir):
        """
        Executes PhyloPhlAn

        Args:
            samples_markers_dir:
            tmp_dir:

        Returns:
            subprocess.CompletedProcess: the finished PhyloPhlAn process
        """
        cmd = f'phylophlan' \
              f' -i {samples_markers_dir} -o . --output_folder {self.output_dir} --nproc {self.nprocs}' \
              f' --strainphlan --{self.phylophlan_mode} --data_folder {tmp_dir} -f {self.phylophlan_configuration}' \
              f' -t n --genome_extension fna --sort' \
              + ''.join(' {} {}'.format(k, v) for k, v in self.FIXED_VALUES.items())

        if self.phylophlan_params:
            cmd += " " + self.phylophlan_params
        if self.mutation_rates:
            cmd += " --mutation_rates"
        info('\tPhyloPhlAn command: {}'.format(cmd))

        r = sb.run(shlex.split(cmd), capture_output=True)
        with open(os.path.join(tmp_dir, "phylophlan_log.stdout"), 'wb') as f:
            f.write(r.stdout)
        with open(os.path.join(tmp_dir, "phylophlan_log.stderr"), 'wb') as f:
            f.write(r.stderr)
        return r


    @staticmethod
    def n_to_gaps(folder):
        """Writes N (and n) in the marker sequences of the per-sample FASTA files as gaps. PhyloPhlAn leaves out a
        sample whose sequence is only gaps, but not one that is only N, on which RAxML stops; N and gaps are both
        missing data for the tree. Returns the number of files changed."""
        n_files = 0
        for f in sorted(os.listdir(folder)):
            path = os.path.join(folder, f)
            if not os.path.isfile(path):
                continue
            opener = bz2.open if f.endswith('.bz2') else open
            with opener(path, 'rt') as fh:
                lines = fh.read().split('\n')
            new = [x if x.startswith('>') else x.replace('N', '-').replace('n', '-') for x in lines]
            if new != lines:
                with opener(path, 'wt') as fh:
                    fh.write('\n'.join(new))
                n_files += 1
        return n_files


    def final_tree(self):
        """The tree PhyloPhlAn built: the refined one ([tree2]) if any, else RAxML's, IQ-TREE's or the plain [tree1]
        output (FastTree, ASTRAL, ASTRID). None if there is none."""
        p = '{}.StrainPhlAn4'.format(self.clade)
        for name in ('RAxML_bestTree.{}_refined.tre', '{}_refined.tre', 'RAxML_bestTree.{}.tre', '{}.tre.treefile',
                     '{}.tre'):
            path = os.path.join(self.output_dir, name.format(p))
            if os.path.isfile(path) and os.path.getsize(path) > 0:
                return path
        return None


    def check_phylophlan_params(self, phylophlan_mode, phylophlan_params):
        """The precision mode, and --phylophlan_params without --fast/--accurate. A mode in --phylophlan_params is used
        if --phylophlan_mode is not given (they must agree if both are); options that set input/output locations or
        the configuration file stop the run; replaced StrainPhlAn defaults are logged."""
        tokens = shlex.split(phylophlan_params or '')
        modes = {t[2:] for t in tokens if t in ('--fast', '--accurate')}
        if len(modes) > 1:
            error('--phylophlan_params has both --fast and --accurate', exit=True)
        if modes and phylophlan_mode is not None and phylophlan_mode not in modes:
            error('--phylophlan_mode {} and --{} in --phylophlan_params disagree'.format(phylophlan_mode,
                                                                                        modes.pop()), exit=True)
        mode = phylophlan_mode or (modes.pop() if modes else 'fast')
        tokens = [t for t in tokens if t not in ('--fast', '--accurate')]
        protected = sorted({t.split('=')[0] for t in tokens} & self.PROTECTED_OPTIONS)
        if protected:
            error('--phylophlan_params cannot set {} (StrainPhlAn sets them; use --phylophlan_configuration for the '
                  'PhyloPhlAn configuration file)'.format(', '.join(protected)), exit=True)
        for i, t in enumerate(tokens):
            k = t.split('=')[0]
            if k in self.FIXED_VALUES:
                v = t.split('=', 1)[1] if '=' in t else (tokens[i + 1] if i + 1 < len(tokens) else '')
                info('PhyloPhlAn option {} {} replaced by --phylophlan_params: {}'.format(k, self.FIXED_VALUES[k], v))
        return mode, ' '.join(shlex.quote(t) for t in tokens) or None


    # ------------------------------------------------------------------------------------------
    # Fallback when RAxML fails on the concatenated alignment (default configuration only).
    # RAxML stops on a sequence that is entirely gaps / N (a sample without data at the kept positions) and on a
    # DNA alignment in which a base never occurs (very few positions left after trimming). First RAxML again
    # without the empty sequences (if at least 4 samples remain), then FastTree on the same alignment.
    # ------------------------------------------------------------------------------------------
    @staticmethod
    def _phylophlan_error(tmp_dir):
        """The last error line ('[e] ...') PhyloPhlAn wrote."""
        path = os.path.join(tmp_dir, 'phylophlan_log.stderr')
        last = None
        if os.path.isfile(path):
            with open(path, errors='replace') as f:
                for line in f:
                    if line.startswith('[e]') and 'error while executing' not in line and 'Command' not in line:
                        last = line[3:].strip()
        return last


    @staticmethod
    def _read_fasta(path):
        seqs, name = {}, None
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    name = line[1:].split()[0]
                    seqs[name] = []
                elif name is not None:
                    seqs[name].append(line)
        return {n: ''.join(s) for n, s in seqs.items()}


    @staticmethod
    def _write_fasta(seqs, path):
        with open(path, 'w') as f:
            for n, s in seqs.items():
                f.write('>{}\n{}\n'.format(n, s))


    def _run_raxml(self, aln):
        """The [tree1] RAxML command of the PhyloPhlAn configuration on aln, with PhyloPhlAn's output names.
        Returns (tree path or None, RAxML's error line)."""
        cfg = configparser.ConfigParser(interpolation=None)
        cfg.read(self.phylophlan_configuration)
        name = '{}.StrainPhlAn4.tre'.format(self.clade)
        for old in glob.glob(os.path.join(self.output_dir, 'RAxML_*.' + name)):
            os.remove(old)  # RAxML refuses to overwrite the files of a run with the same name
        cmd = [cfg['tree1']['program_name'], *shlex.split(cfg['tree1']['params']), '-T', str(min(self.nprocs, 20)),
               '-w', os.path.abspath(self.output_dir), '-s', os.path.abspath(aln), '-n', name]
        r = sb.run(cmd, capture_output=True)
        tree = os.path.join(self.output_dir, 'RAxML_bestTree.' + name)
        if r.returncode == 0 and os.path.isfile(tree) and os.path.getsize(tree) > 0:
            return tree, None
        out = r.stdout.decode(errors='replace').splitlines()
        reason = next((x.strip() for x in out if 'ERROR' in x or 'equal to zero' in x), None)
        for old in glob.glob(os.path.join(self.output_dir, 'RAxML_*.' + name)):
            os.remove(old)
        return None, reason or 'exit status {}'.format(r.returncode)


    def fallback_tree(self, tmp_dir):
        """A tree when PhyloPhlAn's RAxML step failed; returns its path. Exits with a message if none can be built."""
        aln = os.path.join(self.output_dir, '{}.StrainPhlAn4_concatenated.aln'.format(self.clade))
        if not os.path.isfile(aln):
            pp_error = self._phylophlan_error(tmp_dir)
            if pp_error and 'No alignments found to concatenate' in pp_error:
                error('No phylogeny can be built: no marker alignment is left after PhyloPhlAn\'s trimming, which '
                      'happens when the samples do not differ at any kept position (e.g. the same strain in all '
                      'samples)', exit=True)
            error('PhyloPhlAn failed before building the tree: {} (log: {}, kept with --debug)'.format(
                pp_error or 'no concatenated alignment', os.path.join(tmp_dir, 'phylophlan_log.stderr')), exit=True)
        seqs = self._read_fasta(aln)
        n_cols = max((len(x) for x in seqs.values()), default=0)
        if n_cols == 0:
            error('No phylogeny can be built: no position is left in the alignment after trimming (the samples do '
                  'not differ at the kept positions)', exit=True)
        tree, reason = self._run_raxml(aln)  # PhyloPhlAn removes RAxML's files: rerun once to get the reason
        if tree is not None:
            self.tree_note = 'RAxML on the concatenated alignment (run again after PhyloPhlAn failed: {})'.format(
                self._phylophlan_error(tmp_dir))
            warning('PhyloPhlAn failed ({}); the tree was built with RAxML directly'.format(
                self._phylophlan_error(tmp_dir)))
            return tree
        warning('RAxML failed on the concatenated alignment ({} samples, {} positions): {}'.format(
            len(seqs), n_cols, reason))

        empty = [n for n, x in seqs.items() if not any(c in 'ACGTacgt' for c in x)]
        kept = {n: x for n, x in seqs.items() if n not in empty}
        if empty:
            warning('{} samples have no base at any of the {} positions and are left out of the tree: {}'.format(
                len(empty), n_cols, ', '.join(empty)))
            if len(kept) >= 4:
                filtered = os.path.join(tmp_dir, 'concatenated_without_empty.aln')
                self._write_fasta(kept, filtered)
                tree, reason2 = self._run_raxml(filtered)
                if tree is not None:
                    self.tree_note = 'RAxML, without {} samples that have no base at the kept positions ({})'.format(
                        len(empty), ', '.join(empty))
                    warning('Tree built with RAxML without those samples: {}'.format(tree))
                    return tree
                warning('RAxML failed again without them: {}'.format(reason2))

        fasttree = shutil.which('FastTreeMP') or shutil.which('FastTree')
        source = kept if empty else seqs
        if fasttree is None or len(source) < 3:
            error('No phylogeny could be built: RAxML failed'
                  + ('' if fasttree else ' and FastTree is not installed')
                  + (' and only {} samples have data'.format(len(source)) if len(source) < 3 else ''), exit=True)
        src = os.path.join(tmp_dir, 'concatenated_for_fasttree.aln')
        self._write_fasta(source, src)
        tree = os.path.join(self.output_dir, '{}.StrainPhlAn4.FastTree.tre'.format(self.clade))
        with open(tree, 'w') as fo:
            r = sb.run([fasttree, '-nt', '-gtr', '-quiet', src], stdout=fo, stderr=sb.PIPE)
        if r.returncode != 0 or os.path.getsize(tree) == 0:
            error('No phylogeny could be built: RAxML and FastTree both failed. FastTree: '
                  + r.stderr.decode(errors='replace')[-500:], exit=True)
        self.tree_note = 'FastTree (-nt -gtr), because RAxML failed: {}'.format(reason) + (
            '; without {} samples that have no base at the kept positions ({})'.format(len(empty), ', '.join(empty))
            if empty else '')
        warning('Tree built with FastTree instead of RAxML: {}'.format(tree))
        return tree


    def __init__(self, args):
        self.samples = args.samples
        self.phylophlan_mode, self.phylophlan_params = self.check_phylophlan_params(args.phylophlan_mode,
                                                                                    args.phylophlan_params)
        self.phylophlan_configuration = args.phylophlan_configuration
        self.mutation_rates = args.mutation_rates
        self.output_dir = args.output_dir
        self.clade = args.clade
        self.treeshrink = args.treeshrink
        self.nprocs = args.nprocs
        self.tree_note = None
