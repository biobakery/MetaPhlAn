#!/usr/bin/env python

#############################################################################
##  This file is part of PASTA and is forked from SATe2.
##  See "LICENSE.txt" for terms and conditions of usage.
#############################################################################

"""
File and path management.
"""

import os
import sys
import re
import tempfile
import shutil
from threading import Lock
#from pasta import get_logger
#_LOG = get_logger(__name__)

_ILLEGAL_FILENAME_PATTERN = re.compile(r'[^-_a-zA-Z0-9.]')
def get_safe_filename(filename):
    return "".join(_ILLEGAL_FILENAME_PATTERN.split(filename))

def quoted_file_path(path):
    if '"' not in path:
        return '"'+ path + '"'
    elif "'" not in path:
        return "'" + path + "'"
    else:
        path = path.replace('"', r'\"')
        return '"' + path + '"'

def open_with_intermediates(filepath, mode):
    """Opens a `filepath` in the `mode`, but will create intermediate directories if they are absent."""
    d = os.path.dirname(filepath)
    if d:
        if not os.path.exists(d):
            os.makedirs(d)
        elif not os.path.isdir(d):
            raise IOError('The file "%s" cannot be created because "%s" exists but is not a directory' % (filepath, d))
    return open(filepath, mode)

class TempFS(object):
    '''A wrapper for creating temporary directories that will safeguard against
    removing directories that were not created by PASTA (there is no evidence
    that this has happened, but rmtree is a dangerous call and it would be a
    horrible bug to have).

    Note that this class does not protect against incorrect order of deletion.
        If you delete the top level tempdirectory then all of the subdirectories
        will be deleted.
    '''

    run_generated_filenames = [
        ".Job.stderr.txt",
        ".Job.stdout.txt",
        "1.fasta",
        "2.fasta",
        "RAxML_bestTree.default",
        "RAxML_info.default",
        "RAxML_log.default",
        "RAxML_parsimonyTree.default",
        "RAxML_result.default",
        "disttbfast",
        "dvtditr",
        "end_pastaiter_timestamp.txt",
        "end_treeinference_timestamp.txt",
        "hat2",
        "hat3",
        "infile.tree",
        "input.aligned",
        "input.fasta",
        "input.dnd",
        "input.phy",
        "input.phy.reduced",
        "last_used.cfg",
        "order",
        "out.fasta",
        "pairlocalalign",
        "partition.txt",
        "partition.txt.reduced",
        "pre",
        "start.tre",
        "start_align_timestamp.txt",
        "start_pastaiter_timestamp.txt",
        "start_treeinference_timestamp.txt",
        "tbfast",
        "trace",
        "log",
        "results",
        "query.fasta"
        ]

    def __init__(self):
        self._directories_created = set()
        self._top_level_temp = None
        self._top_level_temp_real = None
        self._directories_created_lock = Lock()

    def _is_already_created(self, real_path):
        self._directories_created_lock.acquire()
        b = real_path in self._directories_created
        self._directories_created_lock.release()
        return b

    def _register_created_dir(self, path):
        self._directories_created.add(os.path.abspath(path))

    def create_subdir(self, dir):
        '''Creates a directory `dir`

        `dir` must be a subdirectory of self.top_level_temp and `dir`
        cannot exist.

        The canonical file path is returned.

        Client code is responsible for deleting directory by a call to
        `remove_dir` using the same `TempFS` instance.

        '''
        rp = os.path.realpath(dir)
        if os.path.exists(rp):
            raise OSError("Path exists: '%s'" % rp)
        if not self._is_in_top_level_temp(rp):
            raise OSError("Subdirectory is not under the top level temp dir: '%s'" % rp)
        self._directories_created_lock.acquire()
        try:
            already_in = rp in self._directories_created
            if already_in:
                raise OSError("Subdirectory is flagged as having already been created: '%s'" % rp)
            else:
                os.makedirs(rp)
                self._register_created_dir(rp)
                return rp
        finally:
            self._directories_created_lock.release()


    def create_top_level_temp(self, parent, prefix='temp'):
        '''Creates (and stores the path to) the top-level temporary
        directory under `parent`

        `parent` must already exist
        `prefix` is passed to tempfile.mkdtemp

        Client code is responsible for deleting directory by a call to
        `remove_dir` using the same `TempFS` instance.

        The canonical file path is returned.

        '''
        assert(self._top_level_temp is None)
        assert(self._top_level_temp_real is None)
        r_parent = os.path.realpath(parent)
        if not os.path.exists(r_parent):
            raise OSError("Path does not exist: '%s'" % r_parent)
        if not os.path.isdir(r_parent):
            raise OSError("Path is not a directory: '%s'" % r_parent)
        self._directories_created_lock.acquire()
        try:
            self._top_level_temp = tempfile.mkdtemp(prefix=prefix, dir=r_parent)
            self._top_level_temp_real = os.path.realpath(self._top_level_temp)
            self._register_created_dir(self._top_level_temp_real)
            return self._top_level_temp_real
        finally:
             self._directories_created_lock.release()


    def create_temp_subdir(self, parent, prefix='temp'):
        '''Creates (and stores the path to) a temporary directory under `parent`

        `parent` should be "below" the top_level_temp.
        `prefix` is passed to tempfile.mkdtemp

        Client code is responsible for deleting directory by a call to
        `remove_dir` using the same `TempFS` instance.

        The canonical file path is returned.

        '''
        r_parent = os.path.realpath(parent)
        if not os.path.exists(r_parent):
            raise OSError("Path does not exist: '%s'" % r_parent)
        if not os.path.isdir(r_parent):
            raise OSError("Path is not a directory: '%s'" % r_parent)
        if not self._is_in_top_level_temp(r_parent):
            raise OSError("Subdirectory is not under the top level temp dir: '%s'" % r_parent)
        self._directories_created_lock.acquire()
        try:
            d = tempfile.mkdtemp(prefix=prefix, dir=r_parent)
            rp = os.path.realpath(d)
            self._register_created_dir(rp)
            return rp
        finally:
            self._directories_created_lock.release()

    def remove_dir(self, real_path):
        '''Recursively removes `real_path` from the filesystem if it is
        listed as one of the directories created by this TempFS object (or raises
        a ValueError if it is not listed).
        '''
        self._directories_created_lock.acquire()
        real_path = os.path.abspath(real_path)
        try:
            if real_path in self._directories_created:
                self._directories_created.remove(real_path)
                if (real_path == self._top_level_temp) or (real_path == self._top_level_temp_real):
                    self._top_level_temp = None
                    self._top_level_temp_real = None
            else:
                raise ValueError("'%s' is not registered as a temporary directory that was created by this process!" % real_path)
        finally:
            self._directories_created_lock.release()
        #_LOG.debug("Cleaning temp dir: '%s'" % real_path)
        if os.path.exists(real_path):
            for path in os.listdir(real_path):
                fpath = os.path.join(real_path, path)
                if os.path.isdir(fpath):
                    try:
                        self.remove_dir(fpath)
                    except ValueError as e:
                        sys.stderr.write("Refused to clean '%s': not created by PASTA\n%s\n" % (fpath,str(e)))
                    except OSError:
                        sys.stderr.write("Error trying to remove '%s'\n" % fpath)
            for fname in self.run_generated_filenames:
                try:
                    os.remove(os.path.join(real_path, fname))
                except OSError:
                    pass
            try:
                os.rmdir(real_path)
            except OSError:
                pass
            return True
        return False

    def get_remaining_directories(self):
        '''Returns a copy of the set of directories that have been created but
        not deleted.
        '''
        self._directories_created_lock.acquire()
        c = set(self._directories_created)
        self._directories_created_lock.release()
        return c

    def get_top_level_temp(self):
        return self._top_level_temp_real
    top_level_temp = property(get_top_level_temp)

    def _is_in_top_level_temp(self, real_path):
        if self._top_level_temp_real is None:
            raise ValueError("_top_level_temp has not been set, yet!")
        in_common = os.path.commonprefix([self._top_level_temp_real, real_path])
        return in_common == self._top_level_temp_real

class PastaProducts(object):
    """
    Handles paths to all (final) output produced by PASTA.
    """

    meta_product_types = {
            "score": ".score.txt",
            "tree": ".tre",
            "run_log" : ".out.txt",
            "err_log" : ".err.txt"
            }

    def __init__(self, pasta_user_settings):
        """
        Configures self based on a fully-populated `PastaUserSettings` config
        object.
        """
        self.pasta_user_settings = pasta_user_settings
        self.job_name = self.pasta_user_settings.commandline.job
        if self.job_name is None:
            self.job_name = "pastajob"
        self._job_file_name = get_safe_filename(self.job_name)

        if not self.pasta_user_settings.sate.output_directory:
            self._output_directory = self.get_input_source_directory()
        else:
            self._output_directory = self.pasta_user_settings.sate.output_directory

        if os.path.exists(self._output_directory) and not os.path.isdir(self._output_directory):
            raise Exception("Requested output directory is not valid because a file of the same name already exists: %s" % self._output_directory)

        # ensure input sources populated
        assert self.pasta_user_settings.input_seq_filepaths
        # alignment input/output names
        self._original_input_files = [os.path.abspath(f) for f in self.pasta_user_settings.input_seq_filepaths]
        #self.output_alignment_suffixes = ["." + get_safe_filename(f) + ".aln" \
        #        for f in self.original_input_files]
        self._output_alignment_suffixes = []
        self._input_fpath_alignment_suffix_map = {}
        self._alignment_suffix_input_fpath_map = {}
        fasta_extension_pattern = re.compile(r"\.fast?a?$")
        for fidx, f in enumerate(self._original_input_files):
            safe_fn = get_safe_filename(os.path.basename(f))
            trunc_fn = fasta_extension_pattern.sub("", safe_fn)
            fn_stem = "marker%03d.%s" % (fidx+1, trunc_fn)
            suffix = "." + fn_stem + ".aln"
            sidx = 1
            while suffix in self._output_alignment_suffixes:
                suffix = "." + fn_stem + "-" + str(sidx) + ".aln"
                sidx += 1
            self._output_alignment_suffixes.append(suffix)
            self._input_fpath_alignment_suffix_map[f] = suffix
            self._alignment_suffix_input_fpath_map[suffix] = f

        # initialize/create attributes, setting to dummy values
        for stream_name in self.meta_product_types:
            self._set_stream(stream_name, None)
        self.alignment_streams = []
        self.other_streams = []
        self.input_fpath_alignment_stream_map = {}

        # dummy output prefix
        self.output_prefix = None

        # create working output streams
        self.setup()

    def _compose_stream_attr(self, stream_name):
        return stream_name + "_stream"

    def _get_stream(self, stream_name):
        return getattr(self, self._compose_stream_attr(stream_name), None)

    def _set_stream(self, stream_name, value):
        setattr(self, self._compose_stream_attr(stream_name), value)

    def setup(self):
        """
        Checks for file name clashes, disambiguates if neccessary, and creates
        output files.
        Called by PastaProducts.__init__ (and nowhere else as of Apr 2012)
        """
        self.create_output_prefix()
        self.create_product_paths()

    def create_product_paths(self):
        """
        Opens file streams for writing so that the files will exist during
            the run (and we do not have to keep checking for an output prefix
            that does not overwrite files).
        Called after create_output_prefix by PastaProducts.__init__ 
            via PastaProducts.setup() (and nowhere else as of Apr 2012)
        """
        assert self.output_prefix
        for stream_name, product_extension in list(self.meta_product_types.items()):
            output_path = self.output_prefix + product_extension
            stream = open_with_intermediates(output_path, "w")
            self._set_stream(stream_name, stream)
            self.other_streams.append(stream)
        for asi, sf in enumerate(self._output_alignment_suffixes):
            output_path = self.output_prefix + sf
            stream = open_with_intermediates(output_path, "w")
            self.alignment_streams.append(stream)
            self.input_fpath_alignment_stream_map[self._alignment_suffix_input_fpath_map[sf]] = stream

    def create_output_prefix(self):
        output_prefix_stem = os.path.join(self._output_directory, self._job_file_name)
        idx = 0
        disambiguator = ""
        while True:
            output_prefix = output_prefix_stem + disambiguator
            if not self.check_for_existing_files(output_prefix):
                break
            idx += 1
            disambiguator = "%d" % idx
        self.output_prefix = output_prefix

    def check_for_existing_files(self, output_prefix):
        for ext in list(self.meta_product_types.values()):
            if os.path.exists(output_prefix + ext):
                return True
        for fn in self._output_alignment_suffixes:
            if os.path.exists(output_prefix + fn):
                return True
        return False

    def get_input_source_directory(self):
        """
        Given a configuration object, returns the directory of the input file(s).
        """
        options = self.pasta_user_settings.commandline
        if options.multilocus:
            # multilocus dataset: assume directory is given as input source
            return os.path.abspath(options.input)
        else:
            # single locus dataset: return directory nanme
            return os.path.dirname(os.path.abspath(options.input))
    
    def get_abs_path_for_iter_output(self, iter_num, out_tag, allow_existing=False):
        """
        Returns an absolute path or None for the file for an iteration temporary
            file for iteration `iter_num` with the specificed `out_tag`
        
        The function will try to create numbered versions of the files to avoid
            overwriting an existing file. But it can return None (if numbering
            exceeds a large # of files).
            
        It is not thread-safe.
        """
        p = "iteration_" + str(iter_num) + '_' + out_tag
        o_path = self.output_prefix  + "_temp_" + p
        if not allow_existing:
            if os.path.exists(o_path):
                n = 1
                while os.path.exists(o_path):
                    t_tag = "_temp%d_" % n
                    if n > 100:
                        #_LOG.warn('File %s exists iteration-specific output skipped!' % o_path)
                        return None # don't create a huge # of numbered files
                    o_path = self.output_prefix + t_tag + p
                    n += 1
        return os.path.abspath(o_path)

    def get_abs_path_for_tag(self, out_tag):
        """
        Returns an absolute path or None for the file with the specificed `out_tag`
        
        The function will try to create numbered versions of the files to avoid
            overwriting an existing file. But it can return None (if numbering
            exceeds a large # of files).
            
        It is not thread-safe.
        """
        o_path = self.output_prefix  + "_temp_" + out_tag
        if os.path.exists(o_path):
            n = 1
            while os.path.exists(o_path):
                t_tag = "_temp%d_" % n
                if n > 100:
                    #_LOG.warn('File %s exists iteration-specific output skipped!' % o_path)
                    return None # don't create a huge # of numbered files
                o_path = self.output_prefix + t_tag + out_tag
                n += 1
        return os.path.abspath(o_path)
