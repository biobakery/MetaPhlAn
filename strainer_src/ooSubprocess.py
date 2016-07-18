#!/usr/bin/env python
# Author: Duy Tin Truong (duytin.truong@unitn.it)
#        at CIBIO, University of Trento, Italy


import subprocess
import os
import multiprocessing
from multiprocessing.pool import ThreadPool
import sys
import cStringIO
from tempfile import NamedTemporaryFile 
import which
import functools
import traceback
import numpy


class ooSubprocessException(Exception):
    pass

class ooSubprocess:

    def __init__(self, tmp_dir='tmp/'):
        self.chain_cmds = []
        self.tmp_dir = tmp_dir
        mkdir(tmp_dir)

    def ex(
            self,
            prog,
            args=[],
            get_output=False,
            get_out_pipe=False,
            out_fn=None,
            in_pipe=None,
            verbose=True,
            **kwargs):

        if not which.is_exe(prog):
            raise ooSubprocessException('Error: cannot find the program %s '\
                        'in the executable path!'%prog)

        if isinstance(args, str):
            args = args.split()

        if not isinstance(args, list):
            args = [args]

        cmd = [prog] + args
        print_cmd = 'ooSubprocess: ' + ' '.join(cmd)
        if verbose and out_fn and (not get_output):
            print_stderr(print_cmd + ' > ' + out_fn)
        elif verbose:
            print_stderr(print_cmd)

        if get_output:
            result = subprocess.check_output(
                                                cmd,
                                                stdin=in_pipe,
                                                **kwargs)
        elif get_out_pipe:
            tmp_file = NamedTemporaryFile(dir=self.tmp_dir)
            p = subprocess.Popen(
                                        cmd, 
                                        stdin=in_pipe, 
                                        stdout=tmp_file,
                                        **kwargs)
            p.wait()
            if in_pipe != None:
                in_pipe.close()
            tmp_file.seek(0)
            result = tmp_file
        elif out_fn:
            ofile = open(out_fn, 'w') if out_fn else None
            result = subprocess.check_call(
                                            cmd, 
                                            stdin=in_pipe, 
                                            stdout=ofile, 
                                            **kwargs)
            ofile.close()
        else:
            result = subprocess.check_call(
                                            cmd, 
                                            stdin=in_pipe, 
                                            **kwargs)
        return result

    def chain(
            self,
            prog,
            args=[],
            stop=False,
            in_pipe=None,
            get_output=False,
            get_out_pipe=False,
            out_fn=None,
            verbose=True,
            **kwargs):

        if not which.is_exe(prog):
            raise ooSubprocessException('Error: cannot find the program %s '\
                        'in the executable path!'%prog)

        if in_pipe is None and self.chain_cmds != []:
            raise ooSubprocessException(
                'The pipeline was not stopped before creating a new one!'\
                'In cache: %s' % (' | '.join(self.chain_cmds)))
        if out_fn and stop == False:
            raise ooSubprocessException(
                'out_fn (output_file_name) is only specified when stop = True!')

        if isinstance(args, str):
            args = args.split()

        if not isinstance(args, list):
            args = [args]
        cmd = [prog] + args

        print_cmd = ' '.join(cmd)
        if out_fn and (not get_output):
            print_cmd += ' > ' + out_fn
        self.chain_cmds.append(print_cmd)

        if stop:
            if in_pipe is None:
                raise ooSubprocessException('No input process to create a pipeline!')

            if verbose:
                print_stderr('ooSubprocess: ' + ' | '.join(self.chain_cmds))

            self.chain_cmds = []
            if get_output:
                result = subprocess.check_output(
                                                    cmd,
                                                    stdin=in_pipe,
                                                    **kwargs)
            elif get_out_pipe:
                tmp_file = NamedTemporaryFile(dir=self.tmp_dir)
                p = subprocess.Popen(
                                    cmd,
                                    stdin=in_pipe,
                                    stdout=tmp_file,
                                    **kwargs)
                return_code = p.wait()
                if return_code != 0:
                    raise ooSubprocessException(
                                    'Failed when executing the command: %s\n'\
                                    'return code: %s'\
                                    %(' | '.join(self.chain_cmds), return_code))
                tmp_file.seek(0)
                if in_pipe != None:
                    in_pipe.close()

                result = tmp_file
            elif out_fn:
                ofile = open(out_fn, 'w')
                result = subprocess.check_call(
                                                cmd,
                                                stdin=in_pipe,
                                                stdout=ofile,
                                                **kwargs)
                ofile.close()
            else:
                result = subprocess.check_call(
                                                cmd,
                                                stdin=in_pipe,
                                                **kwargs)
        else:
            tmp_file = NamedTemporaryFile(dir=self.tmp_dir)
            p = subprocess.Popen(
                                cmd,
                                stdin=in_pipe,
                                stdout=tmp_file,
                                **kwargs)
            return_code = p.wait()
            if return_code != 0:
                    raise ooSubprocessException(
                                    'Failed when executing the command: %s\n'\
                                    'return code: %s'\
                                    %(' | '.join(self.chain_cmds), return_code))
            tmp_file.seek(0)
            if in_pipe != None:
                in_pipe.close()

            result = tmp_file
        return result

    def ftmp(self, ifn):
        return os.path.join(self.tmp_dir, os.path.basename(ifn))


def fdir(dir, ifn):
    return os.path.join(dir, os.path.basename(ifn))


def mkdir(dir):
    if not os.path.exists(dir):
        try:
            os.makedirs(dir)
        except OSError as e:
            if e.errno != 17:
                raise
            pass
    elif not os.path.isdir(dir):
        raise ooSubprocessException('Error: %s is not a directory!' % dir)


def replace_ext(ifn, old_ext, new_ext):
    # if not os.path.isfile(ifn):
    #    print_stderr('Error: file %s does not exist!'%(ifn))
    #    exit(1)
    if ifn[len(ifn) - len(old_ext):] != old_ext:
        # print_stderr('Error: the old file extension %s does not match!'%old_ext)
        # exit(1)
        new_ifn = ifn + new_ext
    else:
        new_ifn = ifn[:len(ifn) - len(old_ext)] + new_ext
    return new_ifn


def splitext(ifn):
    basename = os.path.basename(ifn)
    if ifn.endswith('.tar.bz2'):
        ext = '.tar.bz2'
    elif ifn.endswith('.tar.gz'):
        ext = '.tar.gz'
    else:
        ext = basename.split('.')[-1]
        if ext != basename:
            ext = '.' + ext
    base = basename[:-len(ext)]
    for t in ['.sam', '.fastq', '.fasta', '.fna']:
        if base.endswith(t):
            ext = t + ext
            base = base[:-len(t)]
    return base, ext


def trace_unhandled_exceptions(f):
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except:
            #traceback.print_exc()
            #raise Exception(''.join(traceback.format_exception(*sys.exc_info())))
            raise Exception(traceback.format_exc())
    return wrapper


def parallelize(func, args, nprocs=1, use_threads=False):
    if nprocs > 1:
        if use_threads:
            pool = ThreadPool(nprocs)
        else:
            pool = multiprocessing.Pool(nprocs)
        results = pool.map(func, args)
        pool.close()
        pool.join()
    else:
        results = serialize(func, args)
    return results


def parallelize_async(func, args, nprocs=1, use_threads=False):
    if nprocs > 1:
        if use_threads:
            pool = ThreadPool(nprocs)
        else:
            pool = multiprocessing.Pool(nprocs)
        app_results = []
        for a in args:
            app_results.append(pool.apply_async(func, [a]))
        pool.close()
        pool.join()
        results = [r.get() for r in app_results]
    else:
        results = serialize(func, args)
    return results


def serialize(func, args):
    results = []
    for arg in args:
        results.append(func(arg))
    return results


def print_stderr(*args):
    sys.stderr.write(' '.join(map(str, args)) + '\n')
    sys.stderr.flush()


def print_stdout(*args):
    sys.stdout.write(' '.join(map(str, args)) + '\n')
    sys.stdout.flush()



