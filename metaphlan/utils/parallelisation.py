__author__ = ('Aitor Blanco Miguez (aitor.blancomiguez@unitn.it), '
              'Duy Tin Truong (duytin.truong@unitn.it), '
              'Francesco Asnicar (f.asnicar@unitn.it), '
              'Moreno Zolfo (moreno.zolfo@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it)')
__version__ = '4.2.0'
__date__ = '14 May 2025'

from typing import Iterable, Callable, Any

try:
    from .util_fun import error
except ImportError:
    from util_fun import error
from multiprocessing import Event, Pool

CHUNKSIZE = 1


def init_terminating(terminating_):
    """Places terminating in the global namespace of the worker subprocesses. 
    This allows the worker function to access `terminating` even though it is not passed as an argument to the function.


    Args:
        terminating_ (Event): the initialized terminating event
    """
    global terminating
    terminating = terminating_

def parallel_execution(arguments):
    """Executes each parallelised call of a function

    Args:
        arguments ((Callable, *Any)): tuple with the function and its arguments

    Returns:
        function: the call to the function
    """
    function, *args = arguments
    if not terminating.is_set():
        try:
            return function(*args)
        except Exception as e:
            terminating.set()
            raise e
    else:
        terminating.set()

def execute_pool_iter(args, nprocs, ordered):
    try:
        terminating = Event()
        with Pool(initializer=init_terminating, initargs=(terminating,), processes=nprocs) as pool:
            if ordered:
                f = pool.imap
            else:
                f = pool.imap_unordered

            for r in f(parallel_execution, args, chunksize=CHUNKSIZE):
                yield r
    except Exception as e:
        error('Parallel execution fails: {}'.format(e), exit=False)
        raise e

def execute_pool(args, nprocs, return_generator=False, ordered=False):
    """Creates a pool for a parallelised function and returns the results of each execution as a list

    Args:
        args (Iterable[tuple]): tuple with the function and its arguments
        nprocs (int): number of procs to use

    Returns:
        list: the list with the results of the parallel executions
    """
    try:
        gen = execute_pool_iter(args, nprocs, ordered)
                
        if return_generator:
            return gen
        else:
            return [_ for _ in gen]
    
    except Exception as e:
            error('Parallel execution fails: {}'.format(e), exit=True)
