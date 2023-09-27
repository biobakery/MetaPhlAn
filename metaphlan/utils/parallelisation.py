__author__ = ('Aitor Blanco Miguez (aitor.blancomiguez@unitn.it), '
              'Duy Tin Truong (duytin.truong@unitn.it), '
              'Francesco Asnicar (f.asnicar@unitn.it), '
              'Moreno Zolfo (moreno.zolfo@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it)')
__version__ = '4.1.0'
__date__ = '23 Aug 2023'


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


def execute_pool(args, nprocs):
    """Creates a pool for a parallelised function and returns the results of each execution as a list

    Args:
        args (Iterable[tuple]): tuple with the function and its arguments
        nprocs (int): number of procs to use

    Returns:
        list: the list with the results of the parallel executions
    """
    args = list(args)
    if nprocs == 1 or len(args) == 1:  # no need to initialize pool
        return [function(*a) for function, *a in args]
    else:
        terminating = Event()
        with Pool(initializer=init_terminating, initargs=(terminating,), processes=nprocs) as pool:
            try:
                return [_ for _ in pool.imap_unordered(parallel_execution, args, chunksize=CHUNKSIZE)]
            except Exception as e:
                error('Parallel execution fails: {}'.format(e), exit=False)
                raise e
