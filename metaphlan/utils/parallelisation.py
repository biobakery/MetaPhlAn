__author__ = ('Aitor Blanco Miguez (aitor.blancomiguez@unitn.it), '
              'Duy Tin Truong (duytin.truong@unitn.it), '
              'Francesco Asnicar (f.asnicar@unitn.it), '
              'Moreno Zolfo (moreno.zolfo@unitn.it), '
              'Francesco Beghini (francesco.beghini@unitn.it)')
__version__ = '3.0'
__date__ = '21 Feb 2020'

try:
    from .util_fun import error
except ImportError:
    from util_fun import error
from multiprocessing import Event, Pool

CHUNKSIZE = 1

"""
Places terminating in the global namespace of the worker subprocesses.
This allows the worker function to access `terminating` even though it is
not passed as an argument to the function.

:param terminating_: the local terminating variable
"""
def init_terminating(terminating_):    
    global terminating
    terminating = terminating_


"""
Executes each parallelised call of a function

:param arguments: the name of the function and the arguments
:returns: the result of the parallelised execution
"""
def parallel_execution(arguments):
    function, *args = arguments
    if not terminating.is_set():
        try:
            return function(*args)
        except Exception as e:
            terminating.set()
            raise
    else:
        terminating.set()


"""
Creates a pool for a parallelised function and returns the results of each execution
as a list

:param args: the name of the function and the arguments
:param nprocs: the number of threads to use in the pool
:returns: the result of the parallelised funtion
"""
def execute_pool(args, nprocs):
    terminating = Event()
    with Pool(initializer=init_terminating, initargs=(terminating,), processes=nprocs) as pool:
        try:
            return [_ for _ in pool.imap_unordered(parallel_execution, args, chunksize=CHUNKSIZE)] 
        except Exception as e:
            error('Parallel execution fails: '+str(e), init_new_line=True, exit=True)
