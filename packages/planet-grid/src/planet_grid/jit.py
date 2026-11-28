"""numba's njit and prange, or pass-throughs when numba is missing (the kernels then run as plain Python, slowly).

Kernels compiled with parallel=True use prange only over iterations that are independent (each writes its own
outputs; reductions are max or integer counts), so their results are bit-identical at any thread count (journal 154).
set_threads sets how many threads they use in the calling thread (the run's 'threads' parameter)."""
try:
    from numba import njit, prange
    import numba as _numba
    HAVE_NUMBA = True
except ImportError:  # pragma: no cover - exercised only without numba
    HAVE_NUMBA = False
    prange = range

    def njit(*args, **kwargs):
        if args and callable(args[0]):
            return args[0]
        return lambda f: f


def set_threads(n):
    """Threads for the parallel kernels (None: numba's default, all cores; at most numba's configured count)."""
    if not HAVE_NUMBA:
        return
    top = _numba.config.NUMBA_NUM_THREADS
    _numba.set_num_threads(top if n is None else max(1, min(int(n), top)))
