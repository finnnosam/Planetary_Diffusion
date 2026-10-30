"""numba's njit, or a pass-through when numba is missing (the kernels then run as plain Python, slowly)."""
try:
    from numba import njit
    HAVE_NUMBA = True
except ImportError:  # pragma: no cover - exercised only without numba
    HAVE_NUMBA = False

    def njit(*args, **kwargs):
        if args and callable(args[0]):
            return args[0]
        return lambda f: f
