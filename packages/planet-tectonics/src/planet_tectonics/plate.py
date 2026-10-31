"""One rigid plate: its rotation, Euler vector and the crust raster it keeps in its own frame."""
import numpy as np



class Plate:
    def __init__(self, pid, size, axis, speed):
        self.id = pid
        self.Q = np.eye(3)
        self.axis = np.asarray(axis, np.float64)/np.linalg.norm(axis)
        self.speed = float(speed)                 # cm/yr at 90 degrees from the pole
        self.mask = np.zeros(size, bool)
        self.cont = np.zeros(size, bool)
        self.thick = np.zeros(size, np.float32)
        self.born = np.zeros(size, np.float32)    # crust formation time, Myr
        self.orog = np.full(size, np.nan, np.float32)  # last orogeny (collision) time, Myr
        self.under = np.zeros(size, np.int16)           # consecutive steps hidden under another plate
        self.sed = np.zeros(size, np.float32)            # sediment, m (on continental crust too with continental_sediment)
        self.slab = np.zeros(size, np.float32)           # underthrust convergence beside this plate's overriding edge, km

    def omega(self, radius):
        """Angular velocity vector, rad/yr."""
        return self.axis*self.speed*1e-2/radius
