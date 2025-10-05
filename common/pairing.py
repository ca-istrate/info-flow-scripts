import numpy as np
import logging


logger = logging.getLogger("script")

class Pairing:
    x: np.ndarray
    y: np.ndarray

    x_time: np.ndarray
    y_time: np.ndarray

    start: int
    stop: int


    def __init__(self, x, y, *, start=0, end=None, shift=0):
        if len(x[:,1]) != len(y[:,1]):
            logger.warning("Warning: The two inputs do not have the same number of rows.")

        if shift < 0:
            logger.debug(f"Shift {shift} is negative. Dropping {-shift} values from the beginning of y.")
            self.x = x[:,1]
            self.y = y[-shift:,1]
            self.x_time = x[:,0]
            self.y_time = y[-shift:,0]
        else:
            logger.debug(f"Shift {shift} is positive. Dropping {shift} values from the beginning of x.")
            self.x = x[shift:,1]
            self.y = y[:,1]
            self.x_time = x[shift:,0]
            self.y_time = y[:,0]

        self.start = start
        self.end = min(len(self.x), len(self.y))

        if end is not None:
            if self.end < end:
                logger.error(f"Max data index ({self.end}) is less than specified value {end}.")
                exit()

            self.end = end

        logger.debug(f"Truncating data to contain only values between indicies {start} and {end}.")
        self.x = self.x[self.start:self.end]
        self.x_time = self.x_time[self.start:self.end]
        self.y = self.y[self.start:self.end]
        self.y_time = self.y_time[self.start:self.end]

    def subset(self, *, start=0, end=None, shift=0):
        x = np.vstack([self.x_time, self.x]).T
        y = np.vstack([self.y_time, self.y]).T
        return Pairing(x, y, start=start, end=end, shift=shift)


    def reversed(self):
        x = np.vstack([self.x_time[::-1], self.x[::-1]]).T
        y = np.vstack([self.y_time[::-1], self.y[::-1]]).T
        return Pairing(x, y)


    @classmethod
    def from_files(cls, x_file, y_file, **kwargs):
        logger.debug(f"Reading data from files x = {x_file} and y = {y_file}.")
        return Pairing(np.genfromtxt(x_file), np.genfromtxt(y_file), **kwargs)

