from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from numpy import concatenate as concat
from numpy import diff, sqrt
import xarray as xr


k = 15
ym = 1
xm = 15

y = np.linspace(-ym, ym, 200)

r2 = 1 - y**2
dr2 = y**2 / (1 - y**2)

plt.figure(3)
plt.subplot(1, 3, 1)
plt.plot(r2, y)
plt.subplot(1, 3, 2)
plt.plot(dr2, y)
plt.xlim(-0.5, 15)
plt.subplot(1, 3, 3)
plt.plot(dr2 + k * r2, y)
plt.xlim(-0.5, 25)
