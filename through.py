import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

table=ascii.read("/Users/dm/Documents/GitHub/QSOSIM10/BOSS_throughput.dat")#, readme="/Users/dm/Documents/Isotropy/estimation_A_gamma/J_A+A_373_757/ReadMe")

size=np.size(table)
wave=np.zeros(size)
thro=np.zeros(size)
for i in xrange(size):
    wave[i]=table[i][0]
    thro[i]=table[i][1]
    print wave[i],thro[i]

plt.plot(wave,thro)
plt.show()
