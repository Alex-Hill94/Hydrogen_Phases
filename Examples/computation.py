import numpy as np
import h5py as h5
from decimal import *
from numpy import matrix
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import *
from Hydrogen_Phases import Compute_H_Fractions

path = 'test_data.hdf5' # Download at https://drive.google.com/file/d/1FqKmgt7HcPJQ-xrThG6wezKak97ftu8_/view?usp=sharing
u = h5.File(path, 'r')

PIDs = u['PIDs'][()]
Density = u['Density'][()]
Mass = u['Mass'][()]
Temperature = u['Temperature'][()]
H_abundance = u['H_abundance'][()]
Redshift = u['Redshift'][()]
SFR = u['SFR'][()]
u.close()


C = Compute_H_Fractions(Redshift = Redshift,
            H_abun = H_abundance,
            Temp = Temperature,
            Density = Density,
            SFR = SFR,
            Mass = Mass,
            PIDs = PIDs)

C.compute()
C.plot_example()
