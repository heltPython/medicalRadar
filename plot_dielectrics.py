#!/usr/bin/env python
#encoding:latin-1

#Move these later:
import os
import numpy
import scipy
from matplotlib import rc
import pylab
from scipy import signal
import processing
from Measurement import Measurement

############################################
#Setting plot configurations:
#rc('text', usetex=True)
from matplotlib import rc
rc('font',size=16)
rc('lines',lw=2)
#################################################


def getdielectric(freqs=(5e8,2e10),tissue='skin'):
    """
    Returns the relatice permittivity and the conductivity of given tissue over the frequency span.
    @param freqs: tuple containing start and end frequencies.
    """
    mu = numpy.pi*4e-7
    epsilon0 = 8.854e-12
    
    file = open('../pnadata/'+tissue,'r')
    lines = file.readlines()
    file.close()
    freqvec = numpy.linspace(freqs[0],freqs[1],len(lines[3:]))
    permittivity = []
    conductivity = []
    for line in lines[3:]:
        vals = line.split()
        if float(vals[1]) >= freqs[0] and float(vals[1]) <= freqs[1]:
            conductivity.append(float(vals[2]))
            permittivity.append(float(vals[3]))

    return freqvec, permittivity, conductivity



#Reflection coefficient between air and skin for various frequencies:
blood = getdielectric(tissue='blood')
muscle = getdielectric(tissue='muscle')
heart = getdielectric(tissue='heart')
lung = getdielectric(tissue='lung')
skin = getdielectric(tissue='skin')
bone = getdielectric(tissue='bone')
fat = getdielectric(tissue='fat')

#Permittivity:
pylab.figure()
pylab.plot(blood[0]*1e-9,blood[1],label='blood')
pylab.plot(heart[0]*1e-9,heart[1],label='heart')
pylab.plot(muscle[0]*1e-9,muscle[1],label='muscle')
pylab.plot(skin[0]*1e-9,skin[1],label='skin')
pylab.plot(lung[0]*1e-9,lung[1],label='lung')
pylab.plot(bone[0]*1e-9,bone[1],label='bone')
pylab.plot(fat[0]*1e-9,fat[1],label='fat')
pylab.ylabel(r"$\varepsilon_r'$(-)")
pylab.xlabel('Frequency (GHz)')
pylab.legend()

#conductivity
pylab.figure()
pylab.plot(blood[0]*1e-9,blood[2],label='blood')
pylab.plot(heart[0]*1e-9,heart[2],label='heart')
pylab.plot(muscle[0]*1e-9,muscle[2],label='muscle')
pylab.plot(skin[0]*1e-9,skin[2],label='skin')
pylab.plot(lung[0]*1e-9,lung[2],label='lung')
pylab.plot(bone[0]*1e-9,bone[2],label='bone')
pylab.plot(fat[0]*1e-9,fat[2],label='fat')
pylab.ylabel(r'$\sigma\; (\Omega m)$')
pylab.xlabel('Frequency (GHz)')
pylab.legend(loc='best')







pylab.show()
