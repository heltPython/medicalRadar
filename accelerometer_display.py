#!/usr/bin/env python

import sys
import os
import numpy
import scipy
from scipy import signal
import pylab

def simple_filt(input_signal,passband,prf):
    """
    A simple filter creator that creates an iirfilter based on the passband given by passband and the prf, and filters the input signal.
    @param input_signal: The signal to be filtered.
    @param passband: The pass band list [low, high]. If low=0, the filter will be a lowpass filter, if high=0, the filter is a highpass filter.
    @param prf: The prf.
    """
    gpass=4
    gstop=40
    if passband[0]==0:
        #Lowpassfilter
        wp=passband[1]/(prf/2.0)
        ws=min(wp*1.5,1)
    elif passband[1]==0:
        #Highpassfilter
        wp=passband[0]/(prf/2.0)
        ws=wp*0.5
        input_signal = signal.detrend(input_signal)
    else:
        #Bandpassfilter
        wp=[passband[0]/(prf/2.0) , passband[1]/(prf/2.0)]
        #ws=[wp[0]-min(wp[0],.25),wp[1]+.25]
        ws=[wp[0]*0.5 , min(wp[1]*1.5,1)]
        input_signal = signal.detrend(input_signal)
        
    b,a=signal.iirdesign(wp,ws,gpass,gstop,ftype='cheby2')#Design of the filter
    heartbeats=signal.lfilter(b,a,input_signal.conj().T)
    #The time reversed signal is run through the same filter to correct phase distortion:
    if len(heartbeats.shape)>1:
        heartbeats=signal.lfilter(b,a,heartbeats[:,::-1])
        heartbeats=heartbeats[:,::-1]
    else:
        heartbeats=signal.lfilter(b,a,heartbeats[::-1])
        heartbeats=heartbeats[::-1]

    print 'Filter order: ',len(b)

    return heartbeats.conj().T


def integrate(dt,data):
    intgtd = numpy.cumsum(data)*dt
    intgtd = signal.detrend(intgtd,type='constant')

    return intgtd


def plot_movement(ddx,ddy,ddz,fs,t,title=''):
    """Plots accelerometer movement and returns mean deviation from mean"""
    ddx = (ddx-numpy.mean(ddx))*g
    ddx = simple_filt(ddx,[0.5,0],fs)
    dx = integrate(dt,ddx)
    #dx = simple_filt(dx,[0.5,0],fs)
    x = integrate(dt,dx)
    #x = simple_filt(x,[0.5,0],fs)
    ddy = (ddy-numpy.mean(ddy))*g
    ddy = simple_filt(ddy,[0.5,0],fs)
    dy = integrate(dt,ddy)
    #dy = simple_filt(dy,[0.5,0],fs)
    y = integrate(dt,dy)
    #y = simple_filt(y,[0.5,0],fs)
    ddz = (ddz-numpy.mean(ddz))*g
    ddz = simple_filt(ddz,[0.5,0],fs)
    dz = integrate(dt,ddz)
    #dz = simple_filt(dz,[0.5,0],fs)
    z = integrate(dt,dz)
    #z = simple_filt(z,[0.5,0],fs)
    
    pylab.figure()
    pylab.subplot(311)
    pylab.plot(t[0:len(ddx)],ddx)
    pylab.xlabel('time(s)')
    pylab.ylabel(r'$\frac{m}{s^2}$')
    pylab.title(title+' Acceleration along x, y and z axes')
    pylab.subplot(312)
    pylab.plot(t[0:len(ddy)],ddy)
    pylab.xlabel('time(s)')
    pylab.ylabel(r'$\frac{m}{s^2}$')
    pylab.subplot(313)
    pylab.plot(t[0:len(ddz)],ddz)
    pylab.xlabel('time(s)')
    pylab.ylabel(r'$\frac{m}{s^2}$')
    
    pylab.figure()
    pylab.subplot(311)
    pylab.plot(t[0:len(dx)],dx)
    pylab.xlabel('time(s)')
    pylab.ylabel(r'$\frac{m}{s}$')
    pylab.title(title+' Velocity along x, y and z axes')
    pylab.subplot(312)
    pylab.plot(t[0:len(dy)],dy)
    pylab.xlabel('time(s)')
    pylab.ylabel(r'$\frac{m}{s}$')
    pylab.subplot(313)
    pylab.plot(t[0:len(dz)],dz)
    pylab.xlabel('time(s)')
    pylab.ylabel(r'$\frac{m}{s}$')
    
    pylab.figure()
    pylab.subplot(311)
    pylab.plot(t[0:len(x)],x*1e3)
    pylab.xlabel('time(s)')
    pylab.ylabel('mm')
    pylab.title(title+' Movement along x, y and z axes')
    pylab.subplot(312)
    pylab.plot(t[0:len(y)],y*1e3)
    pylab.xlabel('time(s)')
    pylab.ylabel('mm')
    pylab.subplot(313)
    pylab.plot(t[0:len(z)],z*1e3)
    pylab.xlabel('time(s)')
    pylab.ylabel('mm')

    return numpy.mean(numpy.mean(abs(x)),numpy.mean(abs(y)),numpy.mean(abs(z)))


    
try:
    Filename = sys.argv[1]
except:
    print 'Provide filename to be read'
    exit(1)

g=9.81 #m/s^2
if os.path.isdir(Filename):
    files = os.listdir(Filename)
    pleximovements = []
    baremovements = []
    for file in files:
        ifile = open(os.path.join(Filename,file),'r')
        lines = ifile.readlines()
        ifile.close()
        t=[]
        ddx=[]
        ddy=[]
        ddz=[]
        for line in lines[1:]:
            values = line.split()
            t.append(float(values[0]))
            ddx.append(float(values[1]))
            ddy.append(float(values[2]))
            ddz.append(float(values[3]))
        t=numpy.array([ms/1000 for ms in t])
        ddx=numpy.array(ddx)
        ddy=numpy.array(ddy)
        ddz=numpy.array(ddz)
                
        dt=t[1]-t[0]
        fs = 1/dt
        
        moves = plot_movement(ddx,ddy,ddz,fs,t,title=file)
        try:
            freq = int(file[-9:-5])
        except:
            freq =int(file[-8:-5])
        if file[7:11]=='bare':
            baremovements.append( (freq,moves) )
        elif file[7:12]=='plexi':
            pleximovements.append( (freq,moves) )

    baremovements.sort(key=lambda tup: tup[0])
    pleximovements.sort(key=lambda tup: tup[0])
    modulationratios = []
    for n in range(len(baremovements)):
        modulationratios.append(pleximovements[n][1]/baremovements[n][1])
        print 'Modulationratio for '+str(baremovements[n][0])+': '+ str(modulationratios[-1])
    print '\nModulationratio for 600Mhz: '+str(numpy.mean(modulationratios[0:4]))
    print 'Modulationratio for 1Ghz: '+str(numpy.mean(modulationratios[4:8]))
    print 'Modulationratio for 2Ghz: '+str(numpy.mean(modulationratios[8:12]))
    print 'Modulationratio for 3Ghz: '+str(numpy.mean(modulationratios[12:16]))
    print 'Modulationratio for 4Ghz: '+str(numpy.mean(modulationratios[16:20]))

else:
    ifile = open(Filename,'r')
    lines = ifile.readlines()
    ifile.close()
    t=[]
    ddx=[]
    ddy=[]
    ddz=[]
    for line in lines[1:]:
        values = line.split()
        t.append(float(values[0]))
        ddx.append(float(values[1]))
        ddy.append(float(values[2]))
        ddz.append(float(values[3]))
    t=numpy.array([ms/1000 for ms in t])
    ddx=numpy.array(ddx)
    ddy=numpy.array(ddy)
    ddz=numpy.array(ddz)

    dt=t[1]-t[0]
    fs = 1/dt
    
    plot_movement(ddx,ddy,ddz,fs,t)





pylab.show()
    
