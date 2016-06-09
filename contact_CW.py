#!/usr/bin/env python
import numpy
import pylab
import os
import scipy
from Measurement import Measurement,UWBMeasurement
import processing
import readData


############################################
#Setting plot configurations:
#rc('text', usetex=True)
from matplotlib import rc
rc('font',size=20)
rc('lines',lw=4)
#################################################


def filter_ECG(ecg, fs):
    """Filters out the 50 Hz component from the ECG"""
    ecgfiltered = processing.simple_filt(ecg,[45.0,55.0],fs,ps='stop')
    return ecgfiltered


def initialize_ecg(filename,fs):
    """Reads an ECG measurement and makes it synchronous with the radar recording.
    @param filename: The full path of the ECG file.
    """
    try:
        ecg = numpy.array(readData.read_ECG(filename))
        p,ecgn=os.path.split(ecgfile)
        ecgn = os.path.splitext(ecgn)[0]

        indexfile = open(os.path.join(p,'startindexes'),'r')
        lines = indexfile.readlines()
        indexfile.close()
        found = False
        for line in lines:
            info = line.split(';')
            if ecgn == info[0]:
                ecg = numpy.array(ecg[int(info[1]):])
                found = True
        if not found:
            print 'No startindex found for '+ecgn
        ecg = processing.simple_filt(ecg,[0.5,0],fs)
    except:
        ecg=numpy.array([])
        ecg=numpy.zeros(1000)
        print 'No ECG file on path \n'+ecgfile

##     if os.path.splitext(os.path.split(filename)[-1])[0] == 'a1_ecg_out':
##         ecg=ecg[:len(ecg)/2]

    ecg = filter_ECG(ecg,fs)
    return ecg



def initialize_auscultation(filename,fs):
    """Comment yo!"""
    from scipy.io import wavfile
    sound_fs,auscultation = wavfile.read(filename)

    return auscultation[:,0]


def ecgpeaks(ecg):
    """
    Returns the indexes of the peaks in an ECG using a simple threshold detector
    """
    ecgmax = max(ecg)
    treshold = 0.4*ecgmax
    diff = numpy.diff(ecg)
    indexes = []
    for n in range(len(ecg)-2):
        if diff[n]>0 and diff[n+1]<0:
            if ecg[n+1]>treshold:
                if n>0:
                    indexes.append(n+1)
    return indexes







def waterfalls(measurements,time_index_window,type='instfreq',normalize=True,flipper=False):
    """
    Plots a waterfall of a slow time segment from each time instance
    """
    plotvectors=[]
    for n in range(len(measurements)):
        if type == 'instfreq':
            radarvec = processing.phase_estimate(measurements[n][0].data[time_index_window[n,0]:time_index_window[n,1]])
            radarvec = radarvec-numpy.mean(radarvec)
            if normalize:
                radarvec=radarvec/(numpy.max(radarvec)-numpy.min(radarvec))
            radarvec = numpy.diff(radarvec,n=1)*(measurements[n][0].timevector[3]-measurements[n][0].timevector[2])
            radarvec = numpy.append(radarvec,radarvec[-1])

            if any([measurements[n][1] == m for m in ['b2','a2','a32']]):
            #if any([measurements[n][1] == m for m in ['a1']]):
                plotvectors.append(scipy.signal.resample(radarvec,500))
            else:
                plotvectors.append(-1*scipy.signal.resample(radarvec,500))
        elif type == 'phase':
            radarvec = processing.phase_estimate(measurements[n][0].data[time_index_window[n,0]:time_index_window[n,1]])
            radarvec = radarvec-numpy.mean(radarvec)
            if normalize:
                radarvec=radarvec/(numpy.max(radarvec)-numpy.min(radarvec))
            if any([measurements[n][1] == m for m in ['b2','a2','a32']]):
            #if any([measurements[n][1] == m for m in ['a1']]):
                plotvectors.append(scipy.signal.resample(radarvec,500))
            else:
                plotvectors.append(-1*scipy.signal.resample(radarvec,500))

    #Do the waterfall:
    pylab.figure()
    pylab.title(type)
    pylab.ylabel('Measurement number')
    pylab.xlabel('Slow time (two heartbeats)')
    pylab.imshow(plotvectors,aspect='auto',interpolation='nearest',cmap=pylab.cm.jet,origin='upper')#,extent=(1,len(plotvectors),slowtimes[1],slowtimes[0]))




def waterfalls_with_pressure(measurements,time_index_window,pressures,pulses,type='instfreq',normalize=True,flipper=False):
    """
    Plots a waterfall of a slow time segment from each time instance
    """
    plotvectors=[]
    for n in range(len(measurements)):
        if type == 'instfreq':
            radarvec = processing.phase_estimate(measurements[n][0].data[time_index_window[n,0]:time_index_window[n,1]])
            radarvec = radarvec-numpy.mean(radarvec)
            if normalize:
                radarvec=radarvec/(numpy.max(radarvec)-numpy.min(radarvec))
            radarvec = numpy.diff(radarvec,n=1)*(measurements[n][0].timevector[3]-measurements[n][0].timevector[2])
            radarvec = numpy.append(radarvec,radarvec[-1])

            #if any([measurements[n][1] == m for m in ['b2','a2','a32']]):
            if any([measurements[n][1] == m for m in ['a32']]):
                plotvectors.append(scipy.signal.resample(radarvec,500))
            else:
                plotvectors.append(-1*scipy.signal.resample(radarvec,500))
        elif type == 'phase':
            radarvec = processing.phase_estimate(measurements[n][0].data[time_index_window[n,0]:time_index_window[n,1]])
            radarvec = radarvec-numpy.mean(radarvec)
            if normalize:
                radarvec=radarvec/(numpy.max(radarvec)-numpy.min(radarvec))
            #if any([measurements[n][1] == m for m in ['b2','a2','a32']]):
            if any([measurements[n][1] == m for m in ['a1']]):
                plotvectors.append(scipy.signal.resample(radarvec,500))
            else:
                plotvectors.append(-1*scipy.signal.resample(radarvec,500))

    #Do the waterfall:
    fig = pylab.figure()
    pylab.subplot(121)
    pylab.title('Blood pressure and pulse')
    diastolic = []
    systolic = []
    for p in pressures:
        if p[1][0] == '':
            systolic.append(0)
        else:
            systolic.append(p[1][0])
        if p[1][1] =='':
            diastolic.append(0)
        else:
            diastolic.append(p[1][1])
    #vertvec = numpy.arange(len(systolic),0,-1)
    vertvec = numpy.arange(0,len(systolic))
    #pylab.barh(vertvec,systolic,hold=True,height=.3,color='r',label='systolic pressure')
    #pylab.barh(vertvec,diastolic,height=.3,hold=True,label='diastolic pressure')
    #pylab.barh(vertvec-.3,pulses,hold=True,height=.3,color='k',label='Heart rate')
    pylab.plot(systolic,vertvec,'r')
    #pylab.plot(diastolic,vertvec,'b')
    pylab.plot(pulses,vertvec,'k')
    pylab.axis('tight')
    pylab.xlabel('mmHg/bpm')
    pylab.ylabel('Measurement number')
    pylab.xlim( (60,190) )
    pylab.ylim( (len(systolic)-1,0) )

    
    
    ax = pylab.subplot(122)
    #pylab.title(type)
    pylab.title('Heartbeat waveform')
    pylab.ylabel('Measurement number')
    pylab.xlabel('Two heartbeats)')
    #pylab.imshow(plotvectors,aspect='auto',interpolation='nearest',cmap=pylab.cm.jet,origin='upper')#,extent=(1,len(plotvectors),slowtimes[1],slowtimes[0]))

    ## quick fix (fjerner xticks):
    ax.imshow(plotvectors,aspect='auto',interpolation='nearest',cmap=pylab.cm.jet,origin='upper')#,extent=(1,len(plotvectors),slowtimes[1],slowtimes[0]))
    pylab.setp(ax.get_xticklabels(),visible=False)




def read_blood_pressure(measname,trykkfile):
    """Reads an ECG measurement and makes it synchronous with the radar recording.
    @param measname: The name of the measurement.
    @param trykkfile: The full path of the csv file containing blood pressures.
    """
    try:
        trykkf = open(trykkfile,'r')
        lines = trykkf.readlines()
        trykkf.close()
        pressure = ('','')
        for line in lines[1:]:
            info = line.split(';')
            if info[0] == measname:
                try:
                    systolic = int(info[1])
                    diastolic = int(info[2])
                    pressure = (systolic,diastolic)
                except:
                    print 'no blood pressure for '+info[0]
                    systolic = ""
                    diastolic= ""
    except:
        pressure = ('','')
        print 'No blood pressure file on path \n'+trykkfile
    return pressure



def pulse_estimator(ecgpeaks,fs=223.0):
    """
    An estimator of the heartbeat rate based on detected QRS complexes and the sampling frequency.
    """
    pulse = ((len(ecgpeaks)-1)*fs*60)/(ecgpeaks[-1]-ecgpeaks[0])
    return pulse







############################################
#Initializing statements:
dataset=1
#dataset=2
#dataset=3
#dataset=4
#dataset=5


flipper = False

if dataset ==1:
    filebase = os.path.join(os.environ['HOME'],'hubra','MELODY','Measurements','pnafiles','contact130412','')

    throughs21 = Measurement(filebase+'through_long','s21')
    through = numpy.mean(throughs21.data,axis=0)

    #measurementnamelist = [b for b in ['b2','b3','b4','a1','a2','a3','a4','a5','a6','a7','a8','a9','a10','a11','a12','a13','a14','a15','a16','a17','a18','a19','a20','a21','a22','a23','a24','a25','a26','a27','a28','a29','a30']]
    measurementnamelist = [b for b in ['b2','b3','b4','a2','a3','a4','a5','a6','a7','a8','a9','a10','a11','a12','a13','a14','a15','a16','a17','a18','a19','a20','a21','a22','a23','a24','a25','a26','a27','a28','a29','a30']]
    flipper=True
elif dataset ==2:
    filebase = os.path.join(os.environ['HOME'],'hubra','MELODY','Measurements','pnafiles','contact160412','')

    throughs21 = Measurement(filebase+'through','s21')
    through = numpy.mean(throughs21.data,axis=0)

    measurementnamelist = [b for b in ['b1','b2','b3','b4','a1','a2','a3','a4','a5','a6','a7','a8','a9','a10','a11','a12','a13','a14','a15','a16','a17','a18','a20','a22','a23','a24','a25','a26','a27','a28','a29','a30','a31','a32','a33','a34','a35','a36','a37','a38','a39','a40','a41','a42','a43','a44','a45',]]
    #measurementnamelist = [b for b in ['a31','a32','a33','a34','a35','a36','a37','a38','a39','a40','a41','a42','a43','a44','a45',]]
    flipper=True
elif dataset==3:
    filebase = os.path.join(os.environ['HOME'],'hubra','MELODY','Measurements','pnafiles','contact_010612','')

    throughs21 = Measurement(filebase+'through','s21')
    through = numpy.mean(throughs21.data,axis=0)

    measurementnamelist = [b for b in ['b1','b2','b3','b4','b5','a1','a2','a3','a5','a6','a7','a8','a9','a10','a11','a12','a13','a14','a15','a16','a17','a18','a19','a20','a21','a22','a23','a24','a25','a27','a28','a29','a30','a31','a32','a33']]
    flipper=True

elif dataset==4:
    #Non-contact measurements
    filebase = os.path.join(os.environ['HOME'],'hubra','MELODY','Measurements','pnafiles','ecgradar250811','')
    measurementnamelist=['15GHz_nobreath.txt']
    throughs21 = Measurement(filebase+measurementnamelist[0],'s21')
    through = numpy.mean(throughs21.data,axis=0)

elif dataset==5:
    filebase = os.path.join(os.environ['HOME'],'hubra','MELODY','Measurements','pnafiles','Contact180712','')

    throughs21 = Measurement(filebase+'b1','s11')
    through = numpy.mean(throughs21.data,axis=0)

    measurementnamelist = [b for b in ['b1','b2','b3','b4','b5','b6','b7','b8','b9','b10']]
    flipper=True
print 'Through sampling frequency: '+str(throughs21.fs)


measurements = []
modulationsizes = []
meanprofiles = []
pressures = []
heartbeat_indexes = []
pulses = []

###########################################
#Actual reading and processing of data:
sounds = []
for n in range(len(measurementnamelist)):

    #ECG reading:
    ecgfile = os.path.join(filebase,os.path.splitext(measurementnamelist[n])[0]+'_ecg_out.txt')
    tmpecg = initialize_ecg(ecgfile,int(throughs21.fs))
    measurements.append( (Measurement(filebase+measurementnamelist[n],'s21'),os.path.splitext(measurementnamelist[n])[0],initialize_ecg(ecgfile,int(throughs21.fs))) )
    pressures.append( (measurements[n][1],read_blood_pressure(measurements[n][1],filebase+'Blodtrykk.csv')) )

    heartbeat_indexes.append(ecgpeaks(measurements[n][2]))
    pulses.append(pulse_estimator(heartbeat_indexes[n],measurements[n][0].fs))



    #Ausculation:
    if dataset == 3:
        soundfile = os.path.join(filebase+'Sound',measurementnamelist[n]+'.wav')
        sounds.append(initialize_auscultation(soundfile,44100))
        #measurements[n][2] = filter_ECG(measurements[n][2],measurements[n][0].fs)
    
    ###############
    #Filtering of radar signal:
    measurements[n][0].data = processing.simple_filt(measurements[n][0].data,[0,111],measurements[n][0].fs)
    
    print measurements[n][1]+' OK!'+' sampling frequency: '+str(measurements[n][0].fs)



##########################################



pulse=0
#Plotting:
for n in range(len(measurements)):
    slowtimes = [measurements[n][0].timevector[heartbeat_indexes[n][3]] , measurements[n][0].timevector[heartbeat_indexes[n][5]+1]]
    slowtimes = [measurements[n][0].timevector[heartbeat_indexes[n][2]-50] , measurements[n][0].timevector[heartbeat_indexes[n][3]+1-20]]
    #slowtimes = [0.5,3]
    pulse = pulse_estimator(heartbeat_indexes[n],measurements[n][0].fs)

    
    #Plot of ECG and radar together:
    if any([measurements[n][1] == m for m in ['b2','a32']]): #a2?
        measurements[n][0].plot_with_ECG(measurements[n][2], slowtimes, type='phaseinstfreq', flipper=not flipper,title='Blood pressure: '+str(pressures[n][1][0])+'/'+str(pressures[n][1][1])+', pulse: '+'%.1f'%(pulse))
    else:
        measurements[n][0].plot_with_ECG(measurements[n][2], slowtimes, type='phaseinstfreq', flipper=flipper,title='Blood pressure: '+str(pressures[n][1][0])+'/'+str(pressures[n][1][1])+', pulse: '+'%.1f'%(pulse))
    pylab.axis('tight')

    #spectrogram:
    NFFT=2**5
    #measurements[n][0].spectrogram(slowtimes,NFFT=NFFT)




## #plot_onefreq(measurements,1e9,slowtimes,type='phase')
arrayofinds = numpy.array([(hbi[2],hbi[4]) for hbi in heartbeat_indexes])






pylab.show()
