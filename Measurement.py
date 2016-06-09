import numpy
import pylab
import os
import scipy
from scipy import signal
import processing

class Measurement:
    """
    A CW single frequency measurement class.
    """

    def __init__(self,datafilename,param):
        #import numpy
        #import pylab
        
        self.param=param
        fn,ext = os.path.splitext(datafilename)
        if ext == '.s2p':
            self.s2preader(datafilename)
        else:
            self.dataread(datafilename)
            self.configread(datafilename+'.conf')
        self.angle = 0
        self.path = datafilename
        tmp,self.filename = os.path.split(self.path)
        self.timevector = numpy.linspace(self.starttime,self.endtime,self.data.shape[0])
        try:
            self.angle = int(datafilename[-2:])
        except:
            #Do nothing
            self.angle=self.angle
        try:
            self.angle = int(datafilename[-3:])
        except:
            #Do nothing
            self.angle=self.angle

    def configread(self,filename):
        """Reads the config file"""

        self.freq = 0.0
        self.number_of_points = 0
        self.starttime = 0.0
        self.endtime = 0.0
        self.fs = 0.0
        
        file = open(filename,'r')
        lines=file.readlines()
        file.close()
        for line in lines:
            if line[0:10] == '#Frequency':
                self.freq = float(line.split(':')[1])
            elif line[0:17] == '#Number of points':
                self.number_of_points = int(line.split(':')[1])
            elif line[0:10] == '#Timestamp':
                timestamp = float(line.split(':')[1])
            elif line[0:11] == '#Start time':
                self.starttime = float(line.split(':')[1])/timestamp
            elif line[0:9] == '#End time':
                self.endtime = float(line.split(':')[1])/timestamp-self.starttime
                self.starttime=0.0
        self.fs = self.number_of_points/(self.endtime-self.starttime)



    def dataread(self,filename):
        """Reads the datafile"""
        file=open(filename,'r')
        lines=file.readlines()
        file.close()
        #Format of values is as following:
        #Re(S11) Im(S11) Re(S21) Im(S21)
        values=[]
        for line in lines[1:]:
            onetime=[float(x) for x in line.split(',')]
            values.append(onetime)
            
        if self.param == 's21':
            re=[x[2] for x in values]
            im=[x[3] for x in values]
            self.data = numpy.array(re)+1j*numpy.array(im)
        elif self.param == 's11':
            re=[x[0] for x in values]
            im=[x[1] for x in values]
            self.data = numpy.array(re)+1j*numpy.array(im)

    def s2preader(self,filename):
        """Reads s2p dataformat, used when saving directly with the pna and not using the gui."""
        file=open(filename,'r')
        lines=file.readlines()
        file.close()
        #Format of values is as following:
        #Time Re(S11) Im(S11) Re(S21) Im(S21) Re(S12) Im(S12) Re(S22) Im(S22)
        values=[]
        for line in lines[9:]:
            onetime=[float(x) for x in line.split()]
            values.append(onetime)
        if self.param=='s21':
            re=[x[3] for x in values]
            im=[x[4] for x in values]
        elif self.param=='s11':
            re=[x[1] for x in values]
            im=[x[2] for x in values]
        self.data = numpy.array(re)+1j*numpy.array(im)
        self.starttime = float(values[0][0])
        self.endtime = float(values[-1][0])
        self.number_of_points = len(values)
        self.fs = self.number_of_points/(self.endtime-self.starttime)
        #######HARDCODED:
        self.freq = 15e9
        ################
        print 'length of measurement: '+str(self.endtime-self.starttime)



    def write(self):
        """Writes the Measurements parameters"""
        print 'Measurement type: '+str(self.param)
        print 'Frequency: '+str(self.freq*1e-9)+' GHz'
        print 'number of points: '+str(self.number_of_points)
        print 'Measurement duration: '+str(self.endtime-self.starttime)+' seconds'
        print 'Sampling frequency: '+str(self.fs)+' Hz'



    def heartbeatfilter(self):
        """Detrending and highpass filtering to extract the heartbeats."""
        re = signal.detrend(numpy.real(self.data))
        im = signal.detrend(numpy.imag(self.data))
        self.data = processing.simple_filt(re+1j*im,[0.5,0],self.fs)
        self.remove_part()



    def heartbeatrcs(self,rcs_theory):
        """Computes the RCS of the physiological movement in the calibrated data"""
        #filters before returning rcs:
        #re = signal.detrend(numpy.real(self.data))
        #im = signal.detrend(numpy.imag(self.data))
        #tmp = processing.simple_filt(re+1j*im,[0.5,0],self.fs)

        self.rcs = numpy.mean(abs(self.data)**2)*rcs_theory
        return self.rcs


    def plot(self,lim='fixed',title='freq'):
        """Plotting of selected data parts"""
        timevector = numpy.linspace(self.starttime,self.endtime,self.data.shape[0]) #self.number_of_points)
        
##         pylab.figure()
##         pylab.subplot(211)
##         pylab.plot(timevector,20*numpy.log10(self.data))
##         pylab.title('Amplitude, frequency '+str(self.freq*1e-9)+'GHz, angle: '+str(self.angle))
##         pylab.xlabel('Time (s)')
##         pylab.ylabel('Amplitude in dB')
##         pylab.subplot(212)
##         pylab.plot(timevector,180.0*numpy.unwrap(numpy.angle(self.data))/numpy.pi)
##         pylab.title('phase')
##         pylab.xlabel('Time (s)')
##         pylab.ylabel('degrees')

        #add polar plot

        pylab.figure()
##        if title=='freq':
##            pylab.title('Real part, frequency '+str(self.freq*1e-9)+'GHz, angle: '+str(self.angle))
##        else:

##        pylab.subplot(311)
##        pylab.polar(numpy.angle(self.data),abs(self.data))
        
        pylab.subplot(211)
        #pylab.title(self.filename+', real part')
        pylab.title('Real part')
        pylab.plot(timevector,numpy.real(self.data))
        pylab.xlabel('Time (s)')
        if lim=='fixed':
            pylab.ylim(-0.6,0.6)
            #pylab.ylim(-0.006,0.006)
            #pylab.ylim(-0.001,0.001)
        else:
            pylab.axis('tight')
        pylab.ylabel('Voltage (V)')
        pylab.subplot(212)
        pylab.plot(timevector,numpy.imag(self.data))
        pylab.title('Imaginary part')
        pylab.xlabel('Time (s)')
        if lim=='fixed':
            pylab.ylim(-0.6,0.6)
            #pylab.ylim(-0.006,0.006)
            #pylab.ylim(-0.001,0.001)
        else:
            pylab.axis('tight')
        pylab.ylabel('Voltage (V)')



    def polarplot(self,slowtimes,title='none'):
        """
        Polar plot of the data vector.
        """
        startind = numpy.where(self.timevector<=slowtimes[0])[0][-1]
        stopind = numpy.where(self.timevector>=slowtimes[1])[0][0]
        center,tmp = processing.circlecenter(self.data[startind:stopind])
        radarvec = self.data[startind:stopind] #- center
        
        pylab.figure()
        pylab.polar(numpy.unwrap(numpy.angle(radarvec)),abs(radarvec))
        if title=='none':
            pylab.title(self.filename)
        else:
            pylab.title(title)



    def remove_part(self):
        """Searches throught the file 'baddata' for name matching self.name. Removes the parts specified there."""
        dirname,name = os.path.split(self.path)

        try:
            file=open(os.path.join(dirname,'baddata'),'r')
            lines=file.readlines()
            file.close()
            for line in lines:
                info = line.split(';')
                if name == info[0]:
                    startsec = float(info[1])
                    stopsec = float(info[2])
                    startind = int(startsec*self.fs)
                    if stopsec == 0:
                        self.data = self.data[0:startind]
                    else:
                        stopind = int(stopsec*self.fs)
                        if startind == 0:
                            self.data = self.data[stopind:]
                        else:
                            self.data = numpy.concatenate((self.data[0:startind],self.data[stopind-1:]))

        except:
            print 'No baddata file found for '+str(self.path)
                #Gives problem with filtering:
                ## self.data[startind:stopind]=float('nan')



    def rotate(self):
        """
        Rotates the modulation to lie along the real axis
        """
        meandata = numpy.mean(self.data)
        inputdata = self.data-meandata
        maxval = 0.3*numpy.mean(abs(inputdata))#Could be a problem. Use abs?

        #Alternative approach:
##         bigvals = inputdata[abs(inputdata)>maxval]
##         posvals = bigvals[numpy.imag(bigvals)>0]
##         negvals = bigvals[numpy.imag(bigvals)<0]    
##         angles = numpy.angle(posvals)
##         angles = numpy.append(angles,numpy.pi+numpy.angle(negvals))
##         meanangle = numpy.mean(angles)

        #Alternative approach:
        #Makes x,y coordinates an sorts in increasing order:
        verdier = [ (x,y) for x,y in zip(numpy.real(inputdata),numpy.imag(inputdata))]
        verdier = sorted(verdier, key=lambda v: v[0])
        verdier = numpy.array(verdier)
        A = numpy.vstack([verdier[:,0], numpy.ones(verdier.shape[0])]).T
        a, b = numpy.linalg.lstsq(A, verdier[:,1])[0]
        #finner vinkelestimat fra a,b:
        meanangle = numpy.arctan(a)


            
        ##############
        
        print 180.0*meanangle/numpy.pi
        
        self.data = numpy.exp(-1j*meanangle)*inputdata







    def plot_with_ECG(self,ecg,slowtimes,type='real',flipper=False,title='none'):
        """
        Plots the radar on top of the ecg for the given slowtimes. 
        @param ecg: The vector containing ECG data. is assumed to have the same fs as self.
        @param slowtimes: Two element list containing start and stop slowtimes for the plot
        """

        #slowtime indexes computed:
        startind = numpy.where(self.timevector<=slowtimes[0])[0][-1]
        stopind = numpy.where(self.timevector>=slowtimes[1])[0][0]

        c=3e8
        yt = None
        yloc = None
        doplot = True

        ylab='arbitrary units'
        #Add more alternatives later:
        if type=='real':
            radarvec = numpy.real(self.data[startind:stopind])
            ylab = 'real value'
        elif type=='abs':
            radarvec = abs(self.data[startind:stopind])
            ylab = 'absolute value'
        elif type=='dbreal':
            #UPdate!
            dyn = 50
            radarvec = numpy.real(self.data[startind:stopind])
            radarvec = processing.makedB(radarvec,dynamicdb = dyn)
            #OBS:
            #manipulate axes to have correct labels: 20log(max) to 20log(max)-50 for both positive and negative.
            ylab = 'real value, dB scale'
##             yloc = [numpy.min(radarvec),0,numpy.max(radarvec)]
##             yt = [numpy.max(radarvec)-dyn,0,numpy.max(radarvec)-dyn]
        elif type=='phaseapprox':
            #The linear demodulation approximation to phase
            radarvec = numpy.real(self.rotate()[startind:stopind])
            ylab = 'linear phase approximation(arbitrary units)'
            ## self.polarplot(1,slowtimes,altvec=radarvec)
        elif type == 'instfreq':
            #Instantaneous frequency:
            #phase = numpy.real(self.rotate(tauind)[startind:stopind])
            radarvec = processing.phase_estimate(self.data[startind:stopind])
            radarvec = radarvec-numpy.mean(radarvec)
            radarvec = numpy.diff(radarvec,n=1) #First order derivative
            radarvec = numpy.append(radarvec,radarvec[-1])/(self.timevector[3]-self.timevector[2])#*c/(2*2*numpy.pi*self.freq)
            ylab = 'Instantaneous frequency'
        elif type == 'phase':
            radarvec = processing.phase_estimate(self.data[startind:stopind])
            radarvec = radarvec-numpy.mean(radarvec)
            ylab = 'phase'
        elif type == 'movement':
            radarvec = processing.phase_estimate(self.data[startind:stopind])
            radarvec = radarvec-numpy.mean(radarvec)
            radarvec = 1e3*radarvec*3e8/(2*2*numpy.pi*(self.stopfreq-self.startfreq))
            ylab = 'movement (mm)'
        elif type == 'phaseamp':
            center,tmp = processing.circlecenter(self.data[startind:stopind])
            radarvec = self.data[startind:stopind]-center
            radarvec2 = abs(radarvec)-numpy.mean(abs(radarvec))
            radarvec = numpy.unwrap(numpy.angle(radarvec))
            radarvec = radarvec-numpy.mean(radarvec)
            ylab='phase'
            ylab2='amplitude'

            pylab.figure()
            maxecg = max(abs(ecg))
            maxradar = max(abs(radarvec))
            maxradar2 = max(abs(radarvec2))
            pylab.subplot(211)
            if title=='none':
                pylab.title(os.path.splitext(self.filename)[0])
            else:
                pylab.title(title)
            if flipper:
                pylab.plot(self.timevector[startind:stopind],-1*radarvec)
            else:
                pylab.plot(self.timevector[startind:stopind],radarvec)
            pylab.plot(self.timevector[startind:stopind],ecg[startind:stopind]*maxradar/maxecg,'r')
            pylab.xlabel('Slow time (s)')
            pylab.ylabel(ylab)
            pylab.axis('tight')
            
            pylab.subplot(212)
            pylab.plot(self.timevector[startind:stopind],radarvec2)
            pylab.plot(self.timevector[startind:stopind],ecg[startind:stopind]*maxradar2/maxecg,'r')
            pylab.xlabel('Slow time (s)')
            pylab.ylabel(ylab2)
            pylab.axis('tight')

            doplot=False
            
        elif type == 'phaseinstfreq':
            center,tmp = processing.circlecenter(self.data[startind:stopind])
            radarvec = self.data[startind:stopind]-center
            radarvec = numpy.unwrap(numpy.angle(radarvec))
            radarvec = radarvec-numpy.mean(radarvec)
            radarvec2 = numpy.diff(radarvec,n=1) #First order derivative
            radarvec2 = numpy.append(radarvec2,radarvec2[-1])/(self.timevector[3]-self.timevector[2])#*c/(2*2*numpy.pi*self.freq)
            ylab=r'$\phi$'
            ylab2=r'$d\phi/dt$'

            pylab.figure()
            maxecg = max(abs(ecg))
            maxradar = max(abs(radarvec))
            maxradar2 = max(abs(radarvec2))
            pylab.subplot(211)
            if title=='none':
                pylab.title(os.path.splitext(self.filename)[0])
            else:
                pylab.title(title)
            if flipper:
                pylab.plot(self.timevector[startind:stopind],-1*radarvec)
            else:
                pylab.plot(self.timevector[startind:stopind],radarvec)
            pylab.plot(self.timevector[startind:stopind],ecg[startind:stopind]*maxradar/maxecg,'r')
            pylab.xlabel('Slow time (s)')
            pylab.ylabel(ylab)
            pylab.axis('tight')
            
            pylab.subplot(212)
            if flipper:
                pylab.plot(self.timevector[startind:stopind],-1*radarvec2)
            else:
                pylab.plot(self.timevector[startind:stopind],radarvec2)
            pylab.plot(self.timevector[startind:stopind],ecg[startind:stopind]*maxradar2/maxecg,'r')
            pylab.xlabel('Slow time (s)')
            pylab.ylabel(ylab2)
            pylab.axis('tight')

            doplot=False

           
            
        if doplot:
            maxecg = max(abs(ecg))
            maxradar = max(abs(radarvec))
            
            pylab.figure()
            if title=='none':
                pylab.title(os.path.splitext(self.filename)[0]+' '+type)
            else:
                pylab.title(title)
            if flipper:
                pylab.plot(self.timevector[startind:stopind],-1*radarvec)
            else:
                pylab.plot(self.timevector[startind:stopind],radarvec)
            pylab.plot(self.timevector[startind:stopind],ecg[startind:stopind]*maxradar/maxecg,'r')
            pylab.xlabel('Slow time (s)')
            pylab.ylabel(ylab)
            pylab.axis('tight')
            if yt:
                pylab.yticks(yloc,yt)




    def plot_with_ECG_and_sound(self,ecg,sound,slowtimes,type='real',flipper=False,title='none'):
        """
        Plots the radar on top of the ecg for the given slowtimes. 
        @param ecg: The vector containing ECG data. is assumed to have the same fs as self.
        @param slowtimes: Two element list containing start and stop slowtimes for the plot
        """

        #slowtime indexes computed:
        startind = numpy.where(self.timevector<=slowtimes[0])[0][-1]
        stopind = numpy.where(self.timevector>=slowtimes[1])[0][0]

        c=3e8
        yt = None
        yloc = None
        doplot = True


        fsfactor = int(44100/self.fs)
        soundstart = int(self.timevector[startind]*44100)
        soundstop = int(self.timevector[stopind]*44100)
        sound_timevector = numpy.linspace(self.timevector[startind],self.timevector[stopind],soundstop-soundstart)


 
        ylab='arbitrary units'
        #Add more alternatives later:
            
        if type == 'phaseinstfreq':
            center,tmp = processing.circlecenter(self.data[startind:stopind])
            radarvec = self.data[startind:stopind]-center
            radarvec = numpy.unwrap(numpy.angle(radarvec))
            radarvec = radarvec-numpy.mean(radarvec)
            radarvec2 = numpy.diff(radarvec,n=1) #First order derivative
            radarvec2 = numpy.append(radarvec2,radarvec2[-1])/(self.timevector[3]-self.timevector[2])#*c/(2*2*numpy.pi*self.freq)
            ylab=r'$\phi$'
            ylab2=r'$d\phi/dt$'

            pylab.figure()
            maxecg = max(abs(ecg))
            maxradar = max(abs(radarvec))
            maxradar2 = max(abs(radarvec2))
            maxsound = max(abs(sound[soundstart:soundstop]))
            pylab.subplot(211)
            if title=='none':
                pylab.title(os.path.splitext(self.filename)[0])
            else:
                pylab.title(title)

            pylab.plot(self.timevector[startind:stopind],ecg[startind:stopind]*maxradar/maxecg,'r')
            pylab.plot(sound_timevector,sound[soundstart:soundstop]*0.4*maxradar/maxsound,'g')
            if flipper:
                pylab.plot(self.timevector[startind:stopind],-1*radarvec)
            else:
                pylab.plot(self.timevector[startind:stopind],radarvec)
           
            pylab.xlabel('Slow time (s)')
            pylab.ylabel(ylab)
            pylab.axis('tight')
            
            pylab.subplot(212)
            pylab.plot(self.timevector[startind:stopind],ecg[startind:stopind]*maxradar2/maxecg,'r')
            pylab.plot(sound_timevector,sound[soundstart:soundstop]*0.4*maxradar2/maxsound,'g')
            if flipper:
                pylab.plot(self.timevector[startind:stopind],-1*radarvec2)
            else:
                pylab.plot(self.timevector[startind:stopind],radarvec2)

            pylab.xlabel('Slow time (s)')
            pylab.ylabel(ylab2)
            pylab.axis('tight')

            doplot=False

    






    def plot_spectrum(self,flim=40,dyn=100,process='none'):
        """
        Plots the frequency spectrum of the data.
        @param flim: A tuple containing the frequency range to plot.
        @param dyn: the dynamic range of the plot in dB.
        @param process: Determines the amount of processing applied to the data before plotting.
        """
        padlen = numpy.max(2**12,len(self.data))
        if process=='none':
            invector=self.data
        elif process=='meanremove':
            invector = self.data - numpy.mean(self.data)
        elif process == 'highpass':
            invector = self.data - numpy.mean(self.data) #Update this if needed!
        spectrum = scipy.fftpack.fftshift(scipy.fft(processing.window(invector,0,len(invector)),padlen))
        freqvec = numpy.linspace(-self.fs/2,self.fs/2,len(spectrum))
        indexes = (freqvec>=-flim)*(freqvec<=flim)
        
        pylab.figure()
        pylab.plot(freqvec[indexes],20*numpy.log10(abs(spectrum[indexes])))
        pylab.ylim( (numpy.max(spectrum)-dyn,numpy.max(spectrum)) )
        pylab.xlabel('Frequency (Hz)')
        pylab.ylabel('Magnitude (dB)')
        

    def spectrogram(self,slowtimes,NFFT=2**5,type='phase'):
        """
        Plots the spectrogram of the data
        """

        #slowtime indexes computed:
        startind = numpy.where(self.timevector<=slowtimes[0])[0][-1]
        stopind = numpy.where(self.timevector>=slowtimes[1])[0][0]

        c=3e8
        yt = None
        yloc = None
        doplot = True

        ylab='arbitrary units'
        #Add more alternatives later:
        if type=='real':
            radarvec = numpy.real(self.data[startind:stopind])
            ylab = 'real value'
        elif type=='raw':
            radarvec = self.data[startind:stopind]
        elif type=='abs':
            radarvec = abs(self.data[startind:stopind])
            ylab = 'absolute value'
        elif type=='phaseapprox':
            #The linear demodulation approximation to phase
            radarvec = numpy.real(self.rotate()[startind:stopind])
            ylab = 'linear phase approximation(arbitrary units)'
            ## self.polarplot(1,slowtimes,altvec=radarvec)
        elif type == 'instfreq':
            #Instantaneous frequency:
            #phase = numpy.real(self.rotate(tauind)[startind:stopind])
            radarvec = processing.phase_estimate(self.data[startind:stopind])
            radarvec = radarvec-numpy.mean(radarvec)
            radarvec = numpy.diff(radarvec,n=1) #First order derivative
            radarvec = numpy.append(radarvec,radarvec[-1])/(self.timevector[3]-self.timevector[2])#*c/(2*2*numpy.pi*self.freq)
            ylab = 'Instantaneous frequency'
        elif type == 'phase':
            radarvec = processing.phase_estimate(self.data[startind:stopind])
            radarvec = radarvec-numpy.mean(radarvec)
            ylab = 'phase'

        if doplot:
            pylab.figure()
            pylab.title('Spectrogram using '+str(NFFT)+' data points')
            Pxx,freqs,t,im = pylab.specgram(radarvec,NFFT=NFFT,Fs=self.fs, noverlap=NFFT/4)
            pylab.ylabel('Frequency (Hz)')
            pylab.xlabel('time (s)')




class UWBMeasurement(Measurement):
    """A sub-class of Measurement."""
    
    def configread(self,filename):
        """Reads the config file"""
        
        self.startfreq = 0.0
        self.stopfreq = 0.0
        self.number_of_points = 0
        self.groups = 0
        self.starttime = 0.0
        self.endtime = 0.0
        self.fs = 0.0
        
        file = open(filename,'r')
        lines=file.readlines()
        file.close()
        for line in lines:
            if line[0:16] == '#Start Frequency':
                self.startfreq = float(line.split(':')[1])
            elif line[0:15] == '#Stop Frequency':
                self.stopfreq = float(line.split(':')[1])
            elif line[0:17] == '#Number of points':
                self.number_of_points = int(line.split(':')[1])
            elif line[0:17] == '#Number of sweeps':
                self.groups = int(line.split(':')[1])
            elif line[0:10] == '#Timestamp':
                timestamp = float(line.split(':')[1])
            elif line[0:11] == '#Start time':
                self.starttime = float(line.split(':')[1])/timestamp
            elif line[0:9] == '#End time':
                self.endtime = float(line.split(':')[1])/timestamp-self.starttime
                self.starttime=0.0
        self.fs =self.groups/(self.endtime-self.starttime)
    
        #modifies self.data:
        self.data = numpy.reshape(self.data,(self.groups,self.number_of_points))


    def s2preader(self,filename):
        """Reads s2p dataformat, used when saving directly with the pna and not using the gui."""
        file=open(filename,'r')
        lines=file.readlines()
        file.close()
        #Format of values is as following:
        #frequency Re(S11) Im(S11) Re(S21) Im(S21) Re(S12) Im(S12) Re(S22) Im(S22)
        values=[]
        for line in lines[9:]:
            onefreq=[float(x) for x in line.split()]
            values.append(onefreq)

        #Parameters are set:
        self.startfreq = values[0][0]
        self.stopfreq = values[-1][0]
        self.number_of_points = len(values)

        #Data is read:
        if self.param=='s21':
            re=[x[3] for x in values]
            im=[x[4] for x in values]
        elif self.param=='s11':
            re=[x[1] for x in values]
            im=[x[2] for x in values]
        elif self.param=='s12':
            re=[x[5] for x in values]
            im=[x[6] for x in values]
        elif self.param=='s22':
            re=[x[7] for x in values]
            im=[x[8] for x in values]
            
        self.data = numpy.array(re)+1j*numpy.array(im)
        self.number_of_points = len(values)

        #A bit shady this one:
        self.starttime = 0.0
        self.endtime = 1.0
        
        #From the single freq function:
        #self.data = numpy.array(re)+1j*numpy.array(im)
        #self.starttime = float(values[0][0])
        #self.endtime = float(values[-1][0])
        #self.number_of_points = len(values)
        #self.fs = self.number_of_points/(self.endtime-self.starttime)



    def make_tau_profiles(self,start_tau,stop_tau,tau_res,window='hamming'):
        """Performs the chirpz transform along frequency, makes the taudata array.
        @param start_tau: The start of the fast time zoom window in nanoseconds.
        @param stop_tau: The end of the fast time zoom window in nanoseconds.
        @param tau_res: The wanted fast time resolution in nanoseconds."""
        #Ts*fs = self.number_of_points.
        import chirpz
        #####################################
        #Experimental:
        #Adds zeros at the start of the frequencies:
        zeros = numpy.zeros([self.groups,int(self.startfreq*self.number_of_points/(self.stopfreq-self.startfreq))])
        self.data = numpy.append(zeros,self.data,axis=1)
        self.startfreq = 0
        self.number_of_points = self.data.shape[1]
        ####################################
        
        self.tauprofiles = chirpz.zoom_ifft(processing.window(self.data,0,self.data.shape[1],type='hamming'),start_tau,stop_tau,tau_res,self.stopfreq-self.startfreq,self.number_of_points,1,self.number_of_points)
        self.start_tau = start_tau
        self.stop_tau = stop_tau
        self.taures = tau_res
        self.tauvec = numpy.linspace(self.start_tau,self.stop_tau,self.tauprofiles.shape[1])
        return self.tauprofiles




    def filter_profiles(self,passband):
        """Filter the tau profiles along slow time"""
        if passband[0]==0:
            self.tauprofiles = processing.simple_filt(self.tauprofiles,passband,self.fs)
        else:
            self.tauprofiles = self.tauprofiles - numpy.mean(self.tauprofiles,axis=0)
            self.tauprofiles = processing.simple_filt(self.tauprofiles,passband,self.fs)





    def rotate(self,tauind,slowtimeinds=[0,-1]):
        """
        Rotates the modulation of data through slow time at one fast time given by tauind to lie along the real axis
        @param tauind: The fast time index.
        @param slowtimeinds: Two element list containing start and stop indexes of the slow time.
        """
        meandata = numpy.mean(self.tauprofiles[:,tauind])
        inputdata = self.tauprofiles[:,tauind]-meandata
        maxval = 0.3*numpy.mean(abs(inputdata))#Could be a problem. Use abs?

        #Alternative approach:
        #Makes x,y coordinates and sorts in increasing order:
        verdier = [ (x,y) for x,y in zip(numpy.real(inputdata),numpy.imag(inputdata))]
        verdier = sorted(verdier, key=lambda v: v[0])
        verdier = numpy.array(verdier)
        A = numpy.vstack([verdier[:,0], numpy.ones(verdier.shape[0])]).T
        a, b = numpy.linalg.lstsq(A, verdier[:,1])[0]
        #finner vinkelestimat fra a,b:
        meanangle = numpy.arctan(a)
            
        ##############
        
        print 180.0*meanangle/numpy.pi
        
        return numpy.exp(-1j*meanangle)*inputdata
        
            


    def waterfall(self,slowtimes,type = 'real',vmin='None',vmax='None',title='none'):
        """Waterfall plot of the tau profiles
        @param slowtimes: Two element list containing start and stop slowtimes
        """
        print 'Dont go chasing waterfalls!'
        pylab.figure()
        if title=='none':
            pylab.title(self.filename)
        else:
            pylab.title(title)
        pylab.xlabel('Fast time (ns)')
        pylab.ylabel('Slow time (s)')

        #slowtime indexes computed:
        startind = numpy.where(self.timevector<=slowtimes[0])[0][-1]
        stopind = numpy.where(self.timevector>=slowtimes[1])[0][0]
        
        if type=='real':
            if vmin != 'None' and vmax != 'None':
                pylab.imshow(numpy.real(self.tauprofiles[startind:stopind,:]),aspect='auto',cmap=pylab.cm.jet,origin='upper',extent=(self.start_tau,self.stop_tau,slowtimes[1],slowtimes[0]),vmin=vmin,vmax=vmax)
            else:
                pylab.imshow(numpy.real(self.tauprofiles[startind:stopind,:]),aspect='auto',cmap=pylab.cm.jet,origin='upper',extent=(self.start_tau,self.stop_tau,slowtimes[1],slowtimes[0]))
        elif type=='abs':
            if vmin != 'None' and vmax != 'None':
                pylab.imshow(abs(self.tauprofiles[startind:stopind,:]),aspect='auto',cmap=pylab.cm.jet,origin='upper',extent=(self.start_tau,self.stop_tau,slowtimes[1],slowtimes[0]),vmin=vmin,vmax=vmax)
            else:
                pylab.imshow(abs(self.tauprofiles[startind:stopind,:]),aspect='auto',cmap=pylab.cm.jet,origin='upper',extent=(self.start_tau,self.stop_tau,slowtimes[1],slowtimes[0]))
        elif type=='logabs':
            if vmin != 'None' and vmax != 'None':
                pylab.imshow(20*numpy.log10(abs(self.tauprofiles[startind:stopind,:])),aspect='auto',cmap=pylab.cm.jet,origin='upper',extent=(self.start_tau,self.stop_tau,slowtimes[1],slowtimes[0]),vmin=vmin,vmax=vmax)
            else:
                pylab.imshow(20*numpy.log10(abs(self.tauprofiles[startind:stopind,:])),aspect='auto',cmap=pylab.cm.jet,origin='upper',extent=(self.start_tau,self.stop_tau,slowtimes[1],slowtimes[0]))
        elif type=='imag':
            if vmin != 'None' and vmax != 'None':
                pylab.imshow(numpy.imag(self.tauprofiles[startind:stopind,:]),aspect='auto',cmap=pylab.cm.jet,origin='upper',extent=(self.start_tau,self.stop_tau,slowtimes[1],slowtimes[0]),vmin=vmin,vmax=vmax)
            else:
                pylab.imshow(numpy.imag(self.tauprofiles[startind:stopind,:]),aspect='auto',cmap=pylab.cm.jet,origin='upper',extent=(self.start_tau,self.stop_tau,slowtimes[1],slowtimes[0]))
        elif type == 'dbreal':
            plotarray = processing.makedB(numpy.real(self.tauprofiles[startind:stopind,:]),40)
            if vmin != 'None' and vmax != 'None':
                pylab.imshow(plotarray,aspect='auto',cmap=pylab.cm.jet,origin='upper',extent=(self.start_tau,self.stop_tau,slowtimes[1],slowtimes[0]),vmin=vmin,vmax=vmax)
            else:
                pylab.imshow(plotarray,aspect='auto',cmap=pylab.cm.jet,origin='upper',extent=(self.start_tau,self.stop_tau,slowtimes[1],slowtimes[0]))


        pylab.colorbar()
                    





    def UWBplot(self):
        """
        Some selected plots of the UWB measurement
        """
        self.movementindex = 723
        dirname,name = os.path.split(self.path)
        timevector = numpy.linspace(self.starttime,self.endtime,self.data.shape[0])
        
        meandata = numpy.mean(self.data,axis=0)
        freq_ax = numpy.linspace(self.startfreq,self.stopfreq,self.number_of_points)

        
        pylab.figure()
        pylab.plot(freq_ax,meandata)
        pylab.title('Mean raw data')

##        pylab.figure()
##        pylab.title('Raw data')
##        for g in range(self.groups):
##            pylab.plot(freq_ax,self.data[g,:])

        padlen = 2**12
        datafft = scipy.ifft(self.data,padlen,axis=1)
        meanfft = numpy.mean(datafft,axis=0)

        tauvec = 1e9*numpy.linspace(0,self.number_of_points/(self.stopfreq-self.startfreq),padlen)

        pylab.figure()
        pylab.title('Range profiles '+str(name))
        if datafft.shape[0]<1000:
            pylab.plot(tauvec,(abs(datafft[0:-1:21,:].T)))
        elif datafft.shape[0]<5000:
            pylab.plot(tauvec,(abs(datafft[0:-1:35,:].T)))
        else:
            pylab.plot(tauvec,(abs(datafft[0:-1:100,:].T)))
        pylab.xlabel('ns')
        

##        pylab.figure()
##        pylab.title('Mean range profile')
##        pylab.plot(20*numpy.log10(abs(meanfft)))

        pylab.figure()
        pylab.title('Heartbeat movemen, '+str(name))
        pylab.plot(timevector,datafft[:,self.movementindex])






    def plot_with_ECG(self,ecg,tau,slowtimes,type='real',flipper=False,title=''):
        """
        Plots the radar on top of the ecg for the given slowtimes. The radar vector at fast time tau is used.
        @param ecg: The vector containing ECG data. is assumed to have the same fs as self.
        @param tau: The fast time in nanoseconds where the slowtime vector is used.
        @param slowtimes: Two element list containing start and stop slowtimes for the plot
        """

        #slowtime indexes computed:
        startind = numpy.where(self.timevector<=slowtimes[0])[0][-1]
        stopind = numpy.where(self.timevector>=slowtimes[1])[0][0]
        tauind = numpy.where(self.tauvec>=tau)[0][0]

        yt = None
        yloc = None
        doplot = True

        ylab='arbitrary units'
        #Add more alternatives later:
        if type=='real':
            radarvec = numpy.real(self.tauprofiles[startind:stopind,tauind])
            ylab = 'real value'
        elif type=='abs':
            radarvec = abs(self.tauprofiles[startind:stopind,tauind])
            ylab = 'absolute value'
        elif type=='dbreal':
            #UPdate!
            dyn = 50
            radarvec = numpy.real(self.tauprofiles[startind:stopind,tauind])
            radarvec = processing.makedB(radarvec,dynamicdb = dyn)
            #OBS:
            #manipulate axes to have correct labels: 20log(max) to 20log(max)-50 for both positive and negative.
            ylab = 'real value, dB scale'
##             yloc = [numpy.min(radarvec),0,numpy.max(radarvec)]
##             yt = [numpy.max(radarvec)-dyn,0,numpy.max(radarvec)-dyn]
        elif type=='phaseapprox':
            #The linear demodulation approximation to phase
            radarvec = numpy.real(self.rotate(tauind)[startind:stopind])
            ylab = 'linear phase approximation(arbitrary units)'
            ## self.polarplot(1,slowtimes,altvec=radarvec)
        elif type == 'instfreq':
            #Instantaneous frequency:
            #phase = numpy.real(self.rotate(tauind)[startind:stopind])
            radarvec = processing.phase_estimate(self.tauprofiles[startind:stopind,tauind])
            radarvec = radarvec-numpy.mean(radarvec)
            radarvec = numpy.diff(radarvec,n=1) #First order derivative
            radarvec = numpy.append(radarvec,radarvec[-1])/(self.timevector[3]-self.timevector[2])#*c/(2*2*numpy.pi*(self.stopfreq+self.startfreq)/2)
            ylab = 'Instantaneous frequency'
        elif type == 'phase':
            radarvec = processing.phase_estimate(self.tauprofiles[startind:stopind,tauind])
            radarvec = radarvec-numpy.mean(radarvec)
            ylab = 'phase'
        elif type == 'movement':
            radarvec = processing.phase_estimate(self.tauprofiles[startind:stopind,tauind])
            radarvec = radarvec-numpy.mean(radarvec)
            radarvec = 1e3*radarvec*3e8/(2*2*numpy.pi*(self.stopfreq-self.startfreq))
            ylab = 'movement (mm)'
        elif type == 'phaseamp':
            center,tmp = processing.circlecenter(self.tauprofiles[startind:stopind,tauind])
            radarvec = self.tauprofiles[startind:stopind,tauind]-center
            radarvec2 = abs(radarvec)-numpy.mean(abs(radarvec))
            radarvec = numpy.unwrap(numpy.angle(radarvec))
            radarvec = radarvec-numpy.mean(radarvec)
            ylab='phase'
            ylab2='amplitude'

            pylab.figure()
            pylab.title(os.path.splitext(self.filename)[0]+' '+type)
            maxecg = max(abs(ecg))
            maxradar = max(abs(radarvec))
            maxradar2 = max(abs(radarvec2))
            pylab.subplot(211)
            pylab.title(os.path.splitext(self.filename)[0])
            if flipper:
                pylab.plot(self.timevector[startind:stopind],-1*radarvec)
            else:
                pylab.plot(self.timevector[startind:stopind],radarvec)
            pylab.plot(self.timevector[startind:stopind],ecg[startind:stopind]*maxradar/maxecg,'r')
            pylab.xlabel('Slow time (s)')
            pylab.ylabel(ylab)
            pylab.axis('tight')
            
            pylab.subplot(212)
            pylab.plot(self.timevector[startind:stopind],radarvec2)
            pylab.plot(self.timevector[startind:stopind],ecg[startind:stopind]*maxradar2/maxecg,'r')
            pylab.xlabel('Slow time (s)')
            pylab.ylabel(ylab2)
            pylab.axis('tight')

            doplot=False

        elif type == 'phaseinstfreq':
            center,tmp = processing.circlecenter(self.tauprofiles[startind:stopind,tauind])
            radarvec = self.tauprofiles[startind:stopind,tauind]-center
            radarvec = numpy.unwrap(numpy.angle(radarvec))
            radarvec = radarvec-numpy.mean(radarvec)
            radarvec2 = numpy.diff(radarvec,n=1)
            radarvec2 = numpy.append(radarvec2,radarvec2[-1])/(self.timevector[3]-self.timevector[2])
            ylab='phase'
            ylab2='instantaneous frequency'

            pylab.figure()
            if title=='':
                pylab.title(os.path.splitext(self.filename)[0]+' '+type)
            else:
                pylab.title(title)
            maxecg = max(abs(ecg))
            maxradar = max(abs(radarvec))
            maxradar2 = max(abs(radarvec2))
            pylab.subplot(211)
            if title=='':
                pylab.title(os.path.splitext(self.filename)[0])
            else:
                pylab.title(title)
            if flipper:
                pylab.plot(self.timevector[startind:stopind],-1*radarvec)
            else:
                pylab.plot(self.timevector[startind:stopind],radarvec)
            pylab.plot(self.timevector[startind:stopind],ecg[startind:stopind]*maxradar/maxecg,'r')
            pylab.xlabel('Slow time (s)')
            pylab.ylabel(ylab)
            pylab.axis('tight')
            
            pylab.subplot(212)
            if flipper:
                pylab.plot(self.timevector[startind:stopind],-1*radarvec2)
            else:
                pylab.plot(self.timevector[startind:stopind],radarvec2)
            pylab.plot(self.timevector[startind:stopind],ecg[startind:stopind]*maxradar2/maxecg,'r')
            pylab.xlabel('Slow time (s)')
            pylab.ylabel(ylab2)
            pylab.axis('tight')

            doplot=False

           
            
        if doplot:
            maxecg = max(abs(ecg))
            maxradar = max(abs(radarvec))
            
            pylab.figure()
            if title=='':
                pylab.title(os.path.splitext(self.filename)[0]+' '+type)
            else:
                pylab.title(title)
            if flipper:
                pylab.plot(self.timevector[startind:stopind],-1*radarvec)
            else:
                pylab.plot(self.timevector[startind:stopind],radarvec)
            pylab.plot(self.timevector[startind:stopind],ecg[startind:stopind]*maxradar/maxecg,'r')
            pylab.xlabel('Slow time (s)')
            pylab.ylabel(ylab)
            pylab.axis('tight')
            if yt:
                pylab.yticks(yloc,yt)



    def polarplot(self,tau,slowtimes,lim='tight',title='',altvec='None'):
        """
        Polar plot of the data vector.
        @param tau: The fast time in ns along which the data is plotted
        @param slowtimes: A two element list containing the start and stop slowtimes
        @param altvec: An alternative vector for plotting. Will be used instead of self.tauprofiles if not None.
        """
        startind = numpy.where(self.timevector<=slowtimes[0])[0][-1]
        stopind = numpy.where(self.timevector>=slowtimes[1])[0][0]
        pylab.figure()
        if altvec == 'None':
            tauind = numpy.where(self.tauvec>=tau)[0][0]
            center,tmp = processing.circlecenter(self.tauprofiles[startind:stopind,tauind])
            radarvec = self.tauprofiles[startind:stopind,tauind] - center
            pylab.polar(numpy.unwrap(numpy.angle(radarvec)),abs(radarvec))
        else:
            pylab.polar(numpy.unwrap(numpy.angle(altvec[startind:stopind])),abs(altvec[startind:stopind]))
        pylab.title(os.path.splitext(self.filename)[0])

        ###################################################################
        #File writing, remove after use:
##         realfilename = 'phaseamp'+self.filename+'.txt'
##         realoutfile = open(realfilename,'w')
##         phase = numpy.unwrap(numpy.angle(radarvec))
##         amp = abs(radarvec)
##         realoutfile.write('phase; amplitude\n')
##         for p,a in zip(phase,amp):
##             realoutfile.write(str(p)+'; '+str(a)+'\n')
##         realoutfile.close()
        ##################################################################





    def plotprofiles(self,taus,slowtimes,type='real',num=-1):
        """
        Plotting of the tau profiles on top of each other
        @param taus: Two element array containing start and stop fast times (ns)
        @param slowtimes: Two element array containing start and stop slow times (s)
        @param type: The display type of the pulse.
        """
        #slowtime indexes computed:
        startind = numpy.where(self.timevector<=slowtimes[0])[0][-1]
        stopind = numpy.where(self.timevector>=slowtimes[1])[0][0]
        if num<=0:
            timeindexes = range(startind,stopind)
        else:
            timeindexes = range(startind,stopind,int((stopind-startind)/num))

        taustart = numpy.where(self.tauvec<=taus[0])[0][-1]
        taustop = numpy.where(self.tauvec>=taus[1])[0][0]
        tauindexes = range(taustart,taustop)
                    

        ylab='arbitrary units'

        pylab.figure()
        pylab.title(os.path.splitext(self.filename)[0])
        
        if type=='real':
            for n in timeindexes:
                pylab.plot(self.tauvec[tauindexes],numpy.real(self.tauprofiles[n,tauindexes]))
            pylab.ylabel('real value')
        elif type=='abs':
            for n in timeindexes:
                pylab.plot(self.tauvec[tauindexes],abs(self.tauprofiles[n,tauindexes]))
            pylab.ylabel('absolute value')
        elif type=='dbreal':
            #UPdate!
            dyn = 50
            for n in timeindexes:
                radarvec = numpy.real(self.tauprofiles[n,taustart:taustop])
                pylab.plot(self.tauvec[tauindexes],processing.makedB(radarvec,dynamicdb = dyn))
            #OBS:
            #manipulate axes to have correct labels: 20log(max) to 20log(max)-50 for  both positive and negative.
            pylab.ylabel('real value, dB scale')
##             yloc = [numpy.min(radarvec),0,numpy.max(radarvec)]
##             yt = [numpy.max(radarvec)-dyn,0,numpy.max(radarvec)-dyn]
        elif type=='phase':
            for n in timeindexes:
                pylab.plot(self.tauvec[tauindexes],numpy.angle(self.tauprofiles[n,tauindexes]))
            pylab.ylabel('phase')

        pylab.axis('tight')
        





    def plotsum_with_ECG(self,ecg,tau,heartindexes,type='real',flipper=False):
        """
        Plots the radar on top of the ecg for the given slowtimes. The radar vector at fast time tau is used.
        @param ecg: The vector containing ECG data. is assumed to have the same fs as self.
        @param tau: The fast time in nanoseconds where the slowtime vector is used.
        @param slowtimes: Two element list containing start and stop slowtimes for the plot
        """
        tauind = numpy.where(self.tauvec>=tau)[0][0]

        #Sums the heartbeats into two heartbeats:
        N = int(2*self.fs) #Number of points used for each heartbeat
        radarsum = numpy.zeros(N,dtype='complex')
        ecgsum = numpy.zeros(N)

        for n in range(len(heartindexes)-1):
            ecgsum += signal.resample(ecg[heartindexes[n]:heartindexes[n+1]],N)
            radarsum += signal.resample(self.tauprofiles[heartindexes[n]:heartindexes[n+1],tauind],N)
            


        yt = None
        yloc = None
        doplot = True

        ylab='arbitrary units'
        #Add more alternatives later:
        if type=='real':
            radarvec = numpy.real(radarsum)
            ylab = 'real value'
        elif type=='abs':
            radarvec = abs(radarsum)
            ylab = 'absolute value'
        elif type=='dbreal':
            #UPdate!
            dyn = 50
            radarvec = numpy.real(radarsum)
            radarvec = processing.makedB(radarvec,dynamicdb = dyn)
            #OBS:
            #manipulate axes to have correct labels: 20log(max) to 20log(max)-50 for both positive and negative.
            ylab = 'real value, dB scale'
##             yloc = [numpy.min(radarvec),0,numpy.max(radarvec)]
##             yt = [numpy.max(radarvec)-dyn,0,numpy.max(radarvec)-dyn]
        elif type == 'instfreq':
            #Instantaneous frequency:
            #phase = numpy.real(self.rotate(tauind)[startind:stopind])
            radarvec = processing.phase_estimate(radarsum)
            radarvec = radarvec-numpy.mean(radarvec)
            radarvec = numpy.diff(radarvec,n=1) #First order derivative
            radarvec = numpy.append(radarvec,radarvec[-1])/(self.timevector[3]-self.timevector[2])#*c/(2*2*numpy.pi*(self.stopfreq+self.startfreq)/2)
            ylab = 'Instantaneous frequency'
        elif type == 'phase':
            radarvec = processing.phase_estimate(radarsum)
            radarvec = radarvec-numpy.mean(radarvec)
            ylab = 'phase'
        elif type == 'movement':
            radarvec = processing.phase_estimate(radarsum)
            radarvec = radarvec-numpy.mean(radarvec)
            radarvec = 1e3*radarvec*3e8/(2*2*numpy.pi*(self.stopfreq-self.startfreq))
            ylab = 'movement (mm)'
        elif type == 'phaseamp':
            center,tmp = processing.circlecenter(radarsum)
            radarvec = radarsum-center
            radarvec2 = abs(radarvec)-numpy.mean(abs(radarvec))
            radarvec = numpy.unwrap(numpy.angle(radarvec))
            radarvec = radarvec-numpy.mean(radarvec)
            ylab='phase'
            ylab2='amplitude'

            pylab.figure()
            pylab.title(os.path.splitext(self.filename)[0]+' '+type)
            maxecg = max(abs(ecgsum))
            maxradar = max(abs(radarvec))
            maxradar2 = max(abs(radarvec2))
            pylab.subplot(211)
            pylab.title(os.path.splitext(self.filename)[0])
            if flipper:
                pylab.plot(-1*radarvec)
            else:
                pylab.plot(radarvec)
            pylab.plot(ecgsum*maxradar/maxecg,'r')
            pylab.xlabel('Slow time (s)')
            pylab.ylabel(ylab)
            pylab.axis('tight')
            
            pylab.subplot(212)
            pylab.plot(radarvec2)
            pylab.plot(ecgsum*maxradar2/maxecg,'r')
            pylab.xlabel('Slow time (s)')
            pylab.ylabel(ylab2)
            pylab.axis('tight')

            doplot=False

           
            
        if doplot:
            maxecg = max(abs(ecgsum))
            maxradar = max(abs(radarvec))
            
            pylab.figure()
            pylab.title(os.path.splitext(self.filename)[0]+' '+type)
            if flipper:
                pylab.plot(-1*radarvec)
            else:
                pylab.plot(radarvec)
            pylab.plot(ecgsum*maxradar/maxecg,'r')
            pylab.xlabel('Slow time (s)')
            pylab.ylabel(ylab)
            pylab.axis('tight')
            if yt:
                pylab.yticks(yloc,yt)

