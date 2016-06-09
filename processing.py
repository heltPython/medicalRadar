import numpy
import scipy
from scipy import signal

def simple_filt(input_signal,passband,prf,ps='pass'):
    """
    A simple filter creator that creates an iirfilter based on the passband given by passband and the prf, and filters the input signal.
    @param input_signal: The signal to be filtered.
    @param passband: The pass band list [low, high]. If low=0, the filter will be a lowpass filter, if high=0, the filter is a highpass filter.
    @param prf: The prf.
    """
    gpass=4
    gstop=40
    if ps=='pass':
        if passband[0]==0:
            #Lowpassfilter
            wp=passband[1]/(prf/2.0)
            ws=min(wp*1.5,1)
        elif passband[1]==0:
            #Highpassfilter
            wp=passband[0]/(prf/2.0)
            ws=wp*0.5
        else:
            #Bandpassfilter
            wp=[passband[0]/(prf/2.0) , passband[1]/(prf/2.0)]
            #ws=[wp[0]-min(wp[0],.25),wp[1]+.25]
            ws=[wp[0]*0.5 , min(wp[1]*1.5,1)]
    elif ps=='stop':
        ws=[passband[0]/(prf/2.0) , passband[1]/(prf/2.0)]
        wp = [ws[0]*0.8 , min(ws[1]*1.25,1)]
        
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

    ## import pylab
##    w,h=signal.freqz(b,a)
    ## pylab.figure()
##     pylab.title('Filter frequency response for '+str(passband)+' Hz filter')
##     pylab.plot(w*prf/(2*numpy.pi),20*numpy.log10(abs(h)))
##     pylab.ylim( (-60,1) )
    #print 'Largest filter freqeuncy response value: '+str(max(h))
    return heartbeats.conj().T






def least_square_fit(Input,degree):
    '''
    Least squares fit each column of input to a polynomial of degree. The output is the fitted polynomial.
    @param Input: The input matrix. The columns of the matrix will be fitted to a polynomial.
    @param degree: The degree of the polynomial. 0 means that the mean of the Input will be returned.
    '''
    K=Input.shape[0] #Number of rows in Input.
    H=numpy.ones((K,degree+1)) #The polynomial model matrix.
    for i in range(1,degree+1):
        H[:,degree-i]=numpy.power(numpy.arange(K,dtype='float')/K,i)
        #H[:,0]=numpy.arange(K,dtype='float')/K
        #x[:,1]=numpy.power(numpy.arange(K,dtype='float')/K,2)
        #x[:,2]=numpy.power(numpy.arange(K,dtype='float')/K,3)
    ## H=numpy.mat(H)
    PH=numpy.dot(numpy.dot(H,numpy.linalg.inv(numpy.dot(H.conj().T,H))),H.conj().T)
    #This takes forever if the Input matrix is large:
    return numpy.dot(PH,Input)





def zero_pad(f, N_profiles, N_freq, Nstart):
    """
    Zero pads the data from DC to minimum frequency.
    @param f: The data matrix.
    @param N_profiles: The number of profiles in the data set.
    @param N_freq: The number of frequency samples per profile.
    @param Nstart: The number of samples from DC to minimum nonzero frequency component.
    """
    #A new data vector sampled with zeros from DC to min_freq is created:
    data=numpy.zeros((N_profiles,N_freq+Nstart),dtype=float)
    data[:,Nstart:]=f
    return data




def window(f,start,stop,type='blackman'):
    """
    runs the data through a hamming window.
    @param f: The data matrix
    @param start: The start index of the hamming window.
    @param stop: The end index of the hamming window.
    """
    h=numpy.zeros(f.shape,dtype=float)

    if len(h.shape)==1:
        if type=='hamming':
            h[start:stop]=signal.hamming(stop-start)
        elif type=='blackman':
            h[start:stop]=signal.blackman(stop-start)
        elif type=='hann':
            h[start:stop]=signal.hann(stop-start)
        elif type=='blackmanharris':
            h[start:stop]=signal.blackmanharris(stop-start)
        elif type=='rectangular' or type=='rect' or type=='boxcar':
            h[start:stop]=signal.boxcar(stop-start)
    else:
        if type=='hamming':
            h[:,start:stop]=signal.hamming(stop-start)
        elif type=='blackman':
            h[:,start:stop]=signal.blackman(stop-start)
        elif type=='hann':
            h[:,start:stop]=signal.hann(stop-start)
        elif type=='blackmanharris':
            h[:,start:stop]=signal.blackmanharris(stop-start)
        elif type=='rectangular' or type=='rect' or type=='boxcar':
            h[:,start:stop]=signal.boxcar(stop-start)
    return numpy.multiply(f,h)



def moving_average(R,size,dim=1):
    """
    A moving average filter of size size applied to the data matrix R along dimension dim.
    @param R: The data matrix.
    @param size: The size of the moving average filter, preferably an odd number.
    @param dim: The dimension along R in wich the filter is applied. Defaults to dim=1.
    """
    averaged=numpy.mat(numpy.zeros(R.shape))
    if dim==1:
        for i in range(R.shape[dim]):
            if i<size/2:
                averaged[:,i]=numpy.mean(R[:,0:size+1],axis=dim)
            elif i>R.shape[dim]-size/2:
                averaged[:,i]=numpy.mean(R[:,R.shape[dim]-size:R.shape[dim]],axis=dim)
            else:
                averaged[:,i]=numpy.mean(R[:,i-size/2:i+size/2+1],axis=dim)

    elif dim==0:
        for i in range(R.shape[dim]):
            if i<size/2:
                averaged[i,:]=numpy.mean(R[0:size+1,:],axis=dim)
            elif i>R.shape[dim]-size/2:
                averaged[i,:]=numpy.mean(R[R.shape[dim]-size:R.shape[dim],:],axis=dim)
            else:
                averaged[i,:]=numpy.mean(R[i-size/2:i+size/2+1,:],axis=dim)

    return averaged



def get_envelope(R,dim=1):
    """
    Returns the complex version of the input signal R.
    @param R: The input data matrix.
    @param dim: The dimension along which the envelope is to be taken. default: dim=1
    """
    if dim==0:
        R=R.T
    if len(R.shape)==1:
        freqs=scipy.fft(R)
        length=len(R)/2
        freqs[length:]=0
        freqs[1:length]=2*freqs[1:length]
        ## freqs[1:length]=freqs[1:length]
        
        env=scipy.ifft(freqs)
    else:
        freqs=scipy.fft(R)
        length=R.shape[dim]/2
        #Something is fishy here:
        freqs[:,length:]=0
        freqs[:,1:length]=2*freqs[0,1:length]
        ## freqs[:,1:length]=freqs[0,1:length]
        
        env=scipy.ifft(freqs)
        if dim==0:
            return env.T

    return env




def calibration_simple(R, clutter, ref, plot_figures=True):
    """
    A simple function for calibrating the input matrix R with the matrices clutter and ref
    @param R: The input data matrix.
    @param clutter: An empty room measurement data matrix.
    @param ref: The reference target data matrix.
    """

    #First, averaging of the clutter and reference measurements:
    clutter=numpy.mean(clutter,axis=0)
    ref=numpy.mean(ref,axis=0)
    
    #Subtracts the empty room measurements:
    R=R-clutter
    ref=ref-clutter

    #Divides the data R with the envelope of the reference data:
    eps=.01
    ref_env=abs(get_envelope(ref))
    R_cal=R/(ref_env+eps)

    if plot_figures:
        import pylab
        import plotFunctions as pf
        pf.plot_raw_data(ref,'Reference average trace after clutter removal')
        pf.plot_raw_data(R,'Target traces after clutter removal')
        pf.plot_raw_data(clutter,'The average clutter trace.')
        pylab.figure()
        pylab.plot(ref.T)
        pylab.plot(ref_env.T)
        pylab.title('The reference trace and its envelope')

    return R_cal

    #Also include multiplying with analytical ref RCS:
    

def calibration2(R, clutter, ref, plot_figures=True):
    """
    A simple function for calibrating the input matrix R with the matrices clutter and ref. Uses method inspired by (Morgan 1994). Morgan used Tx and Rx antennaes facing directly towards each other as ref.
    @param R: The input data matrix.
    @param clutter: An empty room measurement data matrix.
    @param ref: The reference target data matrix.
    """

    #First, averaging of the clutter and reference measurements:
    clutter=numpy.mean(clutter,axis=0)
    ref=numpy.mean(ref,axis=0)
    
    #Subtracts the empty room measurements:
    R=R-clutter
    ref=ref-clutter

    #Divides the data R with the envelope of the reference data:
    eps=.01
    ref_env=abs(get_envelope(ref))
    R_cal=numpy.multiply(R,numpy.conj(ref))/(numpy.power(ref_env,2)+eps)
    

    if plot_figures:
        import pylab
        import plotFunctions as pf
        pf.plot_raw_data(ref,'Reference average trace after clutter removal')
        pf.plot_raw_data(R,'Target traces after clutter removal')
        pf.plot_raw_data(clutter,'The average clutter trace.')
        pylab.figure()
        pylab.plot(ref.T)
        pylab.plot(ref_env.T)
        pylab.title('The reference trace and its envelope')

    return R_cal



def heartbeat_computations(R,index,prf,bandpass=[.5,5]):
    """
    A function that does a collection of computations on the data matrix R, and returns a vector containing the heartbeats.
    @param R: The input data matrix.
    @param index: The fast time index of R where the person is located.
    @param prf: The PRF.
    """
    #Selects the vector containing the heartbeats:
    R_h = numpy.mat(R[:,index]).H
    ## R_h = R[:,index].conj().T
    R_tilde = R_h-least_square_fit(R_h,1)
    #Proever detrend i stedet for:
    ## R_h=signal.detrend(R,axis=0,type='linear')
    ## R_h=numpy.mat(R_h)
    ## R_tilde=R_h[:,index].conj().T
    #bandpass=[0.5,5.]
    R_heart = simple_filt(R_tilde,bandpass,prf)

    return R_heart






def wienerfilt(V,V_ref):
    """
    A Wiener filter for spike suppression in the calibration procedure
    """
    #First, the Wiener matrix is created from the matrices V and V_ref:
    eps=.0001
    G=V/(V_ref+eps)
    snr=abs(numpy.mean(G,axis=0)/(numpy.std(G,axis=0)+eps))
    V_ref=numpy.mean(V_ref,axis=0)
    ref_env=abs(get_envelope(V_ref))
    eps = .04*max(ref_env)
    ## GW=numpy.multiply(V,numpy.conj(V_ref))/(numpy.power(ref_env,2)+numpy.power(1/snr,2))
    ## GW=V/(ref_env+1/(snr+eps))
    GW = V/(ref_env+eps)
    
    return GW





def calibration_wiener(R, clutter, ref, plot_figures=True):
    """
    A function for calibrating the input matrix R with the matrices clutter and ref, using wiener filter for spike suppression. 
    @param R: The input data matrix.
    @param clutter: An empty room measurement data matrix.
    @param ref: The reference target data matrix.
    """
    #Remember: The matrices must be of the same size for the wiener filter to work.
    #First, averaging of the clutter and reference measurements:
    clutter=numpy.mean(clutter,axis=0)
    #Subtracts the empty room measurements:
    R=R-clutter
    ref=ref-clutter

    ###################Experimentation:
##     BW=3E9
##     N_freq=1000
##     padlength=2**13
##     padfactor=float(padlength)/N_freq
##     f_R=scipy.fft(R,padlength)
##     import pylab
##     f_R=window(f_R,555,600) #fucker opp av uviss grunn.
##     pylab.figure()
##     pylab.plot(20*numpy.log10(abs(f_R[100,:].T+.000001)))
##     R=scipy.ifft(f_R)
##     R=R[:,0:1000]
##     print R.shape[1]
    ## print R.shape
##     BW=3e9
##     N_freq=1000
##     padlength=2**13
##     #The range indices have to be adjusted because of the zeropadding in fft:
##     padfactor=float(padlength)/N_freq
##     R=pre_cal(R,ns_to_samples(23.479-1,BW,padfactor),ns_to_samples(23.479+11,BW,padfactor),padlength)
##     ref=pre_cal(ref,ns_to_samples(23.479-1,BW,padfactor),ns_to_samples(23.479+1,BW,padfactor),padlength)
##     R=R[:,0:N_freq]
##     ref=ref[:,0:N_freq]
##     print R.shape
    ###################

    #Does the calibration:
    R_cal=wienerfilt(R,ref)

    return R_cal


def calibration_wiener_gated(R, clutter, ref, target_tau, BW, padlength, min_f_index, max_f_index=1000, plot_figures=True):
    """
    A function for calibrating the input matrix R with the matrices clutter and ref,
    using wiener filter for spike suppression. First the data is software gated to
    keep only the data in a fast time window around the target. 
    @param R: The input data matrix.
    @param clutter: An empty room measurement data matrix.
    @param ref: The reference target data matrix.
    """
    from numpy import fft
    #############################
    #Litt hardkoding, fy!
    N_freq=R.shape[1]
    padfactor=float(padlength)/N_freq
    
    #############################
    
    #Remember: The matrices must be of the same size for the wiener filter to work.
    #First, averaging of the clutter and reference measurements:
    clutter=numpy.mean(clutter,axis=0)
    #Subtracts the empty room measurements:
    R=R-clutter
    ref=ref-clutter

    import pylab
    ## pylab.figure()
##     pylab.title('Raw data, clutter subtracted.')
##     pylab.plot(numpy.mean(R,axis=0).T)



    #############
    #Bruk dette:
    target_index=ns_to_samples(target_tau, BW, padfactor)
    window_radius=ns_to_samples(10, BW, padfactor)#A radius of 7ns around target
    
    R_f=fft.rfft(R,padlength)

    #Plotting foer gating:
    l=R_f.shape[1]
    tauvec=numpy.arange(l,dtype=float)/l*samples_to_ns(l/2,BW,padfactor)

    R_f=window(R_f, int(target_index-window_radius), int(target_index+window_radius))


    ref_f=fft.rfft(ref,padlength)
    ref_f=window(ref_f, int(target_index-window_radius), int(target_index+window_radius))

    R=fft.irfft(R_f)
    ref=fft.irfft(ref_f)
   #############

   


    ################
    #Hardkoding, skam deg!

    #Cheesy:
    toleranse=0
    R[:,0:(min_f_index-toleranse)] = 0
    R[:,(max_f_index+toleranse):R.shape[1]] = 0

    ref[:,0:(min_f_index-toleranse)] = 0
    ref[:,(max_f_index+toleranse):ref.shape[1]] = 0
    ######################


    #print type(R)
    #Does the calibration:
    R_cal=wienerfilt(R,ref)
    #print type(R)
    #print R_cal

    pylab.figure()
    pylab.title('Calibrated raw data.')
    pylab.plot(numpy.mean(R_cal[:,0:N_freq],axis=0))
    return R_cal




def adjust_calibration(V, V_c, min_f_index, index):
    """
    Outputs a scalar factor such that (V-V_c)/(V-V_c) is equal to 1 at range index index after ifft.
    """
    V_r=calibration_wiener(V, V_c, V)
    padlength=2**14
    norm_ifft_factor=float(padlength)/(float(V.shape[1]-min_f_index)/2)
    V_rifft=scipy.ifft(V_r,padlength)*norm_ifft_factor

    norm_factor=1./(numpy.mean(V_rifft[:,index],axis=0))
    return norm_factor



def ns_to_samples(ns,BW,padfactor):
    """
    A simple function for converting from ns to sample number in the fast time domain.
    Not a good function, assumes N_freq=1/Ts
    """
    #NB: Given N_freq=1/sweeptime
    index=padfactor*ns*float(BW)/1E9
    print 'You are using an obsolete function, switch to ns2samples!'
    return int(index)



def samples_to_ns(index,BW,padfactor,fs=1e6,Ts=1e-3,N_freq=1e3):
    """
    A simple function for converting from ns to sample number in the fast time domain.
    NB: Actual function is tau=index*fs*Ts/(BW*N_freq*padfactor)?
    """
    ns=float(index)*1E9/(BW*padfactor)
    print 'You are using an obsolete function, switch to samples2ns!'
    return ns


def samples2ns(index,maxindex,tau_utvetydig):
    """
    Use this function.
    """
    tau = float(index)/maxindex*tau_utvetydig
    return tau

def ns2samples(tau,maxindex,tau_utvetydig):
    """
    Use this function.
    """
    index = float(tau)/tau_utvetydig*maxindex/2
    return int(index)
    


def my_cal(R,background,cal,min_tau,max_tau,tau_res,BW,N_freq,norm_fft_factor):
    """
    A function for calibrating the data in matrix R, using background matrix background
    and calibration matrix cal. The data matrix is also transformed to the fast time-slow time domain
    using a chirpz-transform.
    """
    import chirpz
    #Subtracts the background:
    average_background=numpy.mean(background,axis=0)
    average_cal=numpy.mean(cal,axis=0)-average_background
    R=R-average_background

    import pylab
    pylab.figure()
    pylab.plot(average_cal.T)
    pylab.title('cal')
    pylab.figure()
    pylab.plot(numpy.mean(R,axis=0).T)
    pylab.title('data')

    #Chirpz transform to the wanted tau-window:
    R_f=chirpz.zoom_fft(R,min_tau,max_tau,tau_res,BW,N_freq)*norm_fft_factor
    cal_f=chirpz.zoom_fft(average_cal,min_tau,max_tau,tau_res,BW,N_freq)*norm_fft_factor

    pylab.figure()
    pylab.plot(10*numpy.log10(abs(cal_f).T))
    pylab.title('cal')
    pylab.figure()
    pylab.plot(10*numpy.log10(abs(numpy.mean(R_f,axis=0)).T))
    pylab.title('data')

    #Squaring of the amplitudes to transform from voltage to power:
    #R_hilb=get_envelope(R_f)
    #cal_hilb=get_envelope(cal_f)
    amp=numpy.power(abs(R_f)/abs(cal_f), 2)
    phase=numpy.exp(1j*numpy.angle(R_f))

    print R_f.shape
    print cal_f.shape
    pylab.figure()
    pylab.plot(10*numpy.log10(numpy.mean(abs(R_f),axis=0)))
    pylab.title('env')
    pylab.figure()
    pylab.plot(10*numpy.log10(abs(cal_f).T))
    pylab.title('cal env')
    
    R_cal=numpy.multiply(amp,phase)

    pylab.figure()
    pylab.plot(10*numpy.log10(abs(numpy.mean(R_cal,axis=0)).T))
    pylab.title('kalibrerte data')
    
    return R_cal




def replace_part(R, Rmean, start, stop):
    """
    Replaces the part of R between index start and index stop with Rmean
    """
    for n in range(R.shape[0]):
        R[n,start:stop]=Rmean[0,start:stop]
    return R



def trace_max(R):
    """
    Traces the max value of the range profile matrix R along the slow time axis.
    @input R: The range profile matrix, preferable zoomed in around the person using chirpz.
    
    """
    ## maxindexes = numpy.zeros(R.shape[0],dtype='float')
    maxindexes = abs(R).argmax(axis=1)#[max(R[n,:]) for n in range(R.shape[0])]
    
    return maxindexes


def myLMS(x,d,p):
    """
    My implementation of a simple LMS adaptive filter algorithm.
    x = input vector
    d = 'noise' vector
    p = filter order
    """
    w = numpy.zeros(p,dtype=x.dtype)
    xk = numpy.zeros(p,dtype=x.dtype)
    y = numpy.ravel(x)
    x=y
    epsilon = numpy.zeros(len(x),dtype=x.dtype)
    for k in xrange(p,len(x)):
        if k==p:
            xk[:]=x[0:k].copy()
        else:
            xk[:]=x[k:k-p:-1].copy()
        mu=0.9/((p+1)*numpy.mean(numpy.power(xk,2)))
        epsilon[k] = d[k] - numpy.dot(xk.conj().T,w) #possibly switch which signal is subtracted from which here
        w=w+2*mu*epsilon[k]*xk
        #Alternative: Normalized lms:
        ## mu=1.0
        # w=w+mu*xk*epsilon[k]/float(numpy.dot(xk.conj().T,xk))

    return epsilon




def makedB(input,dynamicdb = 100):
    """
    Converts the input containing both positive and negative values to dB scale in both positive and negative values.
    @param input: A numpy array of arbitraty size and dimensions.
    """
    input=numpy.real(input)
    lim = 20*numpy.log10(numpy.max(abs(input)))-dynamicdb #Every value below this will be set to lim.
    output = numpy.zeros(input.shape,dtype='float')
    #positive values:
    tmp = 20*numpy.log10(input[input>0])
    tmp[tmp<lim] = lim
    output[input>0] = tmp-lim #Adds the positive lim so that minimum value is set to zero

    #negative values:
    tmp = 20*numpy.log10(-1*input[input<0])
    tmp[tmp<lim] = lim
    output[input<0] = -1*(tmp-lim)
    return output






def circlecenter(invector):
    """
    herp derp
    """
    XY = numpy.array([numpy.real(invector),numpy.imag(invector)])
    n = XY.shape[1]
    centroid = numpy.mean(XY,axis=1)

    Mxx = 0
    Myy = 0
    Mxy = 0
    Mxz = 0
    Myz = 0
    Mzz = 0

    for i in range(n):
        Xi = XY[0,i] - centroid[0]  #  centering data
        Yi = XY[1,i] - centroid[1]  #  centering data
        Zi = Xi*Xi + Yi*Yi
        Mxy = Mxy + Xi*Yi
        Mxx = Mxx + Xi*Xi
        Myy = Myy + Yi*Yi
        Mxz = Mxz + Xi*Zi
        Myz = Myz + Yi*Zi
        Mzz = Mzz + Zi*Zi

    Mxx = Mxx/n
    Myy = Myy/n
    Mxy = Mxy/n
    Mxz = Mxz/n
    Myz = Myz/n
    Mzz = Mzz/n

    Mz = Mxx + Myy
    Cov_xy = Mxx*Myy - Mxy*Mxy
    A3 = 4*Mz
    A2 = -3*Mz*Mz - Mzz
    A1 = Mzz*Mz + 4*Cov_xy*Mz - Mxz*Mxz - Myz*Myz - Mz*Mz*Mz
    A0 = Mxz*Mxz*Myy + Myz*Myz*Mxx - Mzz*Cov_xy - 2*Mxz*Myz*Mxy + Mz*Mz*Cov_xy
    A22 = A2 + A2
    A33 = A3 + A3 + A3
    
    xnew = 0
    ynew = 1e20
    epsilon = 1e-12
    IterMax = 20

    # Newton's method starting at x=0
    for iter in range(IterMax):
        yold = ynew
        ynew = A0 + xnew*(A1 + xnew*(A2 + xnew*A3))
        if abs(ynew) > abs(yold):
            print 'Newton-Taubin goes wrong direction: |ynew| > |yold|'
            xnew = 0
            break
    
        Dy = A1 + xnew*(A22 + xnew*A33)
        xold = xnew
        xnew = xold - ynew/Dy
        if (abs((xnew-xold)/xnew) < epsilon):
            break
        if (iter >= IterMax):
            print'Newton-Taubin will not converge'
            xnew = 0
        
        if (xnew<0.):
            print 'Newton-Taubin negative root: '+str(xnew)
            xnew = 0;

    DET = xnew*xnew - xnew*Mz + Cov_xy
    Center = numpy.array([Mxz*(Myy-xnew)-Myz*Mxy , Myz*(Mxx-xnew)-Mxz*Mxy])/DET/2
    
    Par = [Center[0]+centroid[0] +1j*(Center[1]+centroid[1]), numpy.sqrt(numpy.dot(Center,Center)+Mz)]

    return Par



def phase_estimate(invector):
    """
    Uses the circlecenter method to estimate the phase of invector
    @param invector: The complex input vector.
    """
    center,tmp = circlecenter(invector)
    centralized = invector-center

    
    return numpy.unwrap(numpy.angle(centralized))

