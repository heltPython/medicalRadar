#!/usr/bin/env python

# This program is public domain
# Authors: Paul Kienzle, Nadav Horesh
# Modifisert av Oyvind Aardal
"""
Chirp z-transform.

CZT: callable (x,axis=-1)->array
   define a chirp-z transform that can be applied to different signals
ZoomFFT: callable (x,axis=-1)->array
   define a Fourier transform on a range of frequencies
ScaledFFT: callable (x,axis=-1)->array
   define a limited frequency FFT

czt: array
   compute the chirp-z transform for a signal
zoomfft: array
   compute the Fourier transform on a range of frequencies
scaledfft: array
   compute a limited frequency FFT for a signal
"""
__all__ = ['czt', 'zoomfft', 'scaledfft']

import math, cmath

import numpy as np
from numpy import pi, arange
from scipy.fftpack import fft, ifft, fftshift

class CZT:
    """
    Chirp-Z Transform.

    Transform to compute the frequency response around a spiral.
    Objects of this class are callables which can compute the
    chirp-z transform on their inputs.  This object precalculates
    constants for the given transform.

    If w does not lie on the unit circle, then the transform will be
    around a spiral with exponentially increasing radius.  Regardless,
    angle will increase linearly.

    The chirp-z transform can be faster than an equivalent fft with
    zero padding.  Try it with your own array sizes to see.  It is
    theoretically faster for large prime fourier transforms, but not
    in practice.

    The chirp-z transform is considerably less precise than the
    equivalent zero-padded FFT, with differences on the order of 1e-11
    from the direct transform rather than the on the order of 1e-15 as
    seen with zero-padding.

    See zoomfft for a friendlier interface to partial fft calculations.
    """
    def __init__(self, n, m=None, w=1, a=1):
        """
        Chirp-Z transform definition.

        Parameters:
        ----------
        n: int
          The size of the signal
        m: int
          The number of points desired.  The default is the length of the input data.
        a: complex
          The starting point in the complex plane.  The default is 1.
        w: complex or float
          If w is complex, it is the ratio between points in each step.
          If w is float, it serves as a frequency scaling factor. for instance
          when assigning w=0.5, the result FT will span half of frequncy range
          (that fft would result) at half of the frequncy step size.

        Returns:
        --------
        CZT:
          callable object f(x,axis=-1) for computing the chirp-z transform on x
        """
        if m is None:
            m = n
        if w is None:
            w = cmath.exp(-1j*pi/m)
        elif type(w) in (float, int):
            w = cmath.exp(-1j*pi/m * w)
        else:
            w = cmath.sqrt(w)
        self.w, self.a = w, a
        self.m, self.n = m, n
        #print m
        #print n#Her sjer det noe koedd?
        k = arange(max(m,n))
        wk2 = w**(k**2)
        nfft = 2**nextpow2(n+m-1)
        self._Awk2 = (a**-k * wk2)[:n]
        self._nfft = nfft
        self._Fwk2 = fft(1/np.hstack((wk2[n-1:0:-1], wk2[:m])), nfft)
        self._wk2 = wk2[:m]
        self._yidx = slice(n-1, n+m-1)

    def __call__(self, x, axis=-1):
        """
        Parameters:
        ----------
        x: array
          The signal to transform.
        axis: int
          Array dimension to operate over.  The default is the final
          dimension.

        Returns:
        -------
          An array of the same dimensions as x, but with the length of the
          transformed axis set to m.  Note that this is a view on a much
          larger array.  To save space, you may want to call it as
          y = czt(x).copy()
        """
        x = np.asarray(x)
        if x.shape[axis] != self.n:
            raise ValueError("CZT defined for length %d, not %d" %
                             (self.n, x.shape[axis]))
        # Calculate transpose coordinates, to allow operation on any given axis
        trnsp = np.arange(x.ndim)
        trnsp[[axis, -1]] = [-1, axis]
        x = x.transpose(*trnsp)
        y = ifft(self._Fwk2 * fft(x*self._Awk2, self._nfft))
        y = y[..., self._yidx] * self._wk2
        return y.transpose(*trnsp)


def nextpow2(n):
    """
    Return the smallest power of two greater than or equal to n.
    """
    return int(math.ceil(math.log(n)/math.log(2)))



def czt(x, m=None, w=1.0, a=1, axis=-1):
    """
    Compute the frequency response around a spiral.

    Parameters:
    ----------
    x: array
      The set of data to transform.
    m: int
      The number of points desired.  The default is the length of the input data.
    a: complex
      The starting point in the complex plane.  The default is 1.
    w: complex or float
      If w is complex, it is the ratio between points in each step.
      If w is float, it is the frequency step scale (relative to the
      normal dft frquency step).
    axis: int
      Array dimension to operate over.  The default is the final
      dimension.

    Returns:
    -------
      An array of the same dimensions as x, but with the length of the
      transformed axis set to m.  Note that this is a view on a much
      larger array.  To save space, you may want to call it as
      y = ascontiguousarray(czt(x))

    See zoomfft for a friendlier interface to partial fft calculations.

    If the transform needs to be repeated, use CZT to construct a
    specialized transform function which can be reused without
    recomputing constants.
    """
    x = np.asarray(x)
    transform = CZT(x.shape[axis], m=m, w=w, a=a)
    return transform(x,axis=axis)




def zoom_fft(x,tau1,tau2,tau_res,BW,N_freq,Ts=1e-3,fs=1e6):
    """
    Performs a FFT on the input matrix x, zoomed in on the slow time window (tau1,tau2). The output is the fft in this region, with resolution specified by tau_res.
    @param x: The input data matrix.
    @param tau1: The start of the fast time zoom window in nanoseconds.
    @param tau2: The end of the fast time zoom window in nanoseconds.
    @param tau_res: The wanted fast time resolution in nanoseconds.
    @param BW: The bandwidth.
    @param N_freq: The number of samples in x.
    """
    c=3e8
    
    #Selects range:
    r1=tau1*c/2e9
    r2=tau2*c/2e9
    range=r2-r1
    delta_r_req=c*tau_res/2e9

    rs=c*N_freq/(2*BW)
    #A test in making the function more versatile:
    rs = c*Ts*fs/(2*BW)
    K=np.floor(range/delta_r_req+.5)
    #rvals_czt=np.arange(0,K,dtype=float)*(r2-r1)/K+r1
    w=np.exp(-1j*2*pi*(r2-r1)/(K*rs))
    a=np.exp(1j*2*pi*r1/rs)

    czt_obj=CZT(N_freq,K,w,a)
    zoomed_signal=czt_obj(x)
    return zoomed_signal


def zoom_ifft(x,tau1,tau2,tau_res,BW,N_freq,Ts=1e-3,fs=1e6):
    """
    Performs a FFT on the input matrix x, zoomed in on the slow time window (tau1,tau2). The output is the fft in this region, with resolution specified by tau_res.
    @param x: The input data matrix.
    @param tau1: The start of the fast time zoom window in nanoseconds.
    @param tau2: The end of the fast time zoom window in nanoseconds.
    @param tau_res: The wanted fast time resolution in nanoseconds.
    @param BW: The bandwidth.
    @param N_freq: The number of samples in x.
    """
    c=3e8
    
    #Selects range:
    r1=tau1*c/2e9
    r2=tau2*c/2e9
    range=r2-r1
    delta_r_req=c*tau_res/2e9

    rs=c*N_freq/(2*BW)
    #A test in making the function more versatile:
    rs = c*Ts*fs/(2*BW)
    K=np.floor(range/delta_r_req+.5)
    #rvals_czt=np.arange(0,K,dtype=float)*(r2-r1)/K+r1
    w=np.exp(1j*2*pi*(r2-r1)/(K*rs))
    a=np.exp(-1j*2*pi*r1/rs)

    czt_obj=CZT(N_freq,K,w,a)
    zoomed_signal=czt_obj(x)
    return zoomed_signal

