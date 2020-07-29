'''
Code for quantile matching bias corrections.
Based on `Li, H., J. Sheffield, and E. F. Wood (2010),
Bias correction of monthly precipitation and temperature fields from
Intergovernmental Panel on Climate Change AR4 models using equidistant quantile
matching, J. Geophys. Res., 115, D10101,
doi:10.1029/2009JD012882.
<http://onlinelibrary.wiley.com/doi/10.1029/2009JD012882/abstract;jsessionid=3F0FC6D4A880755443F5DF5697396006.f04t02>`_

Note on general notation:
_oc refers to observations in current day
_mc refers to model data in current day
_mp refers to model data in projections / predictions / forecasts
_amp refers to adjusted model prediction
'''
from __future__ import division
from __future__ import print_function

#from future import standard_library
#standard_library.install_aliases()

from six.moves.builtins import str
from six.moves.builtins import zip
from six.moves.builtins import range
from six.moves.builtins import object

from six.moves import cPickle as pickle #cPickle in python2, pickle in python3

import decimal
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
from scipy.integrate import cumtrapz
import scipy.interpolate
from scipy.constants import pi
import pandas as pd



def get_pdf(x, mu, sigma):
    '''
    Calculate pdf given x values, mean mu and standard deviation sigma

    >>> x = np.arange(20)
    >>> mu = 10.
    >>> sigma = 2.
    >>> pdf = get_pdf(x, mu, sigma)
    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2e}'.format(x)})
    >>> print(pdf)
    [7.43e-07 7.99e-06 6.69e-05 4.36e-04 2.22e-03 8.76e-03 2.70e-02 6.48e-02
     1.21e-01 1.76e-01 1.99e-01 1.76e-01 1.21e-01 6.48e-02 2.70e-02 8.76e-03
     2.22e-03 4.36e-04 6.69e-05 7.99e-06]
    >>> np.set_printoptions()
    '''
    z = (x - mu)/(sigma)
    pdf = np.exp((-z**2)/2.) / (sigma*np.sqrt(2.*pi))
    return pdf

def get_cdf(x, mu, sigma):
    '''
    Calculate cdf (cumulative distribution function),
    given x values, mean mu and standard deviation sigma

    >>> x = np.arange(20)
    >>> mu = 10.
    >>> sigma = 2.
    >>> cdf = get_cdf(x, mu, sigma)
    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2e}'.format(x)})
    >>> print(cdf)
    [2.87e-07 3.40e-06 3.17e-05 2.33e-04 1.35e-03 6.21e-03 2.28e-02 6.68e-02
     1.59e-01 3.09e-01 5.00e-01 6.91e-01 8.41e-01 9.33e-01 9.77e-01 9.94e-01
     9.99e-01 1.00e+00 1.00e+00 1.00e+00]
    >>> np.set_printoptions()
    '''
    z = (x - mu) / (sigma)
    cdf = 0.5+0.5*erf(z/np.sqrt(2.))
    return cdf

def pdf2cdf(x, pdf):
    '''
    Calculate cdf (cumulative distribution function,
    given x values and a pdf of y values, where pdf and cdf
    are considered to be continous, not discrete.

    >>> x = np.arange(20)
    >>> mu = 10.
    >>> sigma = 2.
    >>> pdf = get_pdf(x, mu, sigma)
    >>> cdf = pdf2cdf(x, pdf)
    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2e}'.format(x)})
    >>> print(cdf)
    [0.00e+00 4.37e-06 4.18e-05 2.93e-04 1.62e-03 7.11e-03 2.50e-02 7.09e-02
     1.64e-01 3.12e-01 5.00e-01 6.88e-01 8.36e-01 9.29e-01 9.75e-01 9.93e-01
     9.98e-01 1.00e+00 1.00e+00 1.00e+00]
    >>> np.set_printoptions()
    '''
    cdf = cumtrapz(pdf, x, initial=0.)
    return cdf

def cdf2pdf(x, cdf):
    '''
    Calculate pdf given x values and a cdf of y values

    >>> x = np.arange(20)
    >>> mu = 10.
    >>> sigma = 2.
    >>> cdf = get_cdf(x, mu, sigma)
    >>> pdf = cdf2pdf(x, cdf)
    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2e}'.format(x)})
    >>> print(pdf)
    [3.11e-06 1.57e-05 1.15e-04 6.59e-04 2.99e-03 1.07e-02 3.03e-02 6.80e-02
     1.21e-01 1.71e-01 1.91e-01 1.71e-01 1.21e-01 6.80e-02 3.03e-02 1.07e-02
     2.99e-03 6.59e-04 1.15e-04 2.83e-05]
    >>> np.set_printoptions()
    '''
    #Check distance between x values is consistent
    unique = np.unique(np.diff(x).round(decimals=10))
    if len(unique) > 1:
        raise ValueError('Difference between "x" values must be consistent')
    delta_x = x[1] - x[0]
    pdf = np.gradient(cdf, delta_x)
    return pdf

def hist2cdf(x, hist):
    '''
    Calculate a pdf given x values and a histogram of y values.
    Returns a new x and cdf. The new x has an extra value of x inserted
    at the beginning to allow the cdf to start from zero (if first element
    of histogram array is not zero already).

    >>> x = [1, 2, 3, 4, 5]
    >>> hist = [2, 0, 0, 1, 1]
    >>> newx, cdf = hist2cdf(x,hist)
    >>> print(newx)
    [0 1 2 3 4 5]
    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
    >>> print(cdf)
    [ 0.00  0.50  0.50  0.50  0.75  1.00]
    >>> np.set_printoptions()
    '''
    #Need to ensure histogram goes to zero at either end
    # to ensure sensible pdf and cdf produced
    if hist[0] != 0:
        x = np.insert(x, 0, x[0]-(x[1]-x[0]))
        hist = np.insert(hist, 0, 0)

    hist = np.asarray(hist, dtype=float) #Ensure histogram contains float data
    pdf = hist / sum(hist)
    #Use np.cumsum for discrete bins, not pdf2cdf which is for continuous data.
    cdf = np.cumsum(pdf)
    #Ensure no cdf values > 1. to account for rounding errors
    mask = cdf > 1.
    cdf[mask] = 1.
    #Ensure no cdf values < 0. to account for rounding errors
    mask = cdf < 0.
    cdf[mask] = 0.

    return x, cdf

def check_monotonic_increasing(xarray, strict=False):
    '''
    Checks that values in xarray are montonicly increasing.

    In particular, it will check for strictly increasing away from the
    edges, and equal or increasing near the edges

    If keyword strict=True, then expects strictly increasing everywhere.

    >>> x = [0, 1, 2, 3, 4]
    >>> print(check_monotonic_increasing(x))
    True
    >>> x = [0, 0, 1, 2, 3, 4, 4]
    >>> print(check_monotonic_increasing(x))
    True
    >>> print(check_monotonic_increasing(x, strict=True))
    False
    >>> x = [0, 1, 1, 2, 3, 4]
    >>> print(check_monotonic_increasing(x))
    False
    '''

    found_start = False

    for i in np.arange(len(xarray)-1):
        if strict:
            if xarray[i] >= xarray[i+1]:
                return False
        else:
            if not found_start:
                #Find where xarray starts to increase
                if xarray[i] > xarray[i+1]:
                    #Values start to decrease
                    return False
                if xarray[i] < xarray[i+1]:
                    found_start = True
            else:
                if xarray[i] > xarray[i+1]:
                    #Values start to decrease
                    return False
                if xarray[i] == xarray[i+1]:
                    if xarray[i] != xarray[-1]:
                        #Values level off, but not at final value
                        return False

    return True

def make_strictly_increasing(yarray, xarray=None):
    '''
    Make y strictly increasing.

    This is done by interpolating between surrounding values
    to remove any constant or decreasing values.

    >>> Y = [0,0,0,3,4,5,4,5,6,7,7,7,7]
    >>> print(make_strictly_increasing(Y))
    [0.0, 1.0, 2.0, 3, 4, 4.5, 4.5, 5, 6, 6.25, 6.5, 6.75, 7]
    '''

    if isinstance(yarray, np.ndarray):
        youtput = yarray.copy()
    else:
        youtput = yarray[:]
    if xarray is None:
        xarray = np.arange(len(yarray))

    if len(xarray) > 1 and len(yarray) > 1:
        for i in np.arange(len(yarray)-1):
            fixvalue = False
            if yarray[i] >= yarray[i+1]:
                fixvalue = True
            if i > 0 and yarray[i] <= youtput[i-1]:
                fixvalue = True
            if fixvalue:
                #Need to correct value of yarray[i]

                #Find next different larger value
                for xnext, ynext in zip(xarray[i+1:], yarray[i+1:]):
                    if ynext > yarray[i]:
                        break
                #Find previous different smaller value
                if i > 0:
                    for xprev, yprev in zip(xarray[i-1::-1], youtput[i-1::-1]):
                        if yprev < yarray[i]:
                            break
                else:
                    #Use current value as no previous values available
                    yprev = yarray[i]
                    xprev = xarray[i]

                #Now correct value by linear interpolation
                #between previous and next values
                ftn = scipy.interpolate.interp1d([xprev, xnext],
                                                 [yprev, ynext],
                                                 kind='linear')
                youtput[i] = float(ftn(xarray[i]))

    return youtput

def make_strictly_monotonic_cdf(x, cdf):
    ''' Ensure CDF is monotonic, so there are no duplicate x values,
    and only a single value of y is kept at either end of the CDF.
    .. Note:: Doesn't check rest of cdf is monotonic.

    >>> x = np.arange(20)
    >>> mu = 10.
    >>> sigma = 2.
    >>> cdf = get_cdf(x, mu, sigma)
    >>> new_x, new_cdf = make_strictly_monotonic_cdf(x,cdf)

    Does not change answers if already OK.

    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
    >>> print(cdf - new_cdf)
    [ 0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
      0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00]

    If duplicated x or y values:

    >>> x = [ 2., 2., 3., 4., 5., 6., 7., 8., 9.]
    >>> cdf = [0., 0.1, 0.3, 0.5, 0.7,0.99, 1., 1., 1.]
    >>> new_x, new_cdf = make_strictly_monotonic_cdf(x, cdf)
    >>> print(new_x)
    [2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
    >>> print(new_cdf)
    [0.0, 0.3, 0.5, 0.7, 0.99, 1.0]
    >>> np.set_printoptions()
    '''

    imax = len(x)-1
    imin = 0
    #Ensure that stop cdf as soon as reach a y value of 1,
    #or x values are the same.
    for i in range(len(x)-1, 0, -1):
        if cdf[i-1] < 1 and (x[i] != x[i-1]):
            imax = i
            break
    #Only keep first value of y at zero, and no duplicate x values at start
    for i in range(len(x)-1):
        if cdf[i+1] > 0 and (x[i] != x[i+1]):
            imin = i
            break

    x_new = x[imin:imax+1]
    cdf_new = cdf[imin:imax+1]

    #Check cdf is monotonicly increasing in middle and fix any problems:
    cdf_new = make_strictly_increasing(cdf_new, x_new)
    cdf_new[0] = min(cdf) #Ensure min cdf value is kept
    cdf_new[-1] = max(cdf) #Ensure max cdf value is kept

    #Check x is monotonic and fix any problems:
    x_new = make_strictly_increasing(x_new, cdf_new)

    return x_new, cdf_new

def cdfm(x, cdf_oc, cdf_mc, cdf_mp, minvalidx=None):
    '''
    Calculate CDFm, defined by:
    :math:`x_{amp} = inversef_{oc}(f_{mc}(x_{mp}))`

    Returns a spline function which can be used to calculate x_amp
    afterwards.

    If minvalidx is set, then the minimum value of x_amp allowed
    is set to this value. This may be useful to stop negative concentrations.

    For example, set up idealistic input cdfs:

    >>> x = np.linspace(0,25,10)
    >>> cdf_oc = get_cdf(x, 10.0, 2.0) #oc Obs Current
    >>> cdf_mc = get_cdf(x, 12.0, 2.5) #mc Model Current
    >>> cdf_mp = get_cdf(x, 15.0, 3.0) #mp Model Projection

    Need to ensure that these idealistic cdfs all start with zero, and
    end with 1, so add these values in and extend x:

    >>> x = np.insert(x, 0, x[0]-(x[1]-x[0]))
    >>> cdf_oc = np.insert(cdf_oc,0,0.)
    >>> cdf_mc = np.insert(cdf_mc,0,0.)
    >>> cdf_mp = np.insert(cdf_mp,0,0.)
    >>> x = np.append(x, x[-1]+(x[1]-x[0]))
    >>> cdf_oc = np.append(cdf_oc,1.)
    >>> cdf_mc = np.append(cdf_mc,1.)
    >>> cdf_mp = np.append(cdf_mp,1.)

    >>> cdfm_spline = cdfm(x, cdf_oc, cdf_mc, cdf_mp)
    >>> print(cdfm_spline) # doctest: +ELLIPSIS
    <scipy.interpolate.fitpack2.InterpolatedUnivariateSpline object at ...>

    Can then calculate adjusted values:

    >>> x_adjusted = cdfm_spline(x)
    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2e}'.format(x)})
    >>> print(x_adjusted)
    [-2.78e+00 9.25e-03 2.05e+00 3.81e+00 6.41e+00 9.20e+00 1.18e+01 1.38e+01
     1.66e+01 1.93e+01 2.20e+01 2.78e+01]
    >>> np.set_printoptions()

    Could also just convert a single value:

    >>> print("{:,.6f}".format(float(cdfm_spline(10))))
    8.083867
    '''

    #Need to ensure all cdfs being used are monotonicly increasing
    if not check_monotonic_increasing(cdf_oc):
        raise ValueError("CDF OC is not monotonicly increasing")
    if not check_monotonic_increasing(cdf_mc):
        raise ValueError("CDF MC is not monotonicly increasing")
    if not check_monotonic_increasing(cdf_mp):
        raise ValueError("CDF MP is not monotonicly increasing")

    x_amp = np.zeros(len(x)) #amp = Adjusted Model Projection

    #Need to also ensure that cdfs don't have duplicate values at
    #the beginning and end as this will make inverting them impossible.
    x_oc, cdf_oc = make_strictly_monotonic_cdf(x, cdf_oc)
    inversef_oc = scipy.interpolate.interp1d(cdf_oc, x_oc, kind='linear')

    for i in np.arange(len(x)):
        f_mc = cdf_mc[i]
        x_amp[i] = inversef_oc(f_mc)
        if minvalidx is not None:
            if x_amp[i] < minvalidx:
                x_amp[i] = minvalidx

    spline = scipy.interpolate.UnivariateSpline(x, x_amp, s=0, k=1)

    return spline

def edcdfm(x, cdf_oc, cdf_mc, cdf_mp, ratio=False, minvalidx=None):
    '''
    Calculate EDCDFm, defined by:
    :math:`x_{amp} = x_{mp} + inversef_{oc}(f_{mp}(x_{mp})) \
          - inversef_{mc}(f_{mp}(x_{mp}))`

    Returns a spline function which can be used to calculate x_amp
    afterwards.

    If ratio is set to True, then instead of using an additive difference
    as above, instead uses ratios:
    x_amp = x_mp * inversef_oc(f_mp(x_mp)) / inversef_mc(f_mp(x_mp))

    If minvalidx is set, then the minimum value of x_amp allowed
    is set to this value. This may be useful to stop negative concentrations.

    For example, set up idealistic input cdfs:

    >>> x = np.linspace(0,25,10)
    >>> cdf_oc = get_cdf(x, 10.0, 2.0) #oc Obs Current
    >>> cdf_mc = get_cdf(x, 12.0, 2.5) #mc Model Current
    >>> cdf_mp = get_cdf(x, 15.0, 3.0) #mp Model Projection

    Need to ensure that these idealistic cdfs all start with zero, and
    end with 1, so add these values in and extend x:

    >>> x = np.insert(x, 0, x[0]-(x[1]-x[0]))
    >>> cdf_oc = np.insert(cdf_oc,0,0.)
    >>> cdf_mc = np.insert(cdf_mc,0,0.)
    >>> cdf_mp = np.insert(cdf_mp,0,0.)
    >>> x = np.append(x, x[-1]+(x[1]-x[0]))
    >>> cdf_oc = np.append(cdf_oc,1.)
    >>> cdf_mc = np.append(cdf_mc,1.)
    >>> cdf_mp = np.append(cdf_mp,1.)

    >>> edcdfm_spline = edcdfm(x, cdf_oc, cdf_mc, cdf_mp)
    >>> print(edcdfm_spline) # doctest: +ELLIPSIS
    <scipy.interpolate.fitpack2.InterpolatedUnivariateSpline object at ...>

    Can then calculate adjusted values:

    >>> x_adjusted = edcdfm_spline(x)
    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
    >>> print(x_adjusted)
    [-2.78  1.77  2.64  5.29  7.99  9.32 12.00 14.32 16.76 19.24 20.24 27.78]
    >>> np.set_printoptions()

    Could also just convert a single value:

    >>> print("{:,.6f}".format(float(edcdfm_spline(10))))
    8.788440
    '''

    #Need to ensure all cdfs being used are monotonicly increasing
    if not check_monotonic_increasing(cdf_oc):
        raise ValueError("CDF OC is not monotonicly increasing")
    if not check_monotonic_increasing(cdf_mc):
        raise ValueError("CDF MC is not monotonicly increasing")
    if not check_monotonic_increasing(cdf_mp):
        raise ValueError("CDF MP is not monotonicly increasing")

    x_amp = np.zeros(len(x)) #amp = Adjusted Model Projection

    #Need to also ensure that cdfs don't have duplicate values at
    #the beginning and end as this will make inverting them impossible.
    x_oc, cdf_oc = make_strictly_monotonic_cdf(x, cdf_oc)
    inversef_oc = scipy.interpolate.interp1d(cdf_oc, x_oc, kind='linear')

    x_mc, cdf_mc = make_strictly_monotonic_cdf(x, cdf_mc)
    inversef_mc = scipy.interpolate.interp1d(cdf_mc, x_mc, kind='linear')

    for i in np.arange(len(x)):
        x_mp = x[i]
        f_mp = cdf_mp[i]
        x_oc = inversef_oc(f_mp)
        x_mc = inversef_mc(f_mp)
        if ratio:
            x_amp[i] = x_mp * x_oc / x_mc
        else:
            x_amp[i] = x_mp + x_oc - x_mc
        if minvalidx is not None:
            if x_amp[i] < minvalidx:
                x_amp[i] = minvalidx

    spline = scipy.interpolate.UnivariateSpline(x, x_amp, s=0, k=1)

    return spline


def cdfm_getcdf(x, cdf_oc, cdf_mc, cdf_mp):
    '''
    Calculate the cumulative distribution function of CDFm, defined by:
    :math:`x_{amp} = inversef_{oc}(f_{mc}(x_{mp}))`

    Returns x values and cumulative distribution at these values for the CDFm

    For example, set up idealistic input cdfs:

    >>> x = np.linspace(0,25,10)
    >>> cdf_oc = get_cdf(x, 10.0, 2.0) #oc Obs Current
    >>> cdf_mc = get_cdf(x, 12.0, 2.5) #mc Model Current
    >>> cdf_mp = get_cdf(x, 15.0, 3.0) #mp Model Projection

    Need to ensure that these idealistic cdfs all start with zero, and
    end with 1, so add these values in and extend x:

    >>> x = np.insert(x, 0, x[0]-(x[1]-x[0]))
    >>> cdf_oc = np.insert(cdf_oc,0,0.)
    >>> cdf_mc = np.insert(cdf_mc,0,0.)
    >>> cdf_mp = np.insert(cdf_mp,0,0.)
    >>> x = np.append(x, x[-1]+(x[1]-x[0]))
    >>> cdf_oc = np.append(cdf_oc,1.)
    >>> cdf_mc = np.append(cdf_mc,1.)
    >>> cdf_mp = np.append(cdf_mp,1.)

    Can now calculate the the adjusted model projection cdf.
    Note this is on a new x axis:

    >>> CDF_amp_x, cdf_CDF_amp_y = cdfm_getcdf(x, cdf_oc, cdf_mc, cdf_mp)

    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2e}'.format(x)})
    >>> print(CDF_amp_x)
    [-2.78e+00 4.44e-16 2.78e+00 5.56e+00 8.33e+00 1.11e+01 1.39e+01 1.67e+01
     1.94e+01 2.22e+01 2.50e+01 2.78e+01]

    >>> print(cdf_CDF_amp_y)
    [0.00e+00 2.86e-07 3.53e-04 9.09e-03 7.12e-02 2.88e-01 7.15e-01 9.33e-01
     9.92e-01 1.00e+00 1.00e+00 1.00e+00]
    >>> np.set_printoptions()
    '''

    spline = cdfm(x, cdf_oc, cdf_mc, cdf_mp)
    cdf_amp_xvalues = spline(x)
    cdf_amp_yvalues = cdf_mp

    #Put onto nice linearly spaced gridded values
    spline_cdf_amp = scipy.interpolate.UnivariateSpline(cdf_amp_xvalues,
                                                        cdf_amp_yvalues,
                                                        s=0, k=1)
    cdf_amp_x = np.linspace(cdf_amp_xvalues.min(),
                            cdf_amp_xvalues.max(),
                            len(x))

    cdf_amp_y = spline_cdf_amp(cdf_amp_x)

    return cdf_amp_x, cdf_amp_y

def edcdfm_getcdf(x, cdf_oc, cdf_mc, cdf_mp):
    '''
    Calculate the cumulative distribution function for EDCDFm, defined by:
    :math:`x_{amp} = x_{mp} + inversef_{oc}(f_{mp}(x_{mp})) \
          - inversef_{mc}(f_{mp}(x_{mp}))`

    Returns x values and cumulative distribution at these values for the CDFm

    For example, set up idealistic input cdfs:

    >>> x = np.linspace(0,25,10)
    >>> cdf_oc = get_cdf(x, 10.0, 2.0) #oc Obs Current
    >>> cdf_mc = get_cdf(x, 12.0, 2.5) #mc Model Current
    >>> cdf_mp = get_cdf(x, 15.0, 3.0) #mp Model Projection

    Need to ensure that these idealistic cdfs all start with zero, and
    end with 1, so add these values in and extend x:

    >>> x = np.insert(x, 0, x[0]-(x[1]-x[0]))
    >>> cdf_oc = np.insert(cdf_oc,0,0.)
    >>> cdf_mc = np.insert(cdf_mc,0,0.)
    >>> cdf_mp = np.insert(cdf_mp,0,0.)
    >>> x = np.append(x, x[-1]+(x[1]-x[0]))
    >>> cdf_oc = np.append(cdf_oc,1.)
    >>> cdf_mc = np.append(cdf_mc,1.)
    >>> cdf_mp = np.append(cdf_mp,1.)

    Can now calculate the the adjusted model projection cdf.
    Note this is on a new x axis:

    >>> EDCDF_amp_x, cdf_EDCDF_amp_y = edcdfm_getcdf(x, cdf_oc, cdf_mc, cdf_mp)

    >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2e}'.format(x)})
    >>> print(EDCDF_amp_x)
    [-2.78e+00 4.44e-16 2.78e+00 5.56e+00 8.33e+00 1.11e+01 1.39e+01 1.67e+01
     1.94e+01 2.22e+01 2.50e+01 2.78e+01]

    >>> print(cdf_EDCDF_amp_y)
    [0.00e+00 1.75e-07 6.45e-05 2.02e-03 3.48e-02 2.70e-01 6.45e-01 9.23e-01
     9.94e-01 1.00e+00 1.00e+00 1.00e+00]
    >>> np.set_printoptions()
    '''

    spline = edcdfm(x, cdf_oc, cdf_mc, cdf_mp)
    cdf_amp_xvalues = spline(x)
    cdf_amp_yvalues = cdf_mp

    #Now put the the adjusted values back onto regular grid
    spline_cdf_amp = scipy.interpolate.UnivariateSpline(cdf_amp_xvalues,
                                                        cdf_amp_yvalues,
                                                        s=0, k=1)
    cdf_amp_x = np.linspace(x.min(), cdf_amp_xvalues.max(), len(x))
    cdf_amp_y = spline_cdf_amp(cdf_amp_x)

    return cdf_amp_x, cdf_amp_y

def create_subplot(x, oc, mc, mp, cdf_amp_x, cdf_amp_y,
                   edcdf_amp_x, edcdf_amp_y, ax, title):
    '''
    Plot subplot with appropriate colours and labels
    For use as a quick example only.
    '''

    ax.plot(x, field_layers=oc)
    ax.plot(x, field_layers=mc)
    ax.plot(x, field_layers=mp)
    ax.plot(cdf_amp_x, field_layers=cdf_amp_y)
    ax.plot(edcdf_amp_x, field_layers=edcdf_amp_y)
    ax.set_ylim(bottom=0.) #Min y should be 0 for both pdf and cdfs.

    # Put a legend to the right of the current axis
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size':10})

    ax.set_title(title)


class BinnedData(object):
    '''
    Class for creating and extracting binned data
    and allowing conversion to CDFs.
    '''

    def __init__(self):
        '''
        Initialise class.

        self.binsize is the size of each bin

        self.ref_binloc is a reference bin location. This will be such that
        this point is always an edge of a bin if this bin is included,
        and all other bins will be built based on this location.

        Bins are always labelled according to their maxmium included value.
        So bin labelled 1, with a binsize of 1 incudes 0<values<=1.

        self.databin_dict is a dictionary where the keywords are bin edge
        values, while the values are the total count corresponding
        to that bin edge value.
        '''
        self.binsize = 1.
        self.binsize_fmtstr = None
        self.ref_binloc = 0.
        self.databin_dict = {}

    def calc_bin(self, value):
        '''
        Calculate which bin a given value should be in

        Bins are always labelled according to their maxmium included value.
        So bin labelled 1, with a binsize of 1 incudes 0<values<=1.

        >>> bd = BinnedData()
        >>> print('{:.1f}'.format(bd.calc_bin(2.6)))
        3.0
        >>> print('{:.1f}'.format(bd.calc_bin(2)))
        2.0
        >>> print('{:.1f}'.format(bd.calc_bin(0.0)))
        0.0
        >>> print('{:.1f}'.format(bd.calc_bin(2.0001)))
        3.0
        >>> print('{:.1f}'.format(bd.calc_bin(1.3)))
        2.0
        '''
        #Note use np.fmod instead of python % for modulo to avoid
        #floating point rounding errors
        if (value-self.ref_binloc)%self.binsize == 0.:
            binvalue = value - \
                       np.fmod((value-self.ref_binloc), self.binsize)
        else:
            binvalue = self.binsize + value - \
                       np.fmod((value-self.ref_binloc), self.binsize)

        return binvalue

    def add_datapt(self, value):
        '''
        Add a single value to the data bin

        >>> bd = BinnedData()
        >>> bd.add_datapt(3.2)
        >>> print(bd.databin_dict)
        {'4.0': 1}
        >>> bd.add_datapt(1.3)
        >>> bd.databin_dict == {'2.0': 1, '4.0': 1}
        True
        '''

        #Firstly figure out which bin it should go in.
        binvalue = self.calc_bin(value)

        #Avoid precision errors in the dictionary keywords by using a consistent
        #format string based on the number of decimal points in the required
        #binsize.
        if self.binsize_fmtstr is None:
            binsize_ndp = abs(decimal.Decimal(str(self.binsize)).as_tuple().exponent)
            if binsize_ndp < 6:
                self.binsize_fmtstr = '{:.' + str(binsize_ndp)+'f}'
            else:
                self.binsize_fmtstr = '{:.3e}'

        #Then add to dictionary of data
        if binvalue not in self.databin_dict:
            self.databin_dict[self.binsize_fmtstr.format(binvalue)] = 1
        else:
            self.databin_dict[self.binsize_fmtstr.format(binvalue)] += 1


    def add_cube(self, cube):
        '''
        Add entire iris cube of data to the data bin

        >>> import iris
        >>> cube = iris.cube.Cube(np.arange(10)/2.)
        >>> cube.data[3:6] = 0.
        >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
        >>> print(cube.data)
        [ 0.00  0.50  1.00  0.00  0.00  0.00  3.00  3.50  4.00  4.50]
        >>> np.set_printoptions()
        >>> bd = BinnedData()
        >>> bd.add_cube(cube)
        >>> bd.databin_dict == {'0.0': 4, '1.0': 2, '3.0': 1,
        ... '4.0': 2, '5.0': 1}
        True
        '''


        #for value in cube.data:
        #    self.add_datapt(value)


        mindata = np.nanmin(cube.data)
        maxdata = np.nanmax(cube.data)
        #Figure out which bin these should go in.
        minbin = self.calc_bin(mindata)
        maxbin = self.calc_bin(maxdata)

        #Now generate all possible bins in this range
        bin_edges = np.arange(minbin-self.binsize, maxbin+self.binsize,
                              self.binsize)

        #Now use pandas cut and value_counts to generate binned data histogram
        #where right handsize value of bin is included and left handside is not.
        out = pd.cut(cube.data.flatten(), bin_edges, labels=bin_edges[1:])
        counts = pd.value_counts(out)
        for index, value in zip(counts.index, counts.values):
            binvalue = str(index)
            #Only store data in bin if actually any counts
            #(values) for this binvalue
            if value > 0:
                if binvalue not in self.databin_dict:
                    self.databin_dict[binvalue] = value
                else:
                    self.databin_dict[binvalue] += value

    def get_histarray(self):
        '''
        Converts databin dictionary into an array of bin edges
        and histogram totals

        >>> bd = BinnedData()
        >>> bd.databin_dict = {'1.0':2, '3.0':1, '4.0':1}
        >>> binedges, hist = bd.get_histarray()
        >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
        >>> print(binedges)
        [ 1.00  2.00  3.00  4.00]
        >>> print(hist)
        [ 2.00  0.00  1.00  1.00]
        >>> np.set_printoptions()
        '''

        #Generate full array of all possible bin edges in range
        minval = min([float(k) for k in self.databin_dict.keys()])
        maxval = max([float(k) for k in self.databin_dict.keys()])
        binedges = np.arange(minval,
                             maxval+self.binsize,
                             self.binsize)

        hist = np.zeros((len(binedges)))
        for i, binedge in enumerate(binedges):
            binedge = str(binedge)
            if binedge in self.databin_dict:
                hist[i] = self.databin_dict[binedge]

        return binedges, hist


    def get_cdf(self):
        '''
        Calculate cdf (cumulative distribution function) from data

        >>> bd = BinnedData()
        >>> bd.databin_dict = {'1.0':2, '3.0':1, '4.0':1}
        >>> binedges, cdf = bd.get_cdf()
        >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
        >>> print(binedges)
        [ 0.00  1.00  2.00  3.00  4.00]
        >>> print(cdf)
        [ 0.00  0.25  0.50  0.75  1.00]
        >>> np.set_printoptions()
        '''
        binedges, hist = self.get_histarray()
        binedges, cdf = hist2cdf(binedges, hist)
        cdf = make_strictly_increasing(cdf, binedges)

        return binedges, cdf


    def get_comparablecdfs(self, bdlist):
        '''
        Return cdfs from self, and list of other binned data objects
        with data on same x axis, ie bincentres same for all.

        >>> bd1 = BinnedData()
        >>> bd1.databin_dict = {'0.0': 2, '3.0':1, '4.0':1}
        >>> bd2 = BinnedData()
        >>> bd2.databin_dict = {'1.0': 1, '2.0':1, '5.0':1}
        >>> bd3 = BinnedData()
        >>> bd3.databin_dict = {'1.0': 2, '2.0':1, '3.0':2}

        >>> x, [cdf1, cdf2, cdf3] = bd1.get_comparablecdfs([bd2,bd3])
        >>> np.set_printoptions(formatter={'float':lambda x: '{:5.2f}'.format(x)})
        >>> print(x)
        [-1.00  0.00  1.00  2.00  3.00  4.00  5.00]
        >>> print(cdf1)
        [ 0.00  0.19  0.38  0.50  0.75  1.00  1.00]
        >>> print(cdf2)
        [ 0.00  0.00  0.33  0.50  0.67  0.83  1.00]
        >>> print(cdf3)
        [ 0.00  0.00  0.40  0.60  1.00  1.00  1.00]
        >>> np.set_printoptions()
        '''

        #Firstly check that they can be compared, ie
        #that the binsize is the same and the reference bin location
        #for self could be contained in the bin edges for the bdlist.
        for bd in bdlist:
            if (bd.ref_binloc - self.ref_binloc)%self.binsize != 0:
                raise ValueError("BinnedData objects cannot be compared")

        #Calculate bincentres and cdfs for all BinnedData objects
        bdlist_all = bdlist[:]
        bdlist_all.insert(0, self)
        bincentres_all = []
        cdf_all = []
        for bd in bdlist_all:
            bincentres, cdf = bd.get_cdf()
            bincentres_all.append(bincentres)
            cdf_all.append(cdf)

        #Produce a bincentres_final array which contains all required
        #bincentres to cover all data
        bincentres_max = max([max(v) for v in bincentres_all])
        bincentres_min = min([min(v) for v in bincentres_all])
        bincentres_final = np.arange(bincentres_min,
                                     bincentres_max+self.binsize,
                                     self.binsize)

        #Move the cdfs for each binned data object along to match the
        #new bincentres. Fill missing data at the beginning of the cdf
        #with zeros, and data at the end with 1s to ensure a true cdf still.
        cdf_final_all = []
        for bincentres, cdf in zip(bincentres_all, cdf_all):
            cdf_final = np.zeros((len(bincentres_final)))
            #Find where bincentre fits into bincentres_final:
            #Need to convert to integers by dividing by binsize
            #to allow sensible check for equality
            bincentres_final_int = (np.rint(bincentres_final/
                                            self.binsize)).astype(int)
            bincentres_min_int = (np.rint(bincentres[0]/
                                          self.binsize)).astype(int)
            indices = np.where(bincentres_final_int == bincentres_min_int)
            imin = indices[0][0]
            imax = imin+len(cdf)-1
            if imax >= len(bincentres_final):
                imax = len(bincentres_final)-1

            #Now fit cdf into cdf_final:
            cdf_final[imin:imax+1] = cdf
            #Zeros already set up at start of cdf_final, but 1s need to
            #be added to the end to make a sensible cdf
            cdf_final[imax:] = 1.

            cdf_final_all.append(cdf_final)

        return bincentres_final, cdf_final_all

    def plotbars(self, binedges=None, yvalues=None, edgecolor=None, label=None):
        '''Produce bar plot of histogram'''

        if binedges is None and yvalues is None:
            binedges, yvalues = self.get_histarray()

        plt.bar(binedges-self.binsize, yvalues, width=self.binsize, color='w',
                edgecolor=edgecolor, label=label)

    def save(self, filename='BinnedData.pkl'):
        '''Save BinnedData object to a pickle file'''
        data = [self.binsize, self.ref_binloc, self.databin_dict]
        with open(filename, "wb") as fp:
            pickle.dump(data, fp)

    def load(self, filename='BinnedData.pkl'):
        '''Restore a previously saved BinnedData object from pickle file'''
        with open(filename, "rb") as fp:
            data = pickle.load(fp)
        self.binsize, self.ref_binloc, self.databin_dict = data



def idealised_example():
    '''
    Basic example of applying quantile matching methods to
    idealistic data and plotting resulting PDF and CDFs.
    '''

    npts = 100
    xmin = 0.
    xmax = 23.
    x = np.linspace(xmin, xmax, npts)

    #Obs PDF: o-c
    mu_oc = 10.0
    sig_oc = 2.0

    #Current model PDF: m-c
    mu_mc = 12. #10.0
    sig_mc = 2.5

    #Projection model PDF: m-p
    mu_mp = 15.0
    sig_mp = 3.0

    #Set up PDFs
    pdf_oc = get_pdf(x, mu_oc, sig_oc) #oc Obs Current
    pdf_mc = get_pdf(x, mu_mc, sig_mc) #mc Model Current
    pdf_mp = get_pdf(x, mu_mp, sig_mp) #mp Model Projection

    #Set up CDFs
    cdf_oc = get_cdf(x, mu_oc, sig_oc) #oc Obs Current
    cdf_mc = get_cdf(x, mu_mc, sig_mc) #mc Model Current
    cdf_mp = get_cdf(x, mu_mp, sig_mp) #mp Model Projection

    #Ensure CDF and PDFs all start from zero
    x = np.insert(x, 0, x[0]-(x[1]-x[0]))
    cdf_oc = np.insert(cdf_oc, 0, 0.)
    cdf_mc = np.insert(cdf_mc, 0, 0.)
    cdf_mp = np.insert(cdf_mp, 0, 0.)
    pdf_oc = np.insert(pdf_oc, 0, 0.)
    pdf_mc = np.insert(pdf_mc, 0, 0.)
    pdf_mp = np.insert(pdf_mp, 0, 0.)

    #CDFs must finish with a 1. and PDFs must finish with zero.
    x = np.append(x, x[-1]+(x[1]-x[0]))
    cdf_oc = np.append(cdf_oc, 1.)
    cdf_mc = np.append(cdf_mc, 1.)
    cdf_mp = np.append(cdf_mp, 1.)
    pdf_oc = np.append(pdf_oc, 0.)
    pdf_mc = np.append(pdf_mc, 0.)
    pdf_mp = np.append(pdf_mp, 0.)

    #Calculate CDFm
    cdfm_amp_x, cdf_cdfm_amp_y = cdfm_getcdf(x, cdf_oc, cdf_mc, cdf_mp)
    #Convert back to a pdf
    pdf_cdfm_amp_y = cdf2pdf(cdfm_amp_x, cdf_cdfm_amp_y)

    #Calculate EDCDFm
    edcdfm_amp_x, cdf_edcdfm_amp_y = edcdfm_getcdf(x, cdf_oc, cdf_mc, cdf_mp)
    #Convert back to a pdf
    pdf_edcdfm_amp_y = cdf2pdf(edcdfm_amp_x, cdf_edcdfm_amp_y)

    #Plot
    plt.figure(figsize=(15, 5))

    ax1 = plt.gcf().add_subplot(1, 2, 1)
    create_subplot(x, pdf_oc, pdf_mc, pdf_mp, cdfm_amp_x, pdf_cdfm_amp_y,
                   edcdfm_amp_x, pdf_edcdfm_amp_y, ax1, 'pdfs')

    ax2 = plt.gcf().add_subplot(1, 2, 2)
    create_subplot(x, cdf_oc, cdf_mc, cdf_mp, cdfm_amp_x, cdf_cdfm_amp_y,
                   edcdfm_amp_x, cdf_edcdfm_amp_y, ax2, 'cdfs')

    plt.show()


def simple_example():
    '''
    Simple example. This code can be run to show example of reading in
    data, adjusting using both cdfm and edcdfm, and finally producing
    histogram, cdfs and time-series to illustrate results.

    Load in sample data (here it is as time-series, but could be gridded)
    First load in observations (OD = Observation Data):

    >>> import adaq_data
    >>> import config
    >>> sample_datadir = config.SAMPLE_DATADIR+'sites_cube_list/'
    >>> od = adaq_data.ADAQData()
    >>> od_scl = od.load_ts(sample_datadir+'aurn_1days.nc')

    Then model data that is directly comparable to observations,
    ie for same dates (MDC = Model Data Current)

    >>> mdc = adaq_data.ADAQData()
    >>> mdc_scl = mdc.load_ts(sample_datadir+'aqum_oper_1days.nc')

    And finally model data for future/forecast, that will be adjusted
    (MD_FC = Model Data ForeCast)

    >>> mdfc = adaq_data.ADAQData()
    >>> mdfc_scl = mdfc.load_ts(sample_datadir+'aqum_casestudy_1days.nc')

    Now extract a single cube for a particular species from each of these
    (od_sc = Observation Data Site Cube etc)

    >>> species = 'O3'
    >>> od_sc = od.extract(short_name=species, singlecube=True)
    >>> mdc_sc = mdc.extract(short_name=species, singlecube=True)
    >>> mdfc_sc = mdfc.extract(short_name=species, singlecube=True)

    This data all needs to be binned, using the BinnedData class:

    >>> bd_od = BinnedData()
    >>> bd_od.add_cube(od_sc)
    >>> bd_mdc = BinnedData()
    >>> bd_mdc.add_cube(mdc_sc)
    >>> bd_mdfc = BinnedData()
    >>> bd_mdfc.add_cube(mdfc_sc)

    Further cubes for different time periods etc could be added to the
    BinnedData objects at this stage, using add_cube again.

    CDFs (Cumulative Density Function) should then be generated from all
    the Binned Data. To make this easier, they should all be using the
    same x values, so get some comparable CDFs:

    >>> x, [cdf_od, cdf_mdc, cdf_mdfc] = bd_od.get_comparablecdfs(
    ... [bd_mdc,bd_mdfc])

    Spline functions should then be calculated for CDFm or EDCDFm.
    Here we shall use CDFm:

    >>> cdfm_spline = cdfm(x, cdf_od, cdf_mdc, cdf_mdfc)

    Now take a copy of the Forecast Model Data and get a pointer to the
    same species cube that we are working with.
    (AMDFC = Adjusted Model Data Forecast,
    amdfc_sc = adjusted model data forecast site cube)

    >>> import copy
    >>> amdfc = copy.deepcopy(mdfc)
    >>> for icube, cube in enumerate(mdfc.sites_cube_list):
    ...     if cube.attributes['short_name'] == species:
    ...         amdfc_sc = amdfc.sites_cube_list[icube]

    The spline function takes values from the input data and returns
    the adjusted value according to the CDFm (or EDCDFm) method.
    The spline function however only works on 1D data so need to flatten
    any cubes:

    >>> amdfc_sc.data.flat = cdfm_spline(mdfc_sc.data.flatten())

    This data has now been adjusted:

    >>> print('{:.3f} {:.3f} {:.3f}'.format(mdfc_sc.data.min(), \
mdfc_sc.data.max(), mdfc_sc.data.std()))
    11.754 92.466 16.910
    >>> print('{:.3f} {:.3f} {:.3f}'.format(amdfc_sc.data.min(), \
amdfc_sc.data.max(), amdfc_sc.data.std()))
    25.261 96.000 16.124

    It has also been adjusted in the original ADAQData class, so can
    extract new timeseries to check:

    >>> mdfc_har = mdfc.extract(short_name=species,abbrev='HAR',
    ... singlecube=True)
    >>> amdfc_har = amdfc.extract(short_name=species,abbrev='HAR',
    ... singlecube=True)
    >>> print('{:.3f} {:.3f}'.format(mdfc_har.data.min(), mdfc_har.data.max()))
    11.754 72.911
    >>> print('{:.3f} {:.3f}'.format(amdfc_har.data.min(), amdfc_har.data.max()))
    25.261 79.911
    '''

    import config
    import adaq_data
    import timeseries_plot
    import copy
    sample_datadir = config.SAMPLE_DATADIR+'sites_cube_list/'
    species = 'O3'

    #Load in sample data (here it is as time-series, but could be gridded)
    od = adaq_data.ADAQData()
    od.load_ts(sample_datadir+'AURN_obs.nc')
    od_sc = od.extract(short_name=species, singlecube=True)
    md1 = adaq_data.ADAQData()
    md1.load_ts(sample_datadir+'aqum_oper.nc')
    md1_sc = md1.extract(short_name=species, singlecube=True)
    md2 = adaq_data.ADAQData()
    md2.load_ts(sample_datadir+'aqum_casestudy.nc')
    md2_sc = md2.extract(short_name=species, singlecube=True)

    #Add species site cubes to Bins
    bd_od = BinnedData()
    bd_od.add_cube(od_sc)
    bd_md1 = BinnedData()
    bd_md1.add_cube(md1_sc)
    bd_md2 = BinnedData()
    bd_md2.add_cube(md2_sc)

    #Get cdfs
    x, [cdf_od, cdf_md1, cdf_md2] = bd_od.get_comparablecdfs([bd_md1, bd_md2])

    #Now calculate splines for cdfm and edcdfm
    cdfm_spline = cdfm(x, cdf_od, cdf_md1, cdf_md2)
    edcdfm_spline = edcdfm(x, cdf_od, cdf_md1, cdf_md2)

    #Apply correction values to MD2, to get adjusted values, AMD2
    amd2_cdfm = copy.deepcopy(md2)
    amd2_edcdfm = copy.deepcopy(md2)

    #Get a pointer to the correct cube in the sites_cube_list
    for icube, cube in enumerate(md2.sites_cube_list):
        if cube.attributes['short_name'] == species:
            amd2_cdfm_sc = amd2_cdfm.sites_cube_list[icube]
            amd2_edcdfm_sc = amd2_edcdfm.sites_cube_list[icube]
            break

    #Can now adjust data
    #Spline only takes 1D data, cube data here is 2D, flatten data first
    amd2_cdfm_sc.data.flat = cdfm_spline(md2_sc.data.flatten())
    amd2_edcdfm_sc.data.flat = edcdfm_spline(md2_sc.data.flatten())

    print('Input mean, std dev, max, min:', \
        md2_sc.data.mean(), md2_sc.data.std(), \
        md2_sc.data.max(), md2_sc.data.min())
    print('Output CDFm mean, std dev, max, min:', \
        amd2_cdfm_sc.data.mean(), amd2_cdfm_sc.data.std(), \
        amd2_cdfm_sc.data.max(), amd2_cdfm_sc.data.min())
    print('Output EDCDFm mean, std dev, max, min:', \
        amd2_edcdfm_sc.data.mean(), amd2_edcdfm_sc.data.std(), \
        amd2_edcdfm_sc.data.max(), amd2_edcdfm_sc.data.min())

    #---- Data now adjusted ---#
    #Can now plot as required.

    #Bin data
    bd_amd2_cdfm = BinnedData()
    bd_amd2_cdfm.add_cube(amd2_cdfm_sc)
    bd_amd2_edcdfm = BinnedData()
    bd_amd2_edcdfm.add_cube(amd2_edcdfm_sc)

    #Convert data back to CDFs.
    x_amd2_cdfm, cdf_amd2_cdfm = bd_amd2_cdfm.get_cdf()
    x_amd2_edcdfm, cdf_amd2_edcdfm = bd_amd2_edcdfm.get_cdf()

    #Plot cdfs of initial data
    plt.figure()
    plt.plot(x, cdf_od, 'k', label='OD')
    plt.plot(x, cdf_md1, label='MD1')
    plt.plot(x, cdf_md2, label='MD2')

    #Finish plotting CDFs with adjusted data.
    plt.plot(x_amd2_cdfm, cdf_amd2_cdfm, label='AMD2_CDFm')
    plt.plot(x_amd2_edcdfm, cdf_amd2_edcdfm, label='AMD2_EDCDFm')
    plt.gca().legend(loc='lower right')
    plt.title('Cumulative Distribution Functions')

    #Get and Plot histograms
    plt.figure()
    x_od, od_hist = bd_od.get_histarray()
    plt.plot(x_od, od_hist, 'k', label='OD')
    x_md1, md1_hist = bd_md1.get_histarray()
    plt.plot(x_md1, md1_hist, label='MD1')
    x_md2, md2_hist = bd_md2.get_histarray()
    plt.plot(x_md2, md2_hist, label='MD2')
    x_amd2_cdfm, amd2_cdfm_hist = bd_amd2_cdfm.get_histarray()
    plt.plot(x_amd2_cdfm, amd2_cdfm_hist, label='AMD2 CDFm')
    x_amd2_edcdfm, amd2_edcdfm_hist = bd_amd2_edcdfm.get_histarray()
    plt.plot(x_amd2_edcdfm, amd2_edcdfm_hist, label='AMD2 EDCDFm')
    plt.title('Histograms')
    plt.gca().legend()

    #Get and Plot some time-series
    od_cube = od.extract(short_name=species,
                         abbrev='HAR', singlecube=True)
    md1_cube = md1.extract(short_name=species,
                           abbrev='HAR', singlecube=True)
    md2_cube = md2.extract(short_name=species,
                           abbrev='HAR', singlecube=True)
    amd2_cdfm_cube = amd2_cdfm.extract(short_name=species,
                                       abbrev='HAR', singlecube=True)
    amd2_edcdfm_cube = amd2_edcdfm.extract(short_name=species,
                                           abbrev='HAR', singlecube=True)

    tsp = timeseries_plot.TimeSeriesPlot()
    tsp.add_line(od_cube, label='OD', colour='k')
    tsp.add_line(md1_cube, label='MD1')
    tsp.add_line(md2_cube, label='MD2')
    tsp.add_line(amd2_cdfm_cube, label='AMD2_CDFm')
    tsp.add_line(amd2_edcdfm_cube, label='AMD2_EDCDFm')
    tsp.plot()

    plt.show()


if __name__ == '__main__':

    #idealised_example()
    #simple_example()

    import doctest
    doctest.testmod()
