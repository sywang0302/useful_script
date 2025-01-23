import numpy as np
from scipy.special import factorial
from numpy import log10, exp, log
# follow this paper's method https://articles.adsabs.harvard.edu/pdf/1991ApJ...374..344K
# N:observed number of counts
#B: background conts in source scale
# note that the input value is not an array
def kraftupper(N, B, CL=0.99):
    # Constants
    srcdim = 100000
    srcstp = 0.001
    
    # Set mnrt array
    mnrt = B + srcstp * np.arange(srcdim + 1)

    # Compute psrc based on the value of N
    if N <= 150:
        psrc = exp(-mnrt) * 10**(N * log10(mnrt) - log10(factorial(N)))
    else:
        factn = 0
        for jj in range(1, N + 1):
            factn += log10(jj)
        psrc = exp(-mnrt + log(10) * (N * log10(mnrt) - factn))

    # Find maximum value and index of psrc
    pmax = np.max(psrc)
    nmax = np.argmax(psrc)

    sss = nmax * srcstp

    # Normalize psrc
    t = np.sum(psrc[1:])  # Sum excluding the first element
    area = (psrc[0] + psrc[srcdim]) / 2 + t
    psrc = psrc / area

    # Initialize nupper, nlower, and area
    nupper = nmax
    nlower = nmax
    area = 0

    # Search for nupper and nlower values based on psrc
    for i in range(2):
        if psrc[nlower] == 0 or nlower == 0:
            nupper = nupper + 1
        else:
            if psrc[nlower - 1] >= psrc[nupper + 1]:
                nlower = nlower - 1
            else:
                nupper = nupper + 1

    # Calculate area after searching for nupper and nlower
    area = (psrc[nupper] + psrc[nlower]) / 2 + psrc[nlower + 1]
    k = 0

    # Main loop to find smin and smax
    for i in range(srcdim - 3):
        if k == 0:
            if nlower == 0:
                nupper = nupper + 1
                area = area + (psrc[nupper - 1] + psrc[nupper]) / 2
            else:
                if psrc[nlower - 1] >= psrc[nupper + 1]:
                    nlower = nlower - 1
                    area = area + (psrc[nlower + 1] + psrc[nlower]) / 2
                else:
                    nupper = nupper + 1
                    area = area + (psrc[nupper - 1] + psrc[nupper]) / 2
        
        if k == 0 and area >= CL:
            smax = nupper * srcstp
            smin = nlower * srcstp
            K = 1
            break

    return sss, smin, smax