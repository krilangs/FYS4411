from numpy import log2, zeros, mean, var, sum, loadtxt, arange, \
                  array, cumsum, floor
from time import time

"""
Program for doing a statistical analysis on data imported from data files
using the blocking method.
Uncomment the desired line to execute the blocking analysis.
"""

def block(x):
    # Preliminaries
    d = log2(len(x))
    if (d - floor(d) != 0):
        print("Warning: Data size = %g, is not a power of 2." % floor(2**d))
        print("Truncating data to %g." % 2**floor(d) )
        x = x[:2**int(floor(d))]
    d = int(floor(d))
    n = 2**d
    s, gamma = zeros(d), zeros(d)
    mu = mean(x)
    t0 = time();

    # Estimate the auto-covariance and variances
    # for each blocking transformation
    for i in arange(0,d):
        n = len(x)
        # Estimate autocovariance of x
        gamma[i] = (n)**(-1)*sum( (x[0:(n-1)]-mu)*(x[1:n]-mu) )
        # Estimate variance of x
        s[i] = var(x)
        # Perform blocking transformation
        x = 0.5*(x[0::2] + x[1::2])

    # Generate the test observator M_k from the theorem
    M = (cumsum( ((gamma/s)**2*2**arange(1,d+1)[::-1])[::-1] )  )[::-1]

    # We need a list of magic numbers
    q =array([6.634897,  9.210340,  11.344867, 13.276704, 15.086272,
              16.811894, 18.475307, 20.090235, 21.665994, 23.209251,
              24.724970, 26.216967, 27.688250, 29.141238, 30.577914,
              31.999927, 33.408664, 34.805306, 36.190869, 37.566235,
              38.932173, 40.289360, 41.638398, 42.979820, 44.314105,
              45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

    # Use magic to determine when we should have stopped blocking
    for k in arange(0,d):
        if(M[k] < q[k]):
            break
    if (k >= d-1):
        print("Warning: Use more data")

    ans = s[k]/2**(d-k)
    print("Runtime: %g sec" % (time()-t0))
    print("Blocking Statistics :")
    print("Average E            Iterations      Std. error")
    print("%8g %20g %15g" % (mu, k, ans**.5))
    return ans

#-------------------------------------------
# Non-interacting brute Metropolis
#x = loadtxt("Data/FinalBruteCycle_0.500000_n_0.200000_Np_1_D_1_H_1.dat")
#x = loadtxt("Data/FinalBruteCycle_0.500000_n_0.200000_Np_1_D_1_H_2.dat")
#x = loadtxt("Data/FinalBruteCycle_0.500000_n_0.200000_Np_1_D_2_H_1.dat")
#x = loadtxt("Data/FinalBruteCycle_0.500000_n_0.200000_Np_1_D_2_H_2.dat")
#x = loadtxt("Data/FinalBruteCycle_0.500000_n_0.200000_Np_2_D_1_H_1.dat")
#x = loadtxt("Data/FinalBruteCycle_0.500000_n_0.200000_Np_2_D_1_H_2.dat")
#x = loadtxt("Data/FinalBruteCycle_0.500000_n_0.200000_Np_2_D_2_H_1.dat")
#x = loadtxt("Data/FinalBruteCycle_0.500000_n_0.200000_Np_2_D_2_H_2.dat")
#x = loadtxt("Data/FinalBruteCycle_0.500000_n_0.200000_Np_2_D_2_H_3.dat")
#x = loadtxt("Data/FinalBruteCycle_0.500000_n_0.200000_Np_2_D_2_H_4.dat")

# Non-interacting Importance sampling
#x = loadtxt("Data/FinalImpCycle_0.500000_n_0.200000_Np_1_D_1_H_1.dat")
#x = loadtxt("Data/FinalImpCycle_0.500000_n_0.200000_Np_1_D_1_H_2.dat")
#x = loadtxt("Data/FinalImpCycle_0.500000_n_0.200000_Np_1_D_2_H_1.dat")
#x = loadtxt("Data/FinalImpCycle_0.500000_n_0.200000_Np_1_D_2_H_2.dat")
#x = loadtxt("Data/FinalImpCycle_0.500000_n_0.200000_Np_2_D_1_H_1.dat")
#x = loadtxt("Data/FinalImpCycle_0.500000_n_0.200000_Np_2_D_1_H_2.dat")
#x = loadtxt("Data/FinalImpCycle_0.500000_n_0.200000_Np_2_D_2_H_1.dat")
#x = loadtxt("Data/FinalImpCycle_0.500000_n_0.200000_Np_2_D_2_H_2.dat")
#x = loadtxt("Data/FinalImpCycle_0.500000_n_0.200000_Np_2_D_2_H_3.dat")
#x = loadtxt("Data/FinalImpCycle_0.500000_n_0.200000_Np_2_D_2_H_4.dat")


# Non-interacting Gibbs sampling
#x = loadtxt("Data/FinalGibCycle_0.700000_n_0.200000_Np_1_D_1_H_1.dat")
#x = loadtxt("Data/FinalGibCycle_0.700000_n_0.200000_Np_1_D_1_H_2.dat")
#x = loadtxt("Data/FinalGibCycle_0.700000_n_0.200000_Np_1_D_2_H_1.dat")
#x = loadtxt("Data/FinalGibCycle_0.700000_n_0.200000_Np_1_D_2_H_2.dat")
#x = loadtxt("Data/FinalGibCycle_0.700000_n_0.200000_Np_2_D_1_H_1.dat")
#x = loadtxt("Data/FinalGibCycle_0.700000_n_0.200000_Np_2_D_1_H_2.dat")
#x = loadtxt("Data/FinalGibCycle_0.700000_n_0.200000_Np_2_D_2_H_1.dat")
#x = loadtxt("Data/FinalGibCycle_0.700000_n_0.200000_Np_2_D_2_H_2.dat")
#x = loadtxt("Data/FinalGibCycle_0.700000_n_0.200000_Np_2_D_2_H_3.dat")
#x = loadtxt("Data/FinalGibCycle_0.700000_n_0.200000_Np_2_D_2_H_4.dat")

#-------------------------------------------
# Interacting Metropolis sampling
#x = loadtxt("Data/BruteRegH1.txt")
#x = loadtxt("Data/BruteRegH2.txt")
#x = loadtxt("Data/BruteRegH3.txt")
#x = loadtxt("Data/BruteRegH4.txt")


# Interacting Importance sampling
#x = loadtxt("Data/ImpRegH1.txt")
#x = loadtxt("Data/ImpRegH2.txt")
#x = loadtxt("Data/ImpRegH3.txt")
#x = loadtxt("Data/ImpRegH4.txt")


# Interacting Gibbs sampling
#x = loadtxt("Data/GibRegH1.txt")
#x = loadtxt("Data/GibRegH2.txt")
#x = loadtxt("Data/GibRegH3.txt")
x = loadtxt("Data/GibRegH4.txt")


# Run statistical analysis with blocking method on data
block(x)
