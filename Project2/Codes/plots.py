import matplotlib.pyplot as plt
import numpy as np

"""
Program that produces many plots for Metropolis, Importance and Gibbs sampling
without and with interaction by reading data from files.
"""

fonts = {"font.size": 15}
plt.rcParams.update(fonts)

def read_file(filename):
    """
    Function to read a data file containing energies and export the data
    to an array to be plotted.
    """
    Energy = []
    with open("Data/"+filename+".dat", "r") as infile:
            lines = infile.readlines()
            for i in range(1,len(lines)):
                line = lines[i]
                vals = line.split()
                Energy.append(float(vals[0]))
    return Energy

def LinReg(equil, E, SGD, label, filename, write=False):
    """
    Function which does a linear regression on given energies from data files
    after the system has reached equilibrium. Used to more easily check when
    the system reaches equilibrium for the interacting system.
    If write=True, export the energies after equilibrium to another data file
    which is then imported by the blocking.py file.
    """
    #import pandas as pd
    #E_data = pd.DataFrame(data=E)
    #print(E_data)
    x = SGD[equil:]
    y = E[equil:]
    model = np.polyfit(x,y,1)
    #print(model)
    predict = np.poly1d(model)

    x_lin_reg = np.linspace(int(x[0]), int(x[-1]), len(x))
    y_lin_reg = predict(x_lin_reg)

    plt.figure()
    plt.title("Linear regression, %s" %label)
    plt.plot(SGD[equil:], E[equil:], "x")
    plt.plot(x_lin_reg, y_lin_reg)
    plt.ylabel("Local energy, $E_L$")
    plt.xlabel("SGD cycle")
    plt.tight_layout()

    if write == True:
        E_write = open("Data/" + filename + ".txt", "w")
        for i in range(len(x_lin_reg)):
            E_write.write("%f \n" %y_lin_reg[i])
        E_write.close()


"""Standard Metropolis"""
# Test different step lengths
E_001 = read_file("BruteCycle_0.010000_n_0.100000_Np_1_D_1_H_2")
E_01 = read_file("BruteCycle_0.100000_n_0.100000_Np_1_D_1_H_2")
E_05 = read_file("BruteCycle_0.500000_n_0.100000_Np_1_D_1_H_2")
E_1 = read_file("BruteCycle_1.000000_n_0.100000_Np_1_D_1_H_2")
E_2 = read_file("BruteCycle_2.000000_n_0.100000_Np_1_D_1_H_2")

# Test different learning rates
E_n_001 = read_file("BruteCycle_0.500000_n_0.010000_Np_1_D_1_H_2")
E_n_005 = read_file("BruteCycle_0.500000_n_0.050000_Np_1_D_1_H_2")
E_n_01 = read_file("BruteCycle_0.500000_n_0.100000_Np_1_D_1_H_2")
E_n_02 = read_file("BruteCycle_0.500000_n_0.200000_Np_1_D_1_H_2")
E_n_05 = read_file("BruteCycle_0.500000_n_0.500000_Np_1_D_1_H_2")
E_n_06 = read_file("BruteCycle_0.500000_n_0.600000_Np_1_D_1_H_2")

SGD_cycle = np.linspace(0, 100, len(E_01))

# Plot
plt.figure(figsize=(8, 5))
plt.title("Energy for different step lengths\n with Metropolis sampling")
plt.plot(SGD_cycle, E_001, label="L=0.01")
plt.ylabel("Local energy, $E_L$")
plt.xlabel("SGD cycle")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Figures/L_001.png")

plt.figure(figsize=(8, 5))
plt.title("Energy for different step lengths\n with Metropolis sampling")
plt.plot(SGD_cycle, E_01, label="L=0.1")
plt.plot(SGD_cycle, E_05, label="L=0.5")
plt.plot(SGD_cycle, E_1, label="L=1.0")
plt.plot(SGD_cycle, E_2, label="L=2.0")
plt.ylabel("Local energy, $E_L$")
plt.xlabel("SGD cycle")
plt.legend(loc="best")
plt.tight_layout()
#plt.savefig("Figures/L_01_20.png")

plt.figure(figsize=(8, 5))
plt.title("Energy for different learning rates with\n Metropolis sampling")
plt.plot(SGD_cycle, E_n_001, "--", label="$\\eta=0.01$")
plt.plot(SGD_cycle, E_n_005, "--", label="$\\eta=0.05$")
plt.plot(SGD_cycle, E_n_01, "--", label="$\\eta=0.1$")
plt.plot(SGD_cycle, E_n_02, "--", label="$\\eta=0.2$")
plt.plot(SGD_cycle, E_n_05, "--", label="$\\eta=0.5$")
plt.plot(SGD_cycle, E_n_06, "--", label="$\\eta=0.6$")
plt.ylabel("Local energy, $E_L$")
plt.xlabel("SGD cycle")
plt.legend(loc="best")
plt.tight_layout()
#plt.savefig("Figures/Met_learnrate.png")


"""Importance Sampling"""
# Test different time steps
EI_0001 = read_file("ImpCycle_0.001000_n_0.100000_Np_1_D_1_H_2")
EI_001 = read_file("ImpCycle_0.010000_n_0.100000_Np_1_D_1_H_2")
EI_005 = read_file("ImpCycle_0.050000_n_0.100000_Np_1_D_1_H_2")
EI_01 = read_file("ImpCycle_0.100000_n_0.100000_Np_1_D_1_H_2")
EI_05 = read_file("ImpCycle_0.500000_n_0.100000_Np_1_D_1_H_2")
EI_1 = read_file("ImpCycle_1.000000_n_0.100000_Np_1_D_1_H_2")
EI_2 = read_file("ImpCycle_2.000000_n_0.100000_Np_1_D_1_H_2")

# Test different learning rates
EI_n_001 = read_file("ImpCycle_0.500000_n_0.010000_Np_1_D_1_H_2")
EI_n_005 = read_file("ImpCycle_0.500000_n_0.050000_Np_1_D_1_H_2")
EI_n_01 = read_file("ImpCycle_0.500000_n_0.100000_Np_1_D_1_H_2")
EI_n_02 = read_file("ImpCycle_0.500000_n_0.200000_Np_1_D_1_H_2")
EI_n_05 = read_file("ImpCycle_0.500000_n_0.500000_Np_1_D_1_H_2")
EI_n_06 = read_file("ImpCycle_0.500000_n_0.600000_Np_1_D_1_H_2")

SGD_cycle = np.linspace(0, 100, len(EI_001))

# Plot
plt.figure(figsize=(8, 5))
plt.title("Energy for different $\\Delta t$\n with Importance sampling")
plt.plot(SGD_cycle, EI_0001, label="$\\Delta$t=0.001")
plt.plot(SGD_cycle, EI_001, label="$\\Delta$t=0.01")
plt.plot(SGD_cycle, EI_005, label="$\\Delta$t=0.05")
plt.plot(SGD_cycle, EI_01, label="$\\Delta$t=0.1")
plt.plot(SGD_cycle, EI_05, label="$\\Delta$t=0.5")
plt.plot(SGD_cycle, EI_1, label="$\\Delta$t=1.0")
plt.plot(SGD_cycle, EI_2, label="$\\Delta$t=2.0")
plt.ylabel("Local energy, $E_L$")
plt.xlabel("SGD cycle")
plt.legend(loc="best")
plt.tight_layout()
#plt.savefig("Figures/dt.png")

plt.figure(figsize=(8, 5))
plt.title("Energy for different learning rates\n with Importance sampling")
plt.plot(SGD_cycle, EI_n_001, "--", label="$\\eta=0.01$")
plt.plot(SGD_cycle, EI_n_005, "--", label="$\\eta=0.05$")
plt.plot(SGD_cycle, EI_n_01, "--", label="$\\eta=0.1$")
plt.plot(SGD_cycle, EI_n_02, "--", label="$\\eta=0.2$")
plt.plot(SGD_cycle, EI_n_05, "--", label="$\\eta=0.5$")
plt.plot(SGD_cycle, EI_n_06, "--", label="$\\eta=0.6$")
plt.ylabel("Local energy, $E_L$")
plt.xlabel("SGD cycle")
plt.legend(loc="best")
plt.tight_layout()
#plt.savefig("Figures/Imp_learnrate.png")


"""Gibbs Sampling"""
# Test different sigma values
EG_05 = read_file("GibCycle_0.500000_n_0.100000_Np_1_D_1_H_2")
EG_06 = read_file("GibCycle_0.600000_n_0.100000_Np_1_D_1_H_2")
EG_07 = read_file("GibCycle_0.700000_n_0.100000_Np_1_D_1_H_2")
EG_08 = read_file("GibCycle_0.800000_n_0.100000_Np_1_D_1_H_2")
EG_09 = read_file("GibCycle_0.900000_n_0.100000_Np_1_D_1_H_2")
EG_10 = read_file("GibCycle_1.000000_n_0.100000_Np_1_D_1_H_2")

# Test different learning rates
EG_n_001 = read_file("GibCycle_0.600000_n_0.010000_Np_1_D_1_H_2")
EG_n_005 = read_file("GibCycle_0.600000_n_0.050000_Np_1_D_1_H_2")
EG_n_01 = read_file("GibCycle_0.600000_n_0.100000_Np_1_D_1_H_2")
EG_n_02 = read_file("GibCycle_0.600000_n_0.200000_Np_1_D_1_H_2")
EG_n_05 = read_file("GibCycle_0.600000_n_0.500000_Np_1_D_1_H_2")
EG_n_06 = read_file("GibCycle_0.600000_n_0.600000_Np_1_D_1_H_2")

EG_n7_001 = read_file("GibCycle_0.700000_n_0.010000_Np_1_D_1_H_2")
EG_n7_005 = read_file("GibCycle_0.700000_n_0.050000_Np_1_D_1_H_2")
EG_n7_01 = read_file("GibCycle_0.700000_n_0.100000_Np_1_D_1_H_2")
EG_n7_02 = read_file("GibCycle_0.700000_n_0.200000_Np_1_D_1_H_2")
EG_n7_05 = read_file("GibCycle_0.700000_n_0.500000_Np_1_D_1_H_2")
EG_n7_06 = read_file("GibCycle_0.700000_n_0.600000_Np_1_D_1_H_2")

SGD_cycle = np.linspace(0, 100, len(EG_05))

# Plot
plt.figure(figsize=(8, 5))
plt.title("Energy for different $\\sigma$ with Gibbs sampling")
plt.plot(SGD_cycle, EG_05, label="$\\sigma$=0.5")
plt.plot(SGD_cycle, EG_06, label="$\\sigma$=0.6")
plt.plot(SGD_cycle, EG_07, label="$\\sigma$=0.7")
plt.plot(SGD_cycle, EG_08, label="$\\sigma$=0.8")
plt.plot(SGD_cycle, EG_09, label="$\\sigma$=0.9")
plt.plot(SGD_cycle, EG_10, label="$\\sigma$=1.0")
plt.ylabel("Local energy, $E_L$")
plt.xlabel("SGD cycle")
plt.legend(loc="best")
plt.tight_layout()
#plt.savefig("Figures/sigma.png")

plt.figure(figsize=(8, 5))
plt.title("Energy for different learning rates with\n Gibbs sampling")
plt.plot(SGD_cycle, EG_n_001, "--", label="$\\eta=0.01$")
plt.plot(SGD_cycle, EG_n_005, "--", label="$\\eta=0.05$")
plt.plot(SGD_cycle, EG_n_01, "--", label="$\\eta=0.1$")
plt.plot(SGD_cycle, EG_n_02, "--", label="$\\eta=0.2$")
plt.plot(SGD_cycle, EG_n_05, "--", label="$\\eta=0.5$")
plt.plot(SGD_cycle, EG_n_06, "--", label="$\\eta=0.6$")
plt.ylabel("Local energy, $E_L$")
plt.xlabel("SGD cycle")
plt.legend(loc="best", ncol=2)
plt.tight_layout()
#plt.savefig("Figures/Gib_learnrate.png")

plt.figure(figsize=(8, 5))
plt.title("Test equilibrium speed - Gibbs sampling")
plt.plot(SGD_cycle, EG_n_01, "--", label="$\\sigma=0.6, \\eta=0.1$")
plt.plot(SGD_cycle, EG_n_02, "--", label="$\\sigma=0.6, \\eta=0.2$")
plt.plot(SGD_cycle, EG_n7_01, "--", label="$\\sigma=0.7, \\eta=0.1$")
plt.plot(SGD_cycle, EG_n7_02, "--", label="$\\sigma=0.7, \\eta=0.2$")
plt.ylabel("Local energy, $E_L$")
plt.xlabel("SGD cycle")
plt.legend(loc="best", ncol=2)
plt.tight_layout()
#plt.savefig("Figures/Gib_learnrate_test.png")


""" Interaction Metropolis"""
# Test different learning rates
E01_int = read_file("BruteCycleInt_0.500000_n_0.100000_Np_2_D_2_H_2")
E02_int = read_file("BruteCycleInt_0.500000_n_0.200000_Np_2_D_2_H_2")
E05_int = read_file("BruteCycleInt_0.500000_n_0.500000_Np_2_D_2_H_2")
E06_int = read_file("BruteCycleInt_0.500000_n_0.600000_Np_2_D_2_H_2")

SGD_cycle = np.linspace(0, 600, len(E01_int))

# Plot
plt.figure(figsize=(9, 8))
plt.subplot(211)
plt.title("Energy for different learning rates with\n interacting Metropolis sampling")
plt.plot(SGD_cycle[3:], E01_int[3:], "+", label="$\\eta=0.1$")
plt.plot(SGD_cycle[3:], E02_int[3:], "x", label="$\\eta=0.2$")
plt.ylabel("Local energy, $E_L$")
plt.legend(loc="best", ncol=2)
plt.tight_layout()
plt.subplot(212)
plt.plot(SGD_cycle[1:], E05_int[1:], "+", label="$\\eta=0.5$")
plt.plot(SGD_cycle[1:], E06_int[1:], "x", label="$\\eta=0.6$")
plt.ylabel("Local energy, $E_L$")
plt.xlabel("SGD cycle")
plt.legend(loc="best", ncol=2)
plt.tight_layout()
#plt.savefig("Figures/Int_Met_learnrate.png")

# Linear regression of data
LinReg(400, E01_int, SGD_cycle, label="$\\eta=0.1$", filename="BruteReg01", write=False)
LinReg(300, E02_int, SGD_cycle, label="$\\eta=0.2$", filename="BruteReg02", write=False)
LinReg(300, E05_int, SGD_cycle, label="$\\eta=0.5$", filename="BruteReg05", write=False)
LinReg(300, E06_int, SGD_cycle, label="$\\eta=0.6$", filename="BruteReg06", write=False)

# Test different number of hidden nodes
E_H1_int = read_file("BruteCycleInt_0.500000_n_0.200000_Np_2_D_2_H_1")
E_H2_int = read_file("BruteCycleInt_0.500000_n_0.200000_Np_2_D_2_H_2")
E_H3_int = read_file("BruteCycleInt_0.500000_n_0.200000_Np_2_D_2_H_3")
E_H4_int = read_file("BruteCycleInt_0.500000_n_0.200000_Np_2_D_2_H_4")

# Plot
plt.figure(figsize=(9, 8))
plt.subplot(211)
plt.title("Energy for different # of hidden nodes with\n interacting Metropolis sampling")
plt.plot(SGD_cycle[3:], E_H1_int[3:], "x", label="$N_H=1$")
plt.plot(SGD_cycle[3:], E_H2_int[3:], "+", color="green", label="$N_H=2$")
plt.ylabel("Local energy, $E_L$")
plt.legend(loc="best", ncol=2)
plt.tight_layout()
plt.subplot(212)
plt.plot(SGD_cycle[3:], E_H3_int[3:], "x", label="$N_H=3$")
plt.plot(SGD_cycle[3:], E_H4_int[3:], "+", color="green", label="$N_H=4$")
plt.ylabel("Local energy, $E_L$")
plt.xlabel("SGD cycle")
plt.legend(loc="best", ncol=2)
plt.tight_layout()
#plt.savefig("Figures/Hidden_Int_Met.png")

# Linear regression of data
LinReg(300, E_H1_int, SGD_cycle, label="$N_H=1$", filename="BruteRegH1", write=True)
LinReg(300, E_H2_int, SGD_cycle, label="$N_H=2$", filename="BruteRegH2", write=True)
LinReg(300, E_H3_int, SGD_cycle, label="$N_H=3$", filename="BruteRegH3", write=True)
LinReg(300, E_H4_int, SGD_cycle, label="$N_H=4$", filename="BruteRegH4", write=True)


""" Interaction Importance Sampling"""
# Test different learning rates
EI01_int = read_file("ImpCycleInt_0.500000_n_0.100000_Np_2_D_2_H_2")
EI02_int = read_file("ImpCycleInt_0.500000_n_0.200000_Np_2_D_2_H_2")
EI05_int = read_file("ImpCycleInt_0.500000_n_0.500000_Np_2_D_2_H_2")
EI06_int = read_file("ImpCycleInt_0.500000_n_0.600000_Np_2_D_2_H_2")

SGD_cycle = np.linspace(0, 600, len(EI01_int))

# Plot
plt.figure(figsize=(9, 8))
plt.subplot(211)
plt.title("Energy for different learning rates with\n interacting Importance sampling")
plt.plot(SGD_cycle[5:], EI01_int[5:], "+", label="$\\eta=0.1$")
plt.plot(SGD_cycle[5:], EI02_int[5:], "x", label="$\\eta=0.2$")
plt.ylabel("Local energy, $E_L$")
plt.legend(loc="best", ncol=2)
plt.tight_layout()
plt.subplot(212)
plt.plot(SGD_cycle[2:], EI05_int[2:], "+", label="$\\eta=0.5$")
plt.plot(SGD_cycle[2:], EI06_int[2:], "x", label="$\\eta=0.6$")
plt.ylabel("Local energy, $E_L$")
plt.xlabel("SGD cycle")
plt.legend(loc="best", ncol=2)
plt.tight_layout()
#plt.savefig("Figures/Int_Imp_learnrate.png")

# Linear regression of data
LinReg(450, EI01_int, SGD_cycle, label="$\\eta=0.1$", filename="ImpReg01", write=False)
LinReg(250, EI02_int, SGD_cycle, label="$\\eta=0.2$", filename="ImpReg02", write=False)
LinReg(200, EI05_int, SGD_cycle, label="$\\eta=0.5$", filename="ImpReg05", write=False)
LinReg(200, EI06_int, SGD_cycle, label="$\\eta=0.6$", filename="ImpReg06", write=False)

# Test different number of hidden nodes
EI_H1_int = read_file("ImpCycleInt_0.500000_n_0.200000_Np_2_D_2_H_1")
EI_H2_int = read_file("ImpCycleInt_0.500000_n_0.200000_Np_2_D_2_H_2")
EI_H3_int = read_file("ImpCycleInt_0.500000_n_0.200000_Np_2_D_2_H_3")
EI_H4_int = read_file("ImpCycleInt_0.500000_n_0.200000_Np_2_D_2_H_4")

# Plot
plt.figure(figsize=(9, 8))
plt.subplot(211)
plt.title("Energy for different # of hidden nodes with\n interacting Importance sampling")
plt.plot(SGD_cycle[3:], EI_H1_int[3:], "x", label="$N_H=1$")
plt.plot(SGD_cycle[3:], EI_H2_int[3:], "+", color="green", label="$N_H=2$")
plt.ylabel("Local energy, $E_L$")
plt.legend(loc="best", ncol=2)
plt.tight_layout()
plt.subplot(212)
plt.plot(SGD_cycle[3:], EI_H3_int[3:], "x", label="$N_H=3$")
plt.plot(SGD_cycle[3:], EI_H4_int[3:], "+", color="green", label="$N_H=4$")
plt.ylabel("Local energy, $E_L$")
plt.xlabel("SGD cycle")
plt.legend(loc="best", ncol=2)
plt.tight_layout()
#plt.savefig("Figures/Hidden_Int_Imp.png")

# Linear regression of data
LinReg(250, EI_H1_int, SGD_cycle, label="$N_H=1$", filename="ImpRegH1", write=True)
LinReg(250, EI_H2_int, SGD_cycle, label="$N_H=2$", filename="ImpRegH2", write=True)
LinReg(250, EI_H3_int, SGD_cycle, label="$N_H=3$", filename="ImpRegH3", write=True)
LinReg(250, EI_H4_int, SGD_cycle, label="$N_H=4$", filename="ImpRegH4", write=True)


""" Interaction Gibbs Sampling"""
# Test different learning rates
EG01_int = read_file("GibCycleInt_0.700000_n_0.100000_Np_2_D_2_H_2")
EG02_int = read_file("GibCycleInt_0.700000_n_0.200000_Np_2_D_2_H_2")
EG05_int = read_file("GibCycleInt_0.700000_n_0.500000_Np_2_D_2_H_2")
EG06_int = read_file("GibCycleInt_0.700000_n_0.600000_Np_2_D_2_H_2")

SGD_cycle = np.linspace(0, 600, len(EG01_int))

# Plot
plt.figure(figsize=(9, 8))
plt.subplot(211)
plt.title("Energy for different learning rates with\n interacting Gibbs sampling")
plt.plot(SGD_cycle[5:], EG01_int[5:], "+", label="$\\eta=0.1$")
plt.plot(SGD_cycle[5:], EG02_int[5:], "x", label="$\\eta=0.2$")
plt.ylabel("Local energy, $E_L$")
plt.legend(loc="best", ncol=2)
plt.tight_layout()
plt.subplot(212)
plt.plot(SGD_cycle[5:], EG05_int[5:], "+", label="$\\eta=0.5$")
plt.plot(SGD_cycle[5:], EG06_int[5:], "x", label="$\\eta=0.6$")
plt.ylabel("Local energy, $E_L$")
plt.xlabel("SGD cycle")
plt.legend(loc="best", ncol=2)
plt.tight_layout()
#plt.savefig("Figures/Int_Gib_learnrate.png")

# Linear regression of data
LinReg(200, EG01_int, SGD_cycle, label="$\\eta=0.1$", filename="GibReg01", write=False)
LinReg(200, EG02_int, SGD_cycle, label="$\\eta=0.2$", filename="GibReg02", write=False)
LinReg(200, EG05_int, SGD_cycle, label="$\\eta=0.5$", filename="GibReg05", write=False)
LinReg(200, EG06_int, SGD_cycle, label="$\\eta=0.6$", filename="GibReg06", write=False)

# Test different number of hidden nodes
EG_H1_int = read_file("GibCycleInt_0.700000_n_0.200000_Np_2_D_2_H_1")
EG_H2_int = read_file("GibCycleInt_0.700000_n_0.200000_Np_2_D_2_H_2")
EG_H3_int = read_file("GibCycleInt_0.700000_n_0.200000_Np_2_D_2_H_3")
EG_H4_int = read_file("GibCycleInt_0.700000_n_0.200000_Np_2_D_2_H_4")

# Plot
plt.figure(figsize=(9, 8))
plt.subplot(211)
plt.title("Energy for different # of hidden nodes with\n interacting Gibbs sampling")
plt.plot(SGD_cycle[5:], EG_H1_int[5:], "x", label="$N_H=1$")
plt.plot(SGD_cycle[5:], EG_H2_int[5:], "+", color="green", label="$N_H=2$")
plt.ylabel("Local energy, $E_L$")
plt.legend(loc="best", ncol=2)
plt.tight_layout()
plt.subplot(212)
plt.plot(SGD_cycle[5:], EG_H3_int[5:], "x", label="$N_H=3$")
plt.plot(SGD_cycle[5:], EG_H4_int[5:], "+", color="green", label="$N_H=4$")
plt.ylabel("Local energy, $E_L$")
plt.xlabel("SGD cycle")
plt.legend(loc="best", ncol=2)
plt.tight_layout()
#plt.savefig("Figures/Hidden_Int_Gib.png")

# Linear regression of data
LinReg(200, EG_H1_int, SGD_cycle, label="$N_H=1$", filename="GibRegH1", write=True)
LinReg(200, EG_H2_int, SGD_cycle, label="$N_H=2$", filename="GibRegH2", write=True)
LinReg(200, EG_H3_int, SGD_cycle, label="$N_H=3$", filename="GibRegH3", write=True)
LinReg(200, EG_H4_int, SGD_cycle, label="$N_H=4$", filename="GibRegH4", write=True)
