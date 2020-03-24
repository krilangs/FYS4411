import matplotlib.pyplot as plt

fonts = {"font.size": 14}
plt.rcParams.update(fonts)

def read_file(filename, GP):
    elem1 = []; elem2 = []
    elem3 = []
    if (GP=="Yes"): # Choose when we have GP-energies in the files to plot also
        with open("Data/"+filename+".txt", "r") as infile:
            lines = infile.readlines()
            for i in range(1,len(lines)):
                line = lines[i]
                vals = line.split()
                elem1.append(float(vals[0]))
                elem2.append(float(vals[1]))
                elem3.append(float(vals[2]))
        return elem1, elem2, elem3
    else:
        with open("Data/"+filename+".txt", "r") as infile:
            lines = infile.readlines()
            for i in range(1,len(lines)):
                line = lines[i]
                vals = line.split()
                elem1.append(float(vals[0]))
                elem2.append(float(vals[1]))
        return elem1, elem2

# Non-interacting system
alpha, E, E_GP = read_file("Check_alpha_non_int_d_1_N_1", "Yes")
r_step, A_r = read_file("Check_accept_non_int_d_1_N_1", "No")
time_step, A_t = read_file("Check_timestep_non_int_d_1_N_1", "No")

plt.figure()
plt.title("Local energy as function of $\\alpha$\n for N=1, d=1")
plt.plot(alpha, E, label="Analytic $E_L$")
plt.plot(alpha, E_GP, label="GP $E_L$")
plt.xlabel("Variational parameter, $\\alpha$")
plt.ylabel("Local energy, $E_L$")
plt.legend(loc="best")
plt.tight_layout()
#plt.savefig("Figures/Check_alphas.png")

plt.figure()
plt.title("Acceptance rate as function of\n MC step length for N=1, d=1")
plt.plot(r_step, A_r)
plt.xlabel("MC step length, $r_{step}$")
plt.ylabel("Acceptance rate/$10^5$")
plt.tight_layout()
#plt.savefig("Figures/Check_accept.png")

plt.figure()
plt.title("Acceptance rate as function of\n time step for N=1, d=1")
plt.semilogx(time_step, A_t)
plt.xlabel("Time step, $\\Delta t$")
plt.ylabel("Acceptance rate/$10^5$")
plt.tight_layout()
#plt.savefig("Figures/Check_timestep.png")

# Interacting system
alpha_int10, E_int10, E_GP_int10 = read_file("Check_alpha_int_d_3_N_10", "Yes")
alpha_int50, E_int50, E_GP_int50 = read_file("Check_alpha_int_d_3_N_50", "Yes")
alpha_int100, E_int100, E_GP_int100 = read_file("Check_alpha_int_d_3_N_100", "Yes")

plt.figure()
plt.title("Local energy as function of $\\alpha$\n for N=10, d=3 interaction")
plt.plot(alpha_int10, E_int10, label="Analytic $E_L$")
plt.plot(alpha_int10, E_GP_int10, label="GP $E_L$")
plt.xlabel("Variational parameter, $\\alpha$")
plt.ylabel("Local energy, $E_L$")
plt.legend(loc="best")
plt.tight_layout()
#plt.savefig("Figures/Check_alphas_int10.png")

plt.figure()
plt.title("Local energy as function of $\\alpha$\n for N=50, d=3 interaction")
plt.plot(alpha_int50, E_int50, label="Analytic $E_L$")
plt.plot(alpha_int50, E_GP_int50, label="GP $E_L$")
plt.xlabel("Variational parameter, $\\alpha$")
plt.ylabel("Local energy, $E_L$")
plt.legend(loc="best")
plt.tight_layout()
#plt.savefig("Figures/Check_alphas_int50.png")

plt.figure()
plt.title("Local energy as function of $\\alpha$\n for N=100, d=3 interaction")
plt.plot(alpha_int100, E_int100, label="Analytic $E_L$")
plt.plot(alpha_int100, E_GP_int100, label="GP $E_L$")
plt.xlabel("Variational parameter, $\\alpha$")
plt.ylabel("Local energy, $E_L$")
plt.legend(loc="best")
plt.tight_layout()
#plt.savefig("Figures/Check_alphas_int100.png")

plt.show()