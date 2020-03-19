import numpy as np
import matplotlib.pyplot as plt


def read_file(filename):
    prob = []; r = []
    with open("Data/"+filename+".dat", "r") as infile:
        lines = infile.readlines()
        for i in range(len(lines)):
            line = lines[i]
            vals = line.split()
            prob.append(float(vals[0]))
            r.append(float(vals[1]))

    return prob, r

# Non-int onebody density
prob1, r1 = read_file("density_non-int_alpha_0.500004_beta_1.000000")
prob2, r2 = read_file("density_non-int_alpha_0.500012_beta_1.000000")
prob3, r3 = read_file("density_non-int_alpha_0.500097_beta_1.000000")

# Int (S)
prob4, r4 = read_file("density_int_S_alpha_0.501750_beta_1.000000")
prob5, r5 = read_file("density_int_S_alpha_0.511570_beta_1.000000")
#prob6, r6 = read_file("")

# Int (E)
prob7, r7 = read_file("density_int_E_alpha_0.501750_beta_2.828430")
#prob8, r8 = read_file("")
#prob9, r9 = read_file("")

#plt.plot(r1, prob1, "+", label="N=10, non-int")
#plt.plot(r2, prob2, "+", label="N=50, non-int")
#plt.plot(r3, prob3, "+", label="N=100, non-int")

#plt.plot(r4, prob4, "+", label="N=10, int (S)")
#plt.plot(r5, prob5, "+", label="N=50, int (S)")

plt.plot(r7, prob7, "+", label="N=10, int (E)")


plt.legend()