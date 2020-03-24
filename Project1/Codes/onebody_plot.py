import matplotlib.pyplot as plt

fonts = {"font.size": 14}
plt.rcParams.update(fonts)

def read_file(filename):
    prob = []; r = []
    with open("Data/"+filename+".dat", "r") as infile:
        lines = infile.readlines()
        for i in range(len(lines)):
            line = lines[i]
            vals = line.split()
            prob.append(float(vals[0]))
            r.append(float(vals[1]))

    for j in range(len(prob)):
        if prob[j] > 1.0:
            del r[j]
    for i in prob[:]:
        if i > 1.0:
            prob.remove(i)

    return prob, r


probEJ100, rEJ100 = read_file("density_E_J_alpha_0.522644_beta_2.828430_N_100")
probE100, rE100 = read_file("density_E_noJ_alpha_0.522644_beta_2.828430_N_100")

probEJ50, rEJ50 = read_file("density_E_J_alpha_0.516005_beta_2.828430_N_50")

probE50, rE50 = read_file("density_E_noJ_alpha_0.516005_beta_2.828430_N_50")

probEJ10, rEJ10 = read_file("density_E_J_alpha_0.503261_beta_2.828430_N_10")
probE10, rE10 = read_file("density_E_noJ_alpha_0.503261_beta_2.828430_N_10")


plt.figure()
plt.title("One-body density for N=10 with \n and without Jastrow factor")
plt.plot(rE10, probE10, "+", label="N=10, no-Jastrow")
plt.plot(rEJ10, probEJ10, "+", label="N=10, Jastrow")
plt.xlabel("r/$a_{ho}$")
plt.ylabel("Probability")
plt.legend(loc="best")
plt.tight_layout()
#plt.savefig("Figures/onebody_N_10.png")

plt.figure()
plt.title("One-body density for N=50 with \n and without Jastrow factor")
plt.plot(rE50, probE50, "+", label="N=50, no-Jastrow")
plt.plot(rEJ50, probEJ50, "+", label="N=50, Jastrow")
plt.xlabel("r/$a_{ho}$")
plt.ylabel("Probability")
plt.legend(loc="best")
plt.tight_layout()
#plt.savefig("Figures/onebody_N_50.png")

plt.figure()
plt.title("One-body density for N=100 with \n and without Jastrow factor")
plt.plot(rE100, probE100, "+", label="N=100, no-Jastrow")
plt.plot(rEJ100, probEJ100, "+", label="N=100, Jastrow")
plt.xlabel("r/$a_{ho}$")
plt.ylabel("Probability")
plt.legend(loc="best")
plt.tight_layout()
#plt.savefig("Figures/onebody_N_100.png")

plt.show()