#pragma once
#include "hamiltonian.h"
#include <vector>

using namespace std;

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega);
    double computeLocalEnergy(bool interaction, double GibbsValue, vector<double> X, vector<double> H, vector<double> a, vector<double> b,
                              vector<vector<double>> w);

    double omega() const;
    void setOmega(const double &omega);

private:
    double m_omega = 0;
};

