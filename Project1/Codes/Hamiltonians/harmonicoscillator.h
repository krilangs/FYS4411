#pragma once
#include "hamiltonian.h"
#include <vector>

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega, double omega_z);
    double computeLocalEnergy(std::vector<Particle*> particles);

    std::vector<double> omega() const;
    void setOmega(const std::vector<double> &omega);

private:
    std::vector<double> m_omega = std::vector<double>();
};

