#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

HarmonicOscillator::HarmonicOscillator(System* system, double omega, double omega_z):
                                       Hamiltonian(system) {
    assert(omega > 0);
    assert(omega_z > 0);
    m_omega.reserve(3);
    m_omega.push_back(omega);
    m_omega.push_back(omega);
    m_omega.push_back(omega_z);
}

double HarmonicOscillator::computeLocalEnergy(std::vector<Particle*> particles) {
// Compute the local energy from the kinetic and potential energies.
    double potentialEnergy = 0;
    double kineticEnergy   = 0;
    int dim = m_system->getNumberOfDimensions();
    int N   = m_system->getNumberOfParticles();

    kineticEnergy = -0.5*m_system->getWaveFunction()->computeDoubleDerivative(particles);   //Analytical double derivative
    //kineticEnergy = -0.5*m_system->getHamiltonian()->computeNumericalDoubleDerivative(particles);   //Numerical double derivative

    for (int k=0; k<N; k++ ){
        for (int d=0; d<dim; d++){
            potentialEnergy += m_omega[d]*m_omega[d]*particles.at(k)->getPosition()[d]*particles.at(k)->getPosition()[d];
        }
    }
    potentialEnergy *= 0.5;

    return kineticEnergy + potentialEnergy;

}

std::vector<double> HarmonicOscillator::omega() const
{
    return m_omega;
}

void HarmonicOscillator::setOmega(const std::vector<double> &omega)
{
    m_omega = omega;
}

