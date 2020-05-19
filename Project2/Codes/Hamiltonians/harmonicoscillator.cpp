#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using namespace std;

HarmonicOscillator::HarmonicOscillator(System* system, double omega):
                                       Hamiltonian(system) {
    assert(omega > 0);
    m_omega = omega;
}

double HarmonicOscillator::computeLocalEnergy(bool interaction, double GibbsValue, vector<double> X, vector<double> H, vector<double> a,
                                              vector<double> b, vector<vector<double>> w) {
    // Compute the local energy from the kinetic and potential energies.
    double potentialEnergy   = 0;
    double kineticEnergy     = 0;
    double interactPotential = 0;

    int dim = m_system->getNumberOfDimensions();
    int M   = m_system->getNumberOfVisibleNodes();

    kineticEnergy = -0.5*m_system->getWaveFunction()->computeDoubleDerivative(GibbsValue, X, H, a, b, w);  // From analytical double derivative

    for (int k=0; k<M; k+=dim){
        for (int d=0; d<dim; d++){
            potentialEnergy += m_omega*m_omega*X[k+d]*X[k+d];
        }
    }
    potentialEnergy *= 0.5;

    if (interaction == true){
        for (int j=0; j<m_system->getNumberOfParticles(); j++){
            for (int i=0; i<j; i++){
                interactPotential += 1/m_system->getDistanceMatrixij(i,j);
            }
        }
    }

    return kineticEnergy + potentialEnergy + interactPotential;

}

double HarmonicOscillator::omega() const
{
    return m_omega;
}

void HarmonicOscillator::setOmega(const double &omega)
{
    m_omega = omega;
}

