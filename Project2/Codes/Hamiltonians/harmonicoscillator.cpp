#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include <cmath>
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
    int P   = m_system->getNumberOfParticles();
    int M   = m_system->getNumberOfVisibleNodes();

    kineticEnergy = -0.5*m_system->getWaveFunction()->computeDoubleDerivative(GibbsValue, X, H, a, b, w);  // From analytical double derivative

    for (int k=0; k<M; k+=dim){
        for (int d=0; d<dim; d++){
            potentialEnergy += m_omega*m_omega*X[k+d]*X[k+d];
        }
    }
    potentialEnergy *= 0.5;

    if (interaction == true){
        double par_i = 0;
        double par_j = 0;

        for (int i=1; i<P; i++){
            for (int j=0; j<i; j++){
                double sum = 0;
                for (int d=0; d<dim; d++){
                    par_i = dim*i + d;
                    par_j = dim*j + d;
                    sum += (X[par_i] - X[par_j])*(X[par_i] - X[par_j]);
                }
                if (sqrt(sum) > 0.00005){
                    interactPotential += 1/sqrt(sum);
                }
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

