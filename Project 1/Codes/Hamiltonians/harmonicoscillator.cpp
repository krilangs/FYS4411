#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

HarmonicOscillator::HarmonicOscillator(System* system, double omega) :
        Hamiltonian(system) {
    assert(omega > 0);
    m_omega.reserve(1);
    m_omega.push_back(omega);
}

double HarmonicOscillator::computeLocalEnergy(std::vector<Particle*> particles) {
    /* Here, you need to compute the kinetic and potential energies. Note that
     * when using numerical differentiation, the computation of the kinetic
     * energy becomes the same for all Hamiltonians, and thus the code for
     * doing this should be moved up to the super-class, Hamiltonian.
     *
     * You may access the wave function currently used through the
     * getWaveFunction method in the m_system object in the super-class, i.e.
     * m_system->getWaveFunction()...
     */

    double potentialEnergy = 0;
    double kineticEnergy   = 0;

    kineticEnergy = -0.5*m_system->getWaveFunction()->computeDoubleDerivative(particles);   //analytical double derivative
    //kineticEnergy = -0.5*m_system->getHamiltonian()->computeNumericalDoubleDerivative(particles);   //numerical double derivative

    for (int k = 0; k < m_system->getNumberOfParticles(); k++ ){
            for (int d = 0; d < m_system->getNumberOfDimensions(); d++){
                potentialEnergy += m_omega[0]*m_omega[0]*particles.at(k)->getPosition()[d]*particles.at(k)->getPosition()[d];

            }

        }
    potentialEnergy *= 0.5;


    return kineticEnergy + potentialEnergy;


}

