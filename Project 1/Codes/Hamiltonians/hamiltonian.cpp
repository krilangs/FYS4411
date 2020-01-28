#include "hamiltonian.h"
#include "../system.h"
#include "../WaveFunctions/simplegaussian.h"
#include "../particle.h"

Hamiltonian::Hamiltonian(System* system) {
    m_system = system;
}

double Hamiltonian::computeNumericalDD(std::vector<class Particle *> particles)
{
    // Function to compute the numerical double derivative of the wavefunction
    double h = 0.001;
    double h_sqr = h*h;
    double wf = 0;
    double backward = 0;
    double forward = 0;
    double present = 0;

    present = m_system->getWaveFunction()->evaluate(particles);
    std::vector<double> r(m_system->getNumberOfDimensions());

    for (int j=0; j < m_system->getNumberOfParticles(); j++){
        for (int d=0; d < m_system->getNumberOfDimensions(); d++){
            r[d] = particles.at(j)->getPosition()[d];
        }
        for (int d=0; d < m_system->getNumberOfDimensions(); d++){
            r[d] -= h;
            particles.at(j)->setPosition(r);
            backward += m_system->getWaveFunction()->evaluate(particles);

            r[d] += 2*h;
            particles.at(j)->setPosition(r);
            forward += m_system->getWaveFunction()->evaluate(particles);

            r[d] -= h;
            particles.at(j)->setPosition(r);
        }
    }
    wf = (backward + forward - 2*present*m_system->getNumberOfDimensions()*m_system->getNumberOfParticles())/h_sqr;
    return wf/(present);
}
