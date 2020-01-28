#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

SimpleGaussian::SimpleGaussian(System* system, double alpha) :
        WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

double SimpleGaussian::evaluate(std::vector<class Particle*> particles) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i].getPosition()
     * function.
     *
     * For the actual expression, use exp(-alpha * r^2), with alpha being the
     * (only) variational parameter.
     */
    double r_sqr = 0;

    for (int i=0; i < m_system->getNumberOfParticles(); i++){
        for (int d=0; d < m_system->getNumberOfDimensions(); d++){
            r_sqr += particles.at(i)->getPosition()[d]*particles.at(i)->getPosition()[d]*m_parameters[d];
        }
    }

    return exp(-r_sqr);
}

double SimpleGaussian::computeDoubleDerivative(std::vector<class Particle*> particles) {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */
    double one = 0;
    double interaction = 0;

    // For a single particle
    for (int i=0; i < m_system->getNumberOfParticles(); i++){
        for (int d=0; d < m_system->getNumberOfDimensions(); d++){
            one += m_parameters[d]*m_parameters[d]*
                    particles.at(i)->getPosition()[d]*
                    particles.at(i)->getPosition()[d];
        }
    }
    one *= 4.0;
    one -= 2*((m_system->getNumberOfDimensions() - 1)*m_parameters[0] + m_parameters[2])
            *m_system->getNumberOfParticles();
    return one+interaction;
}
