#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include <algorithm>
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
    double alpha = m_parameters[0];
    int dim = m_system->getNumberOfDimensions();
    int N = m_system->getNumberOfParticles();

    for (int i=0; i < N; i++){
        for (int d=0; d < dim; d++){
            r_sqr += particles.at(i)->getPosition()[d]*particles.at(i)->getPosition()[d]*alpha;
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
    double alpha = m_parameters[0];
    int dim = m_system->getNumberOfDimensions();
    int N = m_system->getNumberOfParticles();

    // For a single particle
    for (int i=0; i < N; i++){
        for (int d=0; d < dim; d++){
            one += alpha*alpha*
                   particles.at(i)->getPosition()[d]*
                   particles.at(i)->getPosition()[d];
        }
    }
    one *= 4.0;
    one -= 2*(dim*alpha)*N;
    return one;
}
