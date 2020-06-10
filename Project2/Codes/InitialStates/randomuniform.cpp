#include "randomuniform.h"
#include <iostream>
#include <cassert>
#include "../Math/random.h"
#include "../particle.h"
#include "../system.h"
#include <cmath>
#include <random>

using std::cout;
using std::endl;

RandomUniform::RandomUniform(System* system, int numberOfDimensions, int numberOfParticles, int numberOfHiddenNodes, int numberOfVisibleNodes,
                             double sigma, vector<double> &m_X, vector<double> &m_H, vector<double> &m_a, vector<double> &m_b,
                             vector<vector<double>> &m_w, double timeStep, int numberOfParameters):
                             InitialState(system) {
    /* The Initial State class is in charge of everything to do with the
     * initialization of the system; this includes among others
     * determining the number of particles and the number of dimensions
     * used. To make sure everything works as intended,
     * this information is passed to the system here.
     */
    assert(numberOfDimensions > 0 && numberOfParticles > 0);
    assert(numberOfHiddenNodes > 0 && numberOfVisibleNodes > 0);

    m_system->setNumberOfDimensions(numberOfDimensions);
    m_system->setNumberOfParticles(numberOfParticles);
    m_system->setTimeStep(timeStep);
    m_system->setSqrtTimeStep(sqrt(timeStep));

    m_system->setSigma(sigma);
    m_system->setSigma_squared(sigma*sigma);
    m_system->setNumberOfHiddenNodes(numberOfHiddenNodes);
    m_system->setNumberOfVisibleNodes(numberOfVisibleNodes);

    setupInitialState(m_X, m_H, m_a, m_b, m_w);
}

void RandomUniform::setupInitialState(vector<double> &m_X, vector<double> &m_H, vector<double> &m_a, vector<double> &m_b,
                                      vector<vector<double>> &m_w) {
    double sigma_0 = 0.5;
    int M = m_system->getNumberOfVisibleNodes();
    int N = m_system->getNumberOfHiddenNodes();

    for (int i=0; i < M; i++) {
        m_X[i] = (Random::nextDouble() - 0.5);
    }

    for (int i=0; i<M; i++){
        m_a[i] = Random::nextGaussian(0, sigma_0);
    }

    for (int i=0; i<N; i++){
        m_b[i] = Random::nextGaussian(0, sigma_0);
    }

    for (int i=0; i<M; i++){
        for (int j=0;j<N; j++){
            m_w[i][j] = Random::nextGaussian(0, sigma_0);
        }
    }
}
