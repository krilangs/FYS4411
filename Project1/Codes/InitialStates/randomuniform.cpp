#include "randomuniform.h"
#include <iostream>
#include <cassert>
#include "../Math/random.h"
#include "../particle.h"
#include "../system.h"
#include <cmath>

using std::cout;
using std::endl;

RandomUniform::RandomUniform(System* system, int numberOfDimensions, int numberOfParticles, double timeStep,
                             double interactionSize, int bins, double bucketSize):
                             InitialState(system) {
    /* The Initial State class is in charge of everything to do with the
     * initialization of the system; this includes determining the number of
     * particles and the number of dimensions used. To make sure everything
     * works as intended, this information is passed to the system here.
     */
    assert(numberOfDimensions > 0 && numberOfParticles > 0);
    m_numberOfDimensions = numberOfDimensions;
    m_numberOfParticles  = numberOfParticles;
    m_system->setNumberOfDimensions(numberOfDimensions);
    m_system->setNumberOfParticles(numberOfParticles);
    m_system->setTimeStep(timeStep);
    m_system->setSqrtTimeStep(sqrt(timeStep));
    m_system->setinteractionSize(interactionSize);
    m_system->setBins(bins);
    m_system->setBucketSize(bucketSize);
    setupInitialState();
}

void RandomUniform::setupInitialState() {
    for (int i=0; i < m_numberOfParticles; i++) {
        std::vector<double> position = std::vector<double>();

        for (int j=0; j < m_numberOfDimensions; j++) {
            /* Place the particles randomly according
             * to a uniform distribution.
             */
            position.push_back(Random::nextDouble()-0.5);
        }
        m_particles.push_back(new Particle());
        m_particles.at(i)->setNumberOfDimensions(m_numberOfDimensions);
        m_particles.at(i)->setPosition(position);
    }
}
