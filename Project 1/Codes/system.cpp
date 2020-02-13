#include "system.h"
#include <cassert>
#include <cmath>
#include <vector>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"

using namespace std;

bool System::metropolisStep() {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */
    int randparticle = Random::nextInt(m_numberOfParticles);

    vector <double> r_old = m_particles.at(randparticle)->getPosition();
    vector <double> r_new(m_numberOfDimensions);

    // Choose new move
    for (int d=0; d < m_numberOfDimensions; d++){
        r_new[d] = r_old[d] + m_stepLength*(Random::nextDouble()-0.5);
    }

    m_particles.at(randparticle)->setPosition(r_new);
    updateDistanceMatrix(m_particles, randparticle);

    double psi_new = m_waveFunction->evaluate(m_particles);

    if (Random::nextDouble() <= psi_new*psi_new/(m_psiOld*m_psiOld)){// Accept new move
        m_psiOld = psi_new;

        getSampler()->setEnergy(getHamiltonian()->computeLocalEnergy(getParticles()));
        return true;
    }
    else{//Don't accept new move
        m_particles.at(randparticle)->setPosition(r_old);
        updateDistanceMatrix(m_particles, randparticle);
        return false;
    }
}

void System::runMetropolisSteps(int numberOfMetropolisSteps) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    getSampler()->setStepNumber(0);
    getSampler()->setAcceptedNumber(0);
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    setDistanceMatrix(computematrixdistance(m_particles));
    m_psiOld = m_waveFunction->evaluate(m_particles);
    getSampler()->setEnergy(getHamiltonian()->computeLocalEnergy(m_particles));

    for (int i=0; i < numberOfMetropolisSteps; i++) {
        bool acceptedStep = metropolisStep();

        /* Here you should sample the energy (and maybe other things using
         * the m_sampler instance of the Sampler class. Make sure, though,
         * to only begin sampling after you have let the system equilibrate
         * for a while. You may handle this using the fraction of steps which
         * are equilibration steps; m_equilibrationFraction.
         */
        m_sampler->sample(acceptedStep);
    }
    m_sampler->computeAverages();
    m_sampler->printOutputToTerminal();
}

bool System::updateDistanceMatrix( std::vector<class Particle*> particles, int randparticle){
    double temp = 0;
    for (int j = 0; j < randparticle; j++){
        temp = 0;
        for (int d = 0; d < m_numberOfDimensions; d++){
            temp += (particles.at(randparticle)->getPosition()[d] - particles.at(j)->getPosition()[d]) *
                    (particles.at(randparticle)->getPosition()[d] - particles.at(j)->getPosition()[d]);
        }
        m_distanceMatrix[randparticle][j] = sqrt(temp);

        m_distanceMatrix[j][randparticle] = m_distanceMatrix[randparticle][j];
    }
    for (int j = randparticle+1; j < m_numberOfParticles; j++){
        temp = 0;
        for (int d = 0; d < m_numberOfDimensions; d++){
            temp += (particles.at(randparticle)->getPosition()[d] - particles.at(j)->getPosition()[d]) *
                    (particles.at(randparticle)->getPosition()[d] - particles.at(j)->getPosition()[d]);
        }
        m_distanceMatrix[randparticle][j] = sqrt(temp);

        m_distanceMatrix[j][randparticle] = m_distanceMatrix[randparticle][j];

    }
    return false;
}

std::vector<vector<double>> System::computematrixdistance(std::vector<class Particle*> particles){

    vector<vector<double>> distancematrix(m_numberOfParticles, vector<double>(m_numberOfParticles));
    double temp=0;
    int j=0;
    while(j < m_numberOfParticles){
        temp = 0;
        for(int i = 0; i < j; i++){

            for(int k=0;k<m_numberOfDimensions;k++){
                temp+=(particles.at(i)->getPosition()[k] - particles.at(j)->getPosition()[k]) *
                      (particles.at(i)->getPosition()[k] - particles.at(j)->getPosition()[k]);
            }
            distancematrix[i][j]=sqrt(temp);
            distancematrix[j][i]=distancematrix[i][j];
        }

        j++;
    }

    return distancematrix;
}

void System::setDistanceMatrix(const std::vector<vector<double> > &distanceMatrix)
{
    m_distanceMatrix = distanceMatrix;
}

std::vector<vector<double> > System::getDistanceMatrix() const
{
    return m_distanceMatrix;
}

double System::getDistanceMatrixij(int i, int j) const
{
    return m_distanceMatrix[i][j];
}

void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
}

double System::getPsiOld() const
{
    return m_psiOld;
}

void System::setPsiOld(double psiOld)
{
    m_psiOld = psiOld;
}
