#include "system.h"
#include <cassert>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"
//#include <omp.h>

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

    if (Random::nextDouble()-0.5 <= psi_new*psi_new/(m_psiOld*m_psiOld)){// Accept new move
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

bool System::metropolisStepImportance(){
    int randparticle = Random::nextInt(m_numberOfParticles);
    vector <double> r_old = m_particles.at(randparticle)->getPosition();
    vector <double> r_new(m_numberOfDimensions);

    vector <vector<double>> QF_old(m_numberOfDimensions,vector<double>(m_numberOfParticles));
    vector <vector<double>> QF_new(m_numberOfDimensions,vector<double>(m_numberOfParticles));
    QF_old=m_QuantumForce;

    // Choose new move
    for (int d=0; d < m_numberOfDimensions; d++){
        r_new[d] = r_old[d] + 0.5*m_QuantumForce[d][randparticle]*m_timeStep + m_sqrtTimeStep*(Random::nextGaussian(0, 1));
    }

    m_particles.at(randparticle)->setPosition(r_new);
    setQuantumForce(m_waveFunction->QuantumForce(m_particles));
    QF_new = m_QuantumForce;
    updateDistanceMatrix(m_particles, randparticle);

    // Green function
    double GreensFunction = 0.;
    for (int d=0; d<m_numberOfDimensions; d++){
        int j = randparticle;
        GreensFunction += 0.5*(QF_old[d][j] + QF_new[d][j])*(0.5*0.5*m_timeStep*(QF_old[d][j] - QF_new[d][j]) - r_new[d] + r_old[d]);
    }
    GreensFunction = exp(GreensFunction);
    double psi_new = m_waveFunction->evaluate(m_particles);

    //Accept new move
    if (Random::nextDouble() <= GreensFunction*psi_new*psi_new/(m_psiOld*m_psiOld)){
        m_psiOld = psi_new;
        getSampler()->setEnergy(getHamiltonian()->computeLocalEnergy(getParticles()));
        return true;
    }
    // Don't accept new move
    else{
        m_particles.at(randparticle)->setPosition(r_old);
        setQuantumForce(QF_old);
        updateDistanceMatrix(m_particles, randparticle);
        return false;
    }
}

void System::runMetropolisSteps(int numberOfMetropolisSteps) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    //m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    getSampler()->setStepNumber(0);
    getSampler()->setAcceptedNumber(0);
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    setDistanceMatrix(computematrixdistance(m_particles));
    m_psiOld = m_waveFunction->evaluate(m_particles);
    getSampler()->setEnergy(getHamiltonian()->computeLocalEnergy(m_particles));
    setQuantumForce(m_waveFunction->QuantumForce(m_particles));

    setHistogram();

    //#pragma omp parallel for ordered schedule(dynamic)
    for (int i=0; i < numberOfMetropolisSteps; i++) {
        bool acceptedStep = metropolisStep();
        //bool acceptedStep = metropolisStepImportance();

        /* Here you should sample the energy (and maybe other things using
         * the m_sampler instance of the Sampler class. Make sure, though,
         * to only begin sampling after you have let the system equilibrate
         * for a while. You may handle this using the fraction of steps which
         * are equilibration steps; m_equilibrationFraction.
         */
        //#pragma omp ordered
        m_sampler->sample(acceptedStep);
        m_sampler->writeToFile();
    }
    m_sampler->computeAverages();
    m_sampler->printOutputToTerminal();
}

double System::gradientDescent(double initialAlpha, string filename, int maxIterations){
//Gradient descent method to find the optimal variational parameter alpha given an initial parameter initialAlpha
    int steepestDescentSteps = int (1e+4);
    double alpha = initialAlpha;
    double beta = getWaveFunction()->getParameters()[2] / getWaveFunction()->getParameters()[0];
    double lambda = -0.001;
    int iterations = 0;
    double energyDerivative = 100;
    double cumulativeAlpha = 0;
    double tol = 1e-10;
    double percentAlphasToSave = 0.3;
    ofstream myFile;
    myFile.open(filename);

    while (iterations < maxIterations && fabs(energyDerivative) > tol){
        vector<double> parameters(3);
        parameters[0] = alpha;
        parameters[1] = alpha;
        parameters[2] = alpha*beta;
        getWaveFunction()->setParameters(parameters);
        runMetropolisSteps(steepestDescentSteps);
        //m_sampler->computeAverages();
        //m_sampler->printOutputToTerminal();
        energyDerivative = findEnergyDerivative();

        // Make sure we accept enough moves (with interaction can get stuck)
        if (double(m_sampler->getAcceptedNumber()) / steepestDescentSteps > 0.64){ //0.64 for brute MC, 0.97 for importance
            alpha += lambda*energyDerivative;
            iterations++;
        }

        cout << " New alpha = "  << alpha <<  endl;
        cout << " Energy derivative = " << energyDerivative << endl;
        cout << " Iterations = " << iterations << endl;

        // Write alpha, mean local energy and st dev to file
        myFile << alpha << "   "  << getSampler()->getEnergy() << "  " <<
                  sqrt(fabs(getSampler()->getCumulativeEnergySquared() - getSampler()->getEnergy()*getSampler()->getEnergy()))/getNumberOfMetropolisSteps() << endl;

        if (double (iterations) / maxIterations > 1-percentAlphasToSave){
            cumulativeAlpha += alpha;
        }
    }
    myFile.close();

    alpha = cumulativeAlpha / (maxIterations*percentAlphasToSave);
    return alpha;
}

void System::oneBodyDensity(){
    //Function to make the histrograms needed to compute the one body density
    vector<int> histogram(m_bins);
    double r2 = 0;
    for (int j = 0; j < getNumberOfParticles(); j++){
        r2 = 0;
        for (int d = 0; d < getNumberOfDimensions(); d++){
            r2 += m_particles.at(j)->getPosition()[d]*m_particles.at(j)->getPosition()[d];
        }
        r2 = sqrt(r2);
        int bucket = (int)floor(r2 / m_bucketSize);
        histogram[bucket] += 1;
    }

    // Update histogram
    for (int k = 0; k < m_bins; k++){
       m_histogram[k] += histogram[k];
    }
}

double System::findEnergyDerivative()
{
    double meanEnergy      = getSampler()->getCumulativeEnergy() / (m_numberOfMetropolisSteps*getEquilibrationFraction());
    double meanWFderiv     = getSampler()->getCumulativeWFderiv() / (m_numberOfMetropolisSteps*getEquilibrationFraction());
    double meanWFderivEloc =  getSampler()->getCumulativeWFderivMultEloc() / (m_numberOfMetropolisSteps*getEquilibrationFraction());
    //printf("E: %d\n", meanEnergy);
    //printf("WF: %d\n", getSampler()->getCumulativeWFderiv());
    //printf("WFD: %d\n", getSampler()->getCumulativeWFderivMultEloc());
    return 2*(meanWFderivEloc - meanEnergy*meanWFderiv);
}

void System::printOneBodyDensity(string filename){
    ofstream myFile;
    myFile.open(filename);
    for (int j = 0; j < m_bins; j++){
        //printf("%d", m_histogram[j]);
        myFile <<(double) m_histogram[j] /(getNumberOfParticles()*getNumberOfMetropolisSteps()*getEquilibrationFraction()*m_bucketSize) << "    " << j*m_bucketSize<< endl;
    }
    cout << "Printed ob density! " << endl;
}

void System::setHistogram()
{
    vector<int> histogram(getBins());
    m_histogram = histogram;
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
        if (m_distanceMatrix[randparticle][j] < getinteractionSize()){
             return true;
        }
        m_distanceMatrix[j][randparticle] = m_distanceMatrix[randparticle][j];
    }
    for (int j = randparticle+1; j < m_numberOfParticles; j++){
        temp = 0;
        for (int d = 0; d < m_numberOfDimensions; d++){
            temp += (particles.at(randparticle)->getPosition()[d] - particles.at(j)->getPosition()[d]) *
                    (particles.at(randparticle)->getPosition()[d] - particles.at(j)->getPosition()[d]);
        }
        m_distanceMatrix[randparticle][j] = sqrt(temp);
        if (m_distanceMatrix[randparticle][j] < getinteractionSize()){
            return true;
        }
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
            for(int k=0; k<m_numberOfDimensions; k++){
                temp += (particles.at(i)->getPosition()[k] - particles.at(j)->getPosition()[k]) *
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

void System::updateQuantumForce(std::vector<std::vector<double> > deltaQuantumForce, bool subtract){
    if (subtract == false){
        for (int d = 0; d < m_numberOfDimensions; d++){
        for(int i=0; i<=m_numberOfParticles; i++){
            m_QuantumForce[d][i] += deltaQuantumForce[d][i];
        }
        }
    }
    else {
        for (int d = 0; d < m_numberOfDimensions; d++){
        for(int i=0; i<=m_numberOfParticles; i++){
            m_QuantumForce[d][i] += deltaQuantumForce[d][i];
        }
        }
    }
}

std::vector<vector<double>> System::getQuantumForce() const
{
    return m_QuantumForce;
}

void System::setQuantumForce(const std::vector<vector<double>> &QuantumForce)
{
    m_QuantumForce = QuantumForce;
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

double System::getTimeStep() const
{
    return m_timeStep;
}

void System::setTimeStep(double timeStep)
{
    m_timeStep = timeStep;
}

double System::getSqrtTimeStep() const
{
    return m_sqrtTimeStep;
}

void System::setSqrtTimeStep(double sqrtTimeStep)
{
    m_sqrtTimeStep = sqrtTimeStep;
}

double System::getinteractionSize() const
{
    return m_interactionSize;
}

void System::setinteractionSize(double interactionSize)
{
    m_interactionSize = interactionSize;
}

int System::getBins() const
{
    return m_bins;
}

void System::setBins(int bins)
{
    m_bins = bins;
}

double System::getBucketSize() const
{
    return m_bucketSize;
}

void System::setBucketSize(double bucketSize)
{
    m_bucketSize = bucketSize;
}

void System::openDataFile(string filename){
    m_sampler->openDataFile(filename);
}
