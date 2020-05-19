#include "system.h"
#include <cassert>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "sampler.h"
#include "particle.h"
#include "Math/random.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"

using namespace std;
ofstream myFile;

bool System::metropolisStep(double GibbsValue, vector<double> &X, vector<double> H, vector<double> a, vector<double> b,
                            vector<vector<double>> w) {
    /* Perform the brute force Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */
    int randparticle = Random::nextInt(getNumberOfVisibleNodes());
    vector<double> X_old(getNumberOfVisibleNodes());
    vector<double> X_new(getNumberOfVisibleNodes());

    X_old = X;
    X_new = X;

    vector <vector<double>> OldDistanceMatrix;
    OldDistanceMatrix = getDistanceMatrix();

    // Choose new move
    X_new[randparticle] = X_old[randparticle] + m_stepLength*(Random::nextDouble() - 0.5);

    setDistanceMatrix(computematrixdistance(X_new));
    double psi_new = m_waveFunction->evaluate(GibbsValue, X_new, H, a, b, w);
    double prob = psi_new*psi_new/(m_psiOld*m_psiOld);

    if (Random::nextDouble()<prob || 1.0<prob){// Accept new move, non-interaction
        m_psiOld = psi_new;
        X = X_new;
        return true;
    }
    else{//Don't accept new move
        X = X_old;
        setDistanceMatrix(OldDistanceMatrix);
        return false;
    }
}

bool System::metropolisStepImportance(double GibbsValue, vector<double> &X, vector<double> H, vector<double> a, vector<double> b,
                                      vector<vector<double>> w){
    // Perform the Importance sampling.
    int randparticle = Random::nextInt(getNumberOfVisibleNodes());
    vector <double> X_old(getNumberOfVisibleNodes());
    vector <double> X_new(getNumberOfVisibleNodes());

    X_old = X;
    X_new = X;

    vector<double> QF_old(m_numberOfVisibleNodes);
    vector<double> QF_new(m_numberOfVisibleNodes);
    QF_old = getQuantumForce();

    vector<vector<double>> OldDistanceMatrix;
    OldDistanceMatrix = getDistanceMatrix();

    // Choose new move
    X_new[randparticle] = X_old[randparticle] + 0.5*m_QuantumForce[randparticle]*m_timeStep + m_sqrtTimeStep*(Random::nextGaussian(0.0,1.0));

    // Update
    setQuantumForce(m_waveFunction->QuantumForce(GibbsValue, X_new, a, b, w));
    QF_new = m_QuantumForce;
    updateDistanceMatrix(X_new, randparticle);

    // Green function
    double GreensFunction = 0.0;
    GreensFunction = 0.5*(QF_old[randparticle] + QF_new[randparticle])
                    *(0.25*m_timeStep*(QF_old[randparticle] - QF_new[randparticle]) - X_new[randparticle] + X_old[randparticle]);
    GreensFunction = exp(GreensFunction);

    double psi_new = m_waveFunction->evaluate(GibbsValue, X_new, H, a, b, w);
    double prob = GreensFunction*psi_new*psi_new/(m_psiOld*m_psiOld);

    if ((Random::nextDouble()<=prob) || (1.0<prob)){//Accept new move
        m_psiOld = psi_new;
        X = X_new;
        return true;
    }

    else{// Don't accept new move
        X = X_old;
        setQuantumForce(QF_old);
        updateDistanceMatrix(X_old, randparticle);
        return false;
    }
}

bool System::GibbsSampling(vector<double> &X, vector<double> H, vector<double> a, vector<double> b, vector<vector<double> > w){
  //NOT WORKING WITH INTERACTION!!!!!
    // Perform Gibbs Sampling.
    int N = getNumberOfHiddenNodes();
    int M = getNumberOfVisibleNodes();

    double sum, sum2, mean, argument;

    for (int j=0; j<N; j++){
        sum = 0;
        for (int i=0; i<M; i++){
            sum += X[i]*w[i][j]/getSigma_squared();
        }

        argument = exp(-b[j] - sum);
        H[j] = 1./(1 + argument);
    }

    for (int i=0; i<M; i++){
        sum2 = 0;
        for (int j=0; j<N; j++){
            sum2 += H[j]*w[i][j];
        }

        mean = a[i] + sum2;
        X[i] = Random::nextGaussian(mean, getSigma());
    }

    return true;
}

void System::runMetropolisSteps(string method, int numberOfMetropolisSteps, vector<double> &G, bool interaction, vector<double> X,
                                vector<double> H, vector<double> a, vector<double> b, vector<vector<double>> w) {
    double GibbsValue = 1.0;
    if (method == "Gibbs") GibbsValue = 0.5;
    m_sampler                 = new Sampler(this);
    m_numberOfMetropolisSteps = numberOfMetropolisSteps;

    getSampler()->setStepNumber(0);
    getSampler()->setAcceptedNumber(0);
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);
    m_sampler->setDimensionOfGradient(m_numberOfParameters);

    // Initial values
    setDistanceMatrix(computematrixdistance(X));
    m_psiOld = m_waveFunction->evaluate(GibbsValue, X, H, a, b, w);
    setQuantumForce(m_waveFunction->QuantumForce(GibbsValue, X, a, b, w));

    vector<double> temp (m_numberOfParameters);
    setCumulativeGradient(temp);
    setCumulativeEnGradient(temp);
    setGradient(temp);
    setEnGradient_average(temp);

    bool acceptedStep;
    // Sample after equilibrium have been reached
    for (int i=0; i<numberOfMetropolisSteps; i++) {
        if (method == "MetropolisBruteForce") acceptedStep = metropolisStep(GibbsValue, X, H, a, b, w);        // Brute force Metropolis
        if (method == "MetropolisImportance") acceptedStep = metropolisStepImportance(GibbsValue, X, H, a, b, w);   // Importance sampling
        if (method == "Gibbs") acceptedStep = GibbsSampling(X, H, a, b, w);       // Gibbs Sampling

        // Sample the energy
        m_sampler->sample(acceptedStep, interaction, GibbsValue, X, H, a, b, w);
        m_sampler->writeToFile();
    }
    m_sampler->computeAverages(G);
    //m_sampler->printOutputToTerminal();
}

void System::updateDistanceMatrix(vector<double> m_X, int randparticle)
{
    double temp = 0;
    double init = computeIndex(randparticle);

    double part = init/m_numberOfDimensions;

    for (int j=0; j<init; j+=m_numberOfDimensions){
        temp = 0;
        for (int d=0; d<m_numberOfDimensions; d++){
            temp += (m_X[j+d] - m_X[init+d])*(m_X[j+d] - m_X[init+d]);
        }
        m_distanceMatrix[part][j/m_numberOfDimensions] = sqrt(temp);
        m_distanceMatrix[j/m_numberOfDimensions][part] = m_distanceMatrix[part][j/m_numberOfDimensions];
    }
    for (int j=init+m_numberOfDimensions; j<m_numberOfVisibleNodes; j+=m_numberOfDimensions){
        temp = 0;
        for (int d=0; d<m_numberOfDimensions; d++){
            temp += (m_X[j+d] - m_X[init+d])*(m_X[j+d] - m_X[init+d]);
        }
        m_distanceMatrix[part][j/m_numberOfDimensions] = sqrt(temp);
        m_distanceMatrix[j/m_numberOfDimensions][part] = m_distanceMatrix[part][j/m_numberOfDimensions];
    }
}

vector<vector<double>> System::computematrixdistance(vector<double>& m_X)
{
    vector<vector<double>> distancematrix(m_numberOfParticles, vector<double>(m_numberOfParticles));
    double temp = 0;
    int j = 0;
    int z;
    int k = 0;

    while (j < m_numberOfVisibleNodes){
        temp = 0;
        z = 0;

        for (int i=0; i<j; i+=m_numberOfDimensions){
            for (int q=0; q<m_numberOfDimensions; q++){
                temp += (m_X[i+q] - m_X[j+q])*(m_X[i+q] - m_X[j+q]);
            }

            distancematrix[z][k] = sqrt(temp);
            distancematrix[k][z] = distancematrix[z][k];
            z++;
        }
        j += m_numberOfDimensions;
        k++;
    }

    return distancematrix;
}

void System::gradientDescent(vector<double> G, vector<double> X, vector<double> &a, vector<double> &b, vector<vector<double>> &w)
{//Gradient descent method to find the optimal variational parameter alpha given an initial parameter initialAlpha
    for (int i=0; i<m_numberOfVisibleNodes; i++){
        a[i] -= m_learningRate*G[i];
    }

    for (int i=0; i<m_numberOfHiddenNodes; i++){
        b[i] -= m_learningRate*G[i + m_numberOfVisibleNodes];
    }

    int z = m_numberOfVisibleNodes + m_numberOfHiddenNodes;

    for (int i=0; i<m_numberOfVisibleNodes; i++){
        for (int j=0; j<m_numberOfHiddenNodes; j++){
            w[i][j] -= m_learningRate*G[z];
            z++;
        }
    }
}

vector<double> System::GradientParameters(double GibbsValue, vector<double> X, vector<double> &a, vector<double> &b,
                                          vector<vector<double>> &w)
{
    vector<double> GradientPsi (m_numberOfParameters);
    vector<double> argument    (m_numberOfHiddenNodes);

    double sum;

    for (int j=0; j<m_numberOfHiddenNodes; j++){
        sum = 0;
        for (int i=0; i<m_numberOfVisibleNodes; i++){
            sum += X[i]*w[i][j]/getSigma_squared();
        }
        argument[j] = b[j] + sum;
    }

    for (int i=0; i<m_numberOfVisibleNodes; i++){
        GradientPsi[i] = (X[i] - a[i])*GibbsValue/getSigma_squared();
    }

    for (int i=m_numberOfVisibleNodes; i<m_numberOfHiddenNodes+m_numberOfVisibleNodes; i++){
        GradientPsi[i] = GibbsValue/(1 + exp(-argument[i-m_numberOfVisibleNodes]));
    }

    int z = m_numberOfHiddenNodes + m_numberOfVisibleNodes;
    for (int i=0; i<m_numberOfVisibleNodes; i++){
        for (int j=0; j<m_numberOfHiddenNodes; j++){
            GradientPsi[z] = GibbsValue*X[i]/(getSigma_squared()*(1 + exp(-argument[j])));
            z++;
        }
    }

    return GradientPsi;
}

double System::findEnergyDerivative()
{
    double meanEnergy      = getSampler()->getCumulativeEnergy()/(m_numberOfMetropolisSteps*getEquilibrationFraction());
    double meanWFderiv     = getSampler()->getCumulativeWFderiv()/(m_numberOfMetropolisSteps*getEquilibrationFraction());
    double meanWFderivEloc = getSampler()->getCumulativeWFderivMultEloc()/(m_numberOfMetropolisSteps*getEquilibrationFraction());

    return 2*(meanWFderivEloc - meanEnergy*meanWFderiv);
}

int System::computeIndex(int index)
{
    int init = index;
    if ((index%m_numberOfDimensions) == 0){
        return index;
    }

    while (index > m_numberOfDimensions){
        index -= m_numberOfDimensions;
    }

    init -= index;

    return init;
}

void System::openDataFile(string filename)
{
    m_sampler->openDataFile(filename);
}

void System::printOut(int cycle)
{
    m_sampler->printOutputToTerminal(cycle);
}

void System::openFile(string filename)
{
    myFile.open(filename);
    myFile << setprecision(12)<<fixed;
    myFile << setw(5)<<fixed;

    myFile << "Energy  " <<  "St. dev  ";
    for( int i = 0; i < m_numberOfVisibleNodes; i++){
        myFile<<"a[" << i << "] " ;
    }

    for(int i = 0; i < m_numberOfHiddenNodes; i++){
        myFile<<"b[" << i << "] " ;
    }

    int z = m_numberOfVisibleNodes + m_numberOfHiddenNodes;
    for (int i = 0; i < m_numberOfVisibleNodes; i++){
        for (int j = 0; j < m_numberOfHiddenNodes; j++){
            z++;
            myFile<<"w[" << i << "]["<< j << "] " ;
        }
    }
    myFile << endl;
}

void System::writeToFile(vector<double> X, vector<double>& a, vector<double>& b, vector<vector<double>>& w)
{
    double energy = getSampler()->getEnergy();
    int MC = getNumberOfMetropolisSteps();
    double ef = getEquilibrationFraction();
    double MC_eq = MC - ef*MC;
    myFile  << energy << "  " << sqrt(getSampler()->getCumulativeEnergySquared() - energy*energy)/sqrt(MC_eq)  << "  ";

    for( int i = 0; i < m_numberOfVisibleNodes; i++){
        myFile << a[i]<<"  " ;
    }

    for(int i = 0; i < m_numberOfHiddenNodes; i++){
        myFile << b[i]<<"  " ;
    }

    for (int i = 0; i < m_numberOfVisibleNodes; i++){
        for (int j = 0; j < m_numberOfHiddenNodes; j++){
            myFile<< w[i][j]<<"  " ;
        }
    }
    myFile << endl;
}

// Under follows all the set and get functions:
void System::setNumberOfParticles(int numberOfParticles)
{
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions)
{
    m_numberOfDimensions = numberOfDimensions;
}

void System::setStepLength(double stepLength)
{
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

double System::getStepLength() const
{
    return m_stepLength;
}

void System::setEquilibrationFraction(double equilibrationFraction)
{
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setHamiltonian(Hamiltonian* hamiltonian)
{
    m_hamiltonian = hamiltonian;
}

void System::setWaveFunction(WaveFunction* waveFunction)
{
    m_waveFunction = waveFunction;
}

void System::setInitialState(InitialState* initialState)
{
    m_initialState = initialState;
}

void System::setTimeStep(double timeStep)
{
    m_timeStep = timeStep;
}

double System::getTimeStep() const
{
    return m_timeStep;
}

void System::setSqrtTimeStep(double sqrtTimeStep)
{
    m_sqrtTimeStep = sqrtTimeStep;
}

double System::getSqrtTimeStep() const
{
    return m_sqrtTimeStep;
}

void System::setPsiOld(double psiOld)
{
    m_psiOld = psiOld;
}

double System::getPsiOld() const
{
    return m_psiOld;
}

int System::getNumberOfParameters() const
{
    return m_numberOfParameters;
}

void System::setNumberOfParameters(int numberOfParameters)
{
    m_numberOfParameters = numberOfParameters;
}

vector<double> System::getCumulativeGradient() const
{
    return m_cumulativeGradient;
}

void System::setCumulativeGradient(const vector<double> &cumulativeGradient)
{
    m_cumulativeGradient = cumulativeGradient;
}

vector<double> System::getCumulativeEnGradient() const
{
    return m_cumulativeEnGradient;
}

void System::setCumulativeEnGradient(const vector<double> &cumulativeEnGradient)
{
    m_cumulativeEnGradient = cumulativeEnGradient;
}

vector<double> System::getGradient() const
{
    return m_Gradient;
}

void System::setGradient(const vector<double> &Gradient)
{
    m_Gradient = Gradient;
}

vector<double> System::getEnGradient_average() const
{
    return m_EnGradient_average;
}

void System::setEnGradient_average(const vector<double> &EnGradient_average)
{
    m_EnGradient_average = EnGradient_average;
}

int System::getNumberOfVisibleNodes() const
{
    return m_numberOfVisibleNodes;
}

void System::setNumberOfVisibleNodes(int numberOfVisibleNodes)
{
    m_numberOfVisibleNodes = numberOfVisibleNodes;
}

int System::getNumberOfHiddenNodes() const
{
    return m_numberOfHiddenNodes;
}

void System::setNumberOfHiddenNodes(int numberOfHiddenNodes)
{
    m_numberOfHiddenNodes = numberOfHiddenNodes;
}

double System::getSigma() const
{
    return m_sigma;
}

void System::setSigma(double sigma)
{
    m_sigma = sigma;
}

double System::getSigma_squared() const
{
    return m_sigma_squared;
}

void System::setSigma_squared(double sigma_squared)
{
    m_sigma_squared = sigma_squared;
}

double System::getLearningRate() const
{
    return m_learningRate;
}

void System::setLearningRate(double learningRate)
{
    m_learningRate = learningRate;
}

void System::setDistanceMatrix(const vector<vector<double> > &distanceMatrix)
{
    m_distanceMatrix = distanceMatrix;
}

vector<vector<double> > System::getDistanceMatrix() const
{
    return m_distanceMatrix;
}

double System::getDistanceMatrixij(int i, int j) const
{
    return m_distanceMatrix[i][j];
}

void System::setQuantumForce(const vector<double> &QuantumForce)
{
    m_QuantumForce = QuantumForce;
}

vector<double> System::getQuantumForce() const
{
    return m_QuantumForce;
}
