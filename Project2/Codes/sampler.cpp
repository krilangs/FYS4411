#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include "sampler.h"
#include "system.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"

using namespace std;
using std::cout;
using std::endl;
std::ofstream ofile;

Sampler::Sampler(System* system) {
    m_system = system;
    m_stepNumber = 0;
}

void Sampler::sample(bool acceptedStep, bool interaction, double GibbsValue, vector<double> X, vector<double> H, vector<double> a,
                     vector<double> b, vector<vector<double>> w) {
    // Make sure the sampling variable(s) are initialized at the first step.
    if (m_stepNumber == 0) {
        m_cumulativeEnergy = 0;
        m_cumulativeEnergySquared = 0;
    }

    m_energy = m_system->getHamiltonian()->computeLocalEnergy(interaction, GibbsValue, X, H, a, b, w);

    if ((double)getStepNumber() >= m_system->getEquilibrationFraction()*getNumberOfMetropolisSteps()){
        // Sample if the system is at equilibrium
        if (acceptedStep == true){
            m_acceptedNumber++;
        }

        m_cumulativeEnergy          += m_energy;
        m_cumulativeEnergySquared   += m_energy*m_energy;

        vector<double> temp(getDimensionOfGradient());
        vector<double> temp2(getDimensionOfGradient());
        vector<double> grad(getDimensionOfGradient());

        grad = m_system->GradientParameters(GibbsValue, X, a, b, w);
        temp = m_system->getCumulativeGradient();
        temp2 = m_system->getCumulativeEnGradient();

        for (int i=0; i<getDimensionOfGradient(); i++){
            temp[i]  += grad[i];
            temp2[i] += m_energy*grad[i];
        }

        m_system->setCumulativeGradient(temp);
        m_system->setCumulativeEnGradient(temp2);
    }

    m_stepNumber++;
}

void Sampler::printOutputToTerminal(int cycle) {
    int     N     = m_system->getNumberOfParticles();
    int     d     = m_system->getNumberOfDimensions();
    int     MC    = m_system->getNumberOfMetropolisSteps();
    int     Np    = m_system->getNumberOfParameters();
    int     M     = m_system->getNumberOfVisibleNodes();
    int     H     = m_system->getNumberOfHiddenNodes();
    double  equil = m_system->getEquilibrationFraction();
    double  dt    = m_system->getTimeStep();
    double  sig   = m_system->getSigma();
    double  MC_eq = MC - equil*MC;
    double  dMC   = m_system->getStepLength();
    double  var   = (m_cumulativeEnergySquared - m_energy*m_energy)/MC_eq;
    double  std   = sqrt(fabs(m_cumulativeEnergySquared - m_energy*m_energy))/sqrt(MC_eq);
    double  A     =  m_acceptedNumber/MC_eq;
    ofile.close();

    if (cycle == 0){
        cout << endl;
        cout << "  -- System info -- " << endl;
        cout << " Number of particles  : " << N << endl;
        cout << " Number of dimensions : " << d << endl;
        cout << " Number of visible nodes : " << M << endl;
        cout << " Number of hidden nodes : " << H << endl;
        cout << " Number of Metropolis steps run : 10^" << log10(MC) << endl;
        cout << " Number of equilibration steps  : 10^" << log10(round(MC*equil)) << endl;
        cout << " Monte Carlo step length : " << dMC << endl;
        cout << " Importance sampling time step : " << dt << endl;
        cout << " Gibbs sampling sigma : " << sig << endl;
        cout << endl;
        cout << "  -- Wave function parameters -- " << endl;
        cout << " Number of parameters : " << Np << endl;
    }
    cout << endl;
    cout << "  -- Results -- " << endl;
    cout << " Energy [a.u.] : " << m_energy << endl;
    cout << " Variance: " << var << endl;
    cout << " St. dev (error): " << std << endl;
    cout << " Acceptance ratio: " << A << endl;
    cout << " Number of cycle: " << cycle << endl;
    cout << endl;
}

void Sampler::computeAverages(vector<double> &G) {
// Compute the averages of the sampled energies.
    double frac = m_system->getNumberOfMetropolisSteps()*(1 - m_system->getEquilibrationFraction());

    m_energy = m_cumulativeEnergy/frac;
    m_cumulativeEnergySquared /= frac;

    vector<double> temp(getDimensionOfGradient());
    vector<double> temp2(getDimensionOfGradient());

    temp = m_system->getCumulativeGradient();
    temp2 = m_system->getCumulativeEnGradient();

    double sum1 = 0;
    double sum2 = 0;

    for (int i=0; i<m_system->getNumberOfParameters(); i++){
        sum1 = temp[i]/frac;
        sum2 = temp2[i]/frac;

        G[i] = 2*(sum2 - m_energy*sum1);
    }
}

void Sampler::openDataFile(string filename){
    if (filename != "0") ofile.open(filename);

    ofile << setprecision(12) << fixed;
    ofile << setw(5) << fixed;
}

void Sampler::writeToFile(){
    if (ofile.is_open()) ofile << m_energy << endl;
}

// Under follows all the set and get functions:
void Sampler::setNumberOfMetropolisSteps(int steps)
{
    m_numberOfMetropolisSteps = steps;
}

int Sampler::getNumberOfMetropolisSteps() const
{
    return m_numberOfMetropolisSteps;
}

void Sampler::setEnergy(double energy)
{
    m_energy = energy;
}

void Sampler::updateEnergy(double dE)
{
    m_energy += dE;
}

void Sampler::setStepNumber(int stepNumber)
{
    m_stepNumber = stepNumber;
}

int Sampler::getStepNumber() const
{
    return m_stepNumber;
}

void Sampler::setAcceptedNumber(int acceptedNumber)
{
    m_acceptedNumber = acceptedNumber;
}

int Sampler::getAcceptedNumber() const
{
    return int(m_acceptedNumber);
}

double Sampler::getCumulativeEnergy() const
{
    return m_cumulativeEnergy;
}

void Sampler::setCumulativeEnergySquared(double cumulativeEnergySquared)
{
    m_cumulativeEnergySquared = cumulativeEnergySquared;
}

double Sampler::getCumulativeEnergySquared() const
{
    return m_cumulativeEnergySquared;
}

int Sampler::getDimensionOfGradient() const
{
    return m_dimensionOfGradient;
}

void Sampler::setDimensionOfGradient(int dimensionOfGradient)
{
    m_dimensionOfGradient = dimensionOfGradient;
}
