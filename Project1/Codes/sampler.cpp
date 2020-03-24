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

using std::cout;
using std::endl;
std::ofstream ofile;

Sampler::Sampler(System* system) {
    m_system = system;
    m_stepNumber = 0;
}

void Sampler::sample(bool acceptedStep) {
    // Make sure the sampling variable(s) are initialized at the first step.
    if (m_stepNumber == 0) {
        m_cumulativeEnergy = 0;
        m_cumulativeEnergySquared = 0;
        m_cumulativeWFderiv = 0;
        m_cumulativeWFderivMultEloc = 0;
    }

    if (acceptedStep==true){
        m_acceptedNumber++;
    }
    // Here we sample the interesting things we want to measure.
    m_energy = m_system->getHamiltonian()->computeLocalEnergy(m_system->getParticles());
    m_WFderiv = 0;
    double beta = m_system->getWaveFunction()->getParameters()[2] / m_system->getWaveFunction()->getParameters()[0];

    for (int i=0; i<m_system->getNumberOfParticles(); i++){
        for (int d=0; d<m_system->getNumberOfDimensions()-1; d++){
            m_WFderiv -= m_system->getParticles().at(i)->getPosition()[d]*m_system->getParticles().at(i)->getPosition()[d];
        }
        int d = m_system->getNumberOfDimensions()-1;
        m_WFderiv -= m_system->getParticles().at(i)->getPosition()[d]*m_system->getParticles().at(i)->getPosition()[d]*beta;
    }
    m_WFderivMultELoc = m_WFderiv*m_energy;
    m_cumulativeEnergy          += m_energy;
    m_cumulativeEnergySquared   += m_energy*m_energy;
    m_cumulativeWFderiv         += m_WFderiv;
    m_cumulativeWFderivMultEloc += m_WFderivMultELoc;

    //m_system->oneBodyDensity();   // Uncomment to do one-body density

    m_stepNumber++;
}

void Sampler::printOutputToTerminal() {
    int     N     = m_system->getNumberOfParticles();
    int     d     = m_system->getNumberOfDimensions();
    int     MC    = m_system->getNumberOfMetropolisSteps();
    int     Np    = m_system->getWaveFunction()->getNumberOfParameters();
    double  equil = m_system->getEquilibrationFraction();
    double  dt    = m_system->getTimeStep();
    double  dMC   = m_system->getStepLength();
    double  var   = (m_cumulativeEnergySquared - m_energy*m_energy)/MC;
    double  std   = sqrt(fabs(m_cumulativeEnergySquared - m_energy*m_energy))/sqrt(MC);
    double  A     =  m_acceptedNumber/m_stepNumber;
    std::vector<double> param = m_system->getWaveFunction()->getParameters();
    ofile.close();

    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << N << endl;
    cout << " Number of dimensions : " << d << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(MC) << endl;
    cout << " Number of equilibration steps  : 10^" << std::log10(std::round(MC*equil)) << endl;
    cout << " Monte Carlo step length : " << dMC << endl;
    cout << " Importance sampling time step : " << dt << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << " Number of parameters : " << Np << endl;
    for (int i=0; i<Np; i++) {
        cout << " Parameter " << i+1 << " : " << param.at(i) << endl;
    }
    cout << endl;
    cout << "  -- Results -- " << endl;
    cout << " Energy : " << m_energy << endl;
    cout << " Variance: " << var << endl;
    cout << " St. dev (error): " << std << endl;
    cout << " Acceptance ratio: " << A << endl;
    cout << endl;
}

void Sampler::computeAverages() {
// Compute the averages of the sampled energies.
    m_energy = m_cumulativeEnergy / (m_system->getNumberOfMetropolisSteps());
    m_cumulativeEnergySquared /= m_system->getNumberOfMetropolisSteps();
    m_cumulativeWFderiv /= m_system->getNumberOfMetropolisSteps();
    m_cumulativeWFderivMultEloc /= m_system->getNumberOfMetropolisSteps();
}

void Sampler::openDataFile(std::string filename){
    if (filename != "0") ofile.open(filename);
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

void Sampler::setCumulativeWF(double cumulativeWF)
{
    m_cumulativeWF = cumulativeWF;
}

double Sampler::getCumulativeWF() const
{
    return m_cumulativeWF;
}

void Sampler::setCumulativeWFderiv(double cumulativeWFderiv)
{
    m_cumulativeWFderiv = cumulativeWFderiv;
}

double Sampler::getCumulativeWFderiv() const
{
    return m_cumulativeWFderiv;
}

void Sampler::setCumulativeWFderivMultEloc(double cumulativeWFderivMultEloc)
{
    m_cumulativeWFderivMultEloc = cumulativeWFderivMultEloc;
}

double Sampler::getCumulativeWFderivMultEloc() const
{
    return m_cumulativeWFderivMultEloc;
}

void Sampler::setWFderiv(double WFderiv)
{
    m_WFderiv = WFderiv;
}

double Sampler::getWFderiv() const
{
    return m_WFderiv;
}

void Sampler::setWFderivMultELoc(double WFderivMultELoc)
{
    m_WFderivMultELoc = WFderivMultELoc;
}

double Sampler::getWFderivMultELoc() const
{
    return m_WFderivMultELoc;
}
