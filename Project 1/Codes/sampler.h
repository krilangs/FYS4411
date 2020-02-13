#pragma once

class Sampler {
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void printOutputToTerminal();
    void computeAverages();
    double getEnergy()          { return m_energy; }
    void setEnergy(double energy);

    int getAcceptedNumber() const;
    void setAcceptedNumber(int acceptedNumber);

    void setStepNumber(int stepNumber);
    int getStepNumber() const;

private:
    int     m_numberOfMetropolisSteps = 0;
    int     m_stepNumber = 0;
    double  m_acceptedNumber          = 0;
    double  m_energy = 0;
    double  m_cumulativeEnergy = 0;
    //double  m_energySquared            = 0;
    //double  m_cumulativeEnergySquared  = 0;

    class System* m_system = nullptr;
};