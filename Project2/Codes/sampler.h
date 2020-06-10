#pragma once
#include <string>
#include <vector>

using namespace std;

class Sampler {
public:
    Sampler(class System* system);
    void sample                             (bool acceptedStep, bool interaction, double GibbsValue, vector<double> X, vector<double> H,
                                             vector<double> a, vector<double> b, vector<vector<double>> w);
    void printOutputToTerminal              (int cycle);
    void openDataFile                       (string filename);
    void writeToFile                        ();
    void computeAverages                    (vector<double>& G);
    void setNumberOfMetropolisSteps         (int steps);
    void setStepNumber                      (int stepNumber);
    void setAcceptedNumber                  (int acceptedNumber);
    void setDimensionOfGradient             (int dimensionOfGradient);
    void setEnergy                          (double energy);
    void updateEnergy                       (double dE);
    void setCumulativeEnergySquared         (double cumulativeEnergySquared);

    int getNumberOfMetropolisSteps()        const;
    int getStepNumber()                     const;
    int getAcceptedNumber()                 const;
    int getDimensionOfGradient() const;

    double getEnergy()                      { return m_energy; }
    double getCumulativeEnergy()            const;
    double getCumulativeEnergySquared()     const;

private:
    int     m_numberOfMetropolisSteps   = 0;
    int     m_stepNumber                = 0;
    int     m_dimensionOfGradient       = 0;
    double  m_acceptedNumber            = 0;
    double  m_energy                    = 0;
    double  m_energySquared             = 0;
    double  m_cumulativeEnergy          = 0;
    double  m_cumulativeEnergySquared   = 0;

    string m_filename;
    class System* m_system = nullptr;
};
