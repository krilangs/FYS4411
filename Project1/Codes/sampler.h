#pragma once
#include <string>

class Sampler {
public:
    Sampler(class System* system);
    void sample                             (bool acceptedStep);
    void printOutputToTerminal              ();
    void openDataFile                       (std::string filename);
    void writeToFile                        ();
    void computeAverages                    ();
    void setNumberOfMetropolisSteps         (int steps);
    void setStepNumber                      (int stepNumber);
    void setAcceptedNumber                  (int acceptedNumber);
    void setEnergy                          (double energy);
    void updateEnergy                       (double dE);
    void setCumulativeEnergySquared         (double cumulativeEnergySquared);
    void setCumulativeWF                    (double cumulativeWF);
    void setCumulativeWFderiv               (double cumulativeWFderiv);
    void setCumulativeWFderivMultEloc       (double cumulativeWFderivMultEloc);
    void setWFderiv                         (double WFderiv);
    void setWFderivMultELoc                 (double WFderivMultELoc);

    int getNumberOfMetropolisSteps()        const;
    int getStepNumber()                     const;
    int getAcceptedNumber()                 const;

    double getEnergy()                      { return m_energy; }
    double getCumulativeEnergy()            const;
    double getCumulativeEnergySquared()     const;
    double getCumulativeWF()                const;
    double getCumulativeWFderiv()           const;
    double getCumulativeWFderivMultEloc()   const;
    double getWFderiv()                     const;
    double getWFderivMultELoc()             const;

private:
    int     m_numberOfMetropolisSteps   = 0;
    int     m_stepNumber                = 0;
    double  m_acceptedNumber            = 0;
    double  m_energy                    = 0;
    double  m_energySquared             = 0;
    double  m_cumulativeEnergy          = 0;
    double  m_cumulativeEnergySquared   = 0;
    double  m_cumulativeWF              = 0;
    double  m_cumulativeWFderiv         = 0;
    double  m_cumulativeWFderivMultEloc = 0;
    double  m_WFderiv                   = 0;
    double  m_WFderivMultELoc           = 0;

    std::string m_filename;
    class System* m_system = nullptr;
};
