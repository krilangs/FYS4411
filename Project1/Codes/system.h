#pragma once
#include <vector>
#include <string>

using namespace std;

class System {
public:
    bool metropolisStep             ();
    bool metropolisStepImportance   ();

    void runMetropolisSteps         (int numberOfMetropolisSteps);
    void setNumberOfParticles       (int numberOfParticles);
    void setNumberOfDimensions      (int numberOfDimensions);
    void setStepLength              (double stepLength);
    void setEquilibrationFraction   (double equilibrationFraction);
    void setHamiltonian             (class Hamiltonian* hamiltonian);
    void setWaveFunction            (class WaveFunction* waveFunction);
    void setInitialState            (class InitialState* initialState);
    void openDataFile               (string filename);

    class WaveFunction*             getWaveFunction()   { return m_waveFunction; }
    class Hamiltonian*              getHamiltonian()    { return m_hamiltonian; }
    class Sampler*                  getSampler()        { return m_sampler; }
    vector<class Particle*>         getParticles()      { return m_particles; }
    vector<vector<double>>          computematrixdistance(vector<class Particle*> particles);
    bool updateDistanceMatrix       (vector<Particle*> particles, int randparticle);

    int getNumberOfParticles()          { return m_numberOfParticles; }
    int getNumberOfDimensions()         { return m_numberOfDimensions; }
    int getNumberOfMetropolisSteps()    { return m_numberOfMetropolisSteps; }
    int getBins()                       const;

    void setTimeStep                    (double timeStep);
    void setSqrtTimeStep                (double sqrtTimeStep);
    void setPsiOld                      (double psiOld);
    void setinteractionSize             (double interactionSize);
    void setDistanceMatrix              (const vector<vector<double> > &distanceMatrix);
    void setQuantumForce                (const vector<vector<double>> &QuantumForce);
    void updateQuantumForce             (vector<vector<double>> deltaQuantumForce, bool subtract);
    void setBucketSize                  (double bucketSize);
    void setBins                        (int bins);
    void setHistogram                   ();
    void oneBodyDensity                 ();
    void printOneBodyDensity            (string filename);

    double getEquilibrationFraction()   { return m_equilibrationFraction; }
    double getStepLength()              const;
    double getTimeStep()                const;
    double getSqrtTimeStep()            const;
    double getPsiOld()                  const;
    double getinteractionSize()         const;
    double getDistanceMatrixij          (int i, int j)      const;
    double computedistance              (int i);
    double gradientDescent              (double initialAlpha, string filename, int maxIterations);
    double findEnergyDerivative();
    double getBucketSize()              const;

    vector<vector<double>>              getDistanceMatrix() const;
    vector<vector<double>>              getQuantumForce()   const;
    vector<int>                         getHistogram()      const;

private:
    int                                 m_numberOfParticles = 0;
    int                                 m_numberOfDimensions = 0;
    int                                 m_numberOfMetropolisSteps = 0;
    int                                 m_bins = 0;
    double                              m_equilibrationFraction = 0;
    double                              m_stepLength = 0.1;
    double                              m_psiOld = 0;
    double                              m_timeStep = 0;
    double                              m_sqrtTimeStep = 0;
    double                              m_interactionSize = 0;
    double                              m_bucketSize = 0;
    class WaveFunction*                 m_waveFunction = nullptr;
    class Hamiltonian*                  m_hamiltonian = nullptr;
    class InitialState*                 m_initialState = nullptr;
    class Sampler*                      m_sampler = nullptr;
    std::vector<class Particle*>        m_particles = std::vector<class Particle*>();
    std::vector<int>                    m_histogram;
    std::vector<std::vector<double>>    m_QuantumForce;
    std::vector<std::vector<double>>    m_distanceMatrix = std::vector<std::vector<double>>();
};

