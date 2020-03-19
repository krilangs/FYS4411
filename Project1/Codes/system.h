#pragma once
#include <vector>
#include <string>

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
    void openDataFile               (std::string filename);
    std::vector<std::vector<double>>    computematrixdistance(std::vector<class Particle*> particles);
    class WaveFunction*             getWaveFunction()   { return m_waveFunction; }
    class Hamiltonian*              getHamiltonian()    { return m_hamiltonian; }
    class Sampler*                  getSampler()        { return m_sampler; }
    std::vector<class Particle*>    getParticles()      { return m_particles; }
    int getNumberOfParticles()          { return m_numberOfParticles; }
    int getNumberOfDimensions()         { return m_numberOfDimensions; }
    int getNumberOfMetropolisSteps()    { return m_numberOfMetropolisSteps; }
    double getEquilibrationFraction()   { return m_equilibrationFraction; }

    double gradientDescent(double initialAlpha, std::string filename, int maxIterations);
    double findEnergyDerivative();

    double getPsiOld() const;
    void setPsiOld(double psiOld);

    double getinteractionSize() const;
    void setinteractionSize(double interactionSize);

    std::vector<std::vector<double> > getDistanceMatrix() const;
    void setDistanceMatrix(const std::vector<std::vector<double> > &distanceMatrix);
    bool updateDistanceMatrix(std::vector<Particle*> particles, int randparticle);
    double getDistanceMatrixij(int i, int j) const;
    double computedistance          (int i);

    std::vector<std::vector<double>> getQuantumForce() const;
    void setQuantumForce(const  std::vector<std::vector<double>> &QuantumForce);
    void updateQuantumForce( std::vector<std::vector<double>> deltaQuantumForce, bool subtract);

    double getTimeStep() const;
    void setTimeStep(double timeStep);
    double getSqrtTimeStep() const;
    void setSqrtTimeStep(double sqrtTimeStep);

    std::vector<int> getHistogram() const;
    void setHistogram();

    void oneBodyDensity();
    void printOneBodyDensity(std::string filename);

    int getBins() const;
    void setBins(int bins);

    double getBucketSize() const;
    void setBucketSize(double bucketSize);

private:
    int                             m_numberOfParticles = 0;
    int                             m_numberOfDimensions = 0;
    int                             m_numberOfMetropolisSteps = 0;
    int                             m_bins = 0;
    double                          m_equilibrationFraction = 0;
    double                          m_stepLength = 0.1;
    double                          m_psiOld = 0;
    double                          m_timeStep = 0;
    double                          m_sqrtTimeStep = 0;
    double                          m_interactionSize = 0;
    double                          m_bucketSize = 0;
    std::vector<int>                m_histogram;
    std::vector<std::vector<double>> m_QuantumForce;
    class WaveFunction*             m_waveFunction = nullptr;
    class Hamiltonian*              m_hamiltonian = nullptr;
    class InitialState*             m_initialState = nullptr;
    class Sampler*                  m_sampler = nullptr;
    std::vector<class Particle*>    m_particles = std::vector<class Particle*>();
    std::vector<std::vector<double>>     m_distanceMatrix = std::vector<std::vector<double>>();
};

