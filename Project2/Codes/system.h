#pragma once
#include <vector>
#include <string>

using namespace std;

class System {
public:
    bool metropolisStep             (double GibbsValue, vector<double> &X, vector<double> H, vector<double> a, vector<double> b,
                                     vector<vector<double>> w);
    bool metropolisStepImportance   (double GibbsValue, vector<double> &X, vector<double> H, vector<double> a, vector<double> b,
                                     vector<vector<double>> w);
    bool GibbsSampling              (vector<double> &X, vector<double> H, vector<double> a, vector<double> b,
                                     vector<vector<double>> w);

    void runMetropolisSteps         (string method, int numberOfMetropolisSteps, vector<double> &G, bool interaction, vector<double> X,
                                     vector<double> H, vector<double> a, vector<double> b, vector<vector<double>> w);
    void setNumberOfParticles       (int numberOfParticles);
    void setNumberOfDimensions      (int numberOfDimensions);
    void setNumberOfVisibleNodes    (int numberOfVisibleNodes);
    void setNumberOfHiddenNodes     (int numberOfHiddenNodes);
    void setNumberOfParameters      (int numberOfParameters);
    void setStepLength              (double stepLength);
    void setEquilibrationFraction   (double equilibrationFraction);
    void setHamiltonian             (class Hamiltonian* hamiltonian);
    void setWaveFunction            (class WaveFunction* waveFunction);
    void setInitialState            (class InitialState* initialState);
    void openDataFile               (string filename);
    void printOut                   (int cycle);
    void writeToFile                (vector<double> X, vector<double>& a, vector<double>& b, vector<vector<double>>& w);
    void openFile                   (string filename);

    class WaveFunction*             getWaveFunction()   { return m_waveFunction; }
    class Hamiltonian*              getHamiltonian()    { return m_hamiltonian; }
    class Sampler*                  getSampler()        { return m_sampler; }
    vector<class Particle*>         getParticles()      { return m_particles; }
    vector<vector<double>>          computematrixdistance(vector<double> &X);

    int getNumberOfParticles()          { return m_numberOfParticles; }
    int getNumberOfDimensions()         { return m_numberOfDimensions; }
    int getNumberOfMetropolisSteps()    { return m_numberOfMetropolisSteps; }
    int getNumberOfVisibleNodes()       const;
    int getNumberOfHiddenNodes()        const;
    int getNumberOfParameters()         const;
    int computeIndex                    (int index);

    void setTimeStep                    (double timeStep);
    void setSqrtTimeStep                (double sqrtTimeStep);
    void setPsiOld                      (double psiOld);
    void setSigma                       (double sigma);
    void setSigma_squared               (double sigma_squared);
    void setLearningRate                (double learningRate);
    void gradientDescent                (vector<double> G, vector<double> X, vector<double> &a, vector<double> &b, vector<vector<double>> &w);
    void setDistanceMatrix              (const vector<vector<double> > &distanceMatrix);
    void updateDistanceMatrix           (vector<double> X, int randparticle);
    void setQuantumForce                (const vector<double> &QuantumForce);
    void setGradient                    (const vector<double> &Gradient);
    void setEnGradient_average          (const vector<double> &EnGradient_average);
    void setCumulativeGradient          (const vector<double> &cumulativeGradient);
    void setCumulativeEnGradient        (const vector<double> &cumulativeEnGradient);

    double getEquilibrationFraction()   { return m_equilibrationFraction; }
    double getStepLength()              const;
    double getTimeStep()                const;
    double getSqrtTimeStep()            const;
    double getPsiOld()                  const;
    double getSigma()                   const;
    double getSigma_squared()           const;
    double getLearningRate()            const;
    double getDistanceMatrixij          (int i, int j)      const;
    double computedistance              (int i);
    double findEnergyDerivative();

    vector<vector<double>>              getDistanceMatrix() const;
    vector<double>                      getQuantumForce()   const;
    vector<double>                      getGradient() const;
    vector<double>                      getEnGradient_average() const;
    vector<double>                      getCumulativeGradient() const;
    vector<double>                      getCumulativeEnGradient() const;
    vector<double>                      GradientParameters(double GibbsValue, vector<double> X, vector<double> &a, vector<double> &b,
                                                           vector<vector<double>> &w);

private:
    int                                 m_numberOfParticles = 0;
    int                                 m_numberOfDimensions = 0;
    int                                 m_numberOfMetropolisSteps = 0;
    int                                 m_numberOfVisibleNodes = 0;
    int                                 m_numberOfHiddenNodes = 0;
    int                                 m_numberOfParameters = 0;
    double                              m_equilibrationFraction = 0;
    double                              m_stepLength = 0.1;
    double                              m_psiOld = 0;
    double                              m_timeStep = 0;
    double                              m_sqrtTimeStep = 0;
    double                              m_sigma = 0;
    double                              m_sigma_squared = 0;
    double                              m_learningRate = 0;
    class WaveFunction*                 m_waveFunction = nullptr;
    class Hamiltonian*                  m_hamiltonian = nullptr;
    class InitialState*                 m_initialState = nullptr;
    class Sampler*                      m_sampler = nullptr;
    vector<class Particle*>             m_particles = vector<class Particle*>();
    vector<double>                      m_cumulativeGradient = vector<double>();
    vector<double>                      m_cumulativeEnGradient = vector<double>();
    vector<double>                      m_Gradient = vector<double>();
    vector<double>                      m_EnGradient_average = vector<double>();
    vector<double>                      m_QuantumForce;
    vector<vector<double>>              m_distanceMatrix = vector<vector<double>>();
};

