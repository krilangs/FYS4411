#pragma once
#include <vector>
#include "../particle.h"

using namespace std;

class InitialState {
public:
    InitialState(class System* system);
    virtual void setupInitialState(vector<double> &m_X, vector<double> &m_H, vector<double> &m_a, vector<double> &m_b,
                                   vector<vector<double>> &m_w) = 0;
    vector<class Particle*> getParticles() { return m_particles; }

protected:
    int m_numberOfDimensions = 0;
    int m_numberOfParticles  = 0;
    class System* m_system = nullptr;
    class Random* m_random = nullptr;
    vector<Particle*> m_particles;
};

