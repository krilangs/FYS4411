#pragma once
#include <vector>
#include "../particle.h"

class InitialState {
public:
    InitialState(class System* system);
    virtual void setupInitialState() = 0;
    std::vector<class Particle*> getParticles() { return m_particles; }

protected:
    int m_numberOfDimensions = 0;
    int m_numberOfParticles = 0;
    class System* m_system = nullptr;
    class Random* m_random = nullptr;
    std::vector<Particle*> m_particles;
};

