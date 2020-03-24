#pragma once
#include "initialstate.h"

class RandomUniform : public InitialState {
public:
    RandomUniform(System* system, int numberOfDimensions, int numberOfParticles, double timeStep, double interactionSize,
                  int bins, double bucketSize, double charLength);
    void setupInitialState();
};

