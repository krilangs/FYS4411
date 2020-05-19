#pragma once
#include "initialstate.h"

using namespace std;

class RandomUniform : public InitialState {
public:
    RandomUniform(System* system, int numberOfDimensions, int numberOfParticles, int numberOfHiddenNodes, int numberOfVisibleNodes,
                  double sigma, vector<double> &m_X, vector<double> &m_H, vector<double> &m_a, vector<double> &m_b,
                  vector<vector<double>> &m_w, double timeStep, int numberOfParameters);
    void setupInitialState(vector<double> &m_X, vector<double> &m_H, vector<double> &m_a, vector<double> &m_b,
                           vector<vector<double>> &m_w);
};
