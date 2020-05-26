#pragma once
#include "wavefunction.h"

using namespace std;

class NeuralQuantumState : public WaveFunction {
public:
    NeuralQuantumState(class System* system);
    double evaluate(double GibbsValue, vector<double> X, vector<double> H, vector<double> a, vector<double> b,
                    vector<vector<double>> w);
    double computeDoubleDerivative(double GibbsValue, vector<double> X, vector<double> H, vector<double> a, vector<double> b,
                                   vector<vector<double>> w);
    vector<double> QuantumForce(vector<double> X, vector<double> a, vector<double> b,
                                vector<vector<double>> w);

protected:
    class WaveFunction* m_wavefunction = nullptr;
};
