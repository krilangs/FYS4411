#pragma once
#include <vector>

using namespace std;

class WaveFunction {
public:
    WaveFunction(class System* system);
    int     getNumberOfParameters() { return m_numberOfParameters; }
    vector<double> getParameters() { return m_parameters; }
    virtual double evaluate(double GibbsValue, vector<double> X, vector<double> H, vector<double> a, vector<double> b,
                            vector<vector<double>> w) = 0;
    virtual double computeDoubleDerivative(double GibbsValue, vector<double> X, vector<double> H, vector<double> a, vector<double> b,
                                           vector<vector<double>> w) = 0;
    virtual vector<double> QuantumForce(double GibbsValue, vector<double> X, vector<double> a, vector<double> b,
                                        vector<vector<double>> w) = 0;

    void setParameters(const vector<double> &parameters);

protected:
    int     m_numberOfParameters = 0;
    vector<double> m_parameters = vector<double>();
    class System* m_system = nullptr;
};

