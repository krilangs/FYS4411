#pragma once
#include <vector>

using namespace std;

class Hamiltonian {
public:
    Hamiltonian(class System* system);
    virtual double computeLocalEnergy(bool interaction, double GibbsValue, vector<double> X, vector<double> H, vector<double> a,
                                      vector<double> b, vector<vector<double>> w) = 0;

protected:
    class System* m_system = nullptr;
};

