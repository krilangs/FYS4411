#include <cmath>
#include <cassert>
#include <algorithm>
#include <vector>
#include "neuralquantumstate.h"
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

using namespace std;

NeuralQuantumState::NeuralQuantumState(System* system):WaveFunction(system) {
    m_numberOfParameters = 3;
    m_parameters.reserve(3);
}

double NeuralQuantumState::evaluate(double GibbsValue, vector<double> X, vector<double> H, vector<double> a, vector<double> b,
                                vector<vector<double>> w) {
    // Implement the NQS wave function.
    double first_sum = 0;
    double prod = 1;
    int M = m_system->getNumberOfVisibleNodes();
    int N = m_system->getNumberOfHiddenNodes();

    for (int i=0; i<M; i++){
        first_sum += (X[i] - a[i])*(X[i] - a[i]);
    }
    first_sum /= 2*m_system->getSigma_squared();
    first_sum = exp(-first_sum*GibbsValue);

    for (int j=0; j<N; j++){
        double second_sum = 0;
        for (int i=0; i<M; i++){
            second_sum += X[i]*w[i][j];
        }
        second_sum /= m_system->getSigma_squared();

        prod *= 1 + exp(b[j] + second_sum);
    }

    if (GibbsValue == 0.5)
        return first_sum*sqrt(prod);
    else
        return first_sum*prod;
}

double NeuralQuantumState::computeDoubleDerivative(double GibbsValue, vector<double> X, vector<double> H, vector<double> a, vector<double> b,
                                               vector<vector<double>> w) {
    /* Computes the double derivative of the wave function analytically.
     * Note that by double derivative, we actually mean the sum of
     * the Laplacians with respect to the coordinates of each particle.
     * This quantity is needed to compute the (local) energy.
     */
    int M = m_system->getNumberOfVisibleNodes();
    int N = m_system->getNumberOfHiddenNodes();

    vector <double> argument(N);

    double firstsum  = 0.0;
    double secondsum = 0.0;
    double kinetic   = 0.0;
    double temp2, temp3, sum;

    for (int j=0; j<N; j++){
        sum = 0;
        for (int i=0; i<M; i++){
            sum += X[i]*w[i][j]/m_system->getSigma_squared();
        }
        argument[j] = exp(-b[j] - sum);
    }

    for (int i=0; i<M; i++){
        temp2 = 0;
        temp3 = 0;
        for (int j=0; j<N; j++){
            double expon = 1.0 + argument[j];
            temp2 += w[i][j]/expon;
            temp3 += w[i][j]*w[i][j]*argument[j]/(expon*expon);
        }
        firstsum = -(X[i] - a[i]) + temp2;
        firstsum /= m_system->getSigma_squared();

        secondsum = temp3/(m_system->getSigma_squared()*m_system->getSigma_squared()) - 1.0/m_system->getSigma_squared();

        kinetic += firstsum*firstsum*GibbsValue*GibbsValue + secondsum*GibbsValue;
    }

    return kinetic;
}

vector<double> NeuralQuantumState::QuantumForce(double GibbsValue, vector<double> X, vector<double> a, vector<double> b,
                                            vector<vector<double>> w) {
    //Function to compute the Quantum Force for the Importance sampling method.
    int M = m_system->getNumberOfVisibleNodes();
    int N = m_system->getNumberOfHiddenNodes();

    vector <double> QuantumForce(M);
    vector <double> argument(N);
    vector <double> temp2(M);

    double sum;

    for (int j=0; j<N; j++){
        sum = 0;
        for (int i=0; i<M; i++){
            sum += X[i]*w[i][j]/m_system->getSigma_squared();
        }
        argument[j] = b[j] + sum;
    }

    for (int i=0; i<M; i++){
        temp2[i] = 0;
        for (int j=0; j<N; j++){
            double temp4 = exp(-argument[j]);
            double expon = 1 + temp4;
            temp2[i] += w[i][j]/(expon);
        }
        QuantumForce[i] = 2*(-(X[i] - a[i]) + temp2[i])*GibbsValue/m_system->getSigma_squared();
    }

    return QuantumForce;
}
