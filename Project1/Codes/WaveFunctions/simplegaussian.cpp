#include <cmath>
#include <cassert>
#include <algorithm>
#include <vector>
#include "simplegaussian.h"
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

SimpleGaussian::SimpleGaussian(System* system, double alpha, double beta):
                               WaveFunction(system) {
    assert(alpha >= 0);
    assert(beta >= 0);
    m_numberOfParameters = 3;
    m_parameters.reserve(3);
    m_parameters.push_back(alpha);
    m_parameters.push_back(alpha);
    m_parameters.push_back(alpha*beta);
}

double SimpleGaussian::evaluate(std::vector<class Particle*> particles) {
// Implement the Gaussian wave function.
    double r_sqr = 0;
    double f     = 1;
    int dim = m_system->getNumberOfDimensions();
    int N   = m_system->getNumberOfParticles();

    for (int j=0; j<N; j++){
        if (N==1) {break;}
        for (int i=0; i<j; i++){
            if (m_system->getDistanceMatrixij(i,j) <= m_system->getinteractionSize()){
                f = 0.0;
                break;
            }
            f *= 1 - m_system->getinteractionSize()/(m_system->getDistanceMatrixij(i,j));
        }
    }
    for (int i=0; i<N; i++){
        for (int d=0; d<dim; d++){
            r_sqr += particles.at(i)->getPosition()[d]*particles.at(i)->getPosition()[d]*m_parameters[d];
        }
    }
    return exp(-r_sqr)*f;
}

double SimpleGaussian::computeDoubleDerivative(std::vector<class Particle*> particles) {
    /* Computes the double derivative of the wave function analytically.
     * Note that by double derivative, we actually mean the sum of
     * the Laplacians with respect to the coordinates of each particle.
     * This quantity is needed to compute the (local) energy.
     */
    double one = 0;
    double interaction = 0;
    double a = m_system->getinteractionSize();
    int dim  = m_system->getNumberOfDimensions();
    int N    = m_system->getNumberOfParticles();

    if (a != 0.0){
    // Interaction terms
    for (int i=0; i<N; i++){
        double r_i_square=0;

        for (int d=0; d<dim-1; d++){
            r_i_square += particles.at(i)->getPosition()[d]*
                          particles.at(i)->getPosition()[d];
        }
        int d = dim-1;
        r_i_square += particles.at(i)->getPosition()[d]*
                      particles.at(i)->getPosition()[d]*
                      m_parameters[2]/(m_parameters[0]);

        double second = 0;
        double third  = 0;
        double fourth = 0;
        double fifth  = 0;
        double temp;

        for (int j=0; j<i; j++) {
            double r_ij = m_system->getDistanceMatrixij(i,j);
            temp    = a/((r_ij-a)*r_ij);
            second += temp;
            double r_ir_j = 0;

            for (int d=0; d<dim-1; d++){
                r_ir_j += particles.at(i)->getPosition()[d]*
                          particles.at(j)->getPosition()[d];
            }
            int d = dim-1;
            r_ir_j += particles.at(i)->getPosition()[d]*
                      particles.at(j)->getPosition()[d]*
                      m_parameters[2]/(m_parameters[0]);
            fourth -= temp*temp;
            fifth  -= 4*m_parameters[0]*(r_i_square - r_ir_j)*temp/(r_ij);
        }
        for (int j=i+1; j<N; j++){
            double r_ij = m_system->getDistanceMatrixij(i,j);
            temp    = a/((r_ij-a)*r_ij);
            second += temp;
            double r_ir_j = 0;

            for (int d=0; d<dim-1; d++){
                r_ir_j += particles.at(i)->getPosition()[d]*
                          particles.at(j)->getPosition()[d];
            }
            int d   = dim-1;
            r_ir_j += particles.at(i)->getPosition()[d]*
                      particles.at(j)->getPosition()[d]*
                      m_parameters[2]/(m_parameters[0]);
            fourth -= temp*temp;
            fifth  -= 4*m_parameters[0]*(r_i_square - r_ir_j)*temp/(r_ij);
        }
        third = second*second;
        interaction += second+third+fourth+fifth;
    }
    }

    // For non-interacting particles
    for (int i=0; i<N; i++){
        for (int d=0; d<dim; d++){
            one += m_parameters[d]*m_parameters[d]*
                   particles.at(i)->getPosition()[d]*
                   particles.at(i)->getPosition()[d];
        }
    }
    one *= 4.0;
    one -= 2*((dim - 1)*m_parameters[0] + m_parameters[2])*N;
    return one + interaction;
}

std::vector<std::vector<double>> SimpleGaussian::QuantumForce(std::vector<class Particle*> particles) {
//Function to compute the Quantum Force for the Importance sampling method.
    double a = m_system->getinteractionSize();
    double constant;
    double r_kj;
    int dim = m_system->getNumberOfDimensions();
    int N   = m_system->getNumberOfParticles();
    std::vector<std::vector<double>> QuantumForce(dim, std::vector<double>(N));

    for (int d=0; d<dim; d++){
        for (int k=0; k<N; k++){
            QuantumForce[d][k] = -2*(m_parameters[d]*particles.at(k)->getPosition()[d]);
            for (int j=0; j<k; j++){
                r_kj = m_system->getDistanceMatrixij(k,j);
                constant = 2*a/(r_kj*r_kj*(r_kj - a));
                QuantumForce[d][k] += (particles.at(k)->getPosition()[d] - particles.at(j)->getPosition()[d])*constant;
            }
        }
    }
    return QuantumForce;
}
