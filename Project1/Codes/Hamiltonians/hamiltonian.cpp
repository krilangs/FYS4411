#include "hamiltonian.h"
#include "../system.h"
#include "../WaveFunctions/simplegaussian.h"
#include "../particle.h"
#include <vector>

Hamiltonian::Hamiltonian(System* system) {
    m_system = system;
}

double Hamiltonian::computeNumericalDoubleDerivative (std::vector<Particle*> particles){
//Function to compute the numerical double derivative of the wavefunction.
    double h=0.001;
    double h_squared=h*h;
    double wf=0;
    double backward=0;
    double forward =0;
    double present=0;

    int dim = m_system->getNumberOfDimensions();
    int N   = m_system->getNumberOfParticles();
    present = m_system->getWaveFunction()->evaluate(particles);
    std::vector <double> r(dim);

    for (int j=0; j<N; j++){
        for (int d=0; d<dim; d++){
            r[d] = particles.at(j)->getPosition()[d];
        }
        for (int d=0; d<dim; d++){
            r[d] -=h;
            particles.at(j)->setPosition(r);
            backward += m_system->getWaveFunction()->evaluate(particles);

            r[d] +=2*h;
            particles.at(j)->setPosition(r);
            forward += m_system->getWaveFunction()->evaluate(particles);

            r[d] -=h;
            particles.at(j)->setPosition(r);
        }
    }
    wf = (backward+forward-2*present*dim*N)/(h_squared);

    return wf/(present);
}
