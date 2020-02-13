#include <iostream>
#include <chrono>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"

using namespace std;


int main() {
    int numberOfDimensions  = 3;
    int numberOfParticles   = 10;
    int numberOfSteps       = (int) 1e6;
    double omega            = 1.0;          // Oscillator frequency.
    double alpha            = 0.5;          // Variational parameter.
    double stepLength       = 0.5;          // Metropolis step length.
    double equilibration    = 0.2;          // Amount of the total steps used
    // for equilibration.

    System* system = new System();
    system->setHamiltonian              (new HarmonicOscillator(system, omega));
    system->setWaveFunction             (new SimpleGaussian(system, alpha));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);

    auto start = chrono::system_clock::now();
    system->runMetropolisSteps          (numberOfSteps);
    auto end = chrono::system_clock::now();
    chrono::duration<double> diff = end-start;
    cout << "Computation time = " << diff.count() << "s\n" << endl;

    return 0;
}
