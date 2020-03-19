#include <iostream>
#include <chrono>
#include <string>
#include <cmath>
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
    int numberOfParticles   = 100;
    int numberOfSteps       = int (1e5);
    double omega            = 1.0;          // Oscillator frequency.
    double alpha            = 0.5;          // Variational parameter.
    double beta             = 2.82843;      // 2.82843 for interaction, 1.0 for non-int.
    double omega_z          = beta;
    double stepLength       = 1.0;     //4.0     // Metropolis step length.
    double timeStep         = 0.001;         // Importance sampling time step
    double equilibration    = 0.2;          // Amount of the total steps used
    double interactionSize  = 0.0043;       // 0.0043 for interaction, 0.0 for non-int

    // Parameters for onebody density histogram
    double bucketSize = 0.01;
    int bins = ceil(4 / bucketSize);

    double GP = numberOfParticles/2*(beta+2)*(alpha+1/(4*alpha));  // Gross-Pitaevski equation 3D

    string filename = "0";
    //string filename = "Data/Alpha_" + to_string(alpha) + "_dim_" + to_string(numberOfDimensions) + "_particles_" + to_string(numberOfParticles) + "_non-int.dat";

    System* system = new System();
    system->setHamiltonian              (new HarmonicOscillator(system, omega, omega_z));
    system->setWaveFunction             (new SimpleGaussian(system, alpha, beta));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles, timeStep, interactionSize, bins, bucketSize));
    system->openDataFile                (filename);
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);

    /*
    // Gradient descent method to find energy minimum

    int maxIterations = 30;
    double initialAlpha = 0.40;
    string minFilename = "Data/Test" + to_string(initialAlpha) +"_N_" + to_string(numberOfParticles) + "_doublecheck.dat";
    alpha = system->gradientDescent(initialAlpha, minFilename, maxIterations);
    cout << "Optimal alpha found by steepest descent: " << alpha << endl;

    vector<double> parameters(3);
    parameters[0] = alpha;
    parameters[1] = alpha;
    parameters[2] = alpha*beta;
    system->getWaveFunction()->setParameters(parameters);
    */
    cout << "Start Metropolis" << endl;
    auto start = chrono::system_clock::now();
    system->runMetropolisSteps          (numberOfSteps);

    //string densityFilename = "Data/density_int_E_alpha_" + to_string(alpha) + "_beta_" + to_string(beta) + ".dat";
    //system->printOneBodyDensity(densityFilename);

    auto end = chrono::system_clock::now();
    chrono::duration<double> diff = end-start;

    cout << "Benchmark energy (GP) = " << GP << endl;
    cout << "Computation time = " << diff.count() << "s\n" << endl;

    return 0;
}
