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
    /* To run the different cases, change the parameters below according to the comments to the right.
     * stepLength:
       - For non-interaction in; 1D: =4.0, 2D: =2.5, 3D: =2.0,
       - For interaction in 3D: =1.5
     * In harmonicoscillator.cpp; in computeLocalEnergy():
       - Choose to do either the analytical or numerical kinetic energy calculation.
     * In system.cpp; in runMetropolisSteps():
       - Choose to do either brute force Metropolis or Importance sampling.
     * To do gradient descent:
       - Uncomment the gradient section below.
     * To write the energies to file for blocking:
       - Change from filename=0 to filename="".
     * In sampler.cpp; in sample():
       - Uncomment m_system->oneBodyDensity(), and below here to run one-body density.
     */
    int numberOfDimensions  = 3;
    int numberOfParticles   = 10;
    int numberOfSteps       = int (1e5);    // Number of MC steps
    double omega            = 1.0;          // Oscillator frequency
    double alpha            = 0.522644;     // Variational parameter
    double beta             = 2.82843;      // 2.82843 for elliptical trap, 1.0 for spherical trap
    double omega_z          = beta;         // Frequency i z-direction
    double stepLength       = 1.5;          // Metropolis step length.
    double timeStep         = 0.001;        // Importance sampling time step
    double equilibration    = 0.2;          // Amount of the total steps used for equilibrium
    double interactionSize  = 0.0043;       // 0.0043 for interaction, 0.0 for non-interaction
    double charLength       = 2.5;          // Characteristic length of box: 1.0 (a_0) for non-interaction, 2.5 for interaction

    // Parameters for one-body density histogram
    double bucketSize = 0.01;
    int bins = int(ceil(4 / bucketSize));

    //for (alpha=0.2; alpha<0.71; alpha += 0.1){ Loop for varying alpha's

    double GP = (numberOfParticles/2.)*(numberOfDimensions/3.)*(beta+2)*(alpha+1/(4*alpha));  // Gross-Pitaevskii equation

    // Choose 0 for not to write to file, used for blocking
    string filename = "0";
    //string filename = "Data/Alpha_" + to_string(alpha) + "_dim_" + to_string(numberOfDimensions) + "_particles_"
    //                   + to_string(numberOfParticles) + ".dat";

    System* system = new System();
    system->setHamiltonian              (new HarmonicOscillator(system, omega, omega_z));
    system->setWaveFunction             (new SimpleGaussian(system, alpha, beta));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles, timeStep,
                                                           interactionSize, bins, bucketSize, charLength));
    system->openDataFile                (filename);
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);

    // Gradient descent method to find energy minimum, uncomment below to run:
/*
    int maxIterations = 200;                 // For N=10: =200, N=50: =50, N=100: =30
    double initialAlpha = 0.40;             // Initial guess, used 0.4 and 0.6
    string minFilename = "Data/find_min_int_" + to_string(initialAlpha) +"_N_" + to_string(numberOfParticles) + "_doublecheck.dat";
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

    // Uncomment to do one-body density and write to file
    string densityFilename = "Data/density_S_J_alpha_"+to_string(alpha)+"_beta_"+to_string(beta)+"_N_"
                                +to_string(numberOfParticles)+".dat";
    //system->printOneBodyDensity(densityFilename);

    auto end = chrono::system_clock::now();
    chrono::duration<double> diff = end-start;

    cout << "Benchmark energy (GP) = " << GP << endl;
    cout << "Computation time = " << diff.count() << "s\n" << endl;
    //}

    return 0;
}
