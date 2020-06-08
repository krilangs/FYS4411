#include <iostream>
#include <chrono>
#include <string>
#include <cmath>
#include <vector>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/neuralquantumstate.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"

using namespace std;

int main() {
    /* The restricted Boltzmann machine (RBM) applied to the quantum many body problem.
     * Uncomment the desired sampling method to use, and set the boolean interaction
     * parameter to be false or true.
     */
    int numberOfDimensions  = 2;    // 1 or 2
    int numberOfParticles   = 2;    // 1 or 2
    int numberOfHiddenNodes = 2;

    int numberOfSteps       = int (1e5);    // Number of MC steps per SGD cycle
    int numberOfCycles      = 600;          // Number of SGD cycles

    string method = "Metropolis";
    //string method = "Importance";
    //string method = "Gibbs";

    bool interaction = true;

    // Variational/RBM parameters
    double learningRate     = 0.2;
    double omega            = 1.0;      // Oscillator frequency
    double sigma            = 1.0;      // Sigma in energy function (Gibbs)
    double stepLength       = 0.5;      // Metropolis step length
    double timeStep         = 0.5;      // Importance sampling time step

    double equilibration    = 0.2;      // Amount of the total steps used for equilibrium

    int numberOfVisible_Nodes = numberOfParticles*numberOfDimensions;  // Number of visible nodes
    int numberOfParameters = numberOfVisible_Nodes + numberOfHiddenNodes + numberOfVisible_Nodes*numberOfHiddenNodes;

    // Vectors for RBM parameters and gradient
    vector<double> X(numberOfVisible_Nodes);
    vector<double> H = vector<double>(numberOfHiddenNodes);
    vector<double> a = vector<double>(numberOfVisible_Nodes);
    vector<double> b = vector<double>(numberOfHiddenNodes);
    vector<vector<double>> w(numberOfVisible_Nodes, vector<double>(numberOfHiddenNodes));
    vector<double> G(numberOfParameters);

    double E_L;
    if (interaction == true)
        E_L = 3.0;   // Analytical energy for 2 interacting electrons in 2D
    else
        E_L = 0.5*numberOfParticles*numberOfDimensions;  // Analytical energy for non-interacting electrons

    System* system = new System();
    system->setHamiltonian              (new HarmonicOscillator(system, omega));
    system->setWaveFunction             (new NeuralQuantumState(system));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles, numberOfHiddenNodes,
                                                           numberOfVisible_Nodes, sigma, X, H, a, b, w, timeStep, numberOfParameters));
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);

    // Data files for SQD cycle and final MC run
    string filename_cycle;
    string filename_final;

    if (method == "Metropolis"){
        if (interaction == true){ // Interaction
            // SGD data file
            filename_cycle = "Data/BruteCycleInt_"+to_string(stepLength)+"_n_"+to_string(learningRate)+"_Np_"+to_string(numberOfParticles)
                             +"_D_"+to_string(numberOfDimensions)+"_H_"+to_string(numberOfHiddenNodes)+".dat";
        }
        else{ // No interaction
        // SGD data file
        filename_cycle = "Data/BruteCycle_"+to_string(stepLength)+"_n_"+to_string(learningRate)+"_Np_"+to_string(numberOfParticles)
                         +"_D_"+to_string(numberOfDimensions)+"_H_"+to_string(numberOfHiddenNodes)+".dat";
        // Instantaneous energy data file (of bigger run after SGD)
        filename_final = "Data/FinalBruteCycle_"+to_string(stepLength)+"_n_"+to_string(learningRate)+"_Np_"+to_string(numberOfParticles)
                         +"_D_"+to_string(numberOfDimensions)+"_H_"+to_string(numberOfHiddenNodes)+".dat";
        }
    }

    if (method == "Importance"){
        if (interaction == true){ // Interaction
            // SGD data file
            filename_cycle = "Data/ImpCycleInt_"+to_string(timeStep)+"_n_"+to_string(learningRate)+"_Np_"+to_string(numberOfParticles)
                             +"_D_"+to_string(numberOfDimensions)+"_H_"+to_string(numberOfHiddenNodes)+".dat";
        }
        else{ // No interaction
        // SGD data file
        filename_cycle = "Data/ImpCycle_"+to_string(timeStep)+"_n_"+to_string(learningRate)+"_Np_"+to_string(numberOfParticles)
                         +"_D_"+to_string(numberOfDimensions)+"_H_"+to_string(numberOfHiddenNodes)+".dat";
        // Instantaneous energy data file (of bigger run after SGD)
        filename_final = "Data/FinalImpCycle_"+to_string(timeStep)+"_n_"+to_string(learningRate)+"_Np_"+to_string(numberOfParticles)
                         +"_D_"+to_string(numberOfDimensions)+"_H_"+to_string(numberOfHiddenNodes)+".dat";
        }
    }

    if (method == "Gibbs"){
        if (interaction == true){ // Interaction
            // SGD data file
            filename_cycle = "Data/GibCycleInt_"+to_string(sigma)+"_n_"+to_string(learningRate)+"_Np_"+to_string(numberOfParticles)
                             +"_D_"+to_string(numberOfDimensions)+"_H_"+to_string(numberOfHiddenNodes)+".dat";
        }
        else{ // No interaction
        // SGD data file
        filename_cycle = "Data/GibCycle_"+to_string(sigma)+"_n_"+to_string(learningRate)+"_Np_"+to_string(numberOfParticles)
                         +"_D_"+to_string(numberOfDimensions)+"_H_"+to_string(numberOfHiddenNodes)+".dat";
        // Instantaneous energy data file (of bigger run after SGD)
        filename_final = "Data/FinalGibCycle_"+to_string(sigma)+"_n_"+to_string(learningRate)+"_Np_"+to_string(numberOfParticles)
                         +"_D_"+to_string(numberOfDimensions)+"_H_"+to_string(numberOfHiddenNodes)+".dat";
        }
    }

    system->setLearningRate         (learningRate);
    system->setNumberOfParameters   (numberOfParameters);
    system->openFile(filename_cycle);
    cout << "Start Metropolis" << endl;
    auto start = chrono::system_clock::now();
    for (int cycles=0; cycles<numberOfCycles; cycles++){
        system->runMetropolisSteps      (method, numberOfSteps, G, interaction, X, H, a, b, w);
        system->gradientDescent         (G, X, a, b, w);
        system->printOut                (cycles);
        system->writeToFile             (X, a, b, w);
    }

    if (interaction == false){
        cout << " Final run" << endl;
        //system->openDataFile            (filename_final);
        int finalNumberOfSteps = 1.5e+6;
        system->runMetropolisSteps      (method, finalNumberOfSteps, G, interaction, X, H, a, b, w);
        system->printOut                (numberOfCycles);
        //system->writeToFile             (X, a, b, w);
    }
    else{
        system->runMetropolisSteps      (method, numberOfSteps, G, interaction, X, H, a, b, w);
        system->printOut                (numberOfCycles);
        system->writeToFile             (X, a, b, w);
    }

    auto end = chrono::system_clock::now();
    chrono::duration<double> diff = end-start;
    cout << " Benchmark energy = " << E_L << " a.u. " << endl;
    cout << " Computation time = " << diff.count() << "s\n" << endl;

    return 0;
}
