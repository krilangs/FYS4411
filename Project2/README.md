FYS4411 - Project 2: The restricted Boltzmann machine applied to the quantum many-body problem

Report-folder: Contains the report files for this project.

Figures-folder: Contains various figures created in this project produced by the Python codes in the Codes-folder.

Data-folder: Contains various data files produced by the programs, and are used to produce the figures and in the blocking program.

Codes-folder: Contains all the codes for this project. All the code was done on a Windows 10 computer with 8GB RAM.
- The C++ codes are done in Qt Creator, and are compiled and run with the following terminal command: 
`gcc -O3 -o main.exe main.cpp sampler.cpp system.cpp particle.cpp Hamiltonians/hamiltonian.cpp Hamiltonians/harmonicoscillator.cpp InitialStates/initialstate.cpp InitialStates/randomuniform.cpp Math/random.cpp WaveFunctions/neuralquantumstate.cpp WaveFunctions/wavefunction.cpp`. To run the desired C++ source code of the VMC tasks for this project, check the Methods section in the report and read the comments in the .cpp files (main comments are in main.cpp) to use the same parameters as in the report.
- The Python (3.7) codes import .txt and .dat files created by the C++ codes, and plots the various figures in the Figures-folder and do a statistical analysis with the blocking method.
