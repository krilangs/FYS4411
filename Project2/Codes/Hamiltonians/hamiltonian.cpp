#include "hamiltonian.h"
#include "../system.h"
#include "../WaveFunctions/neuralquantumstate.h"
#include "../particle.h"
#include <vector>

Hamiltonian::Hamiltonian(System* system) {
    m_system = system;
}

