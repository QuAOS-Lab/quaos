#include "Measurements/pch.h"

#include "Measurements/Aquire/Aquire.h"

namespace Measurements {

Aquire::Aquire(PauliSum const& hamiltonian, std::vector<std::complex<double>> const& psi, AquireConfig const& config)
    : hamiltonian(hamiltonian), config(config) {}

MeasureResult Aquire::Measure(uint64_t shots) { return {}; }

void Aquire::InitializeDiagnosticCircuit() {}

void Aquire::UpdateCovarianceGraph() {}

void Aquire::AllocateMeasurements(std::uint64_t shots) {}

}  // namespace Measurements
