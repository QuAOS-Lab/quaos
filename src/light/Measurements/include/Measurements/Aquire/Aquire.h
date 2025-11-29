#pragma once

#include "Measurements/IMeasurer.h"
#include "Measurements/Aquire/AquireConfig.h"

namespace Measurements {

struct MEASUREMENTS_API Pauli {
    std::uint64_t x;
    std::uint64_t z;
    double phase;
};

using PauliString = std::vector<Pauli>;
using PauliSum = std::vector<PauliString>;

struct AquireConfig;

class MEASUREMENTS_API Aquire : public IMeasurer {
public:
    Aquire(PauliSum const& hamiltonian, std::vector<std::complex<double>> const& psi, AquireConfig const& config);

    // IMeasurer
    MeasureResult Measure(uint64_t shots) override;

private:
    void InitializeDiagnosticCircuit();

    void UpdateCovarianceGraph();

    void AllocateMeasurements(std::uint64_t shots);

    AquireConfig config;
    PauliSum hamiltonian;
};

}  // namespace Measurements