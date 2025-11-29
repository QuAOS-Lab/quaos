#include "Bindings/pch.h"

#include "Bindings/Measurements/Aquire/AquireAdaptor.h"

namespace Bindings {

namespace {

Measurements::PauliSum GetHamiltonian(std::string const& hamiltonian) { return {}; }

std::vector<std::complex<double>> GetPsi(std::vector<std::complex<double>> psi) { return {}; }

Measurements::AquireConfig GetAquireConfig(std::map<std::string, nb::object> const& config) { return {}; }

}  // namespace

AquireAdaptor::AquireAdaptor(std::string hamiltonian,
                             std::vector<std::complex<double>> psi,
                             std::map<std::string, nb::object> config)
    : native_(GetHamiltonian(hamiltonian), GetPsi(psi), GetAquireConfig(config)) {}

nb::dict AquireAdaptor::measure(std::uint64_t samples) { return {}; }

}  // namespace Bindings