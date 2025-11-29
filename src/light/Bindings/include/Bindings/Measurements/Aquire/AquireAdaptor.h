#pragma once

namespace Bindings {

class AquireAdaptor {
public:
    AquireAdaptor(std::string hamiltonian,
                  std::vector<std::complex<double>> psi,
                  std::map<std::string, nb::object> config);

    ~AquireAdaptor() = default;
    AquireAdaptor(const AquireAdaptor&) = default;
    AquireAdaptor(AquireAdaptor&&) = default;

    nb::dict measure(std::uint64_t samples);

private:
    Measurements::Aquire native_;
};

}  // namespace Bindings