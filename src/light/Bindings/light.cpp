#include "Bindings/pch.h"

using ExternalAquireConfig = std::map<std::string, nb::object>;

NB_MODULE(light, m) {
    nb::class_<Bindings::AquireAdaptor>(m, "Aquire")
        .def(nb::init<std::string const&, std::vector<std::complex<double>> const&, ExternalAquireConfig const&>(),
             nb::arg("hamiltonian"), nb::arg("psi"), nb::arg("config"),
             nb::sig("def __init__(self, hamiltonian: str, psi: list[complex], config: dict)"),
             "Initialize AcquireAdaptor with a Hamiltonian, psi vector, and config dictionary")
        .def("measure", &Bindings::AquireAdaptor::measure, nb::arg("samples"),
             nb::sig("def measure(self, samples: int) -> dict"),
             "Run measurement with the given number of samples and return a dictionary");
}