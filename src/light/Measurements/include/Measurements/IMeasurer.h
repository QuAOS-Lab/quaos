#pragma once

namespace Measurements {

struct MEASUREMENTS_API MeasureResult {
    std::vector<std::complex<double>> measures;
};

class MEASUREMENTS_API IMeasurer {
public:
    virtual ~IMeasurer() = default;
    virtual MeasureResult Measure(uint64_t shots) = 0;
};

}  // namespace Measurements