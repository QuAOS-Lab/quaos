#pragma once

namespace Measurements {

class INoiseModel {
public:
    virtual ~INoiseModel() = default;

    virtual double ComputeNoise() = 0;
};

}  // namespace Measurements