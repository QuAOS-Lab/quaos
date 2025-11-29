#pragma once

#include "Measurements/Models/INoiseModel.h"

namespace Measurements {

class MEASUREMENTS_API BaseNoiseModel : public INoiseModel {
public:
    BaseNoiseModel();

    // INoiseModel
    double ComputeNoise() override;
};

}  // namespace Measurements