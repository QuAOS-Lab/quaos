#include "Measurements/pch.h"

#include "Measurements/Models/BaseNoiseModel.h"

namespace Measurements {

BaseNoiseModel::BaseNoiseModel() {}

double BaseNoiseModel::ComputeNoise() { return 0; }

}  // namespace Measurements