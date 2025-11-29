#include "Measurements/pch.h"

#include "Measurements/Models/BaseErrorModel.h"

namespace Measurements {

BaseErrorModel::BaseErrorModel() {}

double BaseErrorModel::ComputeError() { return 0; }

}  // namespace Measurements