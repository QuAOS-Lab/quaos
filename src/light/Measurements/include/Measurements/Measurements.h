#pragma once

#include "Measurements/MeasurementsExport.h"
#include "Measurements/Aquire/AquireConfig.h"
#include "Measurements/Aquire/Aquire.h"
#include "Measurements/Models/BaseErrorModel.h"
#include "Measurements/Models/BaseNoiseModel.h"

namespace Measurements {

constexpr auto VERSION = "1.0.0";

MEASUREMENTS_API const char* GetVersion();

}  // namespace Measurements

#undef MEASUREMENTS_API