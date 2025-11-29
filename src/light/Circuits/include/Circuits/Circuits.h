#pragma once

#include "Circuits/CircuitsExport.h"

namespace Circuits {

constexpr auto VERSION = "1.0.0";

CIRCUITS_API const char* GetVersion();

}  // namespace Circuits

#undef CIRCUITS_API