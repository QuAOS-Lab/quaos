#pragma once

#include "Operators/OperatorsExport.h"

namespace Operators {

constexpr auto VERSION = "1.0.0";

OPERATORS_API const char* GetVersion();

}  // namespace Operators

#undef OPERATORS_API