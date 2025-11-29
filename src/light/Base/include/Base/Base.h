#pragma once

#include "Base/BaseExport.h"
#include "Base/Commons.h"

namespace Base {

constexpr auto VERSION = "1.0.0";

BASE_API const char* GetVersion();

}  // namespace Base

#undef BASE_API
