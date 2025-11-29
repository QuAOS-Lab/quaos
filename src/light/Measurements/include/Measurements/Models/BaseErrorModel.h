#pragma once

#include "Measurements/Models/IErrorModel.h"

namespace Measurements {

class MEASUREMENTS_API BaseErrorModel : public IErrorModel {
public:
    BaseErrorModel();

    // IErrorModel
    double ComputeError() override;
};

}  // namespace Measurements