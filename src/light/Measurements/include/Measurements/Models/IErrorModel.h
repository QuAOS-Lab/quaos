#pragma once

namespace Measurements {

class IErrorModel {
public:
    virtual ~IErrorModel() = default;

    virtual double ComputeError() = 0;
};

}  // namespace Measurements