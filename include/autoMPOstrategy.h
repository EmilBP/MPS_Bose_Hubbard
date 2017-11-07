#ifndef AUTOMPOSTRATEGY_H
#define AUTOMPOSTRATEGY_H

#include "itensor/all.h"

using namespace itensor;

class autoMPOstrategy {
private:

public:
    virtual AutoMPO updateMPO(double control) = 0;
    virtual AutoMPO derivative_control(double control) = 0;
};

#endif
