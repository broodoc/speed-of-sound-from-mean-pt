// Copyright 2026 Caleb Broodo
// all TF1 function bodies: Levy/Tsallis/Hagedorn/TBW/etc + knee math)
#pragma once

#include <cmath>
#include <limits>

#include "TMath.h"

class Models {
public:
    enum class Function {
        LevyTsallis,
        Tsallis,
        Hagedorn,
        TBW,
        BW,
        FokkerPlanck,
        PowerLaw,
        ExpDecay,
        RightKneeFunction,
        RightKneeFunctionApproximation
    };

    Models();
    Models(Function f);
    virtual ~Models(){};

    double operator()(double* x, double* par) const;

private:
    Function type_f = Function::LevyTsallis;  // default to something sane
};


