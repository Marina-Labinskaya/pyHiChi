#pragma once
#include <cmath>

#include "FP.h"
#include "Dimension.h"

namespace pfc {
    inline FP formfactorTSC(FP x)
    {
        FP xa = fabs(x);
#ifdef __MIC
        return (xa <= (FP)0.5) * ((FP)0.75 - x * x) + (x > (FP)0.5) * ((FP)0.5 * ((FP)1.5 - xa) * ((FP)1.5 - xa));
#else
        if (xa <= (FP)0.5)
            return (FP)0.75 - x * x;
        else
            return (FP)0.5 * ((FP)1.5 - xa) * ((FP)1.5 - xa);
#endif
    }

    inline FP formfactorPCS(FP x)
    {
        FP xa = fabs(x);
#ifdef __MIC
        return (xa <= (FP)1.0) * (((FP)4.0 - 6 * xa * xa + 3 * xa * xa * xa) / (FP)6.0)
            + (xa > (FP)1.0) * (((FP)2.0 - xa) * ((FP)2.0 - xa) * ((FP)2.0 - xa) / (FP)6.0);
#else
        if (xa <= (FP)1.0)
            return ((FP)4.0 - 6 * xa * xa + 3 * xa * xa * xa) / (FP)6.0;
        else
            return ((FP)2.0 - xa) * ((FP)2.0 - xa) * ((FP)2.0 - xa) / (FP)6.0;
#endif
    }

    inline void formfactorFourthOrder(FP coeff, FP c[])
    {
        c[0] = (coeff + (FP)1) * coeff * (coeff - (FP)1) * (coeff - (FP)2) / (FP)24;
        c[1] = -(coeff + (FP)2) * coeff * (coeff - (FP)1) * (coeff - (FP)2) / (FP)6;
        c[2] = (coeff + (FP)2) * (coeff + (FP)1) * (coeff - (FP)1) * (coeff - (FP)2) / (FP)4;
        c[3] = -(coeff + (FP)2) * (coeff + (FP)1) * coeff * (coeff - (FP)2) / (FP)6;
        c[4] = (coeff + (FP)2) * (coeff + (FP)1) * coeff * (coeff - (FP)1) / (FP)24;
    }

    template<pfc::Dimension dimension, class T_FP>
    class FormFactor {
    public:
        void operator()(T_FP coords) {};

    };

    template<pfc::Dimension dimension, class T_FP>
    class FormFactorCIC : public FormFactor<dimension, T_FP> {
    public:
        void operator()(T_FP coords) {
            for (int i = 0; i < dimension; ++i) {
                for (int j = 0; j < 2; ++j) {
                    c[i][j] = 1;
                }
            }
        }
    public:
        FP c[dimension][2];
    };

    template<pfc::Dimension dimension, class T_FP>
    class FormFactorTSC : public FormFactor<dimension, T_FP> {
    public:
        void operator()(T_FP coords) {
            for (int i = 0; i < dimension; ++i) {
                for (int j = 0; j < 3; ++j) {
                    c[i][j] = formfactorTSC(FP(j - 1) - coords[i]);
                }
            }
        }
    public:
        FP c[dimension][3];
    };
}