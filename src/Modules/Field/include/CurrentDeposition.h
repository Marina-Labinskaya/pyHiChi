#pragma once

#include "Grid.h"
#include "FormFactor.h"
#include "ParticleArray.h"
#include "Particle.h"
#include <iostream>

namespace pfc
{
    template<class TGrid, class DerivedClass>
    class CurrentDeposition
    {
    public:
        enum class ZeroizeJ {
            NOT_USE_ZEROIZEJ, USE_ZEROIZEJ
        };

        template<class T_Particle>
        void operator()(TGrid* grid, const T_Particle& particle, double dt,
            CurrentDeposition::ZeroizeJ UsingZeroizeJ = CurrentDeposition::ZeroizeJ::NOT_USE_ZEROIZEJ) {
            FormFactorCIC<Dimension::Three, FP3> formFactor;
            if (UsingZeroizeJ == ZeroizeJ::USE_ZEROIZEJ)
                grid->zeroizeJ();
            static_cast<DerivedClass*>(this)->depositOneParticle(grid, &particle, dt, formFactor);
        }

        template<class T_ParticleArray>
        void operator()(TGrid* grid, T_ParticleArray* particleArray, double dt) {
            typedef typename T_ParticleArray::ParticleProxyType ParticleProxyType;
            FormFactorCIC<Dimension::Three, FP3> formFactor;
            grid->zeroizeJ();

            for (int i = 0; i < particleArray->size(); i++) {
                ParticleProxyType particle = (*particleArray)[i];
                static_cast<DerivedClass*>(this)->depositOneParticle(grid, &particle, dt, formFactor);
            }
        }


        template<class T_Particle>
        void depositOneParticle(TGrid* grid, T_Particle* particle, double dt, FormFactorCIC<Dimension::Three, FP3>& formFactor) {
            static_assert(false, "ERROR: CurrentDeposition::depositOnePaticle shouldn't be called");
        }
    };

    template<class TGrid>
    class FirstOrderCurrentDeposition : public CurrentDeposition<TGrid, FirstOrderCurrentDeposition<TGrid>>
    {
    public:

        void depositComponentCurrent(ScalarField<FP> & field, const Int3 & idx,
            const FP3 & internalCoords, const FP & fieldBeforeDeposition, FormFactorCIC<Dimension::Three, FP3>& formFactor)
        {
            field(idx.x, idx.y, idx.z) += formFactor.c[0][0] * formFactor.c[1][0] * formFactor.c[2][0]
                * fieldBeforeDeposition;
            field(idx.x + 1, idx.y, idx.z) += formFactor.c[0][1] * formFactor.c[1][0] * formFactor.c[2][0]
                * fieldBeforeDeposition;
            field(idx.x, idx.y + 1, idx.z) += formFactor.c[0][0] * formFactor.c[1][1] * formFactor.c[2][0]
                * fieldBeforeDeposition;
            field(idx.x, idx.y, idx.z + 1) += formFactor.c[0][0] * formFactor.c[1][0] * formFactor.c[2][1]
                * fieldBeforeDeposition;
            field(idx.x + 1, idx.y + 1, idx.z) += formFactor.c[0][1] * formFactor.c[1][1] * formFactor.c[2][0]
                * fieldBeforeDeposition;
            field(idx.x + 1, idx.y, idx.z + 1) += formFactor.c[0][1] * formFactor.c[1][0] * formFactor.c[2][1]
                * fieldBeforeDeposition;
            field(idx.x, idx.y + 1, idx.z + 1) += formFactor.c[0][0] * formFactor.c[1][1] * formFactor.c[2][1]
                * fieldBeforeDeposition;
            field(idx.x + 1, idx.y + 1, idx.z + 1) += formFactor.c[0][1] * formFactor.c[1][1] * formFactor.c[2][1]
                * fieldBeforeDeposition;
        }

        template<class T_Particle>
        void depositOneParticle(TGrid* grid, T_Particle* particle, double dt, FormFactorCIC<Dimension::Three, FP3>& formFactor)
        {
            Int3 idxJx, idxJy, idxJz;
            FP3 internalCoordsJx, internalCoordsJy, internalCoordsJz;
            FP3 particlePosition = particle->getPosition(); // - (particle->getVelocity() * dt / 2.0);

            FP3 current = (particle->getVelocity() * particle->getCharge() * particle->getWeight()) /
                grid->steps.volume();

            grid->getIndexEJx(particlePosition, idxJx, internalCoordsJx);
            grid->getIndexEJy(particlePosition, idxJy, internalCoordsJy);
            grid->getIndexEJz(particlePosition, idxJz, internalCoordsJz);

            formFactor(internalCoordsJx);
            depositComponentCurrent(grid->Jx, idxJx, internalCoordsJx, current.x, formFactor);

            formFactor(internalCoordsJy);
            depositComponentCurrent(grid->Jy, idxJy, internalCoordsJy, current.y, formFactor);

            formFactor(internalCoordsJz);
            depositComponentCurrent(grid->Jz, idxJz, internalCoordsJz, current.z, formFactor);
        }
    };

    typedef FirstOrderCurrentDeposition<YeeGrid> FirstOrderCurrentDepositionYee;
}