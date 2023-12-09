#pragma once

#include "Constants.h"
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

        CurrentDeposition(double _dt) : dt(_dt) {}

        template<class T_Particle>
        void operator()(TGrid* grid, const T_Particle& particle,
            CurrentDeposition::ZeroizeJ UsingZeroizeJ = CurrentDeposition::ZeroizeJ::NOT_USE_ZEROIZEJ) {
            if (UsingZeroizeJ == ZeroizeJ::USE_ZEROIZEJ)
                grid->zeroizeJ();
            static_cast<DerivedClass*>(this)->depositOneParticle(grid, &particle);
        }

        template<class T_ParticleArray>
        void operator()(TGrid* grid, T_ParticleArray* particleArray) {
            typedef typename T_ParticleArray::ParticleProxyType ParticleProxyType;
            grid->zeroizeJ();

            for (int i = 0; i < particleArray->size(); i++) {
                ParticleProxyType particle = (*particleArray)[i];
                static_cast<DerivedClass*>(this)->depositOneParticle(grid, &particle);
            }
        }


        template<class T_Particle>
        void depositOneParticle(TGrid* grid, T_Particle* particle) {
            static_assert(false, "ERROR: CurrentDeposition::depositOnePaticle shouldn't be called");
        }

    protected:
        double dt;
    };

    template<class TGrid>
    class CurrentDepositionCIC : public CurrentDeposition<TGrid, CurrentDepositionCIC<TGrid>>
    {
    public:
        CurrentDepositionCIC(double _dt) : CurrentDeposition(_dt) {}

        void depositComponentCurrent(ScalarField<FP> & field, const Int3 & idx,
            const FP & fieldBeforeDeposition, FormFactorCIC& formFactor)
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
        void depositOneParticle(TGrid* grid, T_Particle* particle)
        {
            Int3 idxJx, idxJy, idxJz;
            FP3 internalCoordsJx, internalCoordsJy, internalCoordsJz;
            FP3 particlePosition = particle->getPosition(); // - (particle->getVelocity() * dt / 2.0);

            FP3 current = (particle->getVelocity() * particle->getCharge() * particle->getWeight()) /
                grid->steps.volume();

            grid->getIndexEJx(particlePosition, idxJx, internalCoordsJx);
            grid->getIndexEJy(particlePosition, idxJy, internalCoordsJy);
            grid->getIndexEJz(particlePosition, idxJz, internalCoordsJz);

            FormFactorCIC formFactor;
            formFactor(internalCoordsJx);
            depositComponentCurrent(grid->Jx, idxJx, current.x, formFactor);

            formFactor(internalCoordsJy);
            depositComponentCurrent(grid->Jy, idxJy, current.y, formFactor);

            formFactor(internalCoordsJz);
            depositComponentCurrent(grid->Jz, idxJz, current.z, formFactor);
        }
    };

    template<class TGrid>
    class CurrentDepositionTSC : public CurrentDeposition<TGrid, CurrentDepositionTSC<TGrid>>
    {
    public:
        CurrentDepositionTSC(double _dt) : CurrentDeposition(_dt) {}

        void depositComponentCurrent(ScalarField<FP> & field, const Int3 & idx,
            const FP & fieldBeforeDeposition, FormFactorTSC& formFactor)
        {
            for (int i = -1; i <= 1; i++) {
                for (int j = -1; j <= 1; j++) {
                    for (int k = -1; k <= 1; k++) {
                        field(idx.x + i, idx.y + j, idx.z + k) += (formFactor.c[0][i + 1] * formFactor.c[1][j + 1] *
                            formFactor.c[2][k + 1]) * fieldBeforeDeposition;
                    }
                }
            }
        }

        template<class T_Particle>
        void depositOneParticle(TGrid* grid, T_Particle* particle)
        {
            Int3 idxJx, idxJy, idxJz;
            FP3 internalCoordsJx, internalCoordsJy, internalCoordsJz;
            FP3 particlePosition = particle->getPosition();// -(particle->getVelocity() * dt / 2.0);
            FP3 current = (constants::c * particle->getCharge() * particle->getWeight() / particle->getGamma())
                * (particle->getP() / grid->steps.volume());

            grid->getIndexEJx(particlePosition, idxJx, internalCoordsJx);
            grid->getIndexEJy(particlePosition, idxJy, internalCoordsJy);
            grid->getIndexEJz(particlePosition, idxJz, internalCoordsJz);

            FormFactorTSC formFactor;
            formFactor(internalCoordsJx);
            depositComponentCurrent(grid->Jx, idxJx, current.x, formFactor);

            formFactor(internalCoordsJy);
            depositComponentCurrent(grid->Jy, idxJy, current.y, formFactor);

            formFactor(internalCoordsJz);
            depositComponentCurrent(grid->Jz, idxJz, current.z, formFactor);
        }
    };
}
