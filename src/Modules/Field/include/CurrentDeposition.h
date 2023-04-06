#pragma once

#include "Grid.h"
#include "ParticleArray.h"
#include "Particle.h"
#include <iostream>

namespace pfc 
{
    enum ZeroizeJ {
        NOT_USE_ZEROIZEJ, USE_ZEROIZEJ
    };

    template<GridTypes gridType>
    class CurrentDeposition
    {
    public:
        template<class T_Particle>
        void operator()(Grid<FP, gridType>* grid, const T_Particle& particle, pfc::ZeroizeJ UsingZeroizeJ = NOT_USE_ZEROIZEJ) {};

        template<class T_ParticleArray>
        void operator()(Grid<FP, gridType>* grid, T_ParticleArray* particleArray) {};
    };

    template<GridTypes gridType>
    class FirstOrderCurrentDeposition : public CurrentDeposition<gridType>
    {
    public:
        template<class T_Particle>
        void operator()(Grid<FP, gridType>* grid, const T_Particle& particle, pfc::ZeroizeJ UsingZeroizeJ = NOT_USE_ZEROIZEJ) {
            if (UsingZeroizeJ == USE_ZEROIZEJ)
                grid->zeroizeJ();
            oneParticleDeposition(grid, &particle);
        }

        template<class T_ParticleArray>
        void operator()(Grid<FP, gridType>* grid, T_ParticleArray* particleArray) {
            typedef typename T_ParticleArray::ParticleProxyType ParticleProxyType;

            grid->zeroizeJ();

            for (int i = 0; i < particleArray->size(); i++) {
                ParticleProxyType particle = (*particleArray)[i];
                oneParticleDeposition(grid, &particle);
            }
        }

        void CurrentDensityAfterDeposition(ScalarField<FP> & field, const Int3 & idx, const FP3 & internalCoords, const FP & fieldBeforeDeposition)
        {
            //std::cout << "field coords: " << idx.x << ", " << idx.y << ", " << idx.z << std::endl;
            field(idx.x, idx.y, idx.z) += ((FP)1 - internalCoords.x) * ((FP)1 - internalCoords.y) * ((FP)1 - internalCoords.z) * fieldBeforeDeposition;
            field(idx.x + 1, idx.y, idx.z) += internalCoords.x * ((FP)1 - internalCoords.y) * ((FP)1 - internalCoords.z) * fieldBeforeDeposition;
            field(idx.x, idx.y + 1, idx.z) += ((FP)1 - internalCoords.x) * internalCoords.y * ((FP)1 - internalCoords.z) * fieldBeforeDeposition;
            field(idx.x, idx.y, idx.z + 1) += ((FP)1 - internalCoords.x) * ((FP)1 - internalCoords.y) * internalCoords.z * fieldBeforeDeposition;
            field(idx.x + 1, idx.y + 1, idx.z) += internalCoords.x * internalCoords.y * ((FP)1 - internalCoords.z) * fieldBeforeDeposition;
            field(idx.x + 1, idx.y, idx.z + 1) += internalCoords.x * ((FP)1 - internalCoords.y) * internalCoords.z * fieldBeforeDeposition;
            field(idx.x, idx.y + 1, idx.z + 1) += ((FP)1 - internalCoords.x) * internalCoords.y * internalCoords.z * fieldBeforeDeposition;
            field(idx.x + 1, idx.y + 1, idx.z + 1) += internalCoords.x * internalCoords.y * internalCoords.z * fieldBeforeDeposition;
        }

    private:

        template<class T_Particle>
        void oneParticleDeposition(Grid<FP, gridType>* grid, T_Particle* particle)
        {

            Int3 idxJx, idxJy, idxJz;
            FP3 internalCoordsJx, internalCoordsJy, internalCoordsJz;
            FP3 particlePosition = particle->getPosition();

            FP3 JBeforeDeposition = (particle->getVelocity() * particle->getCharge() * particle->getWeight()) / (grid->steps.x * grid->steps.y * grid->steps.z);

            grid->getIndexEJxCoords(particlePosition, idxJx, internalCoordsJx);
            grid->getIndexEJyCoords(particlePosition, idxJy, internalCoordsJy);
            grid->getIndexEJzCoords(particlePosition, idxJz, internalCoordsJz);

            CurrentDensityAfterDeposition(grid->Jx, idxJx, internalCoordsJx, JBeforeDeposition.x);
            CurrentDensityAfterDeposition(grid->Jy, idxJy, internalCoordsJy, JBeforeDeposition.y);
            CurrentDensityAfterDeposition(grid->Jz, idxJz, internalCoordsJz, JBeforeDeposition.z);
        }
    };
}