#pragma once

#include "Dimension.h"
#include "Grid.h"
#include "Particle.h"
#include "Species.h"
#include <cmath>

namespace pfc {
    class ParticleSolver {
    public:
        template <class T_ParticleArray, GridTypes gridType>
        void updateParticlePosition(Grid<FP, gridType> * grid, T_ParticleArray* particleArray) {};

    };

    
    class PeriodicalParticleSolver : ParticleSolver {
    public:
        template <class T_ParticleArray, GridTypes gridType>
        void updateParticlePosition(Grid<FP, gridType> * grid, T_ParticleArray* particleArray);
    };

    template<class T_ParticleArray, GridTypes gridType>
    void PeriodicalParticleSolver::updateParticlePosition(Grid<FP, gridType> * grid, T_ParticleArray* particleArray)
    {
        typedef typename T_ParticleArray::ParticleProxyType ParticleProxyType;
        typedef typename VectorTypeHelper<Dimension::Three, Real>::Type PositionType;
        for (int i = 0; i < particleArray->size(); i++)
        {
            ParticleProxyType particle = (*particleArray)[i];
            PositionType particlePosition = particle.getPosition();
            FP3 minCoords = grid->origin + grid->steps * grid->getNumExternalLeftCells();
            FP3 maxCoords = minCoords + grid->numInternalCells * grid->steps;

            if (particlePosition.x < minCoords.x)
                particlePosition.x += std::ceil((abs(particlePosition.x - minCoords.x) / (grid->numInternalCells.x * grid->steps.x))) * (grid->numInternalCells.x * grid->steps.x);
            if (particlePosition.x > maxCoords.x)
                particlePosition.x -= std::floor((abs(particlePosition.x - minCoords.x) / (grid->numInternalCells.x * grid->steps.x))) * (grid->numInternalCells.x * grid->steps.x);

            if (particlePosition.y < minCoords.y)
                particlePosition.y += std::ceil((abs(particlePosition.y - minCoords.y) / (grid->numInternalCells.y * grid->steps.y))) * (grid->numInternalCells.y * grid->steps.y);
            if (particlePosition.y > maxCoords.y)
                particlePosition.y -= std::floor((abs(particlePosition.y - minCoords.y) / (grid->numInternalCells.y * grid->steps.y))) * (grid->numInternalCells.y * grid->steps.y);

            if (particlePosition.z < minCoords.z)
                particlePosition.z += std::ceil((abs(particlePosition.z - minCoords.z) / (grid->numInternalCells.z * grid->steps.z))) * (grid->numInternalCells.z * grid->steps.z);
            if (particlePosition.z > maxCoords.z)
                particlePosition.z -= std::floor((abs(particlePosition.z - minCoords.z) / (grid->numInternalCells.z * grid->steps.z))) * (grid->numInternalCells.z * grid->steps.z);

            particle.setPosition(particlePosition);
        }
    }
}