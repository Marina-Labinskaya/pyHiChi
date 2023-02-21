#pragma once

#include "Dimension.h"
#include "Grid.h"
#include "Particle.h"
#include "Species.h"
#include <cmath>

namespace pfc {
    class ParticleBoundaryConditions {
    public:
        template <class T_ParticleArray, GridTypes gridType>
        void updateParticlePosition(Grid<FP, gridType> * grid, T_ParticleArray* particleArray) {};

    };

    
    class PeriodicalParticleBoundaryConditions : ParticleBoundaryConditions {
    public:
        template <class T_ParticleArray, GridTypes gridType>
        void updateParticlePosition(Grid<FP, gridType> * grid, T_ParticleArray* particleArray);
    };

    template<class T_ParticleArray, GridTypes gridType>
    void PeriodicalParticleBoundaryConditions::updateParticlePosition(Grid<FP, gridType> * grid, T_ParticleArray* particleArray)
    {
        typedef typename T_ParticleArray::ParticleProxyType ParticleProxyType;
        typedef typename VectorTypeHelper<Dimension::Three, Real>::Type PositionType;
        for (int i = 0; i < particleArray->size(); i++)
        {
            ParticleProxyType particle = (*particleArray)[i];
            PositionType particlePosition = particle.getPosition();
            FP3 minCoords = grid->origin + grid->steps * grid->getNumExternalLeftCells();
            FP3 maxCoords = minCoords + grid->numInternalCells * grid->steps;

            while (particlePosition.x < minCoords.x)
                particlePosition.x = maxCoords.x - (minCoords.x - particlePosition.x);
            while (particlePosition.x > maxCoords.x)
                particlePosition.x = minCoords.x + (particlePosition.x - maxCoords.x);
            
            while (particlePosition.y < minCoords.y)
                particlePosition.y = maxCoords.y - (minCoords.y - particlePosition.y);
            while (particlePosition.y > maxCoords.y)
                particlePosition.y = minCoords.y + (particlePosition.y - maxCoords.y);
            
            while (particlePosition.z < minCoords.z)
                particlePosition.z = maxCoords.z - (minCoords.z - particlePosition.z);
            while (particlePosition.z > maxCoords.z)
                particlePosition.z = minCoords.z + (particlePosition.z - maxCoords.z);

            particle.setPosition(particlePosition);
        }
    }
}
