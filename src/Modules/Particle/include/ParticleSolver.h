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
        FP3 minCoords = grid->origin + grid->steps * grid->getNumExternalLeftCells();
        FP3 maxCoords = minCoords + grid->numInternalCells * grid->steps;
        FP3 AreaEnd = grid->origin + grid->numCells * grid->steps;
        for (int i = 0; i < particleArray->size(); i++)
        {
            ParticleProxyType particle = (*particleArray)[i];
            PositionType particlePosition = particle.getPosition();
            PositionType newParticlePosition = particle.getPosition();

            if (particlePosition.x < minCoords.x) {
                newParticlePosition.x = AreaEnd.x - (minCoords.x - particlePosition.x);
            }
            else if ((particlePosition.x > maxCoords.x) && ((particlePosition.x - maxCoords.x) > FP(0.5) * grid->steps.x)) {
                newParticlePosition.x = grid->origin.x + (particlePosition.x - maxCoords.x);
            }

            if (particlePosition.y < minCoords.y) {
                newParticlePosition.y = AreaEnd.y - (minCoords.y - particlePosition.y);
            }
            else if ((particlePosition.y > maxCoords.y) && ((particlePosition.y - maxCoords.y) > FP(0.5) * grid->steps.y)) {
                newParticlePosition.y = grid->origin.y + (particlePosition.y - maxCoords.y);
            }

            if (particlePosition.z < minCoords.z) {
                newParticlePosition.z = AreaEnd.z - (minCoords.z - particlePosition.z);
            }
            else if ((particlePosition.z > maxCoords.z) && ((particlePosition.z - maxCoords.z) > FP(0.5) * grid->steps.z)) {
                newParticlePosition.z = grid->origin.z + (particlePosition.z - maxCoords.z);
            }            

            /*while (particlePosition.x < minCoords.x)
                newParticlePosition.x = maxCoords.x - (minCoords.x - particlePosition.x);
            while (particlePosition.x > maxCoords.x)
                newParticlePosition.x = minCoords.x + (particlePosition.x - maxCoords.x);
            
            while (particlePosition.y < minCoords.y)
                newParticlePosition.y = maxCoords.y - (minCoords.y - particlePosition.y);
            while (particlePosition.y > maxCoords.y)
                newParticlePosition.y = minCoords.y + (particlePosition.y - maxCoords.y);
            
            while (particlePosition.z < minCoords.z)
                newParticlePosition.z = maxCoords.z - (minCoords.z - particlePosition.z);
            while (particlePosition.z > maxCoords.z)
                newParticlePosition.z = minCoords.z + (particlePosition.z - maxCoords.z);*/


            particle.setPosition(newParticlePosition);
        }
    }
}
