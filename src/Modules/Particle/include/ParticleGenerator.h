#pragma once

#include "Dimension.h"
#include "Particle.h"
#include "VectorsProxy.h"
#include <random>

namespace pfc {
    class ParticleGenerator {
        typedef typename VectorTypeHelper<Three, Real>::Type PositionType;
        typedef typename VectorTypeHelper<Three, Real>::Type MomentumType;
        typedef ParticleTypes TypeIndexType;
        typedef Real WeightType;
     public:

        template<class T_ParticleArray>
        void operator()(T_ParticleArray* particleArray, FP(*f)(FP, FP, FP), FP T0,
            const FP3& minCoords, const FP3& maxCoords, int counter = 1, const MomentumType& momentum = FP3(0, 0, 0), WeightType weight = 1,
            TypeIndexType typeIndex = ParticleTypes::Electron);
    private:
        FP3 GetParticleRandomPosition(FP(*f)(FP, FP, FP), const FP3& minCoord, const FP3& maxCoord) {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> rd_posX(minCoord.x, maxCoord.x);
            std::uniform_real_distribution<> rd_posY(minCoord.y, maxCoord.y);
            std::uniform_real_distribution<> rd_posZ(minCoord.z, maxCoord.z);
            FP3 pos(rd_posX(gen), rd_posY(gen), rd_posZ(gen));
            return pos;
        }
        MomentumType GetParticleRandomMomentum(const MomentumType& momentum, FP T0) {
            FP alpha = 1.5 / (std::pow(T0 / (constants::electronMass * constants::c * constants::c) + 1, 2) - 1);
            FP sigma = std::sqrt(0.5 / alpha) * constants::electronMass * constants::c;
            std::random_device rd;
            std::mt19937 gen(rd());
            std::normal_distribution<> rd_momentum(0.0, 1.0);
            MomentumType newMomentum(momentum.x + sigma*rd_momentum(gen), momentum.y + sigma*rd_momentum(gen), momentum.z + sigma*rd_momentum(gen));
            return newMomentum;
        }
    };

    template<class T_ParticleArray>
    void ParticleGenerator::operator()(T_ParticleArray* particleArray, FP(*density)(FP, FP, FP),
        FP T0 /*temperature is a function of coordinates in general*/,
        const FP3& minCoords, const FP3& maxCoords, int averageCount, const MomentumType& momentum, WeightType weight,
        TypeIndexType typeIndex) {
        FP3 center = (minCoords + maxCoords) * 0.5;
        // number of particle in cell according to density
        FP expectedParticleNum = density(center.x, center.y, center.z) * (maxCoords - minCoords).volume() / weight;

        int particleNum = int(expectedParticleNum);
        // random shift of expectedParticleNum by 1 to eliminate the effect of rounding
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> uniformDist(0.0, 1.0);
        FP randomNumber = uniformDist(gen);
        if (randomNumber < expectedParticleNum - (FP)particleNum)
            ++particleNum;

        for (int i = 0; i < particleNum; ++i) {
            PositionType ParticlePosition(GetParticleRandomPosition(density, minCoords, maxCoords));
            Particle3d NewParticle(ParticlePosition, GetParticleRandomMomentum(momentum, T0), weight, typeIndex);
            particleArray->pushBack(NewParticle);
        }
    }
}
