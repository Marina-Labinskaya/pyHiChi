#pragma once

#include "ParticleArray.h"
#include <algorithm>

template <>
inline void std::iter_swap<iteratorPArray<ParticleArray3d, ParticleArray3d::ParticleProxyType>,
                      iteratorPArray<ParticleArray3d, ParticleArray3d::ParticleProxyType>>
    (iteratorPArray<ParticleArray3d, ParticleArray3d::ParticleProxyType> it1,
    iteratorPArray<ParticleArray3d, ParticleArray3d::ParticleProxyType> it2) {
    pfc::Particle<pfc::Three> temp1 = pfc::Particle<pfc::Three>(*it1);
    pfc::Particle<pfc::Three> temp2 = pfc::Particle<pfc::Three>(*it2);
    it1.IteratorOnRawParticle(temp2);
    it2.IteratorOnRawParticle(temp1);
}

namespace pfc {
    template <typename ArrayType>
    class ParticleSorting {
        public:
        static void sortByX(ArrayType& particleArray) {
            std::sort(particleArray.begin(), particleArray.end(),
                    [](const typename ArrayType::ParticleType& first, const typename ArrayType::ParticleType& second) 
                    { return first.getPosition().x < second.getPosition().x; });
        }
        static void sortByY(ArrayType& particleArray) {
            std::sort(particleArray.begin(), particleArray.end(),
                    [](const typename ArrayType::ParticleType& first, const typename ArrayType::ParticleType& second)
                    { return first.getPosition().y < second.getPosition().y; });
        }
        static void sortByZ(ArrayType& particleArray) {
            std::sort(particleArray.begin(), particleArray.end(),
                    [](const typename ArrayType::ParticleType& first, const typename ArrayType::ParticleType& second)
                    { return first.getPosition().z < second.getPosition().z; });
        }
    };
}
