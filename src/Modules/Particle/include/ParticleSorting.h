#pragma once

#include "ParticleArray.h"
#include <algorithm>
#include <vector>

namespace pfc {
    template <class TGrid, class T_ParticleArray>
    class ParticleSortingOMP {
    public:
        ParticleSorting(TGrid* _grid, T_ParticleArray& particleArray): 
            nn(_grid->numCells.x * _grid->numCells.y * _grid->numCells.z) {
            this->grid = _grid;
            h.resize(n_threads);
            for (int i = 0; i < n_threads; ++i) {
                h[i].resize(nn + 1);
            }

            int portion_size = size / n_threads;
            int tail = size - n_threads * portion_size;
            std::vector<int> portions(n_threads);
            for (int i = 0; i < n_threads; ++i) {
                portions[i] = (i < tail) ? portion_size + 1: portion_size;
            }

            for (int i = 0; i < n_threads; ++i) {
                sortParts[i] = ParticleSorting<TGrid, T_ParticleArray>();
            }
        const int n_threads = omp_get_max_threads();
        TGrid* grid;
        const int nn;
        std::vector<int> portions;
        std::vector<ParticleSorting<TGrid, T_ParticleArray>> sortParts(n_threads);
        
        static void sortOMP(TGrid* _grid, const typename T_ParticleArray::iterator& begin,  const T_ParticleArray::iterator& end) {
            size_t size = end - begin;
                std::vector<ParticleSorting<TGrid, T_ParticleArray>> sortParts(n_threads);
    #pragma omp parallel for
                for (int64_t i = 0; i < n_threads; ++i) {
                    if (i < tail) {
                        portions[i] += 1;
                        n_skip[i] = i * portions[0];
                    }
                    else {
                        n_skip[i] = tail * portions[0] + (i - tail) * portions[tail];
                    }

                    rngs[i].skipahead(n_skip[i]);
                    rngs[i].generate(res + n_skip[i], portions[i]);
                }
            }
            else {
                generate(res, n_gen);
            }
            
            
#pragma omp parallel for
            for (int i = 0; i < particleArray->size(); i++) {
                int thread_num = omp_get_thread_num();
                ParticleProxyType particle = (*particleArray)[i];
                Int3 baseGridIdx = grid->getBaseIndex(
                    particle.getPosition() - (particle.getVelocity() * halfDt));
                if (Depositor[thread_num].getBaseGridIdx() != baseGridIdx) {
                    Depositor[thread_num].addCurrents();
                    Depositor[thread_num].setNewGridCell(baseGridIdx);
                }
                static_cast<DerivedClass*>(this)->depositOneParticle(grid, &particle, Depositor[thread_num]);
            }
            
        }
    };

    template <class TGrid, class ArrayType>
    class ParticleSorting {
    public:
        ParticleSorting(TGrid* _grid, ArrayType& particleArray): 
            nn(_grid->numCells.x * _grid->numCells.y * _grid->numCells.z) {
            this->grid = _grid;
            h.resize(nn + 1);
            int idx = 0;
            for (int i = 0; i < particleArray.size(); ++i) {
                FP3 position = particleArray[i].getPosition();
                int linIdx = linIndex(this->grid->getBaseIndex(position));
                while (linIdx >= idx) {
                    h[idx] = i;
                    idx++;
                    
                }
            }
            for (int i = idx; i < h.size(); ++i) {
                h[i] = particleArray.size();
            }
        }

        void sortL(ArrayType& particleArray) {
            for (int n = 0; n < nn - 1; ++n) {
                sortDiap(h, particleArray, n, this->h[n], this->h[n + 1]);
            }   
        }
        
    private:
        TGrid* grid;
        const int nn;
        std::vector<int> h;


        void moveForward(std::vector<int>& h, ArrayType& particleArray, int u, int n, int m,
            const typename ArrayType::ParticleType& p) {
            int v = h[n + 1] - 1;
            while (m != n) {
                v = h[n + 1] - 1;
                particleArray.changeParticle(particleArray[v], u);
                h[n + 1] = v;
                u = v;
                n = n + 1;
            }
            particleArray.changeParticle(p, v);
        }

        void moveBack(std::vector<int>& h, ArrayType& particleArray, int u, int n, int m,
            const typename ArrayType::ParticleType& p) {
            int v = h[n];
            while (m != n) {
                v = h[n];
                particleArray.changeParticle(particleArray[v], u);
                h[n] = v + 1;
                u = v;
                n = n - 1;
            }
            particleArray.changeParticle(p, v);
        }

        void sortDiap(std::vector<int>& h, ArrayType& particleArray, int n, int u, int s) {
            while (u != s) {
                typename ArrayType::ParticleType p = typename ArrayType::ParticleType(particleArray[u]);
                int m = linIndex(this->grid->getBaseIndex(p.getPosition()));
                if (m == n) {
                    u += 1;
                } 
                else if (m < n) {
                    moveBack(h, particleArray, u, n, m, p);
                    u += 1;
                }
                else {
                    moveForward(h, particleArray, u, n, m, p);
                    s -= 1;
                }
            }
        }
        
        inline int linIndex(const Int3& idx) {
            return (idx.x * this->grid->numCells.y * this->grid->numCells.z +
                    idx.y * this->grid->numCells.z + idx.z);
        }
        inline Int3 cubeIdx(int idx) {
            int x = idx / (this->grid->numCells.y * this->grid->numCells.z);
            idx -= x * (this->grid->numCells.y * this->grid->numCells.z);
            int y = idx / this->grid->numCells.z;
            idx -= y * this->grid->numCells.z;
            return Int3(x, y, idx);
        }
    };
}
