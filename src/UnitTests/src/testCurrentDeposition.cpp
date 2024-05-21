#include "TestingUtility.h"

#include "CurrentDeposition.h"
#include "CurrentBoundaries.h"
#include "Constants.h"
#include "FDTD.h"
#include "FieldGenerator.h"
#include "FieldValue.h"
#include "Particle.h"
#include "ParticleGenerator.h"
#include "ParticleSolver.h"
#include "Pusher.h"
#include <fstream>
#include <iostream>
#include <chrono>
#include <algorithm>

using namespace pfc;

TEST(FormfactorTest, TSCFormFactorTest) {
    int nx = 5, ny = 5, nz = 5;
    FP L = 1.0;
    FP dx = L / nx, dy = dx, dz = dx;
    FP dt = 1e-8;

    CurrentDepositionTSC<YeeGrid> currentDeposition(dt);

    Int3 numInternalCells = Int3(nx, ny, nz);
    FP3 minCoords = FP3(0.0, 0.0, 0.0);
    FP3 steps = FP3(dx, dy, dz);
    Int3 globalGridDims = numInternalCells;

    YeeGrid grid(numInternalCells, minCoords, steps, globalGridDims);

    for (int i = 0; i < grid.numCells.x; ++i)
        for (int j = 0; j < grid.numCells.y; ++j)
            for (int k = 0; k < grid.numCells.z; ++k) {
                grid.Ex(i, j, k) = 0.0;
                grid.Ey(i, j, k) = 0.0;
                grid.Ez(i, j, k) = 0.0;
                grid.Bx(i, j, k) = 0.0;
                grid.By(i, j, k) = 0.0;
                grid.Bz(i, j, k) = 0.0;
            }

    Particle3d particle(FP3(0.25, 0.25, 0.25), FP3(0,0,0));
    particle.setVelocity(FP3(1000, 0, 0));

    FP3 expect_current = (constants::c * particle.getCharge() * particle.getWeight() / particle.getGamma())
        * (particle.getP() / grid.steps.volume());
    Int3 idxJx, idxJy, idxJz;
    FP3 internalCoordsJx, internalCoordsJy, internalCoordsJz;
    grid.getIndexEJx(particle.getPosition(), idxJx, internalCoordsJx);
    grid.getIndexEJy(particle.getPosition(), idxJy, internalCoordsJy);
    grid.getIndexEJz(particle.getPosition(), idxJz, internalCoordsJz);
    FP c[3] = {0.03125, 0.6875, 0.28125};
    FP expect_Jx[3][3][3], expect_Jy[3][3][3], expect_Jz[3][3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                expect_Jx[i][j][k] = 0.0;
                expect_Jy[i][j][k] = 0.0;
                expect_Jz[i][j][k] = 0.0;
            }
        }
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                expect_Jx[i][j][k] = c[0] * c[1] * c[2] * expect_current.x;
            }
        }
    }
    currentDeposition(&grid, particle);
    double delta = 1e-4;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                idxJx = idxJx + Int3(i - 1, j - 1, k - 1);
                idxJy = idxJy + Int3(i - 1, j - 1, k - 1);
                idxJz = idxJz + Int3(i - 1, j - 1, k - 1);
                EXPECT_NEAR(expect_Jx[i][j][k], grid.getJxTSC(grid.JxPosition(idxJx.x, idxJx.y, idxJx.z)), delta);
                EXPECT_NEAR(expect_Jy[i][j][k], grid.getJyTSC(grid.JyPosition(idxJy.x, idxJy.y, idxJy.z)), delta);
                EXPECT_NEAR(expect_Jz[i][j][k], grid.getJzTSC(grid.JzPosition(idxJz.x, idxJz.y, idxJz.z)), delta);
            }
        }
    }
}


TEST(FormfactorTest, TSCSumIsUnit) {

    FP stepx = 0.1, stepy = 0.2, stepz = 0.1;    
    Int3 numInternalCells = Int3(10, 7, 12);
    FP3 minCoords = FP3(0.0, 0.0, 0.0);
    FP3 steps = FP3(stepx, stepy, stepz);
    Int3 globalGridDims = numInternalCells;

    Int3 idxJx, idxJy, idxJz;
    FP3 internalCoordsX, internalCoordsY, internalCoordsZ;
    
    YeeGrid grid(numInternalCells, minCoords, steps, globalGridDims);
    FP x = (FP)0.3 * stepx;
    FP y = (FP)0.1 * stepy;
    FP z = (FP)0.4 * stepz;

    Particle3d::PositionType position(x, y, z);
    grid.getIndexEJx(position, idxJx, internalCoordsX);
    grid.getIndexEJy(position, idxJy, internalCoordsY);
    grid.getIndexEJz(position, idxJz, internalCoordsZ);
    Int3 middleIdxX = truncate(internalCoordsX + FP3(0.5, 0.5, 0.5));
    Int3 middleIdxY = truncate(internalCoordsY + FP3(0.5, 0.5, 0.5));
    Int3 middleIdxZ = truncate(internalCoordsZ + FP3(0.5, 0.5, 0.5));
    FormFactorTSC formFactor;
    double sumX = 0, sumY = 0, sumZ = 0;

    formFactor(internalCoordsX - FP3(middleIdxX));
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                sumX += formFactor.c[0][i] * formFactor.c[1][j] * formFactor.c[2][k];
    }

    formFactor(internalCoordsY - FP3(middleIdxY));
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                sumY += formFactor.c[0][i] * formFactor.c[1][j] * formFactor.c[2][k];
    }

    formFactor(internalCoordsZ - FP3(middleIdxZ));
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                sumZ += formFactor.c[0][i] * formFactor.c[1][j] * formFactor.c[2][k];
    }
    FP maxError = (FP)1e-14;
    EXPECT_NEAR(sumX, 1.0, maxError);
    EXPECT_NEAR(sumY, 1.0, maxError);
    EXPECT_NEAR(sumZ, 1.0, maxError);
}

static FP testCurrentDepositionTestPlasmaOscillationTestParticleDensityFunc(FP x, FP y, FP z) {
    FP T0 = 0.01 * constants::electronMass * constants::c * constants::c;
    FP L = 1.0;
    FP dx = L / 64.0; //remember about nx=64
    FP D = T0 / (8 * constants::pi * constants::electronCharge * constants::electronCharge * 0.25 * dx * dx);
    return D * (1 + 0.05 * std::sin(2 * constants::pi * x / L));
}

static FP testCurrentDepositionTestPlasmaOscillationTestInitialTemperatureFunc(FP x, FP y, FP z) {
    FP T0 = 0.01 * constants::electronMass * constants::c * constants::c;
    return T0;
}

static FP3 testCurrentDepositionTestPlasmaOscillationTestInitialMomentumFunc(FP x, FP y, FP z) {
    return FP3(0.0, 0.0, 0.0);
}

TEST(CurrentDepositionTest, PlasmaOscillationTest) {
    int nx = 64, ny = 8, nz = 8;
    FP L = 1.0;
    FP dx = L / nx, dy = dx, dz = dx;
    int Nip = 256;
    int Np = 100;// nx / (std::sqrt(2) * constants::pi);  // number of periods
    int N = Nip * Np;
    FP T0 = 0.01 * constants::electronMass * constants::c * constants::c;
    FP L_Debay = dx * 0.5;
    FP D = T0 / (8 * constants::pi * constants::electronCharge * L_Debay * constants::electronCharge * L_Debay);
    FP wp = sqrt(4.0 * constants::pi * constants::electronCharge * constants::electronCharge * D
        / constants::electronMass);
    FP dt = 2 * (constants::pi / wp) / Nip;
    FP Nc = 100;
    Particle3d::WeightType w = D * dx * dx * dx / Nc;
    FP A = 0.05;
    FP3 p0 = FP3(0.0, 0.0, 0.0);

    BorisPusher scalarPusher;

    CurrentDepositionTSC<YeeGrid> currentDeposition(dt);

    Int3 numInternalCells = Int3(nx, ny, nz);
    FP3 minCoords = FP3(0.0, 0.0, 0.0);
    FP3 steps = FP3(dx, dy, dz);
    Int3 globalGridDims = numInternalCells;

    YeeGrid grid(numInternalCells, minCoords, steps, globalGridDims);

    for (int i = 0; i < grid.numCells.x; ++i)
        for (int j = 0; j < grid.numCells.y; ++j)
            for (int k = 0; k < grid.numCells.z; ++k) {
                grid.Ex(i, j, k) = -2 * L * D * A * constants::electronCharge 
                    * std::cos(2 * constants::pi * grid.getBaseCoords(i, j, k).x / L);
                grid.Ey(i, j, k) = 0.0;
                grid.Ez(i, j, k) = 0.0;
                grid.Bx(i, j, k) = 0.0;
                grid.By(i, j, k) = 0.0;
                grid.Bz(i, j, k) = 0.0;
            }

    ParticleArray3d particleArray;

    PeriodicalParticleBoundaryConditions particleSolver;
    FDTD fdtd(&grid, dt);
    PeriodicalCurrentBoundaryConditions<YeeGridType> periodicalCurrentBoundary(&fdtd);
    PeriodicalFieldGeneratorYee generator(&fdtd);
    fdtd.setFieldGenerator(&generator);

    ParticleGenerator particleGenerator;
    particleGenerator(&particleArray, &grid,
        testCurrentDepositionTestPlasmaOscillationTestParticleDensityFunc,
        testCurrentDepositionTestPlasmaOscillationTestInitialTemperatureFunc,
        testCurrentDepositionTestPlasmaOscillationTestInitialMomentumFunc,
        w, ParticleTypes::Electron
    );
    std::cout << "n= " << particleArray.size() << std::endl;

    std::ofstream fout("OscillationTestEx.txt");
    std::ofstream fout2("OscillationTestElectronDensity.txt");
    std::ofstream fout_energy("FieldEnergy.txt");
    
    /*std::ofstream fout_sorted_particles("sorted_particles.txt");
    for (int i = 0; i < particleArray.size(); ++i) {
        fout_sorted_particles << particleArray[i].getPosition().x << " " << particleArray[i].getPosition().y
            << " " << particleArray[i].getPosition().z << std::endl;
    }
    fout_sorted_particles.close();*/

    double current_time = 0.0;
    double fdtd_time = 0.0;
    double interpolate_time = 0.0;
    double pusher_time = 0.0;

    //-------------------------------------------initial-field-energy--------------------------------------------------
    double energy = 0.0, energy_particles = 0.0;

    for (int i = 0; i < N; ++i) {
        if (i % 16 == 1) {
            //-------------------------------current-energy--------------------------------------------------------------
            energy = 0.0; energy_particles = 0.0;
            Int3 numExternalLeftCells = grid.getNumExternalLeftCells();
            Int3 size = grid.numInternalCells + numExternalLeftCells;
            for (int i = numExternalLeftCells.x; i < size.x; i++)
                for (int j = numExternalLeftCells.y; j < size.y; j++)
                    for (int k = numExternalLeftCells.z; k < size.z; k++)
                        energy += (grid.Ex(i, j, k) * grid.Ex(i, j, k)
                            + grid.Ey(i, j, k) * grid.Ey(i, j, k)
                            + grid.Ez(i, j, k) * grid.Ez(i, j, k)
                            + grid.Bx(i, j, k) * grid.Bx(i, j, k)
                            + grid.By(i, j, k) * grid.By(i, j, k)
                            + grid.Bz(i, j, k) * grid.Bz(i, j, k));
            energy = energy * grid.steps.volume() * 1e-7 / ((FP)8 * constants::pi);
            fout_energy << i << " " << energy << std::endl;
            for (int j = 0; j < grid.numInternalCells.x; ++j) {
                Int3 idx; FP3 internalCoords;
                fout << i << " " << minCoords.x + j * grid.steps.x << " " << grid.Ex(j
                    + grid.getNumExternalLeftCells().x, grid.getNumExternalLeftCells().y,
                    grid.getNumExternalLeftCells().z) << std::endl;
            }

            for (int k = 0; k < grid.numInternalCells.x; ++k) {
                int electronCount = 0;
                FP minCellCoords = minCoords.x + k * grid.steps.x;
                FP maxCellCoords = minCoords.x + (k + 1) * grid.steps.x;

                for (int j = 0; j < particleArray.size(); ++j) {
                    if ((particleArray[j].getPosition().x >= minCellCoords)
                        && (particleArray[j].getPosition().x <= maxCellCoords))
                        electronCount++;
                }
                // to sinchronize the output with Picador
                FP electronDensity = electronCount * w;  // density through plane
                fout2 << i << " " << minCellCoords << " " << electronDensity << std::endl;
            }
        }
        auto interpolate_start = chrono::steady_clock::now();
        //interpolate
        std::vector<ValueField> fields;
        for (int i = 0; i < particleArray.size(); i++) {
            FP3 E, B;
            grid.getFieldsTSC(particleArray[i].getPosition(), E, B);
            fields.push_back(ValueField(E, B));
        }

        auto interpolate_end = chrono::steady_clock::now();
        double interpolate_local_time = (chrono::duration_cast<chrono::nanoseconds>(interpolate_end 
            - interpolate_start).count()) / 1e9;
        interpolate_time += interpolate_local_time;

        auto pusher_start = chrono::steady_clock::now();
        //pusher
        scalarPusher(&particleArray, fields, dt);

        auto pusher_end = chrono::steady_clock::now();
        double pusher_local_time = (chrono::duration_cast<chrono::nanoseconds>(pusher_end 
            - pusher_start).count()) / 1e9;
        pusher_time += pusher_local_time;

        //periodical particle position
        particleSolver.updateParticlePosition(&grid, &particleArray);
        // current deposition
        fdtd.updateDims();

        auto start = chrono::steady_clock::now();
        if (i % 100 == 0) {
            std::sort(particleArray.begin(), particleArray.end(),
                [](const Particle3d& first, const Particle3d& second) { return first.getPosition().x < second.getPosition().x; });
        }
        currentDeposition(&grid, &particleArray);
        auto end = chrono::steady_clock::now();
        double current_local_time = (chrono::duration_cast<chrono::nanoseconds>(end - start).count()) / 1e9;
        current_time += current_local_time;
        periodicalCurrentBoundary.updateCurrentBoundaries();

        auto fdtd_start = chrono::steady_clock::now();

        //fdtd
        fdtd.updateFields();

        auto fdtd_end = chrono::steady_clock::now();
        double fdtd_local_time = (chrono::duration_cast<chrono::nanoseconds>(fdtd_end - fdtd_start).count()) / 1e9;
        fdtd_time += fdtd_local_time;

        

    }
    fout.close();
    fout2.close();

    //-------------------------final-energy------------------------------------
    energy = 0.0;
    for (int i = grid.numInternalCells.x; i < grid.numCells.x; i++)
        for (int j = grid.numInternalCells.y; j < grid.numCells.y; j++)
            for (int k = grid.numInternalCells.z; k < grid.numCells.z; k++) {
                energy += grid.Ex(i, j, k) * grid.Ex(i, j, k);
                energy += grid.Ey(i, j, k) * grid.Ey(i, j, k);
                energy += grid.Ez(i, j, k) * grid.Ez(i, j, k);
                energy += grid.Bx(i, j, k) * grid.Bx(i, j, k);
                energy += grid.By(i, j, k) * grid.By(i, j, k);
                energy += grid.Bz(i, j, k) * grid.Bz(i, j, k);
            }
    energy = energy * dx * dy * dz * 1e-7 / ((FP)8 * constants::pi);
    std::cout << "in end: " << energy << std::endl;

    std::cout << "fdtd time= " << fdtd_time << std::endl;
    std::cout << "interpolate_time= " << interpolate_time << std::endl;
    std::cout << "pusher_time= " << pusher_time << std::endl;
    std::cout << "current_time= " << current_time << std::endl;
    ASSERT_NO_THROW(true);
}

//static FP OneParticleTestGetParticleVelocity(FP t) {
//    return std::sin(t) * std::sin(t);
//}

//TEST(CurrentDepositionTest, OneParticleTest) {
//    std::ofstream fout("testOneParticle.txt");
//    int nx = 100, ny = 100, nz = 1;
//    FP L = 1.0;
//    FP dx = L / nx, dy = L / ny, dz = L / nz;
//    FP dt = 1e-5;
//    double T = 1e-11; //sec
//
//    BorisPusher scalarPusher;
//    Int3 numInternalCells = Int3(nx, ny, nz);
//    FP3 minCoords = FP3(0.0, 0.0, 0.0);
//    FP3 steps = FP3(dx, dy, dz);
//    Int3 globalGridDims = numInternalCells;
//
//    YeeGrid grid(numInternalCells, minCoords, steps, globalGridDims);
//
//    for (int i = 0; i < grid.numInternalCells.x; ++i)
//        for (int j = 0; j < grid.numInternalCells.y; ++j)
//            for (int k = 0; k < grid.numInternalCells.z; ++k) {
//                grid.Ex(i, j, k) = 1e-10;
//                grid.Ey(i, j, k) = 1e-10;
//                grid.Ez(i, j, k) = 0.0;
//                grid.Bx(i, j, k) = 0.0;
//                grid.By(i, j, k) = 0.0;
//                grid.Bz(i, j, k) = 0.0;
//            }
//
//    ParticleArray3d particles;
//    Particle3d particle;
//    particle.setVelocity(FP3(0, 0, 0));
//    particle.setPosition(FP3(L / 2.0 + dx / 2.0, L / 2.0 + dy / 2.0, L / 2 + dz / 2.0));
//    particles.pushBack(particle);
//    double alpha = constants::pi / 5.0;
//
//
//    PeriodicalParticleBoundaryConditions particleSolver;
//    FDTD fdtd(&grid, dt);
//    std::cout << "dt= " << fdtd.dt << std::endl;
//    dt = fdtd.dt;
//    CurrentDepositionTSC<YeeGrid> currentDeposition(dt);
//    PeriodicalCurrentBoundaryConditions<YeeGridType> periodicalCurrentBoundary(&fdtd);
//    PeriodicalFieldGeneratorYee generator(&fdtd);
//    fdtd.setFieldGenerator(&generator);
//
//    for (double t = 0; t < T; t += dt) {
//        if (particles[0].getVelocity().norm() < 0.3 * constants::c) {
//            particles[0].setVelocity(FP3(OneParticleTestGetParticleVelocity(t) * std::cos(alpha),
//                OneParticleTestGetParticleVelocity(t) * std::sin(alpha),
//                0));
//        }
//        else {
//            particles[0].setVelocity(FP3(0, 0, 0));
//        }
//        
//        //fdtd
//        fdtd.updateFields();
//        //interpolate
//        std::vector<ValueField> fields;
//        FP3 E, B;
//        grid.getFieldsTSC(particles[0].getPosition(), E, B);
//        fields.push_back(ValueField(E, B));
//        //pusher
//        scalarPusher(&particles, fields, dt);
//        //periodical particle position
//        particleSolver.updateParticlePosition(&grid, &particles);
//        // current deposition
//        fdtd.updateDims();
//        currentDeposition(&grid, &particles);
//        periodicalCurrentBoundary.updateCurrentBoundaries();
//    }
//
//    for (int i = 0; i < nx; ++i) {
//        for (int j = 0; j < ny; ++j) {
//            FP3 coords = FP3(i + grid.getNumExternalLeftCells().x, j + grid.getNumExternalLeftCells().y, grid.getNumExternalLeftCells().z);
//            FP norm2 = grid.Ex(coords) * grid.Ex(coords) + grid.Ey(coords) * grid.Ey(coords) + grid.Ez(coords) * grid.Ez(coords);
//            fout << std::sqrt(norm2) << " ";
//        }
//        fout << std::endl;
//    }
//    fout.close();
//    ASSERT_NO_THROW(true);
//}