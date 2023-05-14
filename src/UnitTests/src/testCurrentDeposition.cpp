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

using namespace pfc;

TEST(CurrentDepositionTest, CurrentDepositionForOneParticleIsRight) {
    FP stepx = 0.1, stepy = 0.2, stepz = 0.1;

    Int3 numInternalCells = Int3(10, 7, 12);
    FP3 minCoords = FP3(0.0, 0.0, 0.0);
    FP3 steps = FP3(stepx, stepy, stepz);
    Int3 globalGridDims = numInternalCells;

    Int3 idxX, idxY, idxZ;
    FP3 internalCoordsX, internalCoordsY, internalCoordsZ;
    Real charge = constants::electronCharge;

    YeeGrid grid(numInternalCells, minCoords, steps, globalGridDims);

    ScalarField<FP> expectJx(grid.numCells), expectJy(grid.numCells), expectJz(grid.numCells);

    expectJx.zeroize();
    expectJy.zeroize();
    expectJz.zeroize();

    FP x = (FP)0.3 * grid.steps.x;
    FP y = (FP)0.1 * grid.steps.y;
    FP z = (FP)0.4 * grid.steps.z;
    Particle3d::PositionType position(x, y, z);
    Particle3d::MomentumType momentum(FP3(15.0, 7.0, 10.0));
    Particle3d particle(position, momentum);

    FP3 JBefore = (particle.getVelocity() * particle.getCharge()) / (steps.x * steps.y * steps.z);

    FirstOrderCurrentDeposition<YeeGridType> currentdeposition;
    currentdeposition(&grid, particle);

    grid.getIndexEJxCoords(position, idxX, internalCoordsX);
    currentdeposition.CurrentDensityAfterDeposition(expectJx, idxX, internalCoordsX, JBefore.x);

    grid.getIndexEJyCoords(position, idxY, internalCoordsY);
    currentdeposition.CurrentDensityAfterDeposition(expectJy, idxY, internalCoordsY, JBefore.y);

    grid.getIndexEJzCoords(position, idxZ, internalCoordsZ);
    currentdeposition.CurrentDensityAfterDeposition(expectJz, idxZ, internalCoordsZ, JBefore.z);

    for (int i = 0; i < grid.numCells.x; ++i)
        for (int j = 0; j < grid.numCells.y; ++j)
            for (int k = 0; k < grid.numCells.z; ++k) {
                EXPECT_NEAR (grid.Jx(i, j, k), expectJx(i, j, k), 0.00001);
                EXPECT_NEAR (grid.Jy(i, j, k), expectJy(i, j, k), 0.00001);
                EXPECT_NEAR (grid.Jz(i, j, k), expectJz(i, j, k), 0.00001);
            }
}

TEST(CurrentDepositionTest, CurrentDepositionForParticleArrayIsRight) {
    FP stepX = 0.1, stepY = 0.2, stepZ = 0.1;
    FP minpX = 1.0, minpY = 1.0, minpZ = 1.0;
    FP maxpX = 10.0, maxpY = 10.0, maxpZ = 10.0;
    FP minX = 0.0, minY = 0.0, minZ = 0.0;
    FP maxX, maxY, maxZ;

    Int3 numInternalCells = Int3(10, 7, 12);
    FP3 minCoords = FP3(0.0, 0.0, 0.0);
    FP3 steps = FP3(stepX, stepY, stepZ);
    Int3 globalGridDims = numInternalCells;

    Int3 idxX, idxY, idxZ;
    FP3 internalCoordsX, internalCoordsY, internalCoordsZ;

    YeeGrid grid(numInternalCells, minCoords, steps, globalGridDims);
    maxX = minX + grid.numInternalCells.x * stepX;
    maxY = minY + grid.numInternalCells.y * stepY;
    maxZ = minZ + grid.numInternalCells.z * stepZ;

    ScalarField<FP> expectJx(grid.numCells), expectJy(grid.numCells), expectJz(grid.numCells);
    expectJx.zeroize();
    expectJy.zeroize();
    expectJz.zeroize();

    ParticleArray3d particleArray;
    int numParticles = 10;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> disX(minX, maxX);
    std::uniform_real_distribution<> disY(minY, maxY);
    std::uniform_real_distribution<> disZ(minZ, maxZ);

    std::uniform_real_distribution<> dispX(minpX, maxpX);
    std::uniform_real_distribution<> dispY(minpY, maxpY);
    std::uniform_real_distribution<> dispZ(minpZ, maxpZ);

    for (int i = 0; i < numParticles; i++)
    {
        Particle3d::PositionType position(disX(gen), disY(gen), disZ(gen));
        Particle3d::MomentumType momentum(dispX(gen), dispY(gen), dispZ(gen));
        Particle3d particle(position, momentum);
        particleArray.pushBack(particle);
    }

    FirstOrderCurrentDeposition<YeeGridType> currentdeposition;
    currentdeposition(&grid, &particleArray);

    for (int i = 0; i < numParticles; i++)
    {
        FP3 JBefore = (particleArray[i].getVelocity() * particleArray[i].getCharge()) / (steps.x * steps.y * steps.z);

        grid.getIndexEJxCoords(particleArray[i].getPosition(), idxX, internalCoordsX);
        currentdeposition.CurrentDensityAfterDeposition(expectJx, idxX, internalCoordsX, JBefore.x);

        grid.getIndexEJyCoords(particleArray[i].getPosition(), idxY, internalCoordsY);
        currentdeposition.CurrentDensityAfterDeposition(expectJy, idxY, internalCoordsY, JBefore.y);

        grid.getIndexEJzCoords(particleArray[i].getPosition(), idxZ, internalCoordsZ);
        currentdeposition.CurrentDensityAfterDeposition(expectJz, idxZ, internalCoordsZ, JBefore.z);
    }

    for (int i = 0; i < grid.numCells.x; ++i)
        for (int j = 0; j < grid.numCells.y; ++j)
            for (int k = 0; k < grid.numCells.z; ++k) {
                EXPECT_NEAR(grid.Jx(i, j, k), expectJx(i, j, k), 0.00001);
                EXPECT_NEAR(grid.Jy(i, j, k), expectJy(i, j, k), 0.00001);
                EXPECT_NEAR(grid.Jz(i, j, k), expectJz(i, j, k), 0.00001);
            }
}

TEST(CurrentDepositionTest, particleCanGoThroughACycle) {
    // create grid
    FP Ex = 0.0, Ey = 2.0, Ez = -1.0;
    FP Bx = 2.0, By = 1.0, Bz = 1.0;
    FP stepX = constants::c, stepY = constants::c, stepZ = constants::c;
    Int3 numInternalCells = Int3(10, 20, 10);
    FP3 minCoords = FP3(0.0, 0.0, 0.0);
    FP3 steps = FP3(stepX, stepY, stepZ);
    Int3 globalGridDims = numInternalCells;
    YeeGrid grid(numInternalCells, minCoords, steps, globalGridDims);

    for (int i = 0; i < grid.numCells.x; ++i)
        for (int j = 0; j < grid.numCells.y; ++j)
            for (int k = 0; k < grid.numCells.z; ++k) {
                grid.Ex(i, j, k) = Ex;
                grid.Ey(i, j, k) = Ey;
                grid.Ez(i, j, k) = Ez;
                grid.Bx(i, j, k) = Bx;
                grid.By(i, j, k) = By;
                grid.Bz(i, j, k) = Bz;
            }

    // field solver
    const FP dt = 0.001;
    RealFieldSolver<YeeGridType> realfieldsolver(&grid, dt, 0.0, 0.5 * dt, 0.5 * dt);
    PeriodicalFieldGenerator<YeeGridType> generator(&realfieldsolver);
    FDTD fdtd(&grid, dt);
    fdtd.setFieldGenerator(&generator);
    fdtd.updateFields();

    ParticleArray3d particleArray;
    int numParticles = 10;

    FP minpX = 1.0, minpY = 1.0, minpZ = 1.0;
    FP maxpX = 2.0, maxpY = 2.0, maxpZ = 2.0;
    FP minX = 0.0, minY = 0.0, minZ = 0.0;
    FP maxX = minX + grid.numInternalCells.x * stepX;
    FP maxY = minY + grid.numInternalCells.y * stepY;
    FP maxZ = minZ + grid.numInternalCells.z * stepZ;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> disX(minX, maxX);
    std::uniform_real_distribution<> disY(minY, maxY);
    std::uniform_real_distribution<> disZ(minZ, maxZ);

    std::uniform_real_distribution<> dispX(minpX, maxpX);
    std::uniform_real_distribution<> dispY(minpY, maxpY);
    std::uniform_real_distribution<> dispZ(minpZ, maxpZ);

    for (int i = 0; i < numParticles; i++)
    {
        Particle3d::PositionType position(disX(gen), disY(gen), disZ(gen));
        Particle3d::MomentumType momentum(dispX(gen), dispY(gen), dispZ(gen));
        Particle3d particle(position, momentum);
        particleArray.pushBack(particle);
    }

    // interpolate
    std::vector<ValueField> fields;
    for (int i = 0; i < numParticles; i++) {
        
        FP3 E, B;
        grid.getFieldsCIC(particleArray[i].getPosition(), E, B);
        fields.push_back(ValueField(E, B));
    }

    // pusher
    BorisPusher scalarPusher;
    FP timeStep = 0.01;
    scalarPusher(&particleArray, fields, timeStep);

    // current deposition
    FirstOrderCurrentDeposition<YeeGridType> currentdeposition;
    currentdeposition(&grid, &particleArray);

    ASSERT_NO_THROW(true);
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
    int Np = 2;  // numer of periods
    int N = Nip * Np;
    FP T0 = 0.01 * constants::electronMass * constants::c * constants::c;
    FP D = T0 / (8 * constants::pi * constants::electronCharge * constants::electronCharge * 0.25 * dx * dx);
    FP wp = sqrt(4.0 * constants::pi * constants::electronCharge * constants::electronCharge * D / constants::electronMass);
    FP dt = 2 * (constants::pi / wp) / Nip;
    FP Nc = 30;
    Particle3d::WeightType w = D * dx * dx * dx / Nc;
    FP A = 0.05;
    FP3 p0 = FP3(0.0, 0.0, 0.0);

    BorisPusher scalarPusher;

    FirstOrderCurrentDeposition<YeeGridType> currentDeposition;

    Int3 numInternalCells = Int3(nx, ny, nz);
    FP3 minCoords = FP3(0.0, 0.0, 0.0);
    FP3 steps = FP3(dx, dy, dz);
    Int3 globalGridDims = numInternalCells;

    YeeGrid grid(numInternalCells, minCoords, steps, globalGridDims);

    for (int i = 0; i < grid.numCells.x; ++i)
        for (int j = 0; j < grid.numCells.y; ++j)
            for (int k = 0; k < grid.numCells.z; ++k) {
                grid.Ex(i, j, k) = -2 * L * D * A * constants::electronCharge * std::cos(2 * constants::pi * grid.ExPosition(i, j, k).x / L);;
                grid.Ey(i, j, k) = 0.0;
                grid.Ez(i, j, k) = 0.0;
                grid.Bx(i, j, k) = 0.0;
                grid.By(i, j, k) = 0.0;
                grid.Bz(i, j, k) = 0.0;
            }

    ParticleArray3d particleArray;

    PeriodicalParticleBoundaryConditions particleSolver;  // maybe, PeriodicalParticleBoundaryConditions is better

    FDTD fdtd(&grid, dt);
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

    //std::cout << "minCoords:" << minCoords << std::endl;
    //std::cout << "maxCoords:" << minCoords + grid.numInternalCells * grid.steps << std::endl;
    //std::cout << "origin: " << grid.origin << std::endl;
    //std::cout << "area end: " << grid.origin + grid.numCells * grid.steps << std::endl;
    //std::cout << "step: " << grid.steps << std::endl;

    std::ofstream fout("OscillationTestEx.txt");
    std::ofstream fout2("OscillationTestElectronDensity.txt");

    for (int i = 0; i <= N; ++i) {
        if (i % 16 == 0)
            for (int j = 0; j < grid.numInternalCells.x; ++j) {
                Int3 idx; FP3 internalCoords;
                fout << i << " " << minCoords.x + j * grid.steps.x << " " << grid.Ex(j + grid.getNumExternalLeftCells().x, grid.getNumExternalLeftCells().y, grid.getNumExternalLeftCells().z) << std::endl;
                //fout << i << " " << minCoords.x + j * grid.steps.x << " " << grid.Ex(j, grid.getNumExternalLeftCells().y, grid.getNumExternalLeftCells().z) << std::endl;
            };
        if (i % 16 == 0) {
            for (int k = 0; k < grid.numInternalCells.x; ++k) {
                int electronCount = 0;
                FP minCellCoords = minCoords.x + k * grid.steps.x;
                FP maxCellCoords = minCoords.x + (k + 1) * grid.steps.x;
                //FP minCellCoords = grid.origin.x + k * grid.steps.x;
                //FP maxCellCoords = grid.origin.x + (k + 1) * grid.steps.x;

                for (int j = 0; j < particleArray.size(); ++j) {
                    if ((particleArray[j].getPosition().x >= minCellCoords) && (particleArray[j].getPosition().x <= maxCellCoords))
                        electronCount++;
                }
                // to sinchronize the output with Picador
                FP electronDensity = electronCount * w;  // density through plane
                if (i % 16 == 0)
                    fout2 << i << " " << minCellCoords << " " << electronDensity << std::endl;
            }
        }
        //fdtd
        fdtd.updateFields();
        //interpolate
        std::vector<ValueField> fields;
        for (int i = 0; i < particleArray.size(); i++) {
            FP3 E, B;
            grid.getFieldsCIC(particleArray[i].getPosition(), E, B);
            fields.push_back(ValueField(E, B));
        }
        //pusher
        scalarPusher(&particleArray, fields, dt);
        //periodical particle position
        particleSolver.updateParticlePosition(&grid, &particleArray, dt);
        // current deposition
        fdtd.updateDims();
        //realfieldsolver.updateInternalDims();
        PeriodicalCurrentBoundaryConditions<YeeGridType> periodicalCurrentBoundary(&fdtd);
        currentDeposition(&grid, &particleArray);
        periodicalCurrentBoundary.updateCurrentBoundaries();
    }
    fout.close();
    fout2.close();
    ASSERT_NO_THROW(true);
}
