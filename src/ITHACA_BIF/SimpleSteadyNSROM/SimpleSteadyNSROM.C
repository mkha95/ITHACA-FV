/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2017 by the ITHACA-FV authors
-------------------------------------------------------------------------------

License
    This file is part of ITHACA-FV

    ITHACA-FV is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ITHACA-FV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

/// \file
/// Source file of the SimpleSteadyNSROM class.


#include "SimpleSteadyNSROM.H"

// * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * //
SimpleSteadyNSROM::SimpleSteadyNSROM(SteadyNSSimple& Foamproblem, Eigen::MatrixXd vel)
    :
    //calling the constructor of the base class
    reducedSimpleSteadyNS(Foamproblem, vel),

    // inizializer list using relevant data from the FOM problem
    P(Foamproblem._p()),
    U(Foamproblem._U()),
    mesh(Foamproblem._mesh()),
    phi(Foamproblem._phi()),
    Folder("./ITHACAoutput/Online"),
    NmodesUproj(Foamproblem.NUmodes),
    NmodesPproj(Foamproblem.NPmodes),
    NmodesSup(Foamproblem.NSUPmodes),
    NmodesNut(Foamproblem.NNutModes)
{

    maxIterOn=2000;

    //appending the modes from the FOM problem
    ULmodes.resize(0);

    for (int i = 0; i < problem->inletIndex.rows(); i++)
    {
        ULmodes.append(problem->liftfield[i]);
    }

    for (int i = 0; i < NmodesUproj; i++)
    {
        ULmodes.append(problem->Umodes.toPtrList()[i]);
    }

    for (int i = 0; i < NmodesSup; i++)
    {
        ULmodes.append(problem->supmodes.toPtrList()[i]);
    }

    // assigning the number of modes to be extracted
    UprojN = NmodesUproj;
    PprojN = NmodesPproj;

    // reading values for the residual jump from the ITHACAdict
    residualJumpLim =
        problem->para->ITHACAdict->lookupOrDefault<float>("residualJumpLim", 1e-5);
    normalizedResidualLim =
        problem->para->ITHACAdict->lookupOrDefault<float>("normalizedResidualLim",
                1e-5);
    residual_jump=1 + residualJumpLim;

    // inizialization of expansion coefficients and setting up boundary conditions
    a = Eigen::VectorXd::Zero(UprojN);
    b = Eigen::VectorXd::Zero(PprojN);
    a(0) = vel_now(0, 0);
}


// REDUCED SIMPLE ALGORITHM solve for current value of the viscosity
void SimpleSteadyNSROM::solveOnline(scalar mu_now)
{
    counter++;

    //changing value of the FOM problem, need to be done because the equation are then projected into the reduced space
    problem->change_viscosity(mu_now);
    if (counter==1)
    {
        mkDir(Folder);
        ITHACAutilities::createSymLink(Folder);
    }

    Eigen::VectorXd uresidualOld = Eigen::VectorXd::Zero(UprojN);
    Eigen::VectorXd presidualOld = Eigen::VectorXd::Zero(PprojN);
    Eigen::VectorXd uresidual = Eigen::VectorXd::Zero(UprojN);
    Eigen::VectorXd presidual = Eigen::VectorXd::Zero(PprojN);
    scalar U_norm_res(1);
    scalar P_norm_res(1);

    Time& runTime = problem->_runTime();
    P.rename("p");
    ULmodes.reconstruct(U, a, "U");
    problem->Pmodes.reconstruct(P, b, "p");
    phi = fvc::interpolate(U) & U.mesh().Sf();
    int iter = 0;
    simpleControl& simple = problem->_simple();



    // eddy modes handling for future expansion with turbulence flow problems
    if (ITHACAutilities::isTurbulent())
    {

        Eigen::MatrixXd nutCoeff;
        nutCoeff.resize(NmodesNut, 1);

        for (int i = 0; i < NmodesNut; i++)
        {
            Eigen::MatrixXd muEval;
            muEval.resize(1, 1);
            muEval(0, 0) = mu_now;
            nutCoeff(i, 0) = problem->rbfSplines[i]->eval(muEval);
        }

        volScalarField& nut = const_cast<volScalarField&>
                              (problem->_mesh().lookupObject<volScalarField>("nut"));
        problem->nutModes.reconstruct(nut, nutCoeff, "nut");
        ITHACAstream::exportSolution(nut, name(counter), Folder);
    }

    PtrList<volVectorField> gradModP;

    for (int i = 0; i < NmodesPproj; i++)
    {
        gradModP.append(fvc::grad(problem->Pmodes[i]));
    }

    projGradModP = ULmodes.project(gradModP, NmodesUproj);

    //Simple loop
    while ((residual_jump > residualJumpLim
            || std::max(U_norm_res, P_norm_res) > normalizedResidualLim)
            && iter < maxIterOn)
    {
        iter++;
        std::cout << "Iteration " << iter << std::endl;
#if OFVER == 6
        simple.loop(runTime);
#else
        simple.loop();
#endif
        volScalarField nueff = problem->turbulence->nuEff();
        fvVectorMatrix UEqn
        (
            fvm::div(phi, U)
            - fvm::laplacian(nueff, U)
            - fvc::div(nueff * dev2(T(fvc::grad(U))))
        );
        UEqn.relax();
        List<Eigen::MatrixXd> RedLinSysU = ULmodes.project(UEqn, UprojN);
        RedLinSysU[1] = RedLinSysU[1] - projGradModP * b;
        a = reducedProblem::solveLinearSys(RedLinSysU, a, uresidual, vel_now);
        ULmodes.reconstruct(U, a, "U");
        volScalarField rAU(1.0 / UEqn.A());
        volVectorField HbyA(constrainHbyA(1.0 / UEqn.A() * UEqn.H(), U, P));
        surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
        adjustPhi(phiHbyA, U, P);
        tmp<volScalarField> rAtU(rAU);

        if (simple.consistent())
        {
            rAtU = 1.0 / (1.0 / rAU - UEqn.H1());
            phiHbyA +=
                fvc::interpolate(rAtU() - rAU) * fvc::snGrad(P) * mesh.magSf();
            HbyA -= (rAU - rAtU()) * fvc::grad(P);
        }

        List<Eigen::MatrixXd> RedLinSysP;

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian(rAtU(), P) == fvc::div(phiHbyA)
            );
            RedLinSysP = problem->Pmodes.project(pEqn, PprojN);
            b = reducedProblem::solveLinearSys(RedLinSysP, b, presidual);
            problem->Pmodes.reconstruct(P, b, "p");

            if (simple.finalNonOrthogonalIter())
            {
                phi = phiHbyA - pEqn.flux();
            }
        }

        P.relax();
        U = HbyA - rAtU() * fvc::grad(P);
        U.correctBoundaryConditions();
        uresidualOld = uresidualOld - uresidual;
        presidualOld = presidualOld - presidual;
        uresidualOld = uresidualOld.cwiseAbs();
        presidualOld = presidualOld.cwiseAbs();
        residual_jump = std::max(uresidualOld.sum(), presidualOld.sum());
        uresidualOld = uresidual;
        presidualOld = presidual;
        uresidual = uresidual.cwiseAbs();
        presidual = presidual.cwiseAbs();
        U_norm_res = uresidual.sum() / (RedLinSysU[1].cwiseAbs()).sum();
        P_norm_res = presidual.sum() / (RedLinSysP[1].cwiseAbs()).sum();

        if (problem->para->debug)
        {
            std::cout << "Residual jump = " << residual_jump << std::endl;
            std::cout << "Normalized residual = " << std::max(U_norm_res,
                      P_norm_res) << std::endl;
        }
    }

    std::cout << "Solution " << counter << " converged in " << iter <<
              " iterations." << std::endl;
    std::cout << "Final normalized residual for velocity: " << U_norm_res <<
              std::endl;
    std::cout << "Final normalized residual for pressure: " << P_norm_res <<
              std::endl;

    // reconstruction of the full order solution
    ULmodes.reconstruct(U, a, "Uaux");
    P.rename("Paux");
    problem->Pmodes.reconstruct(P, b, "Paux");
    uRecFields.append(U);
    pRecFields.append(P);

    // exporting the solution into the Online folder
    ITHACAstream::exportSolution(U, name(counter), Folder);
    ITHACAstream::exportSolution(P, name(counter), Folder);
    runTime.setTime(runTime.startTime(), 0);
}
