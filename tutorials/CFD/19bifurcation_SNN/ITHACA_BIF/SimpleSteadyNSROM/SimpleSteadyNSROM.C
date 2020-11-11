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
#include "SimpleSteadyNSROM.H"

SimpleSteadyNSROM::SimpleSteadyNSROM(SteadyNSSimple& Foamproblem, Eigen::MatrixXd vel)
    :
    reducedSimpleSteadyNS(Foamproblem, vel),
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

    ULmodes.resize(0);
    Info<<"checkpoinnt1"<<endl;
    setOnlineVelocity(vel);
    maxIterOn=100;

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


    Info<<"checkpoinnt2"<<endl;


    UprojN = NmodesUproj + NmodesSup;
    PprojN = NmodesPproj;
    Info<<"checkpoinnt3"<<endl;
    residualJumpLim =
        problem->para->ITHACAdict->lookupOrDefault<float>("residualJumpLim", 1e-5);
    normalizedResidualLim =
        problem->para->ITHACAdict->lookupOrDefault<float>("normalizedResidualLim",
                1e-5);

    //problem->restart();

    residual_jump=1 + residualJumpLim;
}



void SimpleSteadyNSROM::solveOnline(scalar mu_now)
{
    counter++;
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
    Eigen::MatrixXd a = Eigen::VectorXd::Zero(UprojN);
    Eigen::MatrixXd b = Eigen::VectorXd::Zero(PprojN);
    a(0) = vel_now(0, 0);
    P.rename("p");
    //ULmodes.reconstruct(U, a, "U");
    Info<<"checkpoinnt4"<<endl;
    //problem->Pmodes.reconstruct(P, b, "p");
    Info<<"checkpoinnt5"<<endl;
    phi = fvc::interpolate(U) & U.mesh().Sf();
    int iter = 0;
    simpleControl& simple = problem->_simple();
    Info<<"checkpoinnt6"<<endl;

    if (ITHACAutilities::isTurbulent())
    {

        Info<<"checkpoinnt7"<<endl;
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
    Info<<"checkpoinnt8"<<endl;

    for (int i = 0; i < NmodesPproj; i++)
    {
        gradModP.append(fvc::grad(problem->Pmodes[i]));
    }

    projGradModP = ULmodes.project(gradModP, NmodesUproj);

    Info<<"checkpoinnt9"<<endl;
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
    ULmodes.reconstruct(U, a, "Uaux");
    P.rename("Paux");
    problem->Pmodes.reconstruct(P, b, "Paux");
    ITHACAstream::exportSolution(U, name(counter), Folder);
    ITHACAstream::exportSolution(P, name(counter), Folder);
    runTime.setTime(runTime.startTime(), 0);
}
