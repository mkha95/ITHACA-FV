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
#include "SteadyNSROM.H"
SteadyNSROM::SteadyNSROM(steadyNS& Foamproblem, Eigen::MatrixXd vel): reducedSteadyNS(Foamproblem)
{


    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();

    // Change initial condition for the lifting function
    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; j++)
        {
            y(j) = vel(j, 0);
        }
    }
    newton_object.BC.resize(N_BC);

    for (int j = 0; j < N_BC; j++)
    {
        newton_object.BC(j) = vel(j, 0);
    }

}

void SteadyNSROM::solveOnline(const scalar & mu)
{

    nu=mu;
    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);
    newton_object.nu = nu;
    Eigen::HybridNonLinearSolver<newton_steadyNS> hnls(newton_object);
    hnls.solve(y);
    Eigen::VectorXd res(y);
    newton_object.operator()(y, res);
    Info << "################## Online solve N° " << count_online_solve <<
         " ##################" << endl;

    if (Pstream::master())
    {
        std::cout << "Solving for the parameter: " << nu << std::endl;
    }

    if (res.norm() < 1e-5 && Pstream::master())
    {
        std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                  hnls.iter << " iterations " << def << std::endl << std::endl;
    }
    else if (Pstream::master())
    {
        std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                  hnls.iter << " iterations " << def << std::endl << std::endl;
    }

    Eigen::MatrixXd tmp_sol(y.rows() + 1, 1);
    tmp_sol<< (count_online_solve + 1),y;
    online_solution.append(tmp_sol);
    reconstruct( "./ITHACAoutput/Online/",count_online_solve);
    count_online_solve += 1;
}


void SteadyNSROM::reconstruct(fileName folder,const int & counter)
 {
     Info<<counter<<endl;
     if (counter==1)
     {
         mkDir(folder);
         ITHACAutilities::createSymLink(folder);
     }
     Eigen::MatrixXd CoeffU= online_solution[counter-1].block(1, 0, Nphi_u, 1);
     Eigen::MatrixXd CoeffP= online_solution[counter-1].bottomRows(Nphi_p);
     volVectorField uRec("uRec", Umodes[0]);
     volScalarField pRec("pRec", problem->Pmodes[0]);
     problem->L_U_SUPmodes.reconstruct(uRec, CoeffU, "uRec");
     problem->Pmodes.reconstruct(pRec, CoeffP, "pRec");
     ITHACAstream::exportSolution(uRec, name(counter), folder);
     ITHACAstream::exportSolution(pRec, name(counter), folder);
}
