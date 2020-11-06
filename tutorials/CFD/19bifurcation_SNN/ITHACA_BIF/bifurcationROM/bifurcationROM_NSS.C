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
#include "ReducedSteadyNS.H"
#include "bifurcationROM_NSS.H"
template <typename T>
BifurcationROM<T>::BifurcationROM(steadyNS& Foamproblem, Eigen::MatrixXd vel): T(Foamproblem)
{


    T::y.resize(Nphi_u + Nphi_p, 1);
    T::y.setZero();

    // Change initial condition for the lifting function
    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < T::N_BC; j++)
        {
            y(j) = vel(j, 0);
        }
    }
    T::newton_object.BC.resize(T::N_BC);

    for (int j = 0; j < N_BC; j++)
    {
        T::newton_object.BC(j) = vel(j, 0);
    }

}

template <typename T>
void BifurcationROM<T>::solveOnline()
{

    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);
    T::newton_object.nu = T::nu;
    Eigen::HybridNonLinearSolver<newton_steadyNS> hnls(newton_object);
    hnls.solve(T::y);
    Eigen::VectorXd res(y);
    T::newton_object.operator()(T::y, res);
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

    T::count_online_solve += 1;
}

template class BifurcationROM<reducedSteadyNS>;
