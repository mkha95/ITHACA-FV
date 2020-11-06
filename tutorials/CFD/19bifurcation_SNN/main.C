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

#include "bifurcationFOM_NSS.H"
#include "bifurcationROM_NSS.H"
#include "SteadyNSSimple.H"


int main(int argc,char * argv[]){
    Bifurcation<steadyNS> bif(argc,argv);
    bif.offlineSolve();
    bif.prepare_POD();

    Eigen::MatrixXd vel_now(1, 1);
    vel_now(0, 0) = 1;

    BifurcationROM<reducedSteadyNS> ridotto(bif, vel_now);
     // Set the inlet velocity
     // used only for penalty approach
     ridotto.tauU = Eigen::MatrixXd::Zero(1, 1);
     ridotto.tauU(0, 0) = 1e-1;
     std::cout<<"party"<<std::endl;
     std::cout<<bif.mu<<std::endl;

     // Perform an online solve for the new values of inlet velocities
     for (label k = 0; k < bif.mu.size(); k++)
     {
         // Set the reduced viscosity
         ridotto.nu = bif.mu(0, k);
         std::cout<<"no problem ye1"<<std::endl;
         ridotto.solveOnline();
         std::cout<<"no problem ye2"<<std::endl;
         Eigen::MatrixXd tmp_sol(ridotto.y.rows() + 1, 1);
         tmp_sol(0) = k + 1;
         tmp_sol.col(0).tail(ridotto.y.rows()) = ridotto.y;
         ridotto.online_solution.append(tmp_sol);
     }

     // Save the online solution
     ITHACAstream::exportMatrix(ridotto.online_solution, "red_coeff", "python",
                                "./ITHACAoutput/red_coeff");
     ITHACAstream::exportMatrix(ridotto.online_solution, "red_coeff", "matlab",
                                "./ITHACAoutput/red_coeff");
     ITHACAstream::exportMatrix(ridotto.online_solution, "red_coeff", "eigen",
                                "./ITHACAoutput/red_coeff");
     // Reconstruct and export the solution
     ridotto.reconstruct(true, "./ITHACAoutput/Reconstruction/");
     return 0;
 }
