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

     // Perform an online solve for the new values of inlet velocities
     for (label k = 0; k < bif.mu.size(); k++)
     {
         // Set the reduced viscosity
         ridotto.solveOnline(bif.mu(0,k));
     }

     return 0;
 }
