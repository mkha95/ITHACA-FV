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
#include "SteadyNSSimple.H"
#include "SimpleSteadyNSROM.H"
#include "SteadyNSROM.H"


int main(int argc,char * argv[])
{
    Bifurcation<SteadyNSSimple> bif(argc,argv);
    bif.offlineSolve();
    bif.prepare_POD();

    Eigen::MatrixXd vel_now(1, 1);
    vel_now(0, 0) = 1;

    SimpleSteadyNSROM reduced(bif,vel_now);
     for (label k = 0; k < bif.mu.size(); k++)
     {
         reduced.solveOnline(bif.mu(0,k));
     }
     Info <<"end of main"<<endl;

//      reducedSimpleSteadyNS reduced(bif);
//      scalar NmodesUproj=bif.NUmodes;
//      scalar NmodesPproj=bif.NPmodes;
// //     reduced.setOnlineVelocity(vel_now);
//
//     for (label k = 0; k < (bif.mu).size(); k++)
//     {
//         scalar mu_now = bif.mu(0, k);
//         reduced.solveOnline_Simple(mu_now, NmodesUproj, NmodesPproj);
//     }

     return 0;
 }
