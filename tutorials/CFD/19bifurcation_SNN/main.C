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
    //creating a bifurcation object for the steadyNS class
    Bifurcation<steadyNS> bif(argc,argv);
    // perform offline solve
    bif.offlineSolve();
    //prepare POD
    bif.prepare_POD();

    //displaying suggested values for modes number
    std::cout<<"suggested value for modesU:   "<<bif.sugUmodes<<std::endl;
    std::cout<<"suggested value for modesp:   "<<bif.sugPmodes<<std::endl;
    std::cout<<"suggested value for modesSup: "<<bif.sugSupmodes<<std::endl;
    if (bif.online=="yes")
    {
        // using the same boundary condition as for the FOM problem
        Eigen::MatrixXd vel_now(1, 1);
        vel_now(0, 0) = 1;
        SteadyNSROM reduced(bif,vel_now);

        // continuation method for the reduced problem
        for (label k = 0; k < bif.mu.size(); k++)
        {
            scalar mu_now=bif.mu(0,k);
            reduced.solveOnline(mu_now);
        }
    Info<<"Total FOM computation time: "<<bif.total_computation_time<<endl;

    // evaluating reconstruction error
    Eigen::MatrixXd errFrobU = ITHACAutilities::errorFrobRel(bif.Ufield,
                               reduced.uRecFields);
    Eigen::MatrixXd errFrobP =  ITHACAutilities::errorFrobRel(bif.Pfield,
                                reduced.pRecFields);
    Eigen::MatrixXd errL2U = ITHACAutilities::errorL2Rel(bif.Ufield,
                             reduced.uRecFields);
    Eigen::MatrixXd errL2P =  ITHACAutilities::errorL2Rel(bif.Pfield,
                              reduced.pRecFields);

    //exporting the errors
    ITHACAstream::exportMatrix(errFrobU, "errFrobU", "matlab",
                               "./ITHACAoutput/ErrorsFrob/");
    ITHACAstream::exportMatrix(errFrobP, "errFrobP", "matlab",
                               "./ITHACAoutput/ErrorsFrob/");
    ITHACAstream::exportMatrix(errL2U, "errL2U", "matlab",
                               "./ITHACAoutput/ErrorsL2/");
    ITHACAstream::exportMatrix(errL2P, "errL2P", "matlab",
                               "./ITHACAoutput/ErrorsL2/");
    }

    return 0;
}
