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
Class
    SteadyNSROM
Description
     A class that provides the utlilities for solving a reduced order model for a Steady navier stokes problem in a fully implicit was
SourceFiles
     SteadyNSROM.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the SimpleSteadyNSROM class.
/// \dir
/// Directory containing the header, source files for the SimpleSteadyNSROM class.
#include "ReducedSteadyNS.H"



#ifndef SteadyNSROM_H
#define SteadyNSROM_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                        Class SteadyNSROM
\*---------------------------------------------------------------------------*/

/// Implementation of class for the reduced model of a <b> Steady navier stokes problem </b> for which the online solve uses a inconsistent fully implicit
/** In this class are implemented the methods for the offline solve and a continuation method for a steady NS problem.
this class inherits from reducedSteadyNS**/
class SteadyNSROM: reducedSteadyNS
{
    //public member
    public:

        //--------------------------------------------------------------------------
        /// constructing a SteadyNSROM starting from the FOM counterpart
        ///
        /// @param[in]      Foamproblem  full order problem which must be an object of the SteadyNSSimple class
        /// @param[in]      vel          eigen matrix that encodes the information for the boundary condition
        SteadyNSROM(steadyNS& Foamproblem, Eigen::MatrixXd vel);

        /// online solve for the current value of the kinematic viscosity
        void solveOnline(const scalar &);

        /// list of pointer of the reduced field for pressure
        PtrList<volScalarField> pRecFields;

        /// list of pointer of the reduced field for velocity
        PtrList<volVectorField> uRecFields;

        // destructor
        ~SteadyNSROM() {};

    // private members
    private:

        //--------------------------------------------------------------------------
        /// private member function used to reconstruct the solution and saving the corresponding outputs
        ///
        /// @param[in]      folder   folder name for reconstruction output="./ITHACAoutput/Online"
        /// @param[in]      counter  counter value corresponding to current viscosity
        void reconstruct(fileName folder,const int & counter);
};


#endif
