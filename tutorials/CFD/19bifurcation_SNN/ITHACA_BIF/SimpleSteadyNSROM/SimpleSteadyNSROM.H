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
     SimpleSteadyNSROM
 Description
     A reduced problem for the stationary NS equations, solved with the "SegregatedROM" approach
 SourceFiles
     reducedSteadyNS.C
 \*---------------------------------------------------------------------------*/

/// \file
/// Header file of the SimpleSteadyNSROM class.
/// \dir
/// Directory containing the header, source files for the SimpleSteadyNSROM class.

#ifndef SimpleSteadyNSROM_H
#define SimpleSteadyNSROM_H
#include "ReducedSimpleSteadyNS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                        Class SimpleSteadyNSROM
\*---------------------------------------------------------------------------*/

/// Implementation of class for the reduced model of a <b> Steady navier stokes problem </b> for which the online solve uses a consistent approach("segregatedROM")
/** In this class are implemented the methods for the offline solve and a continuation method for a steady NS problem.
this class inherits from reducedSimpleSteadyNS**/
class SimpleSteadyNSROM: public reducedSimpleSteadyNS
{

    // private members
    private:

        /// number of modes for the projection of the velocity field
        scalar UprojN;

        /// number of modes for the projection of the velocity field
        scalar PprojN;

        // references to the fields of the FOM problem
        /// reference for the pressure field of SteadyNSSimple
        volScalarField& P;


        /// reference for the velocity field of SteadyNSSimple
        volVectorField& U;


        /// reference for the mesh of SteadyNSSimple
        fvMesh& mesh;

        /// reference for the flux field of SimpleSteadyNSROM
        surfaceScalarField& phi;

        //other members
        /// eigen matrix to store the coefficients for the modes expansion of the velocity field
        Eigen::MatrixXd b;

        /// eigen matrix to store the coefficients for the modes expansion of the velocity field
        Eigen::MatrixXd a;

        /// name for the online results output folder:"./ITHACAoutput/Online"
        word Folder;

        /// number of velocity modes
        label NmodesUproj;

        /// number of pressure modes
        label NmodesPproj;

        /// number of supremizer modes
        label NmodesSup;

        /// number of modes for Eddy viscosity
        label NmodesNut;

        /// data structure for handling the pressure modes
        volScalarModes Pmodes;

        /// limit value for the jump in the residuals
        float residualJumpLim;

        /// normalized residual jump used in the SIMPLE algorithm
        float normalizedResidualLim;

        /// residual jump
        scalar residual_jump;

    // public members
    public:

        //--------------------------------------------------------------------------
        /// constructing a SimpleSteadyNSROM starting from the FOM counterpart
        ///
        /// @param[in]      Foamproblem  full order problem which must be an object of the SteadyNSSimple class
        /// @param[in]      vel          eigen matrix that encodes the information for the boundary condition
        SimpleSteadyNSROM(SteadyNSSimple& Foamproblem, Eigen::MatrixXd vel);

        /// list of pointer of the reduced field for pressure
        PtrList<volScalarField> pRecFields;

        /// list of pointer of the reduced field for velocity
        PtrList<volVectorField> uRecFields;


        /// online solve for the current value of the kinematic viscosity
        void solveOnline(scalar );

        // destructor
        ~SimpleSteadyNSROM() {};
};



#endif
