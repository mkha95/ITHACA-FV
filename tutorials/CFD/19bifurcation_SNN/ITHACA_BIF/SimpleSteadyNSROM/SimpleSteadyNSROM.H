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
     reducedSteadyNS
 Description
     A reduced problem for the stationary NS equations
 SourceFiles
     reducedSteadyNS.C
 \*---------------------------------------------------------------------------*/

#include "ReducedSimpleSteadyNS.H"

class SimpleSteadyNSROM: public reducedSimpleSteadyNS
{
    scalar UprojN;
    scalar PprojN;
    volScalarField& P;
    volVectorField& U;
    fvMesh& mesh;
    Eigen::MatrixXd a;
    Eigen::MatrixXd b;
    surfaceScalarField& phi;
    word Folder;
    label NmodesUproj;
    label NmodesPproj;
    label NmodesSup;
    label NmodesNut;
    volScalarModes Pmodes;
    float residualJumpLim;
    float normalizedResidualLim;
    scalar residual_jump;
    public:
    SimpleSteadyNSROM(SteadyNSSimple& Foamproblem, Eigen::MatrixXd vel);
    void solveOnline(scalar );
    ~SimpleSteadyNSROM() {};
};
