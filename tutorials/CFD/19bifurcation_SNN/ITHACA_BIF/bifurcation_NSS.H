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
Description
    Example of steady NS Reduction Problem
SourceFiles
    03steadyNS.C
\*---------------------------------------------------------------------------*/
#ifndef bifurcation_NSS_H
#define bifurcation_NSS_H
#include "ITHACAstream.C"
#include "ITHACAPOD.H"
#include "ReducedSteadyNS.H"
#include "forces.H"
#include "IOmanip.H"
#include <memory>
#include <functional>

template <typename T>
class Bifurcation:public T
{
private:
    IOdictionary* bif_params;
    //mu range for the bifurcation diagram
    scalar mu_inf;
    scalar mu_sup;
    label N_mu;
    // id of the cell corresponding to the position we want to sample
    label id_cell;
    // physical_quantity for retriving the bifurcation diagram
    word field;
    scalar x_cord;
    scalar y_cord;
    label component;
    using T::mu_range;
    //Eigen::MatrixXd & mu=T::mu;
    scalar sampling_velocity(label id_cell,label component);
    scalar sampling_pressure(label id_cell);
    std::function<scalar(void)> sampling;
    void print_section_name(const std::string &) const;
    Eigen::VectorXd sampled_field;
    word lift_keyword;
    template <typename G> void computeLift(G &Lfield,G & 	liftfield,G & 	omfield );
public:
    using T::mu;
    void lift_solve();
    void print_mu() {
        Info<<mu<<endl;
    };
    Bifurcation()   {
        T();
    }
    Bifurcation(int argc, char* argv[]);
    ~Bifurcation()=default;
    void offlineSolve();
    void prepare_POD();
    //void print_bifurcation() const;
};


#endif
