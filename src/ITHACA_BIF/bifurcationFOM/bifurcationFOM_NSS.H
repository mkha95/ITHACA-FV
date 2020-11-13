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
    Bifurcation<T>
Description
    a template class for bifurcation at FOM level is implemented
SourceFiles
    bifurcationFOM_NSS.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the bifurcationFOM_NSS class.
/// \dir
/// Directory containing the header and source files for the Bifurcation<T> class.

#ifndef bifurcation_NSS_H
#define bifurcation_NSS_H
#include "ITHACAstream.C"
#include "ITHACAPOD.H"
#include "ReducedSteadyNS.H"
#include "forces.H"
#include "IOmanip.H"
#include <memory>
#include <functional>




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                        Class template Bifurcation<T>
\*---------------------------------------------------------------------------*/

/// Implementation of class for inspectig bifurcation phenomena in the context of parametrized full order <b> Steady navier stokes problem </b>
/** In this class are implemented the methods for the offline solve and a continuation method for a steady NS problem.
and the for the generation of the reduced matrices for subsequent online solve, this class is inherits publically from either
steadyNS or SteadyNSSimple, depending on the the reduced-model one wants to use afterwards*/


template <typename T>
class Bifurcation:public T
{
    // public members
public:

    // cunstructor and destructor
    ///cunstructs the class starting from the infomation provided in the system folder
    ///information concerning the bifurcation class are taken from system/BIFURCATIONdict
    Bifurcation(int argc, char* argv[]);

    // defalult destructor
    ~Bifurcation()=default;

    // member function for continuation method and POD
    /// offline solve method which performs a FOM solve within the range provided for the continuation
    /// method. Both the number of snapshots and the interval [mu_inf, mu_sup] are read from BIFURCATIONdict
    void offlineSolve();

    /// method which prepare the POD extracting the modes from the snapshots of the offline solve.
    /// all the data concering the output can be retrived in the ITHACAoutput/POD folder
    void prepare_POD();

    /// lifter solve function which is used to enrich the modes extracted from the POD phase
    /// with a lift function that enables the satisfaction of boundary conditions
    /// in the BIFURCATIONdict one can chose among lift potentatial or snapshot which
    /// the potential setting consists in solving a potential problem and using the solution as a lifting function
    /// the snapshot setting consists in taking the lifting directly from the first snapshot
    void lift_solve();


    // other members
    /// a vector used for storing the all the computation times required to solve each sample in the offlineSolve
    Eigen::VectorXd computation_times;

    /// variable used to store the complessive offline time. this is obtained summing up the elements in computation_times
    scalar total_computation_time;

    /// suggested values for the number of modes of each field
    int sugUmodes,sugPmodes, sugSupmodes;

    /// importing the mu variable from the base class reductionProblem: row matrix containing the parameters
    using T::mu;

    ///user defined variable read from BIFURCATIONdict to decide whether to perform the online stage or not
    word online;


    // Private members
private:
    //importing relevant member from base class and making them private within this context
    /// Matrix that contains informations about the inlet boundaries.
    using T::inletIndex;

    ///  mesh
    using T::_mesh;

    /// Range of the parameter spaces
    using T::mu_range;

    /// unique_ptr to read data from the BIFURCATIONdict dictionary
    std::unique_ptr<IOdictionary> bif_params;

    ///mu range for the bifurcation diagram, it is not required that mu_inf<mu_sup
    scalar mu_inf;
    scalar mu_sup;

    /// number of snapshot for the continuation method, also used for the snapshot matrix of the POD phase
    label N_mu;

    /// id of the cell corresponding to the position of the domain where to sample the scalar field that will be stored in sampled_field
    label id_cell;

    /// physical_quantity for retriving the bifurcation diagram
    word field;

    /// (x,y) coordinates of the sampling point
    scalar x_cord;
    scalar y_cord;

    /// variable used when field="velocity" in order to chose the componet to be sampled
    /// 0=x,1=y
    label component;

    /// variable set from BIFURCATIONdict to decide whether or not to impose an parabolic condition for the inlet profile
    bool parabolic_inlet;

    /// coefficients used in the set_parabolic_inlet member function for defining the quadratic profile
    scalarList parabolic_inlet_coefs;

    /// gives  the scalar field velocity in the (x,y) chosen position
    scalar sampling_velocity(label id_cell,label component);

    /// gives  the pressure in the (x,y) chosen position
    scalar sampling_pressure(label id_cell);

    /// general polymorphic function wrapper which will be used to sample the chosen scalar field
    std::function<scalar(void)> sampling;

    /// routine to print on the standard output the given section name
    void print_section_name(const std::string &) const;

    /// eigen matrix that stores the values of the sampled field returned by sampling
    Eigen::VectorXd sampled_field;

    /// variable used to decide the strategy for computing the lift function
    /// value read from BIFURCATIONdict anc can either be potential or snapshot
    word lift_keyword;


    //--------------------------------------------------------------------------
    /// Homogenize the snapshot matrix, it works with PtrList of volVectorField and volScalarField
    ///
    /// @param[in]      Lfield     The list of snapshots to be homogenized.
    /// @param[in]      liftfield  The list containing.
    /// @param[out]     omfield    The homogenized snapshot matrix.
    ///
    /// @tparam     G          { description }
    ///
    template <typename G> void computeLift(G &Lfield,G & 	liftfield,G & 	omfield );

    /// member function to set the parabolic inlet profile starting from a 3-dimensional list that encodes the coefficients
    void set_parabolic_inlet(const scalarList &);

    /// eigen vectors that containes the cumulative eigenvalues corresponding to SVD decomposition
    Eigen::MatrixXd cumUeig, cumPeig, cumSupeig;

    /// tolerance to be used to suggest the number of modes to extract from the POD
    scalar cum_eig_tol=0.99999;

    //--------------------------------------------------------------------------
    /// @brief     filling an eigen matrix, provided a folder name
    /// @param[in]  mat      reference to the eigen matrix to be filled
    /// @param[in]  folder   name of the folder from which the matrix is filled
    ///
    ///
    // void read_matrix(Eigen::MatrixXd & mat,std::string folder);
    void read_matrix(Eigen::MatrixXd & mat,std::string folder);

};



#endif





// Doxigen Documentation (Main Page)

/*! \mainpage In real Time Highly Advanced Computational Applications for Finite Volumes
 *
 * \section About About ITHACA-FV
 *
 * ITHACA-FV (In real Time Highly Advanced Computational Applications for
Finite Volumes) is C++ library based on the finite volume solver
[<b> OpenFOAM </b>](https://openfoam.org/). It consists of the implementation of several reduced order
modeling techniques for parametrized problems. ITHACA-FV is developed and maintained by
[<b> Dr. Giovanni Stabile </b>](http://people.sissa.it/~gstabile/) under the supervision of [<b> Prof. Gianluigi Rozza </b>](http://people.sissa.it/~grozza/) at [<b> SISSA mathLab </b>](https://mathlab.sissa.it/) an applied math laboratory of the
[<b> International School for advanced studies </b>](https://www.sissa.it/).

 *
 * \section structure ITHACA-FV directory structure
 * ITHACA-FV comprises of four main directories:

- \b src: the core ITHACA-FV library.
- \b applications: some utilities.
- \b tutorials: test-cases that demonstrate a range of ITHACA-FV functionalities.
- \b docs: documentation.
 *
 * \section step1 How to cite IHACA-FV.
 *
Most of the theoretical aspects behind ITHACA-FV are deeply explained in [<b> Stabile2017CAIM </b>](https://doi.org/10.1515/caim-2017-0011) and [<b> Stabile2017CAF </b>](https://doi.org/10.1016/j.compfluid.2018.01.035).
For this reason, if you use this software, please consider citing the mentioned works, reported in the following bibtex entries:


@Article{Stabile2017CAIM,
Title                    = {{POD-Galerkin reduced order methods for CFD using Finite Volume Discretisation: vortex shedding around a circular cylinder}},
Author                   = {Stabile, Giovanni and Hijazi, Saddam and Mola, Andrea and Lorenzi, Stefano and Rozza, Gianluigi},
Journal                  = {Communications in Applied and Industrial Mathematics},
Year                     = {(2017)},
Volume                   = {8},
Number                   = {1},
pages                    = {210-236},
Doi                      = {10.1515/caim-2017-0011}}

@Article{Stabile2017CAF,
Title                    = {{Finite volume POD-Galerkin stabilised reduced order methods for the parametrised incompressible Navier-Stokes equations}},
Author                   = {Stabile, Giovanni and Rozza, Gianluigi},
Journal                  = {Computers & Fluids},
Year                     = {2018},
Doi                      = {10.1016/j.compfluid.2018.01.035}}

And cite the ITHACA-FV [website](https://mathlab.sissa.it/ithaca-fv).

 */
