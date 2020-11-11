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
#include "steadyNS.H"
#include "SteadyNSSimple.H"
#include <algorithm>
#include <map>



template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

template <typename T>
Bifurcation<T>::Bifurcation(int argc, char* argv[])
    :T( argc, argv)
{
    print_section_name("READING BIFURCATION PARAMETERS");
    bif_params = make_unique<IOdictionary>
    (IOobject
     (
         "BIFURCATIONdict",
         T::_runTime().system(),
         T::_mesh(),
         IOobject::MUST_READ,
         IOobject::AUTO_WRITE
     )
    );
    Info<<"\n";
    mu_inf = bif_params->lookupOrDefault<scalar>("mu_inf", 2);
    mu_sup = bif_params->lookupOrDefault<scalar>("mu_sup", 0.2);
    N_mu = bif_params->lookupOrDefault<scalar>("N_mu", 100);
    Info<<"RANGE:                           mu_inf="<<mu_inf<<"         "<<"mu_sup="<<mu_sup<<endl;
    Info<<"Number of snapshot for POD:        N_mu="<<N_mu<<endl;
    x_cord = bif_params->lookupOrDefault<scalar>("x_cord", 14);
    y_cord = bif_params->lookupOrDefault<scalar>("y_cord", 1.25);
    field  =bif_params->lookupOrDefault<word>("field","velocity");
    word coord  =bif_params->lookupOrDefault<word>("component","x");
    lift_keyword = bif_params->lookupOrDefault<word>("lift","snapshot");
    T::inletIndex.resize(1,2);
    T::inletIndex(0,0)=T::_mesh().boundaryMesh().findPatchID("inlet");
    T::inletIndex(0,1)=0;
    parabolic_inlet = bif_params->lookupOrDefault<bool>("parabolic_inlet",true);
    if (parabolic_inlet)
    {
        scalarList coefs(3);
        parabolic_inlet_coefs = bif_params->lookupOrDefault<scalarList>("par_coefs",coefs);
        set_parabolic_inlet(parabolic_inlet_coefs);
    }

    id_cell=T::_mesh().findCell(point(x_cord,y_cord,0));
    if(id_cell<0)
    {
        Info<<"Error:         Coordinate not belonging to the mesh"<<endl;
        Foam::FatalError.exit();
    }
    Info<<"==========================================================================="<<endl;
    Info<<"Coordinate of the sampling point:            (x="<<x_cord<<";y="<<y_cord<<")"<<endl;
    std::vector<std::string> available_fields{"velocity","pressure"};
    if(find(available_fields.begin(),available_fields.end(),field)==available_fields.end())
    {
        Info<<"Error:        given field not available, possible assignment:\n-pressure\n-velocity "<<endl;
        Foam::FatalError.exit();
    }
    else
        Info<<"Field to be sampled:      "<<field<<endl;

    if(field=="velocity")
    {
        std::map<std::string, int> map_coord_component{{"x",0},{"y",1},{"z",2}};
        component=map_coord_component[coord];
        sampling=[this]() {
            return this->sampling_velocity(id_cell,component);
        };
    }
    else
        sampling=[this]() {
        return this->sampling_pressure(id_cell);
    };

    T::Pnumber=1;
    T::Tnumber=N_mu;
    mu.resize(1,N_mu);
    mu_range.resize(1,2);
    mu_range(0,1)=mu_sup;
    mu_range(0,0)=mu_inf;
    std::cout<<"mu_range:"<<std::endl<<mu_range<<std::endl;
    T::genEquiPar();
    sampled_field.resize(N_mu);
};

template <typename T>
void Bifurcation<T>::offlineSolve()
{
    print_section_name("OFFLINE SOLVE");
    List<scalar> mu_now(1);
    // if the offline solution is already performed read the fields
    if (T::offline)
    {
        ITHACAstream::read_fields(T::Ufield,"U", "./ITHACAoutput/Offline/");
        ITHACAstream::read_fields(T::Pfield,"p", "./ITHACAoutput/Offline/");
    }

    else
    {
        for (label i = 0; i < mu.cols(); i++)
        {
            mu_now[0] = mu(0, i);
            T::change_viscosity(mu(0, i));
            T::truthSolve(mu_now);
            sampled_field(i)=sampling();
        }
        ITHACAstream::exportVector(sampled_field,"sampled_field","eigen","./ITHACAoutput/Offline");
    }

};

template <typename T>
void Bifurcation<T>::prepare_POD(void)
{
    T::solvesupremizer();
    ITHACAparameters* para = ITHACAparameters::getInstance(this->_mesh(),
                             this->_runTime());
    T::NUmodesOut = para->ITHACAdict->lookupOrDefault<scalar>("NmodesUout", 15);
    T::NPmodesOut = para->ITHACAdict->lookupOrDefault<scalar>("NmodesPout", 15);
    T::NSUPmodesOut = para->ITHACAdict->lookupOrDefault<scalar>("NmodesSUPout", 15);
    T::NUmodes = para->ITHACAdict->lookupOrDefault<scalar>("NmodesUproj", 10);
    T::NPmodes = para->ITHACAdict->lookupOrDefault<scalar>("NmodesPproj", 10);
    T::NSUPmodes = para->ITHACAdict->lookupOrDefault<scalar>("NmodesSUPproj", 10);
    if (T::bcMethod == "lift")
    {
        lift_solve();
        //ITHACAutilities::normalizeFields(T::liftfield);
        // Homogenize the snapshots
        //T::computeLift(T::Ufield, T::liftfield, T::Uomfield);
        computeLift(T::Ufield, T::liftfield, T::Uomfield);
        // Perform POD on the velocity snapshots
        ITHACAPOD::getModes(T::Uomfield, T::Umodes, T::_U().name(),
                            T::podex, 0, 0, T::NUmodesOut);
    }

    else
    {
        ITHACAPOD::getModes(T::Ufield, T::Umodes, T::_U().name(),
                            T::podex, 0, 0, T::NUmodesOut);
    }

    // Perform POD on velocity pressure and supremizers and store the first modes
    ITHACAPOD::getModes(T::Pfield, T::Pmodes, T::_p().name(),
                        T::podex, 0, 0,
                        T::NPmodesOut);
    ITHACAPOD::getModes(T::supfield, T::supmodes, T::_U().name(),
                        T::podex,
                        T::supex, 1, T::NSUPmodesOut);
    T::projectSUP("./Matrices", T::NUmodes, T::NPmodes, T::NSUPmodes);

}


template <typename T>
scalar Bifurcation<T>::sampling_pressure(label id_cell)
{
    return T::_p().internalField()[id_cell];
}


template <typename T>
scalar Bifurcation<T>::sampling_velocity(label id_cell,label component)
{
    return T::_U().internalField()[id_cell].component(component);
}


template <typename T>
void Bifurcation<T>::print_section_name(const std::string & title) const
{
    Info<<"==========================================================================="<<endl;
    Info<<"             "<<title<<endl;
    Info<<"==========================================================================="<<endl;
}



// Method to compute the lifting function
template <typename T>
void Bifurcation<T>::lift_solve()
{
    if (lift_keyword=="potential")
    {
        print_section_name("SOLVING LIFTING PROBLEM");
        surfaceScalarField& phi =T::_phi();
        volScalarField& p = T::_p0();
        volVectorField& U = T::_U0();
        IOMRFZoneList& MRF = T::_MRF();
        Time& runTime = T::_runTime();
        fvMesh& mesh = T::_mesh();
        label BCind = T::inletIndex(0, 0);
        volVectorField Ulift("Ulift",U);
        instantList Times = runTime.times();
        runTime.setTime(Times[1], 1);
        pisoControl potentialFlow(mesh, "potentialFlow");
        Vector<double> v0(0, 0, 0);

        for (label j = 0; j < U.boundaryField().size(); j++)
        {
            if (j!=BCind)
            {
                T::assignBC(Ulift, j, v0);
            }
        }
        T::assignIF(Ulift, v0);
        phi = linearInterpolate(Ulift) & mesh.Sf();

        Info << "Constructing velocity potential field Phi\n" << endl;

        volScalarField Phi
        (
            IOobject
            (
                "Phi",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Phi", dimLength * dimVelocity, 0),
            p.boundaryField().types()
        );
        label PhiRefCell = 0;
        scalar PhiRefValue = 0;
        setRefCell
        (
            Phi,
            potentialFlow.dict(),
            PhiRefCell,
            PhiRefValue
        );
        mesh.setFluxRequired(Phi.name());
        runTime.functionObjects().start();
        MRF.makeRelative(phi);
        adjustPhi(phi, Ulift, p);

        while (potentialFlow.correctNonOrthogonal())
        {
            fvScalarMatrix PhiEqn
            (
                fvm::laplacian(dimensionedScalar("1", dimless, 1), Phi)
                ==
                fvc::div(phi)
            );
            PhiEqn.setReference(PhiRefCell, PhiRefValue);
            PhiEqn.solve();

            if (potentialFlow.finalNonOrthogonalIter())
            {
                phi -= PhiEqn.flux();
            }
        }

        MRF.makeAbsolute(phi);
        Info << "Continuity error = "
             << mag(fvc::div(phi))().weightedAverage(mesh.V()).value()
             << endl;
        Ulift = fvc::reconstruct(phi);
        Ulift.correctBoundaryConditions();
        Info << "Interpolated velocity error = "
             << (sqrt(sum(sqr((fvc::interpolate(U) & mesh.Sf()) - phi)))
                 / sum(mesh.magSf())).value()
             << endl;
        T::liftfield.append(Ulift);
    }
    else if (lift_keyword=="snapshot")
    {
        bool check_lift=ITHACAutilities::check_folder("Lift");
        if (check_lift)
        {
            Info<<"Lifting function already exist: reading from folder \"./Lift/\""<<endl;
            ITHACAstream::read_fields(this->liftfield, T::_U(), "./Lift/");
        }
        else
        {
            //in questo caso non abbiamo la cartella lift quindi la devo creare e devo anche prendere
            //uno snapshot da usare come lift function
            Info<<"\"Lift\" folder not found, creating \"./Lift/\" from first snapshot"<<endl;
            ITHACAutilities::createSymLink("./Lift");
            volVectorField& Ulift(T::Ufield[0]);
            T::liftfield.append(Ulift);
            ITHACAstream::exportSolution(Ulift, "1", "./Lift/");
        }
    }
}

template <typename T>
template <typename G>
void Bifurcation<T>::computeLift(G& Lfield, G& liftfield, G& omfield)
{

    for (label k = 0; k < T::inletIndex.rows(); k++)
    {

        for (label j = 0; j < Lfield.size(); j++)
        {
            volVectorField C("U", Lfield[j] - liftfield[k]);
            if (k == 0)
                omfield.append(C);
            else
                omfield.set(j, C);
        }
    }
}

template <typename T>
void Bifurcation<T>::set_parabolic_inlet(const scalarList & coef)
{
            label BC_ind = T::inletIndex(0, 0);
            const polyPatch& pp = _mesh().boundaryMesh()[BC_ind];
            forAll(T::_U().boundaryField()[BC_ind], faceI)
            {
               scalar y = pp.faceCentres()[faceI].y();
               T::_U().boundaryFieldRef()[BC_ind][faceI] = vector(coef[0]*y*y+coef[1]*y+coef[2],0,0);
            }
            T::_U().write();
}



template class Bifurcation<steadyNS>;
template class Bifurcation<SteadyNSSimple>;
