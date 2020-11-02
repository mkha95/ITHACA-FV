#include "ReducedSteadyNS.H"
class bifurcationROM: public reducedSteadyNS
{
    public:
    bifurcationROM(steadyNS& Foamproblem, Eigen::MatrixXd vel): reducedSteadyNS(Foamproblem)
    {
        std::cout<<"constructor"<<std::endl;
//             if (problem->bcMethod == "lift")
//     {
//         vel_now = setOnlineVelocity(vel);
//     }
//     else if (problem->bcMethod == "penalty")
//     {
//         vel_now = vel;
//     }

     y.resize(Nphi_u + Nphi_p, 1);
     std::cout<<"y \n"<<y<<std::endl;
     y.setZero();
     std::cout<<"N_BC:    "<<N_BC<<std::endl;
     std::cout<<"vel \n"<<vel<<std::endl;
     std::cout<<"y \n"<<y<<std::endl;


     // Change initial condition for the lifting function
     if (problem->bcMethod == "lift")
     {
         for (int j = 0; j < N_BC; j++)
         {
             std::cout<<"end constructor1"<<std::endl;
             y(j) = vel(j, 0);
         }
     }
        std::cout<<"end constructor2"<<std::endl;
        newton_object.BC.resize(N_BC);

     for (int j = 0; j < N_BC; j++)
     {
         newton_object.BC(j) = vel(j, 0);
     }


        std::cout<<"end constructor3"<<std::endl;
    }

    void solveOnline()
    {
     Color::Modifier red(Color::FG_RED);
     Color::Modifier green(Color::FG_GREEN);
     Color::Modifier def(Color::FG_DEFAULT);
     newton_object.nu = nu;
     Eigen::HybridNonLinearSolver<newton_steadyNS> hnls(newton_object);

     hnls.solve(y);
     Eigen::VectorXd res(y);
     newton_object.operator()(y, res);
     Info << "################## Online solve N° " << count_online_solve <<
          " ##################" << endl;

     if (Pstream::master())
     {
         std::cout << "Solving for the parameter: " << nu << std::endl;
     }

     if (res.norm() < 1e-5 && Pstream::master())
     {
         std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                   hnls.iter << " iterations " << def << std::endl << std::endl;
     }
     else if (Pstream::master())
     {
         std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                   hnls.iter << " iterations " << def << std::endl << std::endl;
     }

     count_online_solve += 1;
 }
};
