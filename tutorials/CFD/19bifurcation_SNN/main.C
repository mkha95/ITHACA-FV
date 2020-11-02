#include "bifurcation_NSS.H"
#include "bifurcationROM_NSS.H"





int main(int argc,char * argv[]){
    Bifurcation<steadyNS> bif(argc,argv);
    bif.offlineSolve();
    bif.prepare_POD();

    Eigen::MatrixXd vel_now(1, 1);
    vel_now(0, 0) = 1;

    bifurcationROM ridotto(bif, vel_now);
     // Set the inlet velocity
     // used only for penalty approach
     ridotto.tauU = Eigen::MatrixXd::Zero(1, 1);
     ridotto.tauU(0, 0) = 1e-1;
     std::cout<<"party"<<std::endl;
     std::cout<<bif.mu<<std::endl;

     // Perform an online solve for the new values of inlet velocities
     for (label k = 0; k < bif.mu.size(); k++)
     {
         // Set the reduced viscosity
         ridotto.nu = bif.mu(0, k);
         std::cout<<"no problem ye1"<<std::endl;
         ridotto.solveOnline();
         std::cout<<"no problem ye2"<<std::endl;
         Eigen::MatrixXd tmp_sol(ridotto.y.rows() + 1, 1);
         tmp_sol(0) = k + 1;
         tmp_sol.col(0).tail(ridotto.y.rows()) = ridotto.y;
         ridotto.online_solution.append(tmp_sol);
     }

     // Save the online solution
     ITHACAstream::exportMatrix(ridotto.online_solution, "red_coeff", "python",
                                "./ITHACAoutput/red_coeff");
     ITHACAstream::exportMatrix(ridotto.online_solution, "red_coeff", "matlab",
                                "./ITHACAoutput/red_coeff");
     ITHACAstream::exportMatrix(ridotto.online_solution, "red_coeff", "eigen",
                                "./ITHACAoutput/red_coeff");
     // Reconstruct and export the solution
     ridotto.reconstruct(true, "./ITHACAoutput/Reconstruction/");
     return 0;
 }
