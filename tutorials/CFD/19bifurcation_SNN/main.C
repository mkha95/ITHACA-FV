#include "bifurcation_NSS.H"





int main(int argc,char * argv[]){
    Bifurcation<steadyNS> bif(argc,argv);
    bif.offlineSolve();
    bif.prepare_POD();
    return 0;
}
