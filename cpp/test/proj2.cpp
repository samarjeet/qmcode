#include"Molecule.h"
#include<iostream>
#include<cmath>
#include<string>

int main(int argc, char* argv[]){
  Molecule mol ;
  mol.readGeometry(argv[1]);
  int status;
  status = mol.readHessian(argv[2]);
  if (status == 0) {
    std::cout << "Error while reading the hessian!\n";
    return 0;
  }
  std::cout << "Eigenvalues of the mass weighted Hessian matrix :\n" << mol.diagEvals << "\n";
  std::cout << "Frequencies :\n" << mol.frequencies << "\n";
  return 0;
}


