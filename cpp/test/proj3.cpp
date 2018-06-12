#include"Molecule.h"
#include<iostream>
#include<cmath>
#include<string>

int main(int argc, char* argv[]){
  Molecule mol ;
  std::string path = "../test/h2o_sto3g/";
  mol.readGeometry(path+"geom.dat");
  mol.readEnuc(path+"enuc.dat");
  mol.readOverlap(path+"s.dat");
  mol.readKinetic(path+"t.dat");
  std::cout << mol.t << "\n" ;
  mol.readNuclearAttraction(path+"v.dat");
  std::cout << mol.v << "\n" ;
  mol.calculateCoreHamiltonian();
  std::cout << "Core Hamiltonian : \n" << mol.ch << "\n" ;
  mol.read2eOverlap(path+"v.dat");
  return 0;
}


