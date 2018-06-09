#include"Molecule.h"
#include<iostream>
#include<cmath>
#include<string>

int main(int argc, char* argv[]){
  Molecule mol ;
  mol.readGeometry(argv[1]);

  // Bonds
  for (int i=0; i < mol.nAtom; ++i){
    for (int j=i+1; j < mol.nAtom; ++j){
      std::cout <<  i  << " " << j << " " << mol.bond(i,j) << "\n"; 
    }
  }

  // Angles
  for (int i=0; i < mol.nAtom; ++i){
    for (int j=i+1; j < mol.nAtom; ++j){
      for (int k=j+1; k < mol.nAtom; ++k){
        if (mol.bond(i,j) < 4.0 && mol.bond(j,k) < 4.0)
          std::cout <<  i  << " " << j << " " << " " << k << " " <<mol.angle(i,j,k)*(180.0/acos(-1.0)) << "\n"; 
      }
    }
  }
  return 0;
}
