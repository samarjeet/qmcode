#include"Molecule.h"
#include<iostream>
#include<cmath>
#include<string>

void printCoords(Molecule m){
  for (auto g: m.geom){
    std::cout << g[0] << " " << 
    g[1] << " " <<
    g[2] << "\n";
  }
}

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
          if (mol.bond(i,j) < 4.0 && mol.bond(j,k) < 4.0 ) {
            std::cout <<  i  << " " << j << " " << " " << k << " "<< 
                mol.angle(i,j,k)*(180.0/acos(-1.0)) << "\n"; 
          }
            
      }
    }
  }

  // Out of plane
  std::cout << "OOP:\n";
  for (int i=0; i < mol.nAtom; ++i){
    for (int j=0; j < mol.nAtom; ++j){
      for (int k=0; k < mol.nAtom; ++k){
        for (int l=k+1; l< mol.nAtom; ++l) {
          if (i!=j && i!=k && i!= l && j!=k && j!= l && mol.bond(k,l) < 4.0 && mol.bond(j,k) < 4.0 && mol.bond(i,k) < 4.0) {
            std::cout <<  i  << " " << j << " " << " " << k << " " << l << " " << 
                mol.oop(i,j,k,l)*(180.0/acos(-1.0)) << "\n"; 
          
          }
            
        }
      }
    }
  }

    

  // Difedrals
  std::cout << "\n Dihedral:\n";
  for (int i=0; i < mol.nAtom; ++i){
    for (int j=i+1; j < mol.nAtom; ++j){
      for (int k=j+1; k < mol.nAtom; ++k){
        for (int l=k+1; l< mol.nAtom; ++l) {
          if (i!=j && i!=k && i!= l && j!=k && j!= l && mol.bond(k,l) < 4.0 && mol.bond(j,k) < 4.0 && mol.bond(i,k) < 4.0) {
            std::cout <<  i  << " " << j << " " << " " << k << " " << l << " " << 
                mol.dihedral(i,j,k,l)*(180.0/acos(-1.0)) << "\n"; 
          
          }
            
        }
      }
    }
  }


  // COM
  std::vector<float> com = mol.com();

  // translate 
  mol.translate(-com[0], -com[1], -com[2]);
  printCoords(mol);
  
  // Moment of interta
  auto tensor = mol.moment();
  std::cout << tensor << "\n";
  std::cout << mol.evecs << "\n";

  // Rotation
  return 0;
}


