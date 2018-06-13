#ifndef MOLECULE_H
#define MOLECULE_H

#include"Eigen/Dense"
#include"Eigen/Eigenvalues"
#include"Eigen/Core"
#include<vector>
#include<string>

//typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXf;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1, 1> Vector;

class Molecule {

public :
  int nAtom;
  int charge;
  std::vector<std::vector<float>> geom;
  std::vector<std::vector<float>> hessian;
  std::vector<int> zvals;

  //Matrix I; 
  Eigen::Matrix3f ten; 
  Eigen::Matrix3f evecs;
  Eigen::MatrixXf hessEigen; //(3*nAtom, 3*nAtom);
  Eigen::MatrixXf diagEvecs;
  Eigen::VectorXf diagEvals;
  Eigen::VectorXf frequencies;
  Molecule(){};
  ~Molecule(){};

  void readGeometry(std::string fileName);
  int readHessian(std::string fileName); // returns 0 if error
  float bond(int i, int j);  
  float angle(int i, int j, int k);  
  float oop(int i, int j, int k, int l);  
  float dihedral(int i, int j, int k, int l);  
  float unitComp(int axis, int i, int j);
  std::vector<float> unit(int i, int j);
  std::vector<float> com();
  void translate(float xsh, float ysh, float zsh);
  //std::vector<std::vector<float>> moment();
  //Matrix moment();
  Eigen::Matrix3f moment();

  // for scf
  float enuc;
  int totalNumOrbitals;
  Eigen::MatrixXf s; // overlap 
  Eigen::MatrixXf t; // kinetic energy
  Eigen::MatrixXf v; // nuclear attraction integrals
  Eigen::MatrixXf ch; // core hamiltonian 
  std::vector<float> eri; // two Electron Repulsion;

  void readEnuc(std::string fileName);
  void readOverlap(std::string fileName);
  void readKinetic(std::string fileName);
  void readNuclearAttraction(std::string fileName);
  void calculateCoreHamiltonian();
  void readTwoERepulsion(std::string fileName);
};

#endif
