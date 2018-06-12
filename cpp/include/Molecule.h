#ifndef MOLECULE_H
#define MOLECULE_H

#include"Eigen/Dense"
#include"Eigen/Eigenvalues"
#include"Eigen/Core"
#include<vector>
#include<string>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1, 1> Vector;

class Molecule {

public :
  int nAtom;
  int charge;
  std::vector<std::vector<float>> geom;
  std::vector<std::vector<float>> hessian;
  std::vector<int> zvals;

  Matrix I; 
  Eigen::Matrix3f t; 
  Eigen::Matrix3f evecs;
  //auto evalues; 

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
};

#endif
