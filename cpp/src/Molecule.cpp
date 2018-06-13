#include"Molecule.h"
#include"constants.h"
#include<fstream>
#include<cmath>
#include<numeric>

void Molecule::readGeometry(std::string fileName){

  std::ifstream input(fileName);
  input >> nAtom ;
  int atomicNum;
  float an,x,y,z;

  std::vector<float> coord;

  for (int i=0; i < nAtom ; ++i)  {
    input >> an >> x >> y >> z;
    atomicNum = (int)an;
    zvals.push_back(atomicNum);

    coord.clear();
    
    coord.push_back(x);
    coord.push_back(y);
    coord.push_back(z);
    
    geom.push_back(coord);
  }
  input.close();
}

float Molecule::bond(int i, int j) {

  return sqrt( 
    pow(geom[i][0] - geom[j][0],2) + 
    pow(geom[i][1] - geom[j][1],2) + 
    pow(geom[i][2] - geom[j][2],2) 
  
  );
}


float Molecule::angle(int i, int j, int k){
  float angle = 0.0;

  std::vector<float> eji, ejk;

  eji.push_back(-(geom[j][0] - geom[i][0] )/bond(i,j));
  eji.push_back(-(geom[j][1] - geom[i][1] )/bond(i,j));
  eji.push_back(-(geom[j][2] - geom[i][2] )/bond(i,j));

  ejk.push_back(-(geom[j][0] - geom[k][0] )/bond(k,j));
  ejk.push_back(-(geom[j][1] - geom[k][1] )/bond(k,j));
  ejk.push_back(-(geom[j][2] - geom[k][2] )/bond(k,j));

  angle = acos(
    eji[0]*ejk[0] + 
    eji[1]*ejk[1] + 
    eji[2]*ejk[2]  
    );
  return angle;
}

float Molecule::unitComp(int axis, int i, int j){
  return -(geom[i][axis] - geom[j][axis])/bond(i,j);
}


std::vector<float> Molecule::unit(int i, int j){
  std::vector<float> u;
  u.push_back(unitComp(0,i,j));
  u.push_back(unitComp(1,i,j));
  u.push_back(unitComp(2,i,j));
  return u;
}

std::vector<float> cross(std::vector<float> v1, std::vector<float> v2){
  std::vector<float> cpd;
  cpd.push_back(v1[1]*v2[2] - v1[2]*v2[1]);
  cpd.push_back(v1[2]*v2[0] - v1[0]*v2[2]);
  cpd.push_back(v1[0]*v2[1] - v1[1]*v2[0]);

  return cpd;
}

float dot(std::vector<float> v1, std::vector<float> v2){
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

float Molecule::oop(int i, int j, int k, int l){
  float a = angle(j,k,l);
  return  asin(dot(cross(unit(k,j), unit(k,l)), unit(k,i))/sin(a));
}

float Molecule::dihedral(int i, int j, int k, int l) {
   
  return acos( dot(cross(unit(i,j), unit(j,k)), cross(unit(j,k), unit(k,l)))/sin(angle(i,j,k))/sin(angle(j,k,l)) );;
}

std::vector<float> Molecule::com(){
  std::vector<float> c;
  c.push_back(0.0);
  c.push_back(0.0);
  c.push_back(0.0);

  int i=0;
  float totalMass = 0.0;
  for (auto g: geom){
   c[0] += g[0] *  mass[zvals[i]];
   c[1] += g[1] *  mass[zvals[i]];
   c[2] += g[2] *  mass[zvals[i]];
   totalMass +=  mass[zvals[i]]; 
   i++;
  }
  c[0]/=totalMass;
  c[1]/=totalMass;
  c[2]/=totalMass;
  return c;
}

void Molecule::translate(float xsh, float ysh, float zsh){
  for (auto& g: geom){
    g[0] += xsh;
    g[1] += ysh;
    g[2] += zsh;
  }
}

Eigen::Matrix3f Molecule::moment(){
  std::vector<std::vector<float>> tensor;
  tensor.push_back({0.0, 0.0, 0.0});
  tensor.push_back({0.0, 0.0, 0.0});
  tensor.push_back({0.0, 0.0, 0.0});

  int i=0;
  for (auto g: geom){
    tensor[0][0] += mass[zvals[i]]*(pow(g[1],2) + pow(g[2],2));
    tensor[1][1] += mass[zvals[i]]*(pow(g[2],2) + pow(g[0],2));
    tensor[2][2] += mass[zvals[i]]*(pow(g[0],2) + pow(g[1],2));

    tensor[0][1] += mass[zvals[i]] * g[0]*g[1];
    tensor[0][2] += mass[zvals[i]] * g[0]*g[2];
    tensor[1][2] += mass[zvals[i]] * g[1]*g[2];

    tensor[1][0] = tensor[0][1];
    tensor[2][0] = tensor[0][2];
    tensor[2][1] = tensor[1][2];

    i++;
  }
  ten(0,0) = tensor[0][0];
  ten(0,1) = tensor[0][1];
  ten(0,2) = tensor[0][2];
  ten(1,0) = tensor[1][0];
  ten(1,1) = tensor[1][1];
  ten(1,2) = tensor[1][2];
  ten(2,0) = tensor[2][0];
  ten(2,1) = tensor[2][1];
  ten(2,2) = tensor[2][2];

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> solver(ten);
  evecs = solver.eigenvectors();
  auto evals = solver.eigenvalues();
  return ten;
}

// Project 2 stuff

int Molecule::readHessian(std::string fileName){
  std::ifstream input(fileName);
  int natomsHessian;
  input >> natomsHessian ;

  if (nAtom != natomsHessian) {
    return 0;
  }

  float val;
  std::vector<float> row;
  hessEigen.resize(3*nAtom, 3*nAtom);
  for (int i=0; i < 3*nAtom ; ++i){
    row.clear();
    for (int j=0; j < 3*nAtom ; ++j){
        input >> val;
        row.push_back(val);
        hessEigen(i,j) = val;
    }
    hessian.push_back(row);
  }

  // mass weighting the matrix
  for (int i=0; i < 3*nAtom ; ++i){
    for (int j=0; j < 3*nAtom ; ++j){
      hessian[i][j] /= sqrt(mass[zvals[i/3]]*mass[zvals[j/3]]);
      hessEigen(i,j) /=sqrt(mass[zvals[i/3]]*mass[zvals[j/3]]);
      //hessian[i][j];
    }
  }

  // Diagonalize the matrix
  // use Eigen library
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> solver(hessEigen);
  diagEvecs.resize(3*nAtom, 3*nAtom);
  diagEvals.resize(3*nAtom);
  frequencies.resize(3*nAtom);

  diagEvecs = solver.eigenvectors();
  diagEvals = solver.eigenvalues();


  for (int i=0; i<3*nAtom; ++i){
    frequencies(i) = sqrt(diagEvals(i))*5140.43206058;
  }
  return 1;
}

void Molecule::readEnuc(std::string fileName){
  std::ifstream input(fileName);
  input >> enuc ;
  input.close();
}

void Molecule::readOverlap(std::string fileName){
  std::ifstream input(fileName);
  totalNumOrbitals = 0;
  for (auto z:zvals) {
    totalNumOrbitals += numOrbitals[z];
  }
  s.resize(totalNumOrbitals, totalNumOrbitals);
  int i, j;
  float overlap;
 
  for (int k=0; k < totalNumOrbitals*(totalNumOrbitals+1)/2; ++k){
    input >> i >> j >> overlap ;
    s(i-1,j-1) = overlap; 
    s(j-1,i-1) = overlap; 
  }
  input.close();
}

void Molecule::readKinetic(std::string fileName){
  std::ifstream input(fileName);
  totalNumOrbitals = 0;
  for (auto z:zvals) {
    totalNumOrbitals += numOrbitals[z];
  }
  t.resize(totalNumOrbitals, totalNumOrbitals);
  int i, j;
  float ke;
 
  for (int k=0; k < totalNumOrbitals*(totalNumOrbitals+1)/2; ++k){
    input >> i >> j >> ke ;
    t(i-1,j-1) = ke; 
    t(j-1,i-1) = ke; 
  }
  input.close();
}

void Molecule::readNuclearAttraction(std::string fileName){
  std::ifstream input(fileName);
  totalNumOrbitals = 0;
  for (auto z:zvals) {
    totalNumOrbitals += numOrbitals[z];
  }
  v.resize(totalNumOrbitals, totalNumOrbitals);
  int i, j;
  float na;
 
  for (int k=0; k < totalNumOrbitals*(totalNumOrbitals+1)/2; ++k){
    input >> i >> j >> na ;
    v(i-1,j-1) = na; 
    v(j-1,i-1) = na; 
  }
  input.close();
}

void Molecule::calculateCoreHamiltonian(){
  ch.resize(totalNumOrbitals, totalNumOrbitals);
  for (int i=0; i < totalNumOrbitals; ++i) {
    for (int j=0; j < totalNumOrbitals; ++j) {
      ch(i,j) = v(i,j) + t(i,j);
    }
  }
}

void Molecule::readTwoERepulsion(std::string fileName){
  std::ifstream input(fileName);
  totalNumOrbitals = 0;
  for (auto z:zvals) {
    totalNumOrbitals += numOrbitals[z];
  }
  input.close();
}


