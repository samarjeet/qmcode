#include"Molecule.h"
#include<fstream>
#include<cmath>


void Molecule::readGeometry(std::string fileName){

  std::ifstream input(fileName);
  input >> nAtom ;
  int atomicNum;
  float x,y,z;

  std::vector<float> coord;

  for (int i=0; i < nAtom ; ++i)  {
    input >> atomicNum >> x >> y >> z;
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

