#include<vector>
#include<string>

class Molecule {

public :
  int nAtom;
  int charge;
  std::vector<std::vector<float>> geom;
  std::vector<int> zvals;

  Molecule(){};
  ~Molecule(){};

  void readGeometry(std::string fileName);
  float bond(int i, int j);  
  float angle(int i, int j, int k);  
};

