import math

class Orbital :
    def __init__(self, Ax, Ay, Az, alpha) :
        self.x0 = Ax;
        self.y0 = Ay;
        self.zo = Az;
        self.alpha = alpha;
        self.N = math.pow(2*alpha/math.pi,0.75);

class Basis :
    """
        n : # of gaussians
        c : contraction coefficients for each gaussian 
        g : array of Orbitals 
    """
    def __init__(self, numGaussians, contractions, orbitals):
        self.orbitals = orbitals;
        self.numGaussians  = numGaussians;
        self.contractions = contractions;

def buildSOrbital(Ax,Ay,Az, alpha):
    return Orbital(Ax,Ay,Az,alpha); 

def main2():
    z = [4, 1, 1, 1];
    al = [[3.0,4.5], [], [], []];

    nuclearRepulsion = buildNuclearRepulsion(z, al);
    basis, n = buildBasis(z,al)
    s = buildOverlap(basis);
    t = buildKinetic(basis);
    nuclearAttraction = buildNuclearAttraction(z, al, basis);
    gabcd = buildElectronRepulsion(basis);
    ho = t + nuclearAttraction;

    minEnergy = scf(ho, gabcd, s, n);
    minEnergy = minEnergy + nuclearRepulsion(z,al);

"""
z is an array of atomic numbers
al is #atoms X 3(x,y,z)
"""
def buildBasis(z, al):
    natoms = len(z);
    n = 0;
    numBasis = 0;

    basis = [];
    for i in range(natoms):
        x0 = al[i][0];
        y0 = al[i][1];
        z0 = al[i][2];

        if z[i] == 1 :
            # Hydrogen
            n = n + 1;
            s = [[18.7311370, 0.0334946], [2.8253937, 0.234769], [0.6401217, 0.81375733]];
            o0 = Orbital(x0,y0,z0,s[0][0])
            o1 = Orbital(x0,y0,z0,s[1][0])
            o2 = Orbital(x0,y0,z0,s[2][0])
            b = Basis(3, s[:1], [o0, o1, o2]);
            basis.append(b);

            b = Basis(1, [1], [Orbital(x0, y0, z0, 0.1612778)]);
            basis.append(b);
    return [basis, n]

def buildOverlap(basis):
  # calculate the overlap of the basis sets

  for ba in  basis :
    for bb in basis :
      #print(ba.numGaussians, bb.numGaussians)
      
def main():
    z = [1, 1, 1];
    al = [[0,0,1], [0, 0, -1], [0,0,3]]
    basis, n  = buildBasis(z, al);

    s = buildOverlap(basis);
if __name__ == '__main__' :
    main();
