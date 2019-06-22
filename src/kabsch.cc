#include"kabsch.h"
#include<iostream>

namespace MDToy{
  Eigen::Matrix3d kabsch(const Molecule &PP, const Molecule &QQ){
    int N;
    if(PP.atoms().size()==QQ.atoms().size() ) {
      N=PP.atoms().size();
    }else{
      throw std::runtime_error("Kabsch algorithm: numbers of atoms mismatch.");
    }

    Molecule P=PP;
    Molecule Q=QQ;

    P.setCenter0();
    Q.setCenter0();

    Eigen::Matrix3d H=Eigen::Matrix3d::Zero();
    for(int i=0; i!=3; ++i){
      for(int j=0; j!=3; ++j){
        for(int k=0; k!=N; ++k){
          if(P.atoms(k).element()!= Q.atoms(k).element()){
            throw std::runtime_error("Kabsch algorithm: types of atoms mismatch.");
          }
          H(i,j)+=P.atoms(k).xyz(i)*Q.atoms(k).xyz(j)*pow(P.phi()(P.atoms(k))*P.weight(k),2);
        }
      }
    }

    Eigen::JacobiSVD<Eigen::Matrix3d> svd(H , Eigen::ComputeFullU | Eigen::ComputeFullV);

    // Find the rotation
    double d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
    d=(d>0?1.0:-1.0);
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
    I(2, 2) = d;
    return svd.matrixV() * I * svd.matrixU().transpose();

  }

}
