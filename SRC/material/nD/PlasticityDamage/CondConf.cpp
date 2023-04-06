/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.4 $
// $Date: 2003-02-14 23:01:24 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/CondConf.cpp,v $

// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the class definition of CondConf.
// The CondConf class is a wrapper class that performs static
// CondConf on a three-dimensional material model to give the 11, 12, and 13
// stress components which can then be integrated over an area to model a
// shear flexible 3D beam.


#include <CondConf.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string.h>
#include <elementAPI.h>
#include <MaterialResponse.h>

Vector CondConf::stress(3);
Matrix CondConf::tangent(3,3);

//      0  1  2  3  4  5
// ND: 11 22 33 12 23 31
// BF: 11 12 31 22 33 23

static int d_out = 0; // Set to 1 to produce a out_cond.txt output file
static double step = 0.0;

void* OPS_CondConf()
{
    int argc = OPS_GetNumRemainingInputArgs() + 2;
    if (argc < 4) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: nDMaterial BeamFiber tag? matTag?" << endln;
	return 0;
    }

    int tags[2];
    int numdata = 2;
    if (OPS_GetIntInput(&numdata, tags) < 0) {
	opserr << "WARNING invalid nDMaterial CondConf tag or matTag" << endln;
	return 0;
    }

    double data[3];
    numdata = 3;
    if (OPS_GetDoubleInput(&numdata, data) < 0) {
        opserr << "WARNING invalid nDMaterial. Want: CondConf tag, matTag, rho_y, rho_z, fy" << endln;
        return 0;
    }

    int tag = tags[0];
    int matTag = tags[1];
    double rho_y = data[0];
    double rho_z = data[1];
    double fy = data[2];

    NDMaterial *threeDMaterial = OPS_getNDMaterial(matTag);
    if (threeDMaterial == 0) {
	opserr << "WARNING nD material does not exist\n";
	opserr << "nD material: " << matTag;
	opserr << "\nBeamFiber nDMaterial: " << tag << endln;
	return 0;
    }

    return new CondConf(tag, *threeDMaterial, rho_y, rho_z, fy);
}

CondConf::CondConf(void)
: NDMaterial(0, ND_TAG_CondConf),
Tstrain22(0.0), Tstrain33(0.0), Tgamma23(0.0),
Cstrain22(0.0), Cstrain33(0.0), Cgamma23(0.0),
theMaterial(0), strain(3), strain_c(3), stress_c(3),
strain_k(3), stress_k(3), strain_c_k(3), stress_c_k(3),
damage(0.0), rho_y(0.0), rho_z(0.0), fy(0.0)
{
	// Nothing to do
}

CondConf::CondConf(int tag, NDMaterial& theMat, double _rho_y, double _rho_z, double _fy)
    : NDMaterial(tag, ND_TAG_CondConf),
    Tstrain22(0.0), Tstrain33(0.0), Tgamma23(0.0),
    Cstrain22(0.0), Cstrain33(0.0), Cgamma23(0.0),
    theMaterial(0), strain(3), strain_c(3), stress_c(3),
    strain_k(3), stress_k(3), strain_c_k(3), stress_c_k(3),
    damage(0.0), rho_y(_rho_y), rho_z(_rho_z), fy(_fy)
{
    // Get a copy of the material
    theMaterial = theMat.getCopy("ThreeDimensional");

    if (theMaterial == 0) {
        opserr << "CondConf::CondConf -- failed to get copy of material\n";
        exit(-1);
    }

    strain_c.Zero();
    stress_c.Zero();
    strain_k.Zero();
    stress_k.Zero();
    strain_c_k.Zero();
    stress_c_k.Zero();

    // Initializing damage output file
    if (d_out == 1) {
        using namespace std;
        ofstream outdata;
        outdata.open("out_cond.txt");
        outdata << "Step eps11 eps22 eps33 gam12 gam23 gam31 sig11 sig22 sig33 tau12 tau23 tau31 Dt Dc D count" << endln;
        outdata.close();
    }
}

CondConf::~CondConf(void) 
{ 
  if (theMaterial != 0)
    delete theMaterial;
} 

NDMaterial*
CondConf::getCopy(void) 
{
  CondConf *theCopy =
    new CondConf(this->getTag(), *theMaterial, rho_y, rho_z, fy);

  theCopy->Tstrain22 = this->Tstrain22;
  theCopy->Tstrain33 = this->Tstrain33;
  theCopy->Tgamma23  = this->Tgamma23;
  theCopy->Cstrain22 = this->Cstrain22;
  theCopy->Cstrain33 = this->Cstrain33;
  theCopy->Cgamma23  = this->Cgamma23;
  
  return theCopy;
}

NDMaterial* 
CondConf::getCopy(const char *type)
{
  if (strcmp(type, "BeamFiber") == 0)
    return this->getCopy();
  else
    return 0;
}

int 
CondConf::getOrder(void) const
{
  return 3;
}

const char*
CondConf::getType(void) const 
{
  return "BeamFiber";
}

int 
CondConf::commitState(void)
{
  Cstrain22 = Tstrain22;
  Cstrain33 = Tstrain33;
  Cgamma23 = Tgamma23;

  // Mantained (fiber ones)
  strain_k = strain;
  stress_k = stress;

  // Condensed out (non fiber ones)
  strain_c_k = strain_c;
  stress_c_k = stress_c;

  return theMaterial->commitState();
}

int 
CondConf::revertToLastCommit(void)
{
  Tstrain22 = Cstrain22;
  Tstrain33 = Cstrain33;
  Tgamma23 = Cgamma23;
  
  return theMaterial->revertToLastCommit();
}

int
CondConf::revertToStart()
{
  this->Tstrain22 = 0.0;
  this->Tstrain33 = 0.0;
  this->Tgamma23  = 0.0;
  this->Cstrain22 = 0.0;
  this->Cstrain33 = 0.0;
  this->Cgamma23  = 0.0;

  return theMaterial->revertToStart();
}

double
CondConf::getRho(void)
{
  return theMaterial->getRho();
}


//receive the strain
int 
CondConf::setTrialStrain(const Vector& strainFromElement)
{
    static const double tol = 1.0e-7;

    strain(0) = strainFromElement(0);
    strain(1) = strainFromElement(1);
    strain(2) = strainFromElement(2);

    //newton loop to solve for out-of-plane strains

    double norm = 1.0;
    static Vector dstrain_m(3);
    static Vector dstrain_c(3);
    static Vector dstress_c(3);
    static Vector threeDstrain(6);
    static Matrix C_cc(3, 3);
    static Matrix C_mc(3, 3);

    // Full strains vector = (eps_11 eps_22 eps_33 gamma_12 gamma_23 gamma_31)
    // Mantained     = 1 4 6 (eps_11 gamma_12 gamma_13) = received
    // Condensed out = 2 3 5 (eps_22 eps_33 eps_23) = must be determined and eventually condensed
    int m[3] = { 0, 3, 5 };
    int c[3] = { 1, 2, 4 };
    int perm[6] = { 0, 3, 5, 1, 2, 4 };

    // Confinement
    double sig_ty; double sig_tz;
    if (strain(0) < 0.0) {
        sig_ty = -rho_y * fy; sig_tz = -rho_z * fy;
    } else {
        sig_ty = 0.0; sig_tz = 0.0;
    }
    stress_c(0) = sig_ty;
    stress_c(1) = sig_tz;
    stress_c(2) = 0.0;

    // Get tangent
    const Matrix& C = theMaterial->getTangent();

    // Other matrices initialization
    C_mc.Zero(); C_cc.Zero();
    int k; int l;
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
        k = m[i]; l = c[j]; C_mc(i, j) = C(k, l);
        k = c[i]; l = c[j]; C_cc(i, j) = C(k, l);
    }

    //opserr << "strain_c_k = " << strain_k(1) << " " << strain_k(2) << " " << strain_k(4) << endln;
    //opserr << "strain_m = " << strain(0) << " " << strain(1) << " " << strain(2) << endln;
    //opserr << "strain_m_k = " << strain_k(0) << " " << strain_k(1) << " " << strain_k(2) << endln;
    //opserr << "stress_c = " << stress_c(0) << " " << stress_c(1) << " " << stress_c(2) << endln;
    //opserr << "stress_c_k = " << stress_c_k(0) << " " << stress_c_k(1) << " " << stress_c_k(2) << endln;

    // Trial strain
    dstrain_m = strain - strain_k;
    dstress_c = stress_c - stress_c_k;
    C_cc.Solve(dstress_c - C_mc * dstrain_m, dstrain_c);

    // Elastic predictor
    strain_c = strain_c_k - dstrain_c;

    int count = 0;
    const int maxCount = 200;

    do {
        //opserr << "Sub iteration n. " << count << endln;
      //set three dimensional strain
        threeDstrain(0) = strain(0);
        threeDstrain(1) = strain_c(0);
        threeDstrain(2) = strain_c(1);
        threeDstrain(3) = strain(1);
        threeDstrain(4) = strain_c(2);
        threeDstrain(5) = strain(2);

        if (theMaterial->setTrialStrain(threeDstrain) < 0) {
            opserr << "CondConf::setTrialStrain - setStrain failed in material with strain " << threeDstrain;
            return -1;
        }

        //three dimensional stress
        const Vector& threeDstress = theMaterial->getStress();

        //three dimensional tangent 
        const Matrix& threeDtangent = theMaterial->getTangent();

        //NDmaterial strain order = 11, 22, 33, 12, 23, 31  
        //CondConf strain order   = 11, 12, 31, 22, 33, 23

        /*
        C_cc(0, 0) = threeDtangent(1, 1);
        C_cc(1, 0) = threeDtangent(2, 1);
        C_cc(2, 0) = threeDtangent(4, 1);

        C_cc(0, 1) = threeDtangent(1, 2);
        C_cc(1, 1) = threeDtangent(2, 2);
        C_cc(2, 1) = threeDtangent(4, 2);

        C_cc(0, 2) = threeDtangent(1, 4);
        C_cc(1, 2) = threeDtangent(2, 4);
        C_cc(2, 2) = threeDtangent(4, 4);
        */

        C_cc.Zero();
        int k; int l;
        for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
            k = c[i]; l = c[j]; C_cc(i, j) = C(k, l);
        }


        // Update dstress_c
        dstress_c(0) = threeDstress(1) - stress_c(0);
        dstress_c(1) = threeDstress(2) - stress_c(1);
        dstress_c(2) = threeDstress(4) - stress_c(2);

        //set norm
        norm = dstress_c.Norm();

        // Solve c-out strains
        if (norm > tol)
            C_cc.Solve(dstress_c, dstrain_c);

        //update out of plane strains
        strain_c -= dstrain_c;

    } while (count++ < maxCount && norm > tol);

    if (count > maxCount - 10) {
        opserr << "norm = " << norm << endln;
        opserr << "count = " << count << endln;
        opserr << "strain = [" << strain(0) << " " << strain(1) << " " << strain(2) << " " << strain_c(0) << " " << strain_c(1) << " " << strain_c(2) << "];" << endln;
        opserr << "stress = [" << stress(0) << " " << stress(1) << " " << stress(2) << " " << stress_c(0) << " " << stress_c(1) << " " << stress_c(2) << "];" << endln;

        opserr << "strain_c = [" << strain_c(0) << " " << strain_c(1) << " " << strain_c(2) << "];" << endln;
        opserr << "stress_c = [" << stress_c(0) << " " << stress_c(1) << " " << stress_c(2) << "];" << endln;

        opserr << "C_cc = [" << C_cc(0,0) << " " << C_cc(1,0) << " " << C_cc(2,0) << ";" << endln;
        opserr << "        " << C_cc(0,1) << " " << C_cc(1,1) << " " << C_cc(2,1) << ";" << endln;
        opserr << "        " << C_cc(0,2) << " " << C_cc(1,2) << " " << C_cc(2,2) << "];" << endln;

    }

    if (d_out == 1) {
        step++;

        // Damage output
        Vector dam(6);
        dam = this->getDamage();
        using namespace std;
        ofstream outdata;
        outdata.open("out_cond.txt", ios::app);
        outdata << step << " ";
        outdata << dam(0) << " " << dam(1) << " " << dam(2) << " " << count << endln;
        outdata.close();
    }

    return 0;
}

const Vector& 
CondConf::getStrain(void)
{
  return strain;
}

const Vector&  
CondConf::getStress()
{
  const Vector &threeDstress = theMaterial->getStress();

  stress(0) = threeDstress(0);
  stress(1) = threeDstress(3);
  stress(2) = threeDstress(5);

  return stress;
}

const Vector&
CondConf::getCommittedStrain(void)
{
    const Vector& threeDstrain = theMaterial->getCommittedStrain();

    strain_k(0) = threeDstrain(0);
    strain_k(1) = threeDstrain(3);
    strain_k(2) = threeDstrain(5);

    return strain_k;
}

const Vector&
CondConf::getCommittedStress()
{
    const Vector& threeDstress = theMaterial->getCommittedStress();

    stress_k(0) = threeDstress(0);
    stress_k(1) = threeDstress(3);
    stress_k(2) = threeDstress(5);

    return stress_k;
}

const Vector& 
CondConf::getStressSensitivity(int gradIndex,
					bool conditional)
{
  const Vector &threeDstress = theMaterial->getStressSensitivity(gradIndex, conditional);

  stress(0) = threeDstress(0);
  stress(1) = threeDstress(3);
  stress(2) = threeDstress(5);

  const Matrix &threeDtangent = theMaterial->getTangent();

  static Matrix dd12(3,3);
  dd12(0,0) = threeDtangent(0,1);
  dd12(1,0) = threeDtangent(3,1);
  dd12(2,0) = threeDtangent(5,1);

  dd12(0,1) = threeDtangent(0,2);
  dd12(1,1) = threeDtangent(3,2);
  dd12(2,1) = threeDtangent(5,2);

  dd12(0,2) = threeDtangent(0,4);
  dd12(1,2) = threeDtangent(3,4);
  dd12(2,2) = threeDtangent(5,4);


  static Matrix C_cc(3,3);
  C_cc(0,0) = threeDtangent(1,1);
  C_cc(1,0) = threeDtangent(2,1);
  C_cc(2,0) = threeDtangent(4,1);
  
  C_cc(0,1) = threeDtangent(1,2);
  C_cc(1,1) = threeDtangent(2,2);
  C_cc(2,1) = threeDtangent(4,2);
  
  C_cc(0,2) = threeDtangent(1,4);
  C_cc(1,2) = threeDtangent(2,4);
  C_cc(2,2) = threeDtangent(4,4);

  
  static Vector sigma2(3);
  sigma2(0) = threeDstress(1);
  sigma2(1) = threeDstress(2);
  sigma2(2) = threeDstress(4);

  static Vector C_ccsigma2(3);
  C_cc.Solve(sigma2,C_ccsigma2);

  stress.addMatrixVector(1.0, dd12, C_ccsigma2, -1.0);

  return stress;
}

const Matrix&  
CondConf::getTangent()
{
  const Matrix &threeDtangent = theMaterial->getTangent();

  static Matrix C_mm(3,3);
  C_mm(0,0) = threeDtangent(0,0);
  C_mm(1,0) = threeDtangent(3,0);
  C_mm(2,0) = threeDtangent(5,0);

  C_mm(0,1) = threeDtangent(0,3);
  C_mm(1,1) = threeDtangent(3,3);
  C_mm(2,1) = threeDtangent(5,3);

  C_mm(0,2) = threeDtangent(0,5);
  C_mm(1,2) = threeDtangent(3,5);
  C_mm(2,2) = threeDtangent(5,5);


  static Matrix dd12(3,3);
  dd12(0,0) = threeDtangent(0,1);
  dd12(1,0) = threeDtangent(3,1);
  dd12(2,0) = threeDtangent(5,1);

  dd12(0,1) = threeDtangent(0,2);
  dd12(1,1) = threeDtangent(3,2);
  dd12(2,1) = threeDtangent(5,2);

  dd12(0,2) = threeDtangent(0,4);
  dd12(1,2) = threeDtangent(3,4);
  dd12(2,2) = threeDtangent(5,4);

  static Matrix dd21(3,3);
  dd21(0,0) = threeDtangent(1,0);
  dd21(1,0) = threeDtangent(2,0);
  dd21(2,0) = threeDtangent(4,0);

  dd21(0,1) = threeDtangent(1,3);
  dd21(1,1) = threeDtangent(2,3);
  dd21(2,1) = threeDtangent(4,3);

  dd21(0,2) = threeDtangent(1,5);
  dd21(1,2) = threeDtangent(2,5);
  dd21(2,2) = threeDtangent(4,5);

  static Matrix C_cc(3,3);
  C_cc(0,0) = threeDtangent(1,1);
  C_cc(1,0) = threeDtangent(2,1);
  C_cc(2,0) = threeDtangent(4,1);

  C_cc(0,1) = threeDtangent(1,2);
  C_cc(1,1) = threeDtangent(2,2);
  C_cc(2,1) = threeDtangent(4,2);

  C_cc(0,2) = threeDtangent(1,4);
  C_cc(1,2) = threeDtangent(2,4);
  C_cc(2,2) = threeDtangent(4,4);


  //int Solve(const Vector &V, Vector &res) const;
  //int Solve(const Matrix &M, Matrix &res) const;
  //CondConf 
  static Matrix C_ccinvdd21(3,3);
  C_cc.Solve(dd21, C_ccinvdd21);

  //this->tangent   = C_mm; 
  //this->tangent  -= (dd12*C_ccinvdd21);
  C_mm.addMatrixProduct(1.0, dd12, C_ccinvdd21, -1.0);
  tangent = C_mm; 

  return tangent;
}

const Matrix&  
CondConf::getInitialTangent()
{
  const Matrix &threeDtangent = theMaterial->getInitialTangent();

  static Matrix C_mm(3,3);
  C_mm(0,0) = threeDtangent(0,0);
  C_mm(1,0) = threeDtangent(3,0);
  C_mm(2,0) = threeDtangent(5,0);

  C_mm(0,1) = threeDtangent(0,3);
  C_mm(1,1) = threeDtangent(3,3);
  C_mm(2,1) = threeDtangent(5,3);

  C_mm(0,2) = threeDtangent(0,5);
  C_mm(1,2) = threeDtangent(3,5);
  C_mm(2,2) = threeDtangent(5,5);


  static Matrix dd12(3,3);
  dd12(0,0) = threeDtangent(0,1);
  dd12(1,0) = threeDtangent(3,1);
  dd12(2,0) = threeDtangent(5,1);

  dd12(0,1) = threeDtangent(0,2);
  dd12(1,1) = threeDtangent(3,2);
  dd12(2,1) = threeDtangent(5,2);

  dd12(0,2) = threeDtangent(0,4);
  dd12(1,2) = threeDtangent(3,4);
  dd12(2,2) = threeDtangent(5,4);

  static Matrix dd21(3,3);
  dd21(0,0) = threeDtangent(1,0);
  dd21(1,0) = threeDtangent(2,0);
  dd21(2,0) = threeDtangent(4,0);

  dd21(0,1) = threeDtangent(1,3);
  dd21(1,1) = threeDtangent(2,3);
  dd21(2,1) = threeDtangent(4,3);

  dd21(0,2) = threeDtangent(1,5);
  dd21(1,2) = threeDtangent(2,5);
  dd21(2,2) = threeDtangent(4,5);

  static Matrix C_cc(3,3);
  C_cc(0,0) = threeDtangent(1,1);
  C_cc(1,0) = threeDtangent(2,1);
  C_cc(2,0) = threeDtangent(4,1);

  C_cc(0,1) = threeDtangent(1,2);
  C_cc(1,1) = threeDtangent(2,2);
  C_cc(2,1) = threeDtangent(4,2);

  C_cc(0,2) = threeDtangent(1,4);
  C_cc(1,2) = threeDtangent(2,4);
  C_cc(2,2) = threeDtangent(4,4);


  //int Solve(const Vector &V, Vector &res) const;
  //int Solve(const Matrix &M, Matrix &res) const;
  //CondConf 
  static Matrix C_ccinvdd21(3,3);
  C_cc.Solve(dd21, C_ccinvdd21);

  //this->tangent   = C_mm; 
  //this->tangent  -= (dd12*C_ccinvdd21);
  C_mm.addMatrixProduct(1.0, dd12, C_ccinvdd21, -1.0);
  tangent = C_mm;

  return tangent;
}

const Vector& CondConf::getDamage(void) {

    return theMaterial->getDamage();
}

/*
Response* CondConf::setResponse(const char** argv, int argc, OPS_Stream& s) {
    
    if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0)
        return new MaterialResponse(this, 1, this->getStress());

    else if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0)
        return new MaterialResponse(this, 2, this->getStrain());

    else if (strcmp(argv[0], "tangent") == 0 || strcmp(argv[0], "Tangent") == 0)
        return new MaterialResponse(this, 3, this->getTangent());

	else if (strcmp(argv[0], "damage") == 0 || strcmp(argv[0], "Damage") == 0)
		return new MaterialResponse(this, 5, this->getDamage());

	else
		return 0;
}

int CondConf::getResponse(int responseID, Information& matInfo) {
    
    switch (responseID) {
    
    case -1:
        return -1;
    case 1:
        if (matInfo.theVector != 0)
            *(matInfo.theVector) = this->getStress();
        return 0;

    case 2:
        if (matInfo.theVector != 0)
            *(matInfo.theVector) = this->getStrain();
        return 0;

    case 3:
        if (matInfo.theMatrix != 0)
            *(matInfo.theMatrix) = this->getTangent();
        return 0;
    
    case 5:
        matInfo.setDouble(this->getDamage());
        return 0;

    }

    return 0;
}
*/

void  
CondConf::Print(OPS_Stream &s, int flag)
{
  s << "CondConf, tag: " << this->getTag() << endln;
  s << "\tWrapped material: "<< theMaterial->getTag() << endln;

  theMaterial->Print(s, flag);
}

int 
CondConf::sendSelf(int commitTag, Channel &theChannel) 
{
  int res = 0;

  // put tag and associated materials class and database tags into an id and send it
  static ID idData(3);
  idData(0) = this->getTag();
  idData(1) = theMaterial->getClassTag();
  int matDbTag = theMaterial->getDbTag();
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    theMaterial->setDbTag(matDbTag);
  }
  idData(2) = matDbTag;

  res = theChannel.sendID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    opserr << "CondConf::sendSelf() - failed to send id data\n";
    return res;
  }

  // put the strains in a vector and send it
  static Vector vecData(3);
  vecData(0) = Cstrain22;
  vecData(1) = Cstrain33;
  vecData(2) = Cgamma23;

  res = theChannel.sendVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) {
    opserr << "CondConf::sendSelf() - failed to send vector data\n";
    return res;
  }

  // now send the materials data
  res = theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) 
    opserr << "CondConf::sendSelf() - failed to send vector material\n";

  return res;
}

int 
CondConf::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;

  // recv an id containing the tag and associated materials class and db tags
  static ID idData(3);
  res = theChannel.sendID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    opserr << "CondConf::sendSelf() - failed to send id data\n";
    return res;
  }

  this->setTag(idData(0));
  int matClassTag = idData(1);
  
  // if the associated material has not yet been created or is of the wrong type
  // create a new material for recvSelf later
  if (theMaterial == 0 || theMaterial->getClassTag() != matClassTag) {
    if (theMaterial != 0)
      delete theMaterial;
    theMaterial = theBroker.getNewNDMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "CondConf::recvSelf() - failed to get a material of type: " << matClassTag << endln;
      return -1;
    }
  }
  theMaterial->setDbTag(idData(2));

  // recv a vector containing strains and set the strains
  static Vector vecData(3);
  res = theChannel.recvVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) {
    opserr << "CondConf::sendSelf() - failed to send vector data\n";
    return res;
  }

  Cstrain22 = vecData(0);
  Cstrain33 = vecData(1);
  Cgamma23  = vecData(2);

  Tstrain22 = Cstrain22;
  Tstrain33 = Cstrain33;
  Tgamma23  = Cgamma23;

  // now receive the materials data
  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) 
    opserr << "CondConf::sendSelf() - failed to send vector material\n";
  
  return res;
}

int
CondConf::setParameter(const char **argv, int argc,
				Parameter &param)
{
  return theMaterial->setParameter(argv, argc, param);
}
