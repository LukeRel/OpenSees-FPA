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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/BeamFiberMaterialEB.cpp,v $

// Written: MHS
// Created: Aug 2001
// Modified: LP - May 2022
//
// Description: This file contains the class definition of BeamFiberMaterialEB.
// The BeamFiberMaterialEB class is a wrapper class that performs static
// condensation on a three-dimensional material model to give only the 11
// stress components which can then be integrated over an area to model an
// Eulero Bernoulli 3D beam.


#include <BeamFiberMaterialEB.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string.h>
#include <elementAPI.h>
#include <MaterialResponse.h>

Vector BeamFiberMaterialEB::stress(3);
Matrix BeamFiberMaterialEB::tangent(3,3);

//      0  1  2  3  4  5
// ND: 11 22 33 12 23 31
// BF: 11 22 33 12 23 31

static int d_out = 0; // Set to 1 to produce a out_J2damage.txt output file
static double step = 0.0;

void* OPS_BeamFiberMaterialEB()
{
    int argc = OPS_GetNumRemainingInputArgs() + 2;
    if (argc < 4) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: nDMaterial BeamFiberEB tag? matTag?" << endln;
	return 0;
    }

    int tags[2];
    int numdata = 2;
    if (OPS_GetIntInput(&numdata, tags) < 0) {
	opserr << "WARNING invalid nDMaterial BeamFiberEB tag or matTag" << endln;
	return 0;
    }

    int tag = tags[0];
    int matTag = tags[1];

    NDMaterial *threeDMaterial = OPS_getNDMaterial(matTag);
    if (threeDMaterial == 0) {
	opserr << "WARNING nD material does not exist\n";
	opserr << "nD material: " << matTag;
	opserr << "\nBeamFiberEB nDMaterial: " << tag << endln;
	return 0;
    }

    return new BeamFiberMaterialEB(tag, *threeDMaterial);
}

BeamFiberMaterialEB::BeamFiberMaterialEB(void)
: NDMaterial(0, ND_TAG_BeamFiberMaterialEB),
Tstrain22(0.0), Tstrain33(0.0), Tgamma12(0.0), Tgamma23(0.0), Tgamma31(0.0),
Cstrain22(0.0), Cstrain33(0.0), Cgamma12(0.0), Cgamma23(0.0), Cgamma31(0.0),
theMaterial(0), strain(3), strain3D(6), stress3D(6),
damage(0.0)
{
	// Nothing to do
}

BeamFiberMaterialEB::BeamFiberMaterialEB(int tag, NDMaterial &theMat)
: NDMaterial(tag, ND_TAG_BeamFiberMaterialEB),
Tstrain22(0.0), Tstrain33(0.0), Tgamma12(0.0), Tgamma23(0.0), Tgamma31(0.0),
Cstrain22(0.0), Cstrain33(0.0), Cgamma12(0.0), Cgamma23(0.0), Cgamma31(0.0),
theMaterial(0), strain(3),
damage(0.0)
{
  // Get a copy of the material
  theMaterial = theMat.getCopy("ThreeDimensional");
  
  if (theMaterial == 0) {
    opserr << "BeamFiberMaterialEB::BeamFiberMaterialEB -- failed to get copy of material\n";
    exit(-1);
  }
  // Initializing damage output file
  if (d_out == 1) {
      using namespace std;
      ofstream outdata;
      outdata.open("out_J2Damage.txt");
      outdata << "Step eps11 eps22 eps33 gam12 gam23 gam31 sig11 sig22 sig33 tau12 tau23 tau31 Dt Dc D" << endln;
      outdata.close();
  }
}

BeamFiberMaterialEB::~BeamFiberMaterialEB(void) 
{ 
  if (theMaterial != 0)
    delete theMaterial;
} 

NDMaterial*
BeamFiberMaterialEB::getCopy(void) 
{
  BeamFiberMaterialEB *theCopy =
    new BeamFiberMaterialEB(this->getTag(), *theMaterial);

  theCopy->Tstrain22 = this->Tstrain22;
  theCopy->Tstrain33 = this->Tstrain33;
  theCopy->Tgamma23  = this->Tgamma12;
  theCopy->Tgamma23  = this->Tgamma23;
  theCopy->Tgamma23  = this->Tgamma31;
  theCopy->Cstrain22 = this->Cstrain22;
  theCopy->Cstrain33 = this->Cstrain33;
  theCopy->Cgamma23  = this->Cgamma12;
  theCopy->Cgamma23  = this->Cgamma23;
  theCopy->Cgamma23  = this->Cgamma31;
  
  return theCopy;
}

NDMaterial* 
BeamFiberMaterialEB::getCopy(const char *type)
{
  if (strcmp(type, "BeamFiber") == 0)
    return this->getCopy();
  else
    return 0;
}

int 
BeamFiberMaterialEB::getOrder(void) const
{
  return 3;
}

const char*
BeamFiberMaterialEB::getType(void) const 
{
  return "BeamFiber";
}

int 
BeamFiberMaterialEB::commitState(void)
{
  Cstrain22 = Tstrain22;
  Cstrain33 = Tstrain33;
  Cgamma12 = Tgamma12;
  Cgamma23 = Tgamma23;
  Cgamma31 = Tgamma31;

  //step += ops_Dt;

  //opserr << "Step = " << step << endln;

  return theMaterial->commitState();
}

int 
BeamFiberMaterialEB::revertToLastCommit(void)
{
  Tstrain22 = Cstrain22;
  Tstrain33 = Cstrain33;
  Tgamma12 = Cgamma12;
  Tgamma23 = Cgamma23;
  Tgamma31 = Cgamma31;
  
  return theMaterial->revertToLastCommit();
}

int
BeamFiberMaterialEB::revertToStart()
{
  this->Tstrain22 = 0.0;
  this->Tstrain33 = 0.0;
  this->Tgamma12  = 0.0;
  this->Tgamma23  = 0.0;
  this->Tgamma31  = 0.0;
  this->Cstrain22 = 0.0;
  this->Cstrain33 = 0.0;
  this->Cgamma12  = 0.0;
  this->Cgamma23  = 0.0;
  this->Cgamma31  = 0.0;

  return theMaterial->revertToStart();
}

double
BeamFiberMaterialEB::getRho(void)
{
  return theMaterial->getRho();
}


//receive the strain
int 
BeamFiberMaterialEB::setTrialStrain(const Vector &strainFromElement)
{
  static const double tolerance = 1.0e-10;

  strain(0) = strainFromElement(0);
  strain(1) = 0.0;// strainFromElement(1);
  strain(2) = 0.0;// strainFromElement(2);

  //newton loop to solve for out-of-plane strains

  double norm;
  static Vector condensedStress(5);
  static Vector strainIncrement(5);
  static Vector threeDstrain(6);
  static Matrix dd22(5,5);

  int count = 0;
  const int maxCount = 2;
  double norm0;

  do {
      //opserr << "Sub iteration n. " << count << endln;
    //set three dimensional strain
    threeDstrain(0) = this->strain(0);
    threeDstrain(1) = this->Tstrain22;
    threeDstrain(2) = this->Tstrain33;
    threeDstrain(3) = this->strain(1); 
    threeDstrain(4) = this->Tgamma23;
    threeDstrain(5) = this->strain(2);

    if (theMaterial->setTrialStrain(threeDstrain) < 0) {
      opserr << "BeamFiberMaterialEB::setTrialStrain - setStrain failed in material with strain " << threeDstrain;
      return -1;   
    }

    //three dimensional stress
    const Vector &threeDstress = theMaterial->getStress();

    //delete LP
    //strain3D = threeDstrain;
    //stress3D = threeDstress;

    //three dimensional tangent 
    const Matrix &threeDtangent = theMaterial->getTangent();

    //NDmaterial strain order          = 11, 22, 33, 12, 23, 31  
    //BeamFiberMaterialEB strain order = 11, 22, 33, 12, 23, 31

    for (int i = 0;i < 5; i++) condensedStress(i) = threeDstress(i+1);

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) dd22(j,i) = threeDtangent(j+1,i+1);
    }

    //set norm
    norm = condensedStress.Norm();
    if (count == 0)
      norm0 = norm;

    //condensation 
    dd22.Solve(condensedStress, strainIncrement);

    //update out of plane strains
    this->Tstrain22 -= strainIncrement(0);
    this->Tstrain33 -= strainIncrement(1);
    this->Tgamma12  -= strainIncrement(2);
    this->Tgamma23  -= strainIncrement(3);
    this->Tgamma31  -= strainIncrement(4);

  } while (count++ < maxCount && norm > tolerance);

  //theMaterial->commitState();

  if (d_out == 1) {
      step++;

      // Damage output
      Vector dam(6);
      dam = this->getDamage();
      using namespace std;
      ofstream outdata;
      outdata.open("out_J2Damage.txt", ios::app);
      outdata << step << " ";
      outdata << strain3D(0) << " " << strain3D(1) << " " << strain3D(2) << " " << strain3D(3) << " " << strain3D(4) << " " << strain3D(5) << " ";
      outdata << stress3D(0) << " " << stress3D(1) << " " << stress3D(2) << " " << stress3D(3) << " " << stress3D(4) << " " << stress3D(5) << " ";
      outdata << dam(0) << " " << dam(1) << " " << dam(2) << endln;
      outdata.close();
  }

  return 0;
}

const Vector& 
BeamFiberMaterialEB::getStrain(void)
{
  return strain;
}

const Vector&  
BeamFiberMaterialEB::getStress()
{
  const Vector &threeDstress = theMaterial->getStress();

  stress(0) = threeDstress(0);
  stress(1) = 0.0;// threeDstress(3);
  stress(2) = 0.0;// threeDstress(5);

  return stress;
}

const Vector& 
BeamFiberMaterialEB::getStressSensitivity(int gradIndex,
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


  static Matrix dd22(3,3);
  dd22(0,0) = threeDtangent(1,1);
  dd22(1,0) = threeDtangent(2,1);
  dd22(2,0) = threeDtangent(4,1);
  
  dd22(0,1) = threeDtangent(1,2);
  dd22(1,1) = threeDtangent(2,2);
  dd22(2,1) = threeDtangent(4,2);
  
  dd22(0,2) = threeDtangent(1,4);
  dd22(1,2) = threeDtangent(2,4);
  dd22(2,2) = threeDtangent(4,4);

  
  static Vector sigma2(3);
  sigma2(0) = threeDstress(1);
  sigma2(1) = threeDstress(2);
  sigma2(2) = threeDstress(4);

  static Vector dd22sigma2(3);
  dd22.Solve(sigma2,dd22sigma2);

  stress.addMatrixVector(1.0, dd12, dd22sigma2, -1.0);

  return stress;
}

const Matrix&  
BeamFiberMaterialEB::getTangent()
{
  const Matrix &threeDtangent = theMaterial->getTangent();

  static Matrix dd11(3,3);
  dd11(0,0) = threeDtangent(0,0);
  dd11(1,0) = threeDtangent(3,0);
  dd11(2,0) = threeDtangent(5,0);

  dd11(0,1) = threeDtangent(0,3);
  dd11(1,1) = threeDtangent(3,3);
  dd11(2,1) = threeDtangent(5,3);

  dd11(0,2) = threeDtangent(0,5);
  dd11(1,2) = threeDtangent(3,5);
  dd11(2,2) = threeDtangent(5,5);


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

  static Matrix dd22(3,3);
  dd22(0,0) = threeDtangent(1,1);
  dd22(1,0) = threeDtangent(2,1);
  dd22(2,0) = threeDtangent(4,1);

  dd22(0,1) = threeDtangent(1,2);
  dd22(1,1) = threeDtangent(2,2);
  dd22(2,1) = threeDtangent(4,2);

  dd22(0,2) = threeDtangent(1,4);
  dd22(1,2) = threeDtangent(2,4);
  dd22(2,2) = threeDtangent(4,4);


  //int Solve(const Vector &V, Vector &res) const;
  //int Solve(const Matrix &M, Matrix &res) const;
  //condensation 
  static Matrix dd22invdd21(3,3);
  dd22.Solve(dd21, dd22invdd21);

  //this->tangent   = dd11; 
  //this->tangent  -= (dd12*dd22invdd21);
  dd11.addMatrixProduct(1.0, dd12, dd22invdd21, -1.0);
  tangent = dd11; 

  return tangent;
}

const Matrix&  
BeamFiberMaterialEB::getInitialTangent()
{
  const Matrix &threeDtangent = theMaterial->getInitialTangent();

  static Matrix dd11(3,3);
  dd11(0,0) = threeDtangent(0,0);
  dd11(1,0) = threeDtangent(3,0);
  dd11(2,0) = threeDtangent(5,0);

  dd11(0,1) = threeDtangent(0,3);
  dd11(1,1) = threeDtangent(3,3);
  dd11(2,1) = threeDtangent(5,3);

  dd11(0,2) = threeDtangent(0,5);
  dd11(1,2) = threeDtangent(3,5);
  dd11(2,2) = threeDtangent(5,5);


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

  static Matrix dd22(3,3);
  dd22(0,0) = threeDtangent(1,1);
  dd22(1,0) = threeDtangent(2,1);
  dd22(2,0) = threeDtangent(4,1);

  dd22(0,1) = threeDtangent(1,2);
  dd22(1,1) = threeDtangent(2,2);
  dd22(2,1) = threeDtangent(4,2);

  dd22(0,2) = threeDtangent(1,4);
  dd22(1,2) = threeDtangent(2,4);
  dd22(2,2) = threeDtangent(4,4);


  //int Solve(const Vector &V, Vector &res) const;
  //int Solve(const Matrix &M, Matrix &res) const;
  //condensation 
  static Matrix dd22invdd21(3,3);
  dd22.Solve(dd21, dd22invdd21);

  //this->tangent   = dd11; 
  //this->tangent  -= (dd12*dd22invdd21);
  dd11.addMatrixProduct(1.0, dd12, dd22invdd21, -1.0);
  tangent = dd11;

  return tangent;
}

const Vector& BeamFiberMaterialEB::getDamage(void) {

    return theMaterial->getDamage();
}

/*
Response* BeamFiberMaterialEB::setResponse(const char** argv, int argc, OPS_Stream& s) {
    
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

int BeamFiberMaterialEB::getResponse(int responseID, Information& matInfo) {
    
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
BeamFiberMaterialEB::Print(OPS_Stream &s, int flag)
{
  s << "BeamFiberMaterialEB, tag: " << this->getTag() << endln;
  s << "\tWrapped material: "<< theMaterial->getTag() << endln;

  theMaterial->Print(s, flag);
}

int 
BeamFiberMaterialEB::sendSelf(int commitTag, Channel &theChannel) 
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
    opserr << "BeamFiberMaterialEB::sendSelf() - failed to send id data\n";
    return res;
  }

  // put the strains in a vector and send it
  static Vector vecData(5);
  vecData(0) = Cstrain22;
  vecData(1) = Cstrain33;
  vecData(2) = Cgamma12;
  vecData(3) = Cgamma23;
  vecData(4) = Cgamma31;

  res = theChannel.sendVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) {
    opserr << "BeamFiberMaterialEB::sendSelf() - failed to send vector data\n";
    return res;
  }

  // now send the materials data
  res = theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) 
    opserr << "BeamFiberMaterialEB::sendSelf() - failed to send vector material\n";

  return res;
}

int 
BeamFiberMaterialEB::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;

  // recv an id containing the tag and associated materials class and db tags
  static ID idData(3);
  res = theChannel.sendID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    opserr << "BeamFiberMaterialEB::sendSelf() - failed to send id data\n";
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
      opserr << "BeamFiberMaterialEB::recvSelf() - failed to get a material of type: " << matClassTag << endln;
      return -1;
    }
  }
  theMaterial->setDbTag(idData(2));

  // recv a vector containing strains and set the strains
  static Vector vecData(5);
  res = theChannel.recvVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) {
    opserr << "BeamFiberMaterialEB::sendSelf() - failed to send vector data\n";
    return res;
  }

  Cstrain22 = vecData(0);
  Cstrain33 = vecData(1);
  Cgamma12  = vecData(2);
  Cgamma23  = vecData(3);
  Cgamma31  = vecData(4);

  Tstrain22 = Cstrain22;
  Tstrain33 = Cstrain33;
  Tgamma12  = Cgamma12;
  Tgamma23  = Cgamma23;
  Tgamma31  = Cgamma31;

  // now receive the materials data
  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) 
    opserr << "BeamFiberMaterialEB::sendSelf() - failed to send vector material\n";
  
  return res;
}

int
BeamFiberMaterialEB::setParameter(const char **argv, int argc,
				Parameter &param)
{
  return theMaterial->setParameter(argv, argc, param);
}
