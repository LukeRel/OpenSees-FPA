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

// $Revision$
// $Date$
// $Source$

// Written: MHS
// Created: 2012
// Revision: 
//
// Description: This file contains the implementation for the
// NDFiberIS3d class. NDFiberIS3d provides the abstraction of an
// n-dimensional fiber with initial strains added to it.
//
// What: "@(#) NDFiberIS3d.h, revA"

#include <stdlib.h>

#include <NDMaterial.h>
#include <NDFiberIS3d.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ID.h>
#include <SectionForceDeformation.h>
#include <Information.h>
#include <Parameter.h>
#include <FiberResponse.h>
#include <elementAPI.h>

Matrix NDFiberIS3d::ks(6,6); 
Vector NDFiberIS3d::fs(6); 
ID NDFiberIS3d::code(6);

static int numNDFiberIS3d = 0;

void* OPS_NDFiberIS3d()
{
    // Number of arguments after "NDFiberIS":
    int numArgs = OPS_GetNumRemainingInputArgs();
    
    // Minimum arguments required
    if(numArgs < 4) {
	opserr<<"insufficient arguments for NDFiberIS3d. Required:\n";
    opserr<<"y position, z position, Area, Material tag.\nOptional: initial strains.\n";
	return 0;
    }

    // Get data [0]=y [1]=z [2]=A [3]=eps0
    int numData = 4;
    double data[4];
    data[3] = 0;
    if(OPS_GetDoubleInput(&numData,&data[0]) < 0) return 0;

    // Get mat tag
    int tag;
    numData = 1;
    if(OPS_GetIntInput(&numData,&tag) < 0) return 0;

    // Check
    opserr << tag << " " << data[0] << " " << data[1] << " " << data[2] << " " << data[3] << "\n";

    // get material
    NDMaterial* theMat = OPS_getNDMaterial(tag);
    if(theMat == 0) {
	opserr<<"invalid NDMaterial tag\n";
	return 0;
    }

    return new NDFiberIS3d(numNDFiberIS3d++,*theMat,data[2],data[0],data[1],data[3]);
}


// constructor:
NDFiberIS3d::NDFiberIS3d(int tag, NDMaterial &theMat, double Area, double yy, double zz, double Eps0):
  Fiber(tag, FIBER_TAG_NDIS3d), theMaterial(0), area(Area), y(yy), z(zz), eps0(Eps0)
{
  theMaterial = theMat.getCopy("BeamFiber");
  
  if (theMaterial == 0) {
    opserr << "NDFiberIS3d::NDFiberIS3d -- failed to get copy of NDMaterial\n";
    exit(-1);
  }
  
  if (code(0) != SECTION_RESPONSE_P) {
    code(0) = SECTION_RESPONSE_P;
    code(1) = SECTION_RESPONSE_MZ;
    code(2) = SECTION_RESPONSE_MY;
    code(3) = SECTION_RESPONSE_VY;
    code(4) = SECTION_RESPONSE_VZ;
    code(5) = SECTION_RESPONSE_T;
  }
}

// constructor for blank object that recvSelf needs to be invoked upon
NDFiberIS3d::NDFiberIS3d(): 
  Fiber(0, FIBER_TAG_NDIS3d), theMaterial(0), area(0), y(0.0), z(0.0), eps0(0.0)
{
  if (code(0) != SECTION_RESPONSE_P) {
    code(0) = SECTION_RESPONSE_P;
    code(1) = SECTION_RESPONSE_MZ;
    code(2) = SECTION_RESPONSE_MY;
    code(3) = SECTION_RESPONSE_VY;
    code(4) = SECTION_RESPONSE_VZ;
    code(5) = SECTION_RESPONSE_T;
  }
}


// Destructor: 
NDFiberIS3d::~NDFiberIS3d ()
{
  if (theMaterial != 0)
    delete theMaterial;
}

int   
NDFiberIS3d::setTrialFiberStrain(const Vector &vs)
{
  static Vector strain(3);
  strain(0) = 0;
  strain(1) = 0;
  strain(2) = 0;
  
  opserr << "NDFiberIS3d::setTrialFiberStrain() -- not implemented" << endln;

  return theMaterial->setTrialStrain(strain);
}

// get fiber stress resultants 
Vector &
NDFiberIS3d::getFiberStressResultants (void) 
{
  fs.Zero();
  
  opserr << "NDFiberIS3d::getFiberStressResultants() -- not implemented" << endln;

  return fs;
}

// get contribution of fiber to section tangent stiffness
Matrix &
NDFiberIS3d::getFiberTangentStiffContr(void) 
{
  ks.Zero();

  opserr << "NDFiberIS3d::getFiberTangentStiffContr() -- not implemented" << endln;

  return ks;
}

Fiber*
NDFiberIS3d::getCopy (void)
{
   // make a copy of the fiber 
  NDFiberIS3d *theCopy = new NDFiberIS3d (this->getTag(), 
				      *theMaterial, area, y, z, eps0);

  return theCopy;
}  

int
NDFiberIS3d::getOrder(void)
{
  return 6;
}

const ID&
NDFiberIS3d::getType(void)
{
  return code;
}

int   
NDFiberIS3d::commitState(void)
{
  return theMaterial->commitState();
}


int   
NDFiberIS3d::revertToLastCommit(void)
{
  return theMaterial->revertToLastCommit();
}


int   
NDFiberIS3d::revertToStart(void)
{
  return theMaterial->revertToStart();
}


int   
NDFiberIS3d::sendSelf(int commitTag, Channel &theChannel)
{
  // 
  // store tag and material info in an ID and send it
  //
  int res = 0;
  
  int dbTag = this->getDbTag();
  
  static ID idData(3);
  
  idData(0) = this->getTag();
  idData(1) = theMaterial->getClassTag();
  
  int matDbTag = theMaterial->getDbTag();
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    if (matDbTag != 0)
      theMaterial->setDbTag(matDbTag);
  }
  
  idData(2) = matDbTag;
  
  res += theChannel.sendID(dbTag, commitTag, idData);
  if (res < 0) {
    opserr << "NDFiberIS3d::sendSelf - failed to send ID data\n";
    return res;
  }    
  
  // 
  // store area and position data in a vector and send it
  //
  static Vector dData(4);
  
  dData(0) = area;
  dData(1) = y;
  dData(2) = z;
  dData(2) = eps0;
  
  res += theChannel.sendVector(dbTag, commitTag, dData);
  if (res < 0) {
    opserr << "NDFiberIS3d::sendSelf - failed to send Vector data\n";
    return res;
  }    

  // now invoke sendSelf on the material
  res += theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) {
    opserr << "NDFiberIS3d::sendSelf - failed to send UniaxialMaterial\n";
      return res;
  }
    
  return res;
}

int   
NDFiberIS3d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  // 
  // get tag and material info from an ID
  //
  
  int res = 0;
  
  int dbTag = this->getDbTag();
  
  static ID idData(3);
    
  res += theChannel.recvID(dbTag, commitTag, idData);
  if (res < 0) {
    opserr << "NDFiberIS3d::recvSelf - failed to receive ID data\n";
    return res;
  }    
  
  this->setTag(idData(0));

  // 
  // get area from a vector received from channel
  //
  
  static Vector dData(4);
  
  res += theChannel.recvVector(dbTag, commitTag, dData);
  if (res < 0) {
      opserr << "NDFiberIS3d::recvSelf - failed to receive Vector data\n";
      return res;
  }
  
  area = dData(0);
  y = dData(1);
  z = dData(2);
  eps0 = dData(3);

  //
  // now we do the material stuff
  //
  
  int matClassTag = idData(1);    
  
    // if we have a material, check it is of correct type
  if (theMaterial != 0) {
    if (matClassTag != theMaterial->getClassTag()) {
      delete theMaterial;
      theMaterial = 0;
    }
    }
  
  // if no material we need to get one,
  // NOTE: not an else if in case deleted in if above
  if (theMaterial == 0) {
    theMaterial = theBroker.getNewNDMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "NDFiberIS3d::recvSelf() - " <<
	  "failed to get a NDMaterial of type " << matClassTag << endln;
      return -1;
    }
  }
  
    // set the materials dbTag and invoke recvSelf on the material
  theMaterial->setDbTag(idData(2));
  
  // now invoke recvSelf on the material
  res += theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) {
    opserr << "NDFiberIS3d::recvSelf() - the material failed in recvSelf()\n";
    return res;
  }    	
  
  return res;
}


void NDFiberIS3d::Print(OPS_Stream &s, int flag)
{
  s << "\nNDFiberIS3d, tag: " << this->getTag() << endln;
  s << "\tArea: " << area << endln; 
  s << "\tLocation (y,z): " << y << " " << z << endln; 
  s << "\tGiven initial strains: " << eps0 << endln;
  s << "\tMaterial, tag: " << theMaterial->getTag() << endln;
}

Response*
NDFiberIS3d::setResponse(const char **argv, int argc, OPS_Stream &s)
{
  if (argc == 0)
    return 0;
  
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0)
    return new FiberResponse(this, 1, Vector(2));
  
  else
    return theMaterial->setResponse(argv, argc, s);
}

int
NDFiberIS3d::getResponse(int responseID, Information &fibInfo)
{
  switch(responseID) {
  case 1:
    return fibInfo.setVector(this->getFiberStressResultants());
    
  default:
    return -1;
  }
}

void 
NDFiberIS3d::getFiberLocation(double &yLoc, double &zLoc)
{
  yLoc = y;
  zLoc = z;
}

int
NDFiberIS3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strcmp(argv[0],"A") == 0)
    return param.addObject(1, this);

  if (strcmp(argv[0],"y") == 0)
    return param.addObject(2, this);

  if (strcmp(argv[0],"z") == 0)
    return param.addObject(3, this);

  if (strcmp(argv[0],"eps0") == 0)
    return param.addObject(4, this);

  else
    return theMaterial->setParameter(argv, argc, param);
}

int
NDFiberIS3d::updateParameter(int parameterID, Information &info)
{
  switch(parameterID) {
  case 1:
    area = info.theDouble;
    return 0;
  case 2:
    y = info.theDouble;
    return 0;
  case 3:
    z = info.theDouble;
    return 0;
  case 4:
    eps0 = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
NDFiberIS3d::activateParameter(int parameterID)
{
  return -1;
}

const Vector&
NDFiberIS3d::getFiberSensitivity(int gradNumber, bool cond)
{
  return Fiber::getFiberSensitivity(gradNumber, cond);
}

int 
NDFiberIS3d::commitSensitivity(const Vector &dedh, int gradNumber, int numGrads)
{
  return -1;
}

