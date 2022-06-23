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
// Modified: 03/22 - L. Parente
// 
// Description: This file contains the class implementation of FiberSection3d.
//
// The class is modified to include additional consistent fiber strains.
// These strains are defined in the Tcl input script and
// they are added to usual fiber deformations after the effects from
// the section deformation are computed:
//    eps_x = a*eps_s + eps_0
// It can be used to emulate prestress, shrinkage, viscosity and other
// long term effects.
// Deformations entity must be set within the Tcl script.

// To be added: gamma parabolic effects in the tangent methods and the other ones

#include <NDFiberSection3d.h>

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixUtil.h>
#include <Fiber.h>
#include <classTags.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <NDMaterial.h>
#include <SectionIntegration.h>
#include <Parameter.h>
#include <elementAPI.h>

// For damage output
#include <iostream>
#include <fstream>
#include <Domain.h>

ID NDFiberSection3d::code(6);

// Settings
static int set_shear = 1; // Set to 1 to produce a parabolic distribution for shear deformations
static int d_out = 1; // Set to 1 to produce a damage.txt output file

void* OPS_NDFiberSection3d()
{
    int numData = OPS_GetNumRemainingInputArgs();
    if(numData < 1) {
	opserr<<"insufficient arguments for NDFiberSection3d\n";
	return 0;
    }

    numData = 1;
    int tag;
    if (OPS_GetIntInput(&numData,&tag) < 0) return 0;

    bool computeCentroid = true;
    if (OPS_GetNumRemainingInputArgs() > 0) {
      const char* opt = OPS_GetString();
      if (strcmp(opt, "-noCentroid") == 0)
	computeCentroid = false;
    }

    int num = 30;
    return new NDFiberSection3d(tag, num, computeCentroid);
}

// constructors:
NDFiberSection3d::NDFiberSection3d(int tag, int num, Fiber** fibers, double a, bool compCentroid) :
    SectionForceDeformation(tag, SEC_TAG_NDFiberSection3d),
    numFibers(num), sizeFibers(num), theMaterials(0), matData(0), gammaData(0),
    Abar(0.0), QyBar(0.0), QzBar(0.0), yBar(0.0), zBar(0.0), computeCentroid(compCentroid),
    alpha(a), sectionIntegr(0), e(6), s(0), ks(0), t0(0),
    parameterID(0), dedh(6), theDomain(0), step(0.0), eps0_creep(0), eps0_creep_k(0), deps0_creep(0)
{
    if (numFibers != 0) {
        theMaterials = new NDMaterial * [numFibers];

        if (theMaterials == 0) {
            opserr << "NDFiberSection3d::NDFiberSection3d -- failed to allocate Material pointers";
            exit(-1);
        }

        matData = new double[numFibers * 5];
        if (matData == 0) {
            opserr << "NDFiberSection3d::NDFiberSection3d -- failed to allocate double array for material data\n";
            exit(-1);
        }
        gammaData = new double[numFibers * 2];       // Gamma parabolic scale factors (LP)
        if (gammaData == 0) {
            opserr << "NDFiberSection3d::NDFiberSection3d -- failed to allocate double array for gamma scale factors\n";
            exit(-1);
        }
        // Creep
        eps0_creep = new double[numFibers];
        eps0_creep_k = new double[numFibers];
        deps0_creep = new double[numFibers];
        t0 = new double[numFibers];
        if (eps0_creep == 0 || eps0_creep_k == 0 || deps0_creep == 0) {
            opserr << "NDFiberSection3d::NDFiberSection3d -- failed to allocate double array for material data\n";
            exit(-1);
        }
        double yLoc, zLoc, Area, eps0, beta;

        // Top and bottom coords determination
        double y1, y2, z1, z2;

        for (int i = 0; i < numFibers; i++) {
            Fiber* theFiber = fibers[i];
            theFiber->getFiberLocation(yLoc, zLoc);
            Area = theFiber->getArea();
            eps0 = theFiber->getEps0();   // LP
            beta = theFiber->getBeta();   // LP

            //if (eps0 > 0) opserr << "Received eps0 = " << eps0 << " on fiber " << i+1;

            Abar += Area;
            QzBar += yLoc * Area;
            QyBar += zLoc * Area;
            matData[5 * i] = yLoc;
            matData[5 * i + 1] = zLoc;
            matData[5 * i + 2] = Area;
            matData[5 * i + 3] = eps0;       // LP
            matData[5 * i + 4] = beta;       // LP
            NDMaterial* theMat = theFiber->getNDMaterial();
            theMaterials[i] = theMat->getCopy("BeamFiber");

            //opserr << matData[i * 4] << " " << matData[i * 4 + 1] << " " << matData[i * 4 + 2] << " " << matData[i * 4 + 3] << "\n";

            // Top and bottom coords determination
            if (i == 0) {
                y1 = yLoc; y2 = yLoc;
                z1 = zLoc; z2 = zLoc;
            }
            y1 = fmin(yLoc, y1); y2 = fmax(yLoc, y2);
            z1 = fmin(zLoc, z1); z2 = fmax(zLoc, z2);

            if (theMaterials[i] == 0) {
                opserr << "NDFiberSection3d::NDFiberSection3d -- failed to get copy of a Material\n";
                exit(-1);
            }

            // Creep
            eps0_creep[i] = 0.0;       // LP
            eps0_creep_k[i] = 0.0;     // LP
            deps0_creep[i] = 0.0;     // LP
            t0[i] = 0.0;

        }

        // Scale factors
        for (int i = 0; i < numFibers; i++) {
            Fiber* theFiber = fibers[i];
            theFiber->getFiberLocation(yLoc, zLoc);
            if (set_shear == 1) {
                gammaData[2 * i] = 6 / pow(y1 - y2, 2) * (-yLoc * yLoc + yLoc * (y1 + y2) - y1 * y2);
                gammaData[2 * i + 1] = 6 / pow(z1 - z2, 2) * (-zLoc * zLoc + zLoc * (z1 + z2) - z1 * z2);
            }
            else {
                gammaData[2 * i] = 1;
                gammaData[2 * i + 1] = 1;
            }
        }

        if (computeCentroid) {
            yBar = QzBar / Abar;
            zBar = QyBar / Abar;
        }
    }

    s = new Vector(sData, 6);
    ks = new Matrix(kData, 6, 6);

    for (int i = 0; i < 6; i++)
        sData[i] = 0.0;

    for (int i = 0; i < 6 * 6; i++)
        kData[i] = 0.0;

    code(0) = SECTION_RESPONSE_P;
    code(1) = SECTION_RESPONSE_MZ;
    code(2) = SECTION_RESPONSE_MY;
    code(3) = SECTION_RESPONSE_VY;
    code(4) = SECTION_RESPONSE_VZ;
    code(5) = SECTION_RESPONSE_T;

    theDomain = new Domain;

    // Initializing damage output file
    if (d_out == 1) {
        using namespace std;
        ofstream outdata;
        outdata.open("out_fibers.txt");
        outdata << "Step y z eps11 gam12 gam13 sig11 tau12 tau13 Dt Dc D" << endln;
        outdata.close();
        outdata.open("out_sections.txt");
        outdata << "Step ex ky gz N My Vz E" << endln;
        outdata.close();
    }
}

NDFiberSection3d::NDFiberSection3d(int tag, int num, double a, bool compCentroid) :
    SectionForceDeformation(tag, SEC_TAG_NDFiberSection3d),
    numFibers(0), sizeFibers(num), theMaterials(0), matData(0), gammaData(0),
    Abar(0.0), QyBar(0.0), QzBar(0.0), yBar(0.0), zBar(0.0), computeCentroid(compCentroid),
    alpha(a), sectionIntegr(0), e(6), s(0), ks(0), t0(0),
    parameterID(0), dedh(6), theDomain(0), step(0.0), eps0_creep(0), eps0_creep_k(0), deps0_creep(0) // LP
{
    if (sizeFibers != 0) {
        theMaterials = new NDMaterial * [sizeFibers];

        if (theMaterials == 0) {
            opserr << "NDFiberSection3d::NDFiberSection3d -- failed to allocate Material pointers";
            exit(-1);
        }

        matData = new double[sizeFibers * 5];
        if (matData == 0) {
            opserr << "NDFiberSection3d::NDFiberSection3d -- failed to allocate double array for material data\n";
            exit(-1);
        }
        gammaData = new double[sizeFibers * 2];
        if (gammaData == 0) {
            opserr << "NDFiberSection3d::NDFiberSection3d -- failed to allocate double array for shear factor data\n";
            exit(-1);
        }
        // Creep
        eps0_creep = new double[numFibers];
        eps0_creep_k = new double[numFibers];
        deps0_creep = new double[numFibers];
        t0 = new double[numFibers];
        if (eps0_creep == 0 || eps0_creep_k == 0 || deps0_creep == 0) {
            opserr << "NDFiberSection3d::NDFiberSection3d -- failed to allocate double array for material data\n";
            exit(-1);
        }

        for (int i = 0; i < sizeFibers; i++) {
            matData[5 * i] = 0.0;
            matData[5 * i + 1] = 0.0;
            matData[5 * i + 2] = 0.0;
            matData[5 * i + 3] = 0.0;
            matData[5 * i + 4] = 0.0;
            theMaterials[i] = 0;
            gammaData[2 * i] = 0.0;
            gammaData[2 * i + 1] = 0.0;
            eps0_creep[i] = 0.0;
            eps0_creep_k[i] = 0.0;
            deps0_creep[i] = 0.0;
            t0[i] = 0.0;
        }
    }

    s = new Vector(sData, 6);
    ks = new Matrix(kData, 6, 6);

    for (int i = 0; i < 6; i++) sData[i] = 0.0;

    for (int i = 0; i < 6 * 6; i++) kData[i] = 0.0;

    code(0) = SECTION_RESPONSE_P;
    code(1) = SECTION_RESPONSE_MZ;
    code(2) = SECTION_RESPONSE_MY;
    code(3) = SECTION_RESPONSE_VY;
    code(4) = SECTION_RESPONSE_VZ;
    code(5) = SECTION_RESPONSE_T;
}

NDFiberSection3d::NDFiberSection3d(int tag, int num, NDMaterial **mats,
				   SectionIntegration &si, double a, bool compCentroid):
  SectionForceDeformation(tag, SEC_TAG_NDFiberSection3d),
  numFibers(num), sizeFibers(num), theMaterials(0), matData(0), gammaData(0),
  Abar(0.0), QyBar(0.0), QzBar(0.0), yBar(0.0), zBar(0.0), computeCentroid(compCentroid),
  alpha(a), sectionIntegr(0), e(6), s(0), ks(0), t0(0),
  parameterID(0), dedh(6), theDomain(0), step(0.0), eps0_creep(0), eps0_creep_k(0), deps0_creep(0)
{
  if (numFibers != 0) {
    theMaterials = new NDMaterial *[numFibers];

    if (theMaterials == 0) {
      opserr << "NDFiberSection3d::NDFiberSection3d -- failed to allocate Material pointers";
      exit(-1);
    }
    matData = new double [numFibers*5];
    if (matData == 0) {
      opserr << "NDFiberSection3d::NDFiberSection3d -- failed to allocate double array for material data\n";
      exit(-1);
    }
    gammaData = new double[numFibers * 2];
    if (gammaData == 0) {
        opserr << "NDFiberSection3d::NDFiberSection3d -- failed to allocate double array for shear factor data\n";
        exit(-1);
    }
  }

  sectionIntegr = si.getCopy();
  if (sectionIntegr == 0) {
    opserr << "Error: NDFiberSection3d::NDFiberSection3d: could not create copy of section integration object" << endln;
    exit(-1);
  }

  static double yLocs[10000];
  static double zLocs[10000];
  sectionIntegr->getFiberLocations(numFibers, yLocs, zLocs);
  
  static double fiberArea[10000];
  sectionIntegr->getFiberWeights(numFibers, fiberArea);

  for (int i = 0; i < numFibers; i++) {

    Abar  += fiberArea[i];
    QzBar += yLocs[i]*fiberArea[i];
    QyBar += zLocs[i]*fiberArea[i];

    theMaterials[i] = mats[i]->getCopy("BeamFiber");
    
    if (theMaterials[i] == 0) {
      opserr << "NDFiberSection3d::NDFiberSection3d -- failed to get copy of a Material\n";
      exit(-1);
    }
  }    

  if (computeCentroid) {
    yBar = QzBar/Abar;  
    zBar = QyBar/Abar;
  }

  s = new Vector(sData, 6);
  ks = new Matrix(kData, 6, 6);
  
  for (int i = 0; i < 6; i++) sData[i] = 0.0;

  for (int i = 0; i < 6*6; i++) kData[i] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_MY;
  code(3) = SECTION_RESPONSE_VY;
  code(4) = SECTION_RESPONSE_VZ;
  code(5) = SECTION_RESPONSE_T;

  theDomain = new Domain;
}

// constructor for blank object that recvSelf needs to be invoked upon
NDFiberSection3d::NDFiberSection3d():
  SectionForceDeformation(0, SEC_TAG_NDFiberSection3d),
  numFibers(0), sizeFibers(0), theMaterials(0), matData(0), gammaData(0),
  Abar(0.0), QyBar(0.0), QzBar(0.0), yBar(0.0), zBar(0.0), computeCentroid(true),
  alpha(1.0), sectionIntegr(0), e(6), s(0), ks(0), t0(0),
  parameterID(0), dedh(6), theDomain(0), step(0.0), eps0_creep(0), eps0_creep_k(0), deps0_creep(0)
{
  s = new Vector(sData, 6);
  ks = new Matrix(kData, 6, 6);

  for (int i = 0; i < 6; i++) sData[i] = 0.0;

  for (int i = 0; i < 6*6; i++) kData[i] = 0.0;
   
  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_MY;
  code(3) = SECTION_RESPONSE_VY;
  code(4) = SECTION_RESPONSE_VZ;
  code(5) = SECTION_RESPONSE_T;

  theDomain = new Domain;
}

int
NDFiberSection3d::addFiber(Fiber &newFiber)
{
  // need to create larger arrays
  if(numFibers == sizeFibers) {
      opserr << "AddFiber started" << endln;

      int newSize = 2*sizeFibers;
      NDMaterial **newArray = new NDMaterial *[newSize]; 

      double *newMatData = new double [5 * newSize];
      if (newArray == 0 || newMatData == 0) {
          opserr << "NDFiberSection3d::addFiber -- failed to allocate Fiber pointers\n";
          return -1;
      }
      double* newGammaData = new double[2 * newSize];
      if (newArray == 0 || newGammaData == 0) {
          opserr << "NDFiberSection3d::addFiber -- failed to allocate Fiber pointers\n";
          return -1;
      }
      
      // copy the old pointers and data
      for (int i = 0; i < numFibers; i++) {
	  newArray[i] = theMaterials[i];
	  newMatData[5*i] = matData[5*i];
	  newMatData[5*i+1] = matData[5*i+1];
	  newMatData[5*i+2] = matData[5*i+2];
      newMatData[5*i+3] = matData[5*i+3];
      newMatData[5*i+4] = matData[5*i+4];
      newGammaData[2*i] = gammaData[5*i];
      newGammaData[2*i+1] = gammaData[2*i+1];
      }

      // initialize new memory
      for (int i = numFibers; i < newSize; i++) {
	  newArray[i] = 0;
	  newMatData[5*i] = 0.0;
	  newMatData[5*i+1] = 0.0;
	  newMatData[5*i+2] = 0.0;
      newMatData[5*i+3] = 0.0;
      newMatData[5*i+4] = 0.0;
      newGammaData[2*i] = 0.0;
      newGammaData[2*i+1] = 0.0;
      }

      sizeFibers = newSize;

      // set new memory
      if (theMaterials != 0) {
	  delete [] theMaterials;
	  delete [] matData;
      delete [] gammaData;
      }

      theMaterials = newArray;
      matData = newMatData;
      gammaData = newGammaData;
  }

  // set the new pointers and data
  double yLoc, zLoc, Area, eps0, beta;
  newFiber.getFiberLocation(yLoc, zLoc);
  Area = newFiber.getArea();
  eps0 = newFiber.getEps0();
  beta = newFiber.getBeta();
  matData[numFibers*5] = yLoc;
  matData[numFibers*5+1] = zLoc;
  matData[numFibers*5+2] = Area;
  matData[numFibers*5+3] = eps0;
  matData[numFibers*5+4] = beta;

  NDMaterial *theMat = newFiber.getNDMaterial();
  theMaterials[numFibers] = theMat->getCopy("BeamFiber");

  if (theMaterials[numFibers] == 0) {
    opserr <<"NDFiberSection3d::addFiber -- failed to get copy of a Material\n";
    return -1;
  }

  numFibers++;

  // Recompute centroid
  if (computeCentroid) {
    Abar  += Area;
    QzBar += yLoc*Area;
    QyBar += zLoc*Area;
    
    yBar = QzBar/Abar;
    zBar = QyBar/Abar;
  }
  
  return 0;
}


// destructor:
NDFiberSection3d::~NDFiberSection3d()
{
  if (theMaterials != 0) {
    for (int i = 0; i < numFibers; i++)
      if (theMaterials[i] != 0)
	delete theMaterials[i];
    delete [] theMaterials;
  }

  if (matData != 0)         delete [] matData;
  if (eps0_creep != 0)      delete [] eps0_creep;
  if (eps0_creep_k != 0)    delete [] eps0_creep_k;
  if (deps0_creep != 0)     delete [] deps0_creep;
  if (t0 != 0)              delete [] t0;
  if (gammaData != 0)       delete [] gammaData;
  if (s != 0)               delete s;
  if (ks != 0)              delete ks;
  if (sectionIntegr != 0)   delete sectionIntegr;
  if (theDomain != 0)       delete theDomain;
}

// a = [1 -y z       0       0  0
//      0  0 0 sqrt(a)       0 -z
//      0  0 0       0 sqrt(a)  y]
int
NDFiberSection3d::setTrialSectionDeformation (const Vector &deforms)
{
  int res = 0;

  e = deforms;

  s->Zero();
  ks->Zero();

  double d0 = deforms(0);
  double d1 = deforms(1);
  double d2 = deforms(2);
  double d3 = deforms(3);
  double d4 = deforms(4);
  double d5 = deforms(5);

  static double yLocs[10000];
  static double zLocs[10000];
  static double fiberArea[10000];
  static double iStrains[10000];
  static double betaAngles[10000];

  // Shear factors arrays
  static double yScale[10000];
  static double zScale[10000];

  //opserr << "y z A eps0:\n";

  if (sectionIntegr != 0) {
    sectionIntegr->getFiberLocations(numFibers, yLocs, zLocs);
    sectionIntegr->getFiberWeights(numFibers, fiberArea);
  }  
  else {
    for (int i = 0; i < numFibers; i++) {
      yLocs[i] = matData[5*i];
      zLocs[i] = matData[5*i+1];
      fiberArea[i] = matData[5*i+2];
      iStrains[i] = matData[5*i+3];
      betaAngles[i] = matData[5*i+4];
      yScale[i] = gammaData[2*i];
      zScale[i] = gammaData[2*i+1];

      //opserr << yLocs[5*i] << " " << zLocs[5*i+1] << " " << fiberArea[5*i+2] << " " << iStrains[5*i+3] << "\n";
    }
  }
  //opserr << "\n";
  // Fiber strains
  static Vector eps(3);
  static Vector stress_0(3);

  // Shear shape factor
  double rootAlpha = 1.0;
  if (alpha != 1.0)
    rootAlpha = sqrt(alpha);

  // Cycle on every fiber
  for (int i = 0; i < numFibers; i++) {
    // Material, position, area
    NDMaterial *theMat = theMaterials[i];
    double y = yLocs[i] - yBar;
    double z = zLocs[i] - zBar;
    double A = fiberArea[i];
    double sgy = yScale[i];
    double sgz = zScale[i];
    //opserr << "sgy = " << gammaData[2 * i] << "; sgz = " << gammaData[2 * i + 1] << "\n";
    
    // Additional fiber strains
    double eps0 = iStrains[i];
    double beta = betaAngles[i];

    // Products
    double y2 = y*y;
    double z2 = z*z;
    double yz = y*z;
    double tmp;

    // Determine material strain and set it
    eps(0) = d0 - y * d1 + z * d2; // Axial strain eps_G
    eps(1) = sgy * rootAlpha * d3 - z * d5; // Shear strain gamma_12
    eps(2) = sgz * rootAlpha * d4 + y * d5; // Shear strain gamma_13

    // Adding initial strains
    double cosBeta = cos(beta);
    double sinBeta = sin(beta);

    // Set fiber strains to evaluate stress and tangent
    if (eps0 != 0) eps(0) += eps0;
    res += theMat->setTrialStrain(eps);
    const Vector& stress = theMat->getStress();
    const Matrix& tangent = theMat->getTangent();

    /* Prestress as external loads - Unused
    stress_0.Zero();
    int ps_set = 0; // 0 = internal load, 1 = external
    if (ps_set == 1) {
        theMat->setTrialStrain(eps);
        const Vector& stress_p = theMat->getStress();

        // Loads to be carried out
        stress_0 = stress_p - stress;
        // stress_p is the total stress;
        // stress   is the stress coming from the fiber deformation only
        // stress_0 is the difference, so it's the contribute of the external strains only
    }*/

    // Section tangent matrix construction - K(x)_s
    double d00 = tangent(0,0)*A;
    double d01 = tangent(0,1)*A;
    double d02 = tangent(0,2)*A;
    double d10 = tangent(1,0)*A;
    double d11 = tangent(1,1)*A;
    double d12 = tangent(1,2)*A;
    double d20 = tangent(2,0)*A;
    double d21 = tangent(2,1)*A;
    double d22 = tangent(2,2)*A;

    Matrix &ksi = *ks;
    Vector &si = *s;

    // Bending terms
    ksi(0,0) += d00;
    ksi(1,1) += y2*d00;
    ksi(2,2) += z2*d00;
    tmp = -y*d00;
    ksi(0,1) += tmp;
    ksi(1,0) += tmp;
    tmp = z*d00;
    ksi(0,2) += tmp;
    ksi(2,0) += tmp;
    tmp = -yz*d00;
    ksi(1,2) += tmp;
    ksi(2,1) += tmp;
    
    // Shear terms
    ksi(3,3) += sgy*sgy*alpha*d11;
    ksi(3,4) += sgy*sgz*alpha*d12;
    ksi(4,3) += sgz*sgz*alpha*d21;
    ksi(4,4) += sgz*sgy*alpha*d22;
    
    // Torsion term
    ksi(5,5) += z2*d11 - yz*(d12+d21) + y2*d22;
    
    // Bending-torsion coupling terms
    tmp = -z*d01 + y*d02;
    ksi(0,5) += tmp;
    ksi(1,5) -= y*tmp;
    ksi(2,5) += z*tmp;
    tmp = -z*d10 + y*d20;
    ksi(5,0) += tmp;
    ksi(5,1) -= y*tmp;
    ksi(5,2) += z*tmp;
    
    // Hit tangent terms with rootAlpha
    d01 *= rootAlpha; d02 *= rootAlpha;
    d10 *= rootAlpha; d11 *= rootAlpha; d12 *= rootAlpha;
    d20 *= rootAlpha; d21 *= rootAlpha; d22 *= rootAlpha;
    
    // Bending-shear coupling terms
    ksi(0,3) += sgy*d01;
    ksi(0,4) += sgz*d02;
    ksi(1,3) -= sgy*y*d01;
    ksi(1,4) -= sgz*y*d02;
    ksi(2,3) += sgy*z*d01;
    ksi(2,4) += sgz*z*d02;
    ksi(3,0) += sgy*d10;
    ksi(4,0) += sgz*d20;
    ksi(3,1) -= sgy*y*d10;
    ksi(4,1) -= sgz*y*d20;
    ksi(3,2) += sgy*z*d10;
    ksi(4,2) += sgz*z*d20;
    
    // Torsion-shear coupling terms
    y2 =  y*d22;
    z2 = -z*d11;
    ksi(5,3) +=  sgy*( z2 + y*d21);
    ksi(5,4) +=  sgz*(-z*d12 + y2);
    ksi(3,5) +=  sgy*( z2 + y*d12);
    ksi(4,5) +=  sgz*(-z*d21 + y2);

    // Section stress vector - sigma(x)_s - [N Mz My Vy Vz T]'
    double sig0 = stress(0)*A;
    double sig1 = stress(1)*A;
    double sig2 = stress(2)*A;

    /* Additional stress - Unused
    if (eps0 != 0) {
        sig0 += stress_0(0) * A;
        sig1 += stress_0(1) * A;
        sig2 += stress_0(2) * A;

    }*/

    // Rotations if necessary
    if ((beta > 1e-4) || (beta < 1e-4)) {
        //opserr << "Received beta is " << beta << " on fiber " << i+1 << endln;
        //opserr << "Its cosine is " << cosBeta << endln;
        //opserr << "Its sine is " << sinBeta << endln;
        // Normal fiber stress rotation
        sig0 = sig0 * cosBeta -sig2 * sinBeta;
        sig2 = sig0 * sinBeta +sig2 * cosBeta;
        //sig1 = sig1; // No rotations needed because beta is about y-axis
    }

    si(0) += sig0;
    si(1) += -y*sig0;
    si(2) += z*sig0;
    si(3) += rootAlpha*sig1;
    si(4) += rootAlpha*sig2;
    si(5) += -z*sig1 + y*sig2;
  }

  if (alpha != 1.0) {

  }

  return res;
}

const Vector&
NDFiberSection3d::getSectionDeformation(void)
{
  return e;
}

const Matrix&
NDFiberSection3d::getInitialTangent(void)
{
  static double kInitial[36];
  static Matrix ki(kInitial, 6, 6);
  ki.Zero();

  static double yLocs[10000];
  static double zLocs[10000];
  static double fiberArea[10000];
  // Shear factors arrays
  static double yScale[10000];
  static double zScale[10000];

  if (sectionIntegr != 0) {
    sectionIntegr->getFiberLocations(numFibers, yLocs, zLocs);
    sectionIntegr->getFiberWeights(numFibers, fiberArea);
  }  
  else {
    for (int i = 0; i < numFibers; i++) {
      yLocs[i] = matData[5*i];
      zLocs[i] = matData[5*i+1];
      fiberArea[i] = matData[5*i+2];
      yScale[i] = gammaData[2 * i];
      zScale[i] = gammaData[2 * i + 1];
    }
  }

  double rootAlpha = 1.0;
  if (alpha != 1.0)
    rootAlpha = sqrt(alpha);

  for (int i = 0; i < numFibers; i++) {
    NDMaterial *theMat = theMaterials[i];
    double y = yLocs[i] - yBar;
    double z = zLocs[i] - zBar;
    double A = fiberArea[i];
    double sgy = yScale[i];
    double sgz = zScale[i];

    double y2 = y*y;
    double z2 = z*z;
    double yz = y*z;
    double tmp;

    const Matrix &tangent = theMat->getInitialTangent();

    double d00 = tangent(0,0)*A;
    double d01 = tangent(0,1)*A;
    double d02 = tangent(0,2)*A;
    double d10 = tangent(1,0)*A;
    double d11 = tangent(1,1)*A;
    double d12 = tangent(1,2)*A;
    double d20 = tangent(2,0)*A;
    double d21 = tangent(2,1)*A;
    double d22 = tangent(2,2)*A;

    // Bending terms
    ki(0,0) += d00;
    ki(1,1) += y2*d00;
    ki(2,2) += z2*d00;
    tmp = -y*d00;
    ki(0,1) += tmp;
    ki(1,0) += tmp;
    tmp = z*d00;
    ki(0,2) += tmp;
    ki(2,0) += tmp;
    tmp = -yz*d00;
    ki(1,2) += tmp;
    ki(2,1) += tmp;
    
    // Shear terms
    ki(3,3) += alpha*d11;
    ki(3,4) += alpha*d12;
    ki(4,3) += alpha*d21;
    ki(4,4) += alpha*d22;
    
    // Torsion term
    ki(5,5) += z2*d11 - yz*(d12+d21) + y2*d22;
    
    // Bending-torsion coupling terms
    tmp = -z*d01 + y*d02;
    ki(0,5) += tmp;
    ki(1,5) -= y*tmp;
    ki(2,5) += z*tmp;
    tmp = -z*d10 + y*d20;
    ki(5,0) += tmp;
    ki(5,1) -= y*tmp;
    ki(5,2) += z*tmp;
    
    // Hit tangent terms with rootAlpha
    d01 *= rootAlpha; d02 *= rootAlpha;
    d10 *= rootAlpha; d11 *= rootAlpha; d12 *= rootAlpha;
    d20 *= rootAlpha; d21 *= rootAlpha; d22 *= rootAlpha;
    
    // Bending-shear coupling terms
    ki(0,3) += d01;
    ki(0,4) += d02;
    ki(1,3) -= y*d01;
    ki(1,4) -= y*d02;
    ki(2,3) += z*d01;
    ki(2,4) += z*d02;
    ki(3,0) += d10;
    ki(4,0) += d20;
    ki(3,1) -= y*d10;
    ki(4,1) -= y*d20;
    ki(3,2) += z*d10;
    ki(4,2) += z*d20;
    
    // Torsion-shear coupling terms
    y2 =  y*d22;
    z2 = -z*d11;
    ki(5,3) +=  z2 + y*d21;
    ki(5,4) += -z*d12 + y2;
    ki(3,5) +=  z2 + y*d12;
    ki(4,5) += -z*d21 + y2;
  }

  if (alpha != 1.0) {

  }

  return ki;
}

const Matrix&
NDFiberSection3d::getSectionTangent(void)
{
  return *ks;
}

const Vector&
NDFiberSection3d::getStressResultant(void)
{
  return *s;
}

SectionForceDeformation*
NDFiberSection3d::getCopy(void)
{
  NDFiberSection3d *theCopy = new NDFiberSection3d ();
  theCopy->setTag(this->getTag());

  theCopy->numFibers = numFibers;
  theCopy->sizeFibers = numFibers;

  if (numFibers != 0) {
    theCopy->theMaterials = new NDMaterial *[numFibers];

    if (theCopy->theMaterials == 0) {
      opserr <<"NDFiberSection3d::getCopy -- failed to allocate Material pointers\n";
      exit(-1);
    }
  
    theCopy->matData = new double [numFibers*5];
    if (theCopy->matData == 0) {
      opserr << "NDFiberSection3d::getCopy -- failed to allocate double array for material data\n";
      exit(-1);
    }
    // Gamma
    theCopy->gammaData = new double[numFibers * 2];
    if (theCopy->gammaData == 0) {
        opserr << "NDFiberSection3d::getCopy -- failed to allocate double array for shear scale factor\n";
        exit(-1);
    }
    // Creep
    theCopy->eps0_creep = new double[numFibers];
    theCopy->eps0_creep_k = new double[numFibers];
    theCopy->deps0_creep = new double[numFibers];
    theCopy->t0 = new double[numFibers];
    if (theCopy->eps0_creep == 0 || theCopy->eps0_creep_k == 0 || theCopy->deps0_creep == 0) {
        opserr << "NDFiberSection3d::getCopy -- failed to allocate double array for material data\n";
        exit(-1);
    }

    for (int i = 0; i < numFibers; i++) {
      theCopy->matData[5*i] = matData[5*i];
      theCopy->matData[5*i+1] = matData[5*i+1];
      theCopy->matData[5*i+2] = matData[5*i+2];
      theCopy->matData[5*i+3] = matData[5*i+3];
      theCopy->matData[5*i+4] = matData[5*i+4];
      theCopy->theMaterials[i] = theMaterials[i]->getCopy("BeamFiber");
      theCopy->gammaData[2*i] = gammaData[2*i];
      theCopy->gammaData[2*i+1] = gammaData[2*i+1];
      theCopy->eps0_creep[i] = eps0_creep[i]; // Creep strains at current step
      theCopy->eps0_creep_k[i] = eps0_creep_k[i]; // Creep strains at previous step
      theCopy->deps0_creep[i] = deps0_creep[i]; // Delta creep strains
      theCopy->t0[i] = t0[i]; // Creep time t0

      if (theCopy->theMaterials[i] == 0) {
	opserr <<"NDFiberSection3d::getCopy -- failed to get copy of a Material";
	exit(-1);
      }
    }  
  }

  theCopy->e = e;
  theCopy->QzBar = QzBar;
  theCopy->QyBar = QyBar;
  theCopy->Abar = Abar;
  theCopy->yBar = yBar;
  theCopy->zBar = zBar;
  theCopy->computeCentroid = computeCentroid;
  theCopy->alpha = alpha;
  theCopy->parameterID = parameterID;

  for (int i = 0; i < 6; i++)
    theCopy->sData[i] = sData[i];

  for (int i = 0; i < 6*6; i++)
    theCopy->kData[i] = kData[i];

  if (sectionIntegr != 0)
    theCopy->sectionIntegr = sectionIntegr->getCopy();
  else
    theCopy->sectionIntegr = 0;

  return theCopy;
}

const ID&
NDFiberSection3d::getType ()
{
  return code;
}

int
NDFiberSection3d::getOrder () const
{
  return 6;
}

int
NDFiberSection3d::commitState(void)
{
  int err = 0;
  for (int i = 0; i < numFibers; i++) {
      err += theMaterials[i]->commitState();
      
      // Load reset
      t0[i] = 1;
  }

  // Time step
  step += ops_Dt;

  // Damage output
  if (d_out == 1) {
      
      //opserr << "Export damage at step " << step*100 << "%" << endln;
      for (int i = 0; i < numFibers; i++) {
          double y = matData[5 * i];
          double z = matData[5 * i + 1];
          const Vector& dam = theMaterials[i]->getDamage();
          const Vector& eps = theMaterials[i]->getStrain();
          const Vector& sig = theMaterials[i]->getStress();

          // Damage output
          using namespace std;
          ofstream outdata;
          outdata.open("out_fibers.txt", ios::app);
          outdata << step << " " << y << " " << z << " ";
          outdata << eps(0) << " " << eps(1) << " " << eps(2) << " ";
          outdata << sig(0) << " " << sig(1) << " " << sig(2) << " ";
          outdata << dam(0) << " " << dam(1) << " " << dam(2) << endln;
          outdata.close();

          //opserr << "At fiber " << i + 1 << " eps_creep = " << deps0_creep[i] << endln;
      }

      // Strains and stresses output
      const Vector& strain = this->getSectionDeformation();
      const Vector& stress = this->getStressResultant();
      double energy = this->getEnergy();
      using namespace std;
      ofstream outdata;
      outdata.open("out_sections.txt", ios::app);
      outdata << step << " " << strain(0) << " " << strain(2) << " " << strain(4) << " " << stress(0) << " " << stress(2) << " " << stress(4) << " " << energy << endln;
      outdata.close();
      
  }

  return err;
}

int
NDFiberSection3d::revertToLastCommit(void)
{
  int err = 0;

  ks->Zero();
  s->Zero();
  
  static double yLocs[10000];
  static double zLocs[10000];
  static double fiberArea[10000];
  static double iStrains[10000];
  static double betaAngle[10000];

  if (sectionIntegr != 0) {
    sectionIntegr->getFiberLocations(numFibers, yLocs, zLocs);
    sectionIntegr->getFiberWeights(numFibers, fiberArea);
  }  
  else {
    for (int i = 0; i < numFibers; i++) {
      yLocs[i] = matData[5*i];
      zLocs[i] = matData[5*i+1];
      fiberArea[i] = matData[5*i+2];
      iStrains[i] = matData[5*i+3];
      betaAngle[i] = matData[5*i+4];
    }
  }

  double rootAlpha = 1.0;
  if (alpha != 1.0)
    rootAlpha = sqrt(alpha);

  for (int i = 0; i < numFibers; i++) {
    NDMaterial *theMat = theMaterials[i];
    double y = yLocs[i] - yBar;
    double z = zLocs[i] - zBar;
    double A = fiberArea[i];

    double y2 = y*y;
    double z2 = z*z;
    double yz = y*z;
    double tmp;

    // invoke revertToLast on the material
    err += theMat->revertToLastCommit();

    // get material stress & tangent for this strain and determine ks and fs
    const Matrix &tangent = theMat->getTangent();
    const Vector &stress = theMat->getStress();

    double d00 = tangent(0,0)*A;
    double d01 = tangent(0,1)*A;
    double d02 = tangent(0,2)*A;
    double d10 = tangent(1,0)*A;
    double d11 = tangent(1,1)*A;
    double d12 = tangent(1,2)*A;
    double d20 = tangent(2,0)*A;
    double d21 = tangent(2,1)*A;
    double d22 = tangent(2,2)*A;

    Matrix &ksi = *ks;
    Vector &si = *s;

    // Bending terms
    ksi(0,0) += d00;
    ksi(1,1) += y2*d00;
    ksi(2,2) += z2*d00;
    tmp = -y*d00;
    ksi(0,1) += tmp;
    ksi(1,0) += tmp;
    tmp = z*d00;
    ksi(0,2) += tmp;
    ksi(2,0) += tmp;
    tmp = -yz*d00;
    ksi(1,2) += tmp;
    ksi(2,1) += tmp;
    
    // Shear terms
    ksi(3,3) += alpha*d11;
    ksi(3,4) += alpha*d12;
    ksi(4,3) += alpha*d21;
    ksi(4,4) += alpha*d22;
    
    // Torsion term
    ksi(5,5) += z2*d11 - yz*(d12+d21) + y2*d22;
    
    // Bending-torsion coupling terms
    tmp = -z*d01 + y*d02;
    ksi(0,5) += tmp;
    ksi(1,5) -= y*tmp;
    ksi(2,5) += z*tmp;
    tmp = -z*d10 + y*d20;
    ksi(5,0) += tmp;
    ksi(5,1) -= y*tmp;
    ksi(5,2) += z*tmp;
    
    // Hit tangent terms with rootAlpha
    d01 *= rootAlpha; d02 *= rootAlpha;
    d10 *= rootAlpha; d11 *= rootAlpha; d12 *= rootAlpha;
    d20 *= rootAlpha; d21 *= rootAlpha; d22 *= rootAlpha;
    
    // Bending-shear coupling terms
    ksi(0,3) += d01;
    ksi(0,4) += d02;
    ksi(1,3) -= y*d01;
    ksi(1,4) -= y*d02;
    ksi(2,3) += z*d01;
    ksi(2,4) += z*d02;
    ksi(3,0) += d10;
    ksi(4,0) += d20;
    ksi(3,1) -= y*d10;
    ksi(4,1) -= y*d20;
    ksi(3,2) += z*d10;
    ksi(4,2) += z*d20;
    
    // Torsion-shear coupling terms
    y2 =  y*d22;
    z2 = -z*d11;
    ksi(5,3) +=  z2 + y*d21;
    ksi(5,4) += -z*d12 + y2;
    ksi(3,5) +=  z2 + y*d12;
    ksi(4,5) += -z*d21 + y2;

    double sig0 = stress(0)*A;
    double sig1 = stress(1)*A;
    double sig2 = stress(2)*A;

    si(0) += sig0;
    si(1) += -y*sig0;
    si(2) += z*sig0;
    si(3) += rootAlpha*sig1;
    si(4) += rootAlpha*sig2;
    si(5) += -z*sig1 + y*sig2;
  }

  if (alpha != 1.0) {

  }

  return err;
}

int
NDFiberSection3d::revertToStart(void)
{
  // revert the fibers to start    
  int err = 0;

  ks->Zero();
  s->Zero();
  
  static double yLocs[10000];
  static double zLocs[10000];
  static double fiberArea[10000];

  if (sectionIntegr != 0) {
    sectionIntegr->getFiberLocations(numFibers, yLocs, zLocs);
    sectionIntegr->getFiberWeights(numFibers, fiberArea);
  }  
  else {
    for (int i = 0; i < numFibers; i++) {
      yLocs[i] = matData[5*i];
      zLocs[i] = matData[5*i+1];
      fiberArea[i] = matData[5*i+2];
    }
  }

  double rootAlpha = 1.0;
  if (alpha != 1.0)
    rootAlpha = sqrt(alpha);

  for (int i = 0; i < numFibers; i++) {
    NDMaterial *theMat = theMaterials[i];

    double y = yLocs[i] - yBar;
    double z = zLocs[i] - zBar;
    double A = fiberArea[i];

    double y2 = y*y;
    double z2 = z*z;
    double yz = y*z;
    double tmp;

    // invoke revertToLast on the material
    err += theMat->revertToStart();

    // get material stress & tangent for this strain and determine ks and fs
    const Matrix &tangent = theMat->getTangent();
    const Vector &stress = theMat->getStress();

    double d00 = tangent(0,0)*A;
    double d01 = tangent(0,1)*A;
    double d02 = tangent(0,2)*A;
    double d10 = tangent(1,0)*A;
    double d11 = tangent(1,1)*A;
    double d12 = tangent(1,2)*A;
    double d20 = tangent(2,0)*A;
    double d21 = tangent(2,1)*A;
    double d22 = tangent(2,2)*A;

    Matrix &ksi = *ks;
    Vector &si = *s;

    // Bending terms
    ksi(0,0) += d00;
    ksi(1,1) += y2*d00;
    ksi(2,2) += z2*d00;
    tmp = -y*d00;
    ksi(0,1) += tmp;
    ksi(1,0) += tmp;
    tmp = z*d00;
    ksi(0,2) += tmp;
    ksi(2,0) += tmp;
    tmp = -yz*d00;
    ksi(1,2) += tmp;
    ksi(2,1) += tmp;
    
    // Shear terms
    ksi(3,3) += alpha*d11;
    ksi(3,4) += alpha*d12;
    ksi(4,3) += alpha*d21;
    ksi(4,4) += alpha*d22;
    
    // Torsion term
    ksi(5,5) += z2*d11 - yz*(d12+d21) + y2*d22;
    
    // Bending-torsion coupling terms
    tmp = -z*d01 + y*d02;
    ksi(0,5) += tmp;
    ksi(1,5) -= y*tmp;
    ksi(2,5) += z*tmp;
    tmp = -z*d10 + y*d20;
    ksi(5,0) += tmp;
    ksi(5,1) -= y*tmp;
    ksi(5,2) += z*tmp;
    
    // Hit tangent terms with rootAlpha
    d01 *= rootAlpha; d02 *= rootAlpha;
    d10 *= rootAlpha; d11 *= rootAlpha; d12 *= rootAlpha;
    d20 *= rootAlpha; d21 *= rootAlpha; d22 *= rootAlpha;
    
    // Bending-shear coupling terms
    ksi(0,3) += d01;
    ksi(0,4) += d02;
    ksi(1,3) -= y*d01;
    ksi(1,4) -= y*d02;
    ksi(2,3) += z*d01;
    ksi(2,4) += z*d02;
    ksi(3,0) += d10;
    ksi(4,0) += d20;
    ksi(3,1) -= y*d10;
    ksi(4,1) -= y*d20;
    ksi(3,2) += z*d10;
    ksi(4,2) += z*d20;
    
    // Torsion-shear coupling terms
    y2 =  y*d22;
    z2 = -z*d11;
    ksi(5,3) +=  z2 + y*d21;
    ksi(5,4) += -z*d12 + y2;
    ksi(3,5) +=  z2 + y*d12;
    ksi(4,5) += -z*d21 + y2;

    double sig0 = stress(0)*A;
    double sig1 = stress(1)*A;
    double sig2 = stress(2)*A;

    si(0) += sig0;
    si(1) += -y*sig0;
    si(2) += z*sig0;
    si(3) += rootAlpha*sig1;
    si(4) += rootAlpha*sig2;
    si(5) += -z*sig1 + y*sig2;
  }

  if (alpha != 1.0) {

  }

  return err;
}

int
NDFiberSection3d::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;

  // create an id to send objects tag and numFibers, 
  //     size 3 so no conflict with matData below if just 1 fiber
  static ID data(3);
  data(0) = this->getTag();
  data(1) = numFibers;
  data(2) = computeCentroid ? 1 : 0; // Now the ID data is really 3    
  int dbTag = this->getDbTag();
  res += theChannel.sendID(dbTag, commitTag, data);
  if (res < 0) {
    opserr <<  "NDFiberSection3d::sendSelf - failed to send ID data\n";
    return res;
  }    

  if (numFibers != 0) {
    
    // create an id containingg classTag and dbTag for each material & send it
    ID materialData(2*numFibers);
    for (int i=0; i<numFibers; i++) {
      NDMaterial *theMat = theMaterials[i];
      materialData(2*i) = theMat->getClassTag();
      int matDbTag = theMat->getDbTag();
      if (matDbTag == 0) {
	matDbTag = theChannel.getDbTag();
	if (matDbTag != 0)
	  theMat->setDbTag(matDbTag);
      }
      materialData(2*i+1) = matDbTag;
    }    
    
    res += theChannel.sendID(dbTag, commitTag, materialData);
    if (res < 0) {
      opserr <<  "NDFiberSection3d::sendSelf - failed to send material data\n";
      return res;
    }    

    // send the fiber data, i.e. area and loc
    Vector fiberData(matData, numFibers*5);
    res += theChannel.sendVector(dbTag, commitTag, fiberData);
    if (res < 0) {
      opserr <<  "NDFiberSection3d::sendSelf - failed to send material data\n";
      return res;
    }    

    // now invoke send(0 on all the materials
    for (int j=0; j<numFibers; j++)
      theMaterials[j]->sendSelf(commitTag, theChannel);

  }

  return res;
}

int
NDFiberSection3d::recvSelf(int commitTag, Channel &theChannel,
			 FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static ID data(3);
  
  int dbTag = this->getDbTag();
  res += theChannel.recvID(dbTag, commitTag, data);
  if (res < 0) {
    opserr <<  "NDFiberSection3d::recvSelf - failed to recv ID data\n";
    return res;
  }    
  this->setTag(data(0));

  // recv data about materials objects, classTag and dbTag
  if (data(1) != 0) {
    ID materialData(2*data(1));
    res += theChannel.recvID(dbTag, commitTag, materialData);
    if (res < 0) {
      opserr <<  "NDFiberSection3d::recvSelf - failed to recv material data\n";
      return res;
    }    

    // if current arrays not of correct size, release old and resize
    if (theMaterials == 0 || numFibers != data(1)) {
      // delete old stuff if outa date
      if (theMaterials != 0) {
	for (int i=0; i<numFibers; i++)
	  delete theMaterials[i];
	delete [] theMaterials;
	if (matData != 0)
	  delete [] matData;
	matData = 0;
	theMaterials = 0;
      }

      // create memory to hold material pointers and fiber data
      numFibers = data(1);
      sizeFibers = data(1);
      if (numFibers != 0) {
	theMaterials = new NDMaterial *[numFibers];
	
	if (theMaterials == 0) {
	  opserr <<"NDFiberSection3d::recvSelf -- failed to allocate Material pointers\n";
	  exit(-1);
	}
	
	for (int j=0; j<numFibers; j++)
	  theMaterials[j] = 0;

	matData = new double [numFibers*5];

	if (matData == 0) {
	  opserr <<"NDFiberSection3d::recvSelf  -- failed to allocate double array for material data\n";
	  exit(-1);
	}
      }
    }

    Vector fiberData(matData, numFibers*5);
    res += theChannel.recvVector(dbTag, commitTag, fiberData);
    if (res < 0) {
      opserr <<  "NDFiberSection3d::recvSelf - failed to recv material data\n";
      return res;
    }    

    int i;
    for (i=0; i<numFibers; i++) {
      int classTag = materialData(2*i);
      int dbTag = materialData(2*i+1);

      // if material pointed to is blank or not of corrcet type, 
      // release old and create a new one
      if (theMaterials[i] == 0)
	theMaterials[i] = theBroker.getNewNDMaterial(classTag);
      else if (theMaterials[i]->getClassTag() != classTag) {
	delete theMaterials[i];
	theMaterials[i] = theBroker.getNewNDMaterial(classTag);      
      }

      if (theMaterials[i] == 0) {
	opserr <<"NDFiberSection3d::recvSelf -- failed to allocate double array for material data\n";
	exit(-1);
      }

      theMaterials[i]->setDbTag(dbTag);
      res += theMaterials[i]->recvSelf(commitTag, theChannel, theBroker);
    }

    QzBar = 0.0;
    QyBar = 0.0;
    Abar  = 0.0;
    double yLoc, zLoc, Area;

    computeCentroid = data(2) ? true : false;
    
    // Recompute centroid
    for (i = 0; computeCentroid && i < numFibers; i++) {
      yLoc = matData[5*i];
      zLoc = matData[5*i+1];
      Area = matData[5*i+2];
      Abar  += Area;
      QzBar += yLoc*Area;
      QyBar += zLoc*Area;
    }

    if (computeCentroid) {
      yBar = QzBar/Abar;
      zBar = QyBar/Abar;
    } else {
      yBar = 0.0;
      zBar = 0.0;      
    }
  }    

  return res;
}

void
NDFiberSection3d::Print(OPS_Stream &s, int flag)
{
  s << "\nNDFiberSection3d, tag: " << this->getTag() << endln;
  s << "\tSection code: " << code;
  s << "\tNumber of Fibers: " << numFibers << endln;
  s << "\tCentroid (y,z): " << yBar << ' ' << zBar << endln;
  s << "\tShape factor, alpha = " << alpha << endln;

  if (flag == 1) {
    for (int i = 0; i < numFibers; i++) {
      s << "\nLocation (y,z) = " << matData[5*i] << ' ' << matData[5*i+1];
      s << "\nArea = " << matData[5*i+2] << endln;
      s << "\neps0 = " << matData[5*i+3] << endln;
      s << "\nbeta = " << matData[5*i+4] << endln;
      theMaterials[i]->Print(s, flag);
    }
  }
}

Response*
NDFiberSection3d::setResponse(const char **argv, int argc,
			      OPS_Stream &output)
{
  Response *theResponse =0;

  if (argc > 2 && strcmp(argv[0],"fiber") == 0) {

    static double yLocs[10000];
    static double zLocs[10000];
    
    if (sectionIntegr != 0) {
      sectionIntegr->getFiberLocations(numFibers, yLocs, zLocs);
    }  
    else {
      for (int i = 0; i < numFibers; i++) {
	yLocs[i] = matData[5*i];
	zLocs[i] = matData[5*i+1];
      }
    }
    
    int key = numFibers;
    int passarg = 2;
    
    if (argc <= 3) {		  // fiber number was input directly
      
      key = atoi(argv[1]);
      
    } else if (argc > 5) {  // find fiber closest to coord. with mat tag
      
      int matTag = atoi(argv[3]);
      double yCoord = atof(argv[1]);
      double zCoord = atof(argv[2]);
      double closestDist = 0;
      double ySearch, zSearch, dy, dz;
      double distance;
      int j;
      // Find first fiber with specified material tag
      for (j = 0; j < numFibers; j++) {
	if (matTag == theMaterials[j]->getTag()) {
	  //ySearch = matData[4*j];
	  //zSearch = matData[4*j+1];
	  ySearch = yLocs[j];
	  zSearch = zLocs[j];	    	  
	  dy = ySearch-yCoord;
	  dz = zSearch-zCoord;
	  closestDist = dy*dy + dz*dz;
	  key = j;
	  break;
	}
      }
      // Search the remaining fibers
      for ( ; j < numFibers; j++) {
	if (matTag == theMaterials[j]->getTag()) {
	  //ySearch = matData[4*j];
	  //zSearch = matData[4*j+1];
	  ySearch = yLocs[j];
	  zSearch = zLocs[j];	    	  	  
	  dy = ySearch-yCoord;
	  dz = zSearch-zCoord;
	  distance = dy*dy + dz*dz;
	  if (distance < closestDist) {
	    closestDist = distance;
	    key = j;
	  }
	}
      }
      passarg = 4;
    }
    
    else {                  // fiber near-to coordinate specified
      
      double yCoord = atof(argv[1]);
      double zCoord = atof(argv[2]);
      double closestDist;
      double ySearch, zSearch, dy, dz;
      double distance;
      
      //ySearch = matData[0];
      //zSearch = matData[1];
      ySearch = yLocs[0];
      zSearch = zLocs[0];	    	  	        
      dy = ySearch-yCoord;
      dz = zSearch-zCoord;
      closestDist = dy*dy + dz*dz;
      key = 0;
      for (int j = 1; j < numFibers; j++) {
	//ySearch = matData[4*j];
	//zSearch = matData[4*j+1];
	ySearch = yLocs[j];
	zSearch = zLocs[j];	    	
	dy = ySearch-yCoord;
	dz = zSearch-zCoord;
	distance = dy*dy + dz*dz;
	if (distance < closestDist) {
	  closestDist = distance;
	  key = j;
	}
      }
      passarg = 3;
    }
    
    if (key < numFibers && key >= 0) {
      output.tag("FiberOutput");
      output.attr("yLoc",matData[5*key]);
      output.attr("zLoc",matData[5*key+1]);
      output.attr("area",matData[5*key+2]);
      output.attr("eps0",matData[5*key+3]);
      output.attr("beta",matData[5*key+4]);
      
      theResponse = theMaterials[key]->setResponse(&argv[passarg], argc-passarg, output);
      
      output.endTag();
    }

  }

  if (theResponse == 0)
      return SectionForceDeformation::setResponse(argv, argc, output);

  return theResponse;
}


int 
NDFiberSection3d::getResponse(int responseID, Information &sectInfo)
{
  // Just call the base class method ... don't need to define
  // this function, but keeping it here just for clarity
  return SectionForceDeformation::getResponse(responseID, sectInfo);
}


void
NDFiberSection3d::zeroLoad(void)
{

}

int
NDFiberSection3d::addLoad(int iFib, int jFib, double eps0, double phi_t_t0, int creep)
{
    // Save previous phi
    phi_k = fmax(phi_t_t0, phi_k);

    for (int i = iFib - 1; i < jFib; i++) {
        // Extract strains at previous step
        const Vector& eps = theMaterials[i]->getStrain();

        // Creep correction only
        if (creep == 1) {
            // Creep eps0
            if (eps(0) < 0) eps0_creep[i] = -eps(0);
            else eps0_creep[i] = 0.0;

        }
        // Save previous eps0
        if (phi_k == phi_t_t0) {
            eps0_creep_k[i] = eps0_creep[i];

            // Update eps0
            matData[5 * i + 3] = eps0 + eps0_creep[i] * phi_t_t0;
        }
        else
            matData[5 * i + 3] = eps0 + eps0_creep_k[i]*phi_k + (eps0_creep[i] - eps0_creep_k[i]) * phi_t_t0;

        /*if (i == 45) {
            opserr << "eps0c = " << eps0_creep[i] << " eps0ck = " << eps0_creep_k[i] << endln;
            opserr << "phi = " << phi_t_t0 << " phik = " << phi_k << endln;
        }*/
    }

    return 0;
}


// AddingSensitivity:BEGIN ////////////////////////////////////
int
NDFiberSection3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int result = -1;

  if (strstr(argv[0],"alpha") != 0)
    return param.addObject(1, this);

  // Check if the parameter belongs to the material (only option for now)
  if (strstr(argv[0],"material") != 0) {
    
    if (argc < 3)
      return 0;

    // Get the tag of the material
    int materialTag = atoi(argv[1]);
    
    // Loop over fibers to find the right material
    for (int i = 0; i < numFibers; i++)
      if (materialTag == theMaterials[i]->getTag()) {
	int ok = theMaterials[i]->setParameter(&argv[2], argc-2, param);
	if (ok != -1)
	  result = ok;
      }
    return result;
  }

  // Check if it belongs to the section integration
  else if (strstr(argv[0],"integration") != 0) {
    if (sectionIntegr != 0)
      return sectionIntegr->setParameter(&argv[1], argc-1, param);
    else
      return -1;
  }

  int ok = 0;
  
  for (int i = 0; i < numFibers; i++) {
    ok = theMaterials[i]->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }

  if (sectionIntegr != 0) {
    ok = sectionIntegr->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }

  return result;
}

int
NDFiberSection3d::updateParameter(int paramID, Information &info)
{
  switch(paramID) {
  case 1:
    alpha = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
NDFiberSection3d::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

const Vector &
NDFiberSection3d::getSectionDeformationSensitivity(int gradIndex)
{
  return dedh;
}

const Vector &
NDFiberSection3d::getStressResultantSensitivity(int gradIndex, bool conditional)
{
  static Vector ds(6);
  
  ds.Zero();
  
  double y, z, A;
  static Vector stress(3);
  static Vector dsigdh(3);
  static Vector sig_dAdh(3);
  static Matrix tangent(3,3);

  static double yLocs[10000];
  static double zLocs[10000];
  static double fiberArea[10000];

  if (sectionIntegr != 0) {
    sectionIntegr->getFiberLocations(numFibers, yLocs, zLocs);
    sectionIntegr->getFiberWeights(numFibers, fiberArea);
  }  
  else {
    for (int i = 0; i < numFibers; i++) {
      yLocs[i] = matData[5*i];
      zLocs[i] = matData[5*i+1];
      fiberArea[i] = matData[5*i+2];
    }
  }

  static double dydh[10000];
  static double dzdh[10000];
  static double areaDeriv[10000];

  if (sectionIntegr != 0) {
    sectionIntegr->getLocationsDeriv(numFibers, dydh, dzdh);  
    sectionIntegr->getWeightsDeriv(numFibers, areaDeriv);
  }
  else {
    for (int i = 0; i < numFibers; i++) {
      dydh[i] = 0.0;
      dzdh[i] = 0.0;
      areaDeriv[i] = 0.0;
    }
  }
  
  double rootAlpha = 1.0;
  if (alpha != 1.0)
    rootAlpha = sqrt(alpha);

  double drootAlphadh = 0.0;
  if (parameterID == 1)
    drootAlphadh = 0.5/rootAlpha;

  for (int i = 0; i < numFibers; i++) {
    y = yLocs[i] - yBar;
    z = zLocs[i] - zBar;
    A = fiberArea[i];
    
    dsigdh = theMaterials[i]->getStressSensitivity(gradIndex,true);

    ds(0) += dsigdh(0)*A;
    ds(1) += -y*dsigdh(0)*A;
    ds(2) +=  z*dsigdh(0)*A;
    ds(3) += rootAlpha*dsigdh(1)*A;
    ds(4) += rootAlpha*dsigdh(2)*A;
    ds(5) += (-z*dsigdh(1)+y*dsigdh(2))*A;

    if (areaDeriv[i] != 0.0 || dydh[i] != 0.0 ||  dzdh[i] != 0.0 || parameterID == 1)
      stress = theMaterials[i]->getStress();

    if (dydh[i] != 0.0 || dzdh[i] != 0.0 || parameterID == 1)
      tangent = theMaterials[i]->getTangent();

    if (areaDeriv[i] != 0.0) {
      sig_dAdh(0) = stress(0)*areaDeriv[i];
      sig_dAdh(1) = stress(1)*areaDeriv[i];
      sig_dAdh(2) = stress(2)*areaDeriv[i];
      
      ds(0) += sig_dAdh(0);
      ds(1) += -y*sig_dAdh(0);
      ds(2) +=  z*sig_dAdh(0);
      ds(3) += rootAlpha*sig_dAdh(1);
      ds(4) += rootAlpha*sig_dAdh(2);
      ds(5) += -z*sig_dAdh(1)+y*sig_dAdh(2);
    }

    if (dydh[i] != 0.0) {
      ds(1) += -dydh[i] * (stress(0)*A);
      ds(5) +=  dydh[i] * (stress(2)*A);
    }

    if (dzdh[i] != 0.0) {
      ds(2) +=  dzdh[i] * (stress(0)*A);
      ds(5) += -dzdh[i] * (stress(1)*A);
    }

    if (parameterID == 1) {
      ds(3) += drootAlphadh * (stress(1)*A);
      ds(4) += drootAlphadh * (stress(2)*A);
    }

    static Matrix as(3,6);
    as(0,0) = 1;
    as(0,1) = -y;
    as(0,2) = z;
    as(1,3) = rootAlpha;
    as(2,4) = rootAlpha;
    as(1,5) = -z;
    as(2,5) = y;
    
    static Matrix dasdh(3,6);
    dasdh(0,1) = -dydh[i];
    dasdh(0,2) = dzdh[i];
    dasdh(1,3) = drootAlphadh;
    dasdh(2,4) = drootAlphadh;
    dasdh(1,5) = -dzdh[i];
    dasdh(2,5) = dydh[i];
    
    static Matrix tmpMatrix(6,6);
    tmpMatrix.addMatrixTripleProduct(0.0, as, tangent, dasdh, 1.0);
    
    ds.addMatrixVector(1.0, tmpMatrix, e, A);
  }

  return ds;
}

const Matrix &
NDFiberSection3d::getInitialTangentSensitivity(int gradIndex)
{
  static Matrix dksdh(6,6);
  
  dksdh.Zero();
  /*
  double y, A, dydh, dAdh, tangent, dtangentdh;

  static double fiberLocs[10000];
  static double fiberArea[10000];

  if (sectionIntegr != 0) {
    sectionIntegr->getFiberLocations(numFibers, fiberLocs);
    sectionIntegr->getFiberWeights(numFibers, fiberArea);
  }  
  else {
    for (int i = 0; i < numFibers; i++) {
      fiberLocs[i] = matData[2*i];
      fiberArea[i] = matData[2*i+1];
    }
  }

  static double locsDeriv[10000];
  static double areaDeriv[10000];

  if (sectionIntegr != 0) {
    sectionIntegr->getLocationsDeriv(numFibers, locsDeriv);  
    sectionIntegr->getWeightsDeriv(numFibers, areaDeriv);
  }
  else {
    for (int i = 0; i < numFibers; i++) {
      locsDeriv[i] = 0.0;
      areaDeriv[i] = 0.0;
    }
  }
  
  for (int i = 0; i < numFibers; i++) {
    y = fiberLocs[i] - yBar;
    A = fiberArea[i];
    dydh = locsDeriv[i];
    dAdh = areaDeriv[i];
    
    tangent = theMaterials[i]->getInitialTangent();
    dtangentdh = theMaterials[i]->getInitialTangentSensitivity(gradIndex);

    dksdh(0,0) += dtangentdh*A + tangent*dAdh;

    dksdh(0,1) += -y*(dtangentdh*A+tangent*dAdh) - dydh*(tangent*A);

    dksdh(1,1) += 2*(y*dydh*tangent*A) + y*y*(dtangentdh*A+tangent*dAdh);
  }

  dksdh(1,0) = dksdh(0,1);
  */
  return dksdh;
}

int
NDFiberSection3d::commitSensitivity(const Vector& defSens,
				    int gradIndex, int numGrads)
{
  double d0 = defSens(0);
  double d1 = defSens(1);
  double d2 = defSens(2);
  double d3 = defSens(3);
  double d4 = defSens(4);
  double d5 = defSens(5);

  dedh = defSens;

  static double yLocs[10000];
  static double zLocs[10000];

  if (sectionIntegr != 0)
    sectionIntegr->getFiberLocations(numFibers, yLocs, zLocs);
  else {
    for (int i = 0; i < numFibers; i++) {
      yLocs[i] = matData[5*i];
      zLocs[i] = matData[5*i+1];
    }
  }

  static double dydh[10000];
  static double dzdh[10000];

  if (sectionIntegr != 0)
    sectionIntegr->getLocationsDeriv(numFibers, dydh, dzdh);  
  else {
    for (int i = 0; i < numFibers; i++) {
      dydh[i] = 0.0;
      dzdh[i] = 0.0;
    }
  }

  double y, z;

  static Vector depsdh(3);

  double rootAlpha = 1.0;
  if (alpha != 1.0)
    rootAlpha = sqrt(alpha);

  double drootAlphadh = 0.0;
  if (parameterID == 1)
    drootAlphadh = 0.5/rootAlpha;

  for (int i = 0; i < numFibers; i++) {
    NDMaterial *theMat = theMaterials[i];
    y = yLocs[i] - yBar;
    z = zLocs[i] - zBar;

    // determine material strain and set it
    depsdh(0) = d0 - y*d1 + z*d2 - dydh[i]*e(1) + dzdh[i]*e(2);
    depsdh(1) = rootAlpha*d3 - z*d5 + drootAlphadh*e(3) - dzdh[i]*e(5);
    depsdh(2) = rootAlpha*d4 + y*d5 + drootAlphadh*e(4) + dydh[i]*e(5);

    theMat->commitSensitivity(depsdh,gradIndex,numGrads);
  }

  return 0;
}

// AddingSensitivity:END ///////////////////////////////////

// L. Parente
double NDFiberSection3d::getEnergy() const
{
    static double fiberArea[10000];

    if (sectionIntegr != 0) {
        sectionIntegr->getFiberWeights(numFibers, fiberArea);
    }
    else {
        for (int i = 0; i < numFibers; i++) {
            fiberArea[i] = matData[5 * i + 2];
        }
    }
    double energy = 0;
    for (int i = 0; i < numFibers; i++)
    {
        double A = fiberArea[i];
        double dE;
        dE = A * theMaterials[i]->getEnergy();
        if (dE == 0) {
            const Vector& s = theMaterials[i]->getStress();
            const Vector& e = theMaterials[i]->getStrain();
            const Vector& sk = theMaterials[i]->getCommittedStress();
            const Vector& ek = theMaterials[i]->getCommittedStrain();
            for (int i = 0; i < 3; i++) dE = A * 0.5 * (s(i) + sk(i)) * (e(i) - ek(i));
        }
        energy += dE;
    }
    return energy;
}
