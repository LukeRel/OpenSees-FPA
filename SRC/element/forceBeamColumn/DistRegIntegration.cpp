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

// $Revision: 1.3 $
// $Date: 2007-10-26 04:50:14 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/DistRegIntegration.cpp,v $

#include <DistRegIntegration.h>
#include <ElementalLoad.h>

#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <math.h>

// Creazione della libreria dinamica
#ifdef _USRDLL
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif

// Funzione DistRegIntegration all'interno della classe DistRegIntegration -> nella classe BeamIntegration
DistRegIntegration::DistRegIntegration(int nIP, int nc, double Lc, BeamIntegration& bi):
  BeamIntegration(BEAM_INTEGRATION_TAG_DistReg),
  nIP(nIP), nc(nc), Lc(Lc), beamInt(0), parameterID(0)
{
  beamInt = bi.getCopy();
  if (beamInt == 0) {
    opserr << "DistRegIntegration::DistRegIntegration -- failed to get copy of BeamIntegration" << endln;
  }
}

DistRegIntegration::DistRegIntegration():
  BeamIntegration(BEAM_INTEGRATION_TAG_DistReg),
    nIP(0), nc(0), Lc(0.0), beamInt(0), parameterID(0)
{

}

DistRegIntegration::~DistRegIntegration()
{
  if (beamInt != 0)
    delete beamInt;
}

void
DistRegIntegration::getSectionLocations(int numSections, double L, double* xi)
{
    // Dati
    int nm = nIP - 2 * nc + 2; // Punti regione centrale
    double betas;
    if (nc == 2) betas = Lc / (1.0 / 2.0) / L; // 2 punti di Gauss
    else         betas = Lc / (1.0 / 6.0) / L; // 3 punti di Gauss
    double betam = 1.0 - 2.0 * betas; // Lunghezza regione centrale

    // Regione centrale
    beamInt->getSectionLocations(nm, betam, xi);
    for (int i = 0; i < nm; i++) {
        xi[nIP-nc-1-i] = xi[nm-1-i]*betam+betas;
    }

    // Regioni caratteristiche
    if (nc == 2) {
        xi[0] = 0.0;
        xi[nIP - 3] = 1.0;
        xi[nIP - 2] = betas;
        xi[nIP - 1] = 1.0 - betas;
    }
    else {
        xi[0] = 0.0;
        xi[nIP - 3] = 1.0;
        xi[nIP - 2] = betas/2;
        xi[nIP - 1] = 1 - betas/2;
    }

    // Verifica
    /*
    opserr <<  "xi: ";
    for (int i = 0; i < nIP; i++) {
        opserr << xi[i] << " ";
    }
    opserr << endln;
    */
}

void
DistRegIntegration::getSectionWeights(int numSections, double L, double *wt)
{
    // Dati
    int nm = nIP - 2 * nc + 2; // Punti regione centrale
    double betas;
    if (nc == 2) betas = Lc / (1.0 / 2.0) / L; // 2 punti di Gauss
    else         betas = Lc / (1.0 / 6.0) / L; // 3 punti di Gauss
    double betam = 1 - 2.0 * betas; // Lunghezza regione centrale

    // Regione centrale
    beamInt->getSectionWeights(nm, betam, wt);
    for (int i = 0; i < nm; i++) {
        wt[nIP - nc - 1 - i] = wt[nm - 1 - i] * betam;
    }

    // Regioni caratteristiche
    if (nc == 2) {
        wt[nIP - 2] = wt[0] + betas / 2;
        wt[nIP - 1] = wt[nIP - nc - 1] + betas / 2;
        wt[0] = betas / 2;
        wt[nIP - 3] = betas / 2;
    }
    else {
        wt[0] = betas / 6;
        wt[1] = wt[2] + wt[0];
        wt[nIP - nc - 2] = wt[nIP - nc - 2] + wt[0];
        wt[nIP - 3] = wt[0];
        wt[nIP - 2] = betas *2/3;
        wt[nIP - 1] = wt[nIP - 2];
    }

    // Verifica
    /*
    opserr << "wt: ";
    for (int i = 0; i < nIP; i++) {
        opserr << wt[i] << " ";
    }
    opserr << endln;
    */
}

BeamIntegration*
DistRegIntegration::getCopy(void)
{
  return new DistRegIntegration(nIP, nc, Lc, *beamInt);
}

int
DistRegIntegration::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(3);

  data(0) = nIP;
  data(1) = nc;
  data(2) = Lc;

  int dbTag = this->getDbTag();

  if (theChannel.sendVector(dbTag, cTag, data) < 0) {
    opserr << "DistRegIntegration::sendSelf() - failed to send Vector data\n";
    return -1;
  }    

  return 0;
}

int
DistRegIntegration::recvSelf(int cTag, Channel &theChannel,
				 FEM_ObjectBroker &theBroker)
{
  static Vector data(3);

  int dbTag = this->getDbTag();

  if (theChannel.recvVector(dbTag, cTag, data) < 0)  {
    opserr << "DistRegIntegration::recvSelf() - failed to receive Vector data\n";
    return -1;
  }
  
  nIP = data(0);
  nc = data(1);
  Lc = data(2);

  return 0;
}

int
DistRegIntegration::setParameter(const char **argv, int argc,
				   Parameter &param)
{
  if (argc < 1)
    return -1;

  if (strcmp(argv[0],"nIP") == 0) {
    param.setValue(nIP);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"nc") == 0) {
    param.setValue(nc);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"Lc") == 0) {
    param.setValue(Lc);
    return param.addObject(3, this);
  }
  return -1;
}

int
DistRegIntegration::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
  case 1:
    nIP = info.theInt;
    return 0;
  case 2:
    nc = info.theInt;
    return 0;
  case 3:
    Lc = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
DistRegIntegration::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

void
DistRegIntegration::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "{\"type\": \"DistReg\", ";
		s << "\"nIP\": " << nIP << ", ";
		s << "\"nc\": " << nc << ", ";
		s << "\"integration\": ";
		beamInt->Print(s, flag);
		s << "}";
	}
	
	else {
		s << "DistReg" << endln;
		s << " nIP = " << nIP;
		s << " nc = " << nc << endln;
		beamInt->Print(s, flag);
	}
}

void 
DistRegIntegration::getLocationsDeriv(int numSections, double L,
					double dLdh, double *dptsdh)
{
    opserr << "getLocationsDeriv was used";
    int numPerHinge = (numSections-2)/2;

  double oneOverL = 1/L;
  double betaI = Lc*oneOverL;
  double betaJ = Lc*oneOverL;

  beamInt->getSectionLocations(numPerHinge, L, dptsdh);

  //opserr << "DistRegIntegration::getLocationsDeriv -- implementation for interior not yet finished" << endln;

  if (parameterID == 1) { // lpI
    for (int i = 0; i < numPerHinge; i++) {
      dptsdh[i] = oneOverL*dptsdh[i];
      dptsdh[numSections-3-i] = 0.0;
    }
  }
  else if (parameterID == 2) { // lpJ
    for (int i = 0; i < numPerHinge; i++) {
      dptsdh[numSections-3-i] = -oneOverL*dptsdh[i];
      dptsdh[i] = 0.0;
    }
  }
  else if (dLdh != 0.0) {
    for (int i = 0; i < numPerHinge; i++) {
      dptsdh[numSections-3-i] = betaJ*oneOverL*dLdh*dptsdh[i];
      dptsdh[i] = -betaI*oneOverL*dLdh*dptsdh[i];
    }
  }
  else {
    for (int i = 0; i < numSections; i++)
      dptsdh[i] = 0.0;
  }

  return;
}

void
DistRegIntegration::getWeightsDeriv(int numSections, double L,
					double dLdh, double *dwtsdh)
{
    opserr << "getWeightsDeriv was used";
  int numPerHinge = (numSections-2)/2;

  double oneOverL = 1/L;
  double betaI = Lc*oneOverL;
  double betaJ = Lc*oneOverL;

  beamInt->getSectionWeights(numPerHinge, L, dwtsdh);

  //opserr << "DistRegIntegration::getWeightsDeriv -- implementation for interior not yet finished" << endln;

  if (parameterID == 1) { // lpI
    for (int i = 0; i < numPerHinge; i++) {
      dwtsdh[i] = oneOverL*dwtsdh[i];
      dwtsdh[numSections-3-i] = 0.0;
    }
  }
  else if (parameterID == 2) { // lpJ
    for (int i = 0; i < numPerHinge; i++) {
      dwtsdh[numSections-3-i] = oneOverL*dwtsdh[i];
      dwtsdh[i] = 0.0;
    }
  }
  else if (dLdh != 0.0) {
    for (int i = 0; i < numPerHinge; i++) {
      dwtsdh[numSections-3-i] = -betaJ*oneOverL*dLdh*dwtsdh[i];
      dwtsdh[i] = -betaI*oneOverL*dLdh*dwtsdh[i];
    }
  }
  else {
    for (int i = 0; i < numSections; i++)
      dwtsdh[i] = 0.0;
  }

  return;
}
