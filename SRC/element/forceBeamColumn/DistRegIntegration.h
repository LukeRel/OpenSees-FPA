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

// $Revision: 1.2 $
// $Date: 2006-09-05 22:57:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/DistRegIntegration.h,v $

#ifndef DistRegIntegration_h
#define DistRegIntegration_h

#include <BeamIntegration.h>

class Matrix;
class ElementalLoad;
class Channel;
class FEM_ObjectBroker;

class DistRegIntegration : public BeamIntegration
{
 public:
  DistRegIntegration(int nIP, int nc, double Lc, BeamIntegration &bi);
  DistRegIntegration();
  ~DistRegIntegration();
  
  void getSectionLocations(int numSections, double L, double *xi);
  void getSectionWeights(int numSections, double L, double *wt);

  BeamIntegration *getCopy(void);

  int sendSelf(int cTag, Channel &theChannel);
  int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);
  int activateParameter(int parameterID);

  void Print(OPS_Stream &s, int flag = 0);

  void getLocationsDeriv(int nIP, double L, double dLdh, double *dptsdh);
  void getWeightsDeriv(int nIP, double L, double dLdh, double *dwtsdh);

 private:
	 // void LobattoPositions(int n, double a, double b, double* x);
	 void getLobattoWeights(int n, double L, double* wt);

	 // Punti di integrazione
	 int nIP; // Totali
	 int nc;  // Regione caratteristica
	 
	 // Lunghezza caratteristica
	 double Lc;	// Zona di localizzazione

  BeamIntegration *beamInt;

  int parameterID;
};

#endif
