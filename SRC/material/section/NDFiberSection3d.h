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
//
// Description: This file contains the class definition for 
// NDFiberSection3d.h. NDFiberSection3d provides the abstraction of a 
// 2d beam section discretized by fibers. The section stiffness and
// stress resultants are obtained by summing fiber contributions.

#ifndef NDFiberSection3d_h
#define NDFiberSection3d_h

#include <SectionForceDeformation.h>
#include <Vector.h>
#include <Matrix.h>
#include <Domain.h>
#include <DomainComponent.h>

class NDMaterial;
class Fiber;
class Response;
class SectionIntegration;

class NDFiberSection3d : public SectionForceDeformation
{
  public:
    NDFiberSection3d(); 
    NDFiberSection3d(int tag, int numFibers, Fiber **fibers, double a = 1.0, bool compCentroid=true);
    NDFiberSection3d(int tag, int numFibers, double a = 1.0, bool compCentroid=true);
    NDFiberSection3d(int tag, int numFibers, NDMaterial **mats,
		     SectionIntegration &si, double a = 1.0, bool compCentroid=true);
    ~NDFiberSection3d();

    const char *getClassType(void) const {return "NDFiberSection3d";};

    int   setTrialSectionDeformation(const Vector &deforms); 
    const Vector &getSectionDeformation(void);

    const Vector &getStressResultant(void);
    const Matrix &getSectionTangent(void);
    const Matrix &getInitialTangent(void);

    int   commitState(void);
    int   revertToLastCommit(void);    
    int   revertToStart(void);
 
    SectionForceDeformation *getCopy(void);
    const ID &getType (void);
    int getOrder (void) const;
    
    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);
	    
    Response *setResponse(const char **argv, int argc, 
			  OPS_Stream &s);
    int getResponse(int responseID, Information &info);

    // Adding fiber loads (L. Parente)
    void zeroLoad(void);
    int addLoad(int iFib, int jFib, double eps0, double phi_t_t0, int creep);

    int addFiber(Fiber &theFiber);

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);
    int activateParameter(int parameterID);
    const Vector& getStressResultantSensitivity(int gradIndex,
						bool conditional);
    const Vector& getSectionDeformationSensitivity(int gradIndex);
    const Matrix& getInitialTangentSensitivity(int gradIndex);
    int commitSensitivity(const Vector& sectionDeformationGradient,
			  int gradIndex, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////

    double getEnergy() const; // L. Parente

  protected:
    
    //  private:
    int numFibers, sizeFibers;       // number of fibers in the section
    NDMaterial **theMaterials;       // array of pointers to materials
    double   *matData;               // data for the materials [yloc and area]
    double   kData[36];              // data for ks matrix 
    double   sData[6];               // data for s vector

    // y and z gamma scale factors to obtain a parabolic distribution
    double   *gammaData;            // y and z scale factors

    double Abar,QyBar, QzBar;
    double yBar;       // Section centroid
    double zBar;       // Section centroid
    bool computeCentroid;
    double alpha;      // Shear shape factor

    SectionIntegration *sectionIntegr;

    static ID code;

    Vector e;          // trial section deformations 
    Vector *s;         // section resisting forces  (axial force, bending moment)
    Matrix *ks;        // section stiffness

// AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
    Vector dedh; // MHS hack
// AddingSensitivity:END ///////////////////////////////////////////

    // For damage output
    Domain* theDomain;
    double step;

    // For creep deformations
    double phi_k = 0.0;     // Previous phi if higher
    double *eps0_creep;     // eps0 saved at t0
    double *eps0_creep_k;   // eps0 saved at previous t0
    double *deps0_creep;    // Delta eps0

    double *t0;     // Save time at which creep is applied
};

#endif
