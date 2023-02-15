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

// $Revision: 1.0 $
// $Date: 2012-05-23 21:11:45 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/Condensation1D.h,v $

// Condensation algorithm to obtain a three component fiber vector e11 g12 g13
// from a uniaxial material. t12 and t13 stresses are simply set to 0.
// From PlateRebarMaterial, by Luca Parente (feb 2023)

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>
#include <ID.h> 
#include <UniaxialMaterial.h>
#include <NDMaterial.h>



class Condensation1D : public NDMaterial {
public:
    Condensation1D();
    Condensation1D(int tag, UniaxialMaterial& uniMat);

    virtual ~Condensation1D();

    //make a clone of this material
    NDMaterial* getCopy();
    NDMaterial* getCopy(const char* type);

    //send back order of strain in vector form
    int getOrder() const;

    //send back order of strain in vector form
    const char* getType() const;

    //swap history variables
    int commitState();

    //revert to last saved state
    int revertToLastCommit();

    //revert to start
    int revertToStart();

    //get the strain 
    int setTrialStrain(const Vector& strainFromElement);

    //send back the strain
    const Vector& getStrain();

    //send back the stress 
    const Vector& getStress();

    //send back the tangent 
    const Matrix& getTangent();
    const Matrix& getInitialTangent();  // AV Not Sure if it works
    //const Vector& getDamage();

    //density
    double getRho();

    //print out data
    void Print(OPS_Stream& s, int flag);

    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

private:
    UniaxialMaterial* theMat;

    Vector strain;
    static Vector stress;
    static Matrix tangent;

};





