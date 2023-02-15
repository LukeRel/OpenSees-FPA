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
// $Date: 2012-05-26 22:03:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/Condensation1D.cpp,v $

// Condensation algorithm to obtain a three component fiber vector e11 g12 g13
// from a uniaxial material. t12 and t13 stresses are simply set to 0.
// From PlateRebar, by Luca Parente (feb 2023)

#include <Condensation1D.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MaterialResponse.h>   //Antonios Vytiniotis used for the recorder
#include <math.h>
#include <elementAPI.h>

//static vector and matrices
Vector  Condensation1D::stress(5);
Matrix  Condensation1D::tangent(5, 5);

void* OPS_Condensation1D()
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 3) {
        opserr << "WARNING insufficient arguments\n";
        opserr << "Want: nDMaterial Condensation tag? matTag? angle?" << endln;
        return 0;
    }

    int tag[2];
    numdata = 2;
    if (OPS_GetIntInput(&numdata, tag) < 0) {
        opserr << "WARNING invalid nDMaterial Condensation tag or matTag" << endln;
        return 0;
    }

    UniaxialMaterial* theMaterial = OPS_getUniaxialMaterial(tag[1]);
    if (theMaterial == 0) {
        opserr << "WARNING uniaxialmaterial does not exist\n";
        opserr << "UniaxialMaterial: " << tag[1];
        opserr << "\nCondensation nDMaterial: " << tag[0] << endln;
        return 0;
    }

    NDMaterial* mat = new Condensation1D(tag[0], *theMaterial);

    if (mat == 0) {
        opserr << "WARNING: failed to create Condensation material\n";
        return 0;
    }

    return mat;
}

//null constructor
Condensation1D::Condensation1D() :
    NDMaterial(0, ND_TAG_Condensation1D),
    strain(5), theMat(0)
{ }

//full constructor
Condensation1D::Condensation1D(int tag, UniaxialMaterial& uniMat) :
    NDMaterial(tag, ND_TAG_Condensation1D),
    strain(3), theMat(0)
{
    theMat = uniMat.getCopy();
}


//destructor
Condensation1D::~Condensation1D()
{
    if (theMat != 0) delete theMat;
}


//make a clone of this material
NDMaterial*
Condensation1D::getCopy()
{
    Condensation1D* clone;   //new instance of this class

    clone = new Condensation1D(this->getTag(), *theMat); //make the copy

    return clone;
}


//make a clone of this material
NDMaterial*
Condensation1D::getCopy(const char* type)
{
    if (strcmp(type, "BeamFiber") == 0)
        return this->getCopy();
    else
        return 0;
}


//send back order of strain in vector form
int
Condensation1D::getOrder() const
{
    return 5;
}


const char*
Condensation1D::getType() const
{
    return "PlateFiber";
}



//swap history variables
int
Condensation1D::commitState()
{
    return theMat->commitState();
}



//revert to last saved state
int
Condensation1D::revertToLastCommit()
{
    return theMat->revertToLastCommit();
}


//revert to start
int
Condensation1D::revertToStart()
{
    strain.Zero();
    return theMat->revertToStart();
}


//mass per unit volume
double
Condensation1D::getRho()
{
    return theMat->getRho();
}


//receive the strain
int
Condensation1D::setTrialStrain(const Vector& strainFromElement)
{
    strain(0) = strainFromElement(0);
    strain(1) = strainFromElement(1);
    strain(2) = strainFromElement(2);

    return theMat->setTrialStrain(strain(0));
}


//send back the strain
const Vector&
Condensation1D::getStrain()
{
    return strain;
}


//send back the stress 
const Vector&
Condensation1D::getStress()
{
    double sig = theMat->getStress();

    stress(0) = sig;
    return stress;
}


//send back the tangent 
const Matrix&
Condensation1D::getTangent()
{
    double tan = theMat->getTangent();

    tangent.Zero();
    tangent(0, 0) = tan;

    return tangent;
}

const Matrix&
Condensation1D::getInitialTangent()
{
    double tan = theMat->getInitialTangent();
    tangent.Zero();
    tangent(0, 0) = tan;

    return tangent;
}
/*
const Vector& Condensation1D::getDamage(void)
{
    Vector dam(3);
    dam(0) = 0.0;
    dam(1) = 0.0;
    dam(2) = 0.0;

    opserr << "D = [" << dam(0) << " " << dam(1) << " " << dam(2) << "]" << endln;

    return dam;
}
*/


//print out data
void
Condensation1D::Print(OPS_Stream& s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "Condensation Material tag: " << this->getTag() << endln;
        s << "using uniaxial material: " << endln;
        theMat->Print(s, flag);
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"Condensation1D\", ";
        s << "\"material\": \"" << theMat->getTag() << "\"}";
    }
}


int
Condensation1D::sendSelf(int commitTag, Channel& theChannel)
{
    int res = 0;

    int dataTag = this->getDbTag();

    int matDbTag;

    static ID idData(3);
    idData(0) = dataTag;
    idData(1) = theMat->getClassTag();
    matDbTag = theMat->getDbTag();
    if (matDbTag == 0) {
        matDbTag = theChannel.getDbTag();
        theMat->setDbTag(matDbTag);
    }
    idData(2) = matDbTag;

    res = theChannel.sendID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "Condensation1D::sendSelf() - failed to send data" << endln;
        return res;
    }

    static Vector vecData(1);
    vecData(0) = 0.0;

    res = theChannel.sendVector(dataTag, commitTag, vecData);
    if (res < 0) {
        opserr << "Condensation1D::sendSelf() - failed to send data" << endln;
        return res;
    }

    // now send the materials data
    res += theMat->sendSelf(commitTag, theChannel);
    if (res < 0)
        opserr << "Condensation1D::sendSelf() - failed to send material1" << endln;

    return res;
}

int
Condensation1D::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
    int res = 0;

    int dataTag = this->getDbTag();

    // recv an id containing the tag and associated materials class and db tags
    static ID idData(3);
    res = theChannel.recvID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "Condensation1D::sendSelf() - failed to receive id data" << endln;
        return res;
    }

    this->setTag(idData(0));
    int matClassTag = idData(1);
    if (theMat->getClassTag() != matClassTag) {
        if (theMat != 0) delete theMat;
        theMat = theBroker.getNewUniaxialMaterial(matClassTag);
        if (theMat == 0) {
            opserr << "Condensation1D::recvSelf() - failed to get a material of type: " << matClassTag << endln;
            return -1;
        }
    }
    theMat->setDbTag(idData(2));

    static Vector vecData(1);
    res = theChannel.recvVector(dataTag, commitTag, vecData);
    if (res < 0) {
        opserr << "Condensation1D::sendSelf() - failed to receive vector data" << endln;
        return res;
    }

    // now receive the materials data
    res = theMat->recvSelf(commitTag, theChannel, theBroker);
    if (res < 0)
        opserr << "Condensation1D::sendSelf() - failed to receive material1" << endln;

    return res;
}

