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
                                                                        
// $Revision: 1.7 $
// $Date: 2007-10-17 22:11:35 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/FiberLoad.cpp,v $
                                                                        

// Written: fmk 
//
// Purpose: This file contains the class implementation of FiberLoad.

#include <FiberLoad.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>

Vector FiberLoad::data(2);

FiberLoad::FiberLoad(int tag, int _fiber, double _eps0, int theSectionTag)
  :SectionalLoad(tag, LOAD_TAG_FiberLoad, theSectionTag),
   fiber(_fiber), eps0(_eps0), parameterID(0)
{

}

FiberLoad::FiberLoad()
  :SectionalLoad(LOAD_TAG_FiberLoad),
   fiber(0.0), eps0(0.0), parameterID(0)
{

}

FiberLoad::~FiberLoad()
{

}

const Vector &
FiberLoad::getData(int &type, double loadFactor)
{
  type = LOAD_TAG_FiberLoad;
  data(0) = fiber;
  data(1) = eps0;
  return data;
}


int 
FiberLoad::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  static Vector vectData(3);
  vectData(0) = eps0;
  vectData(1) = fiber;
  vectData(2) = this->getTag();

  int result = theChannel.sendVector(dbTag, commitTag, vectData);
  if (result < 0) {
    opserr << "FiberLoad::sendSelf - failed to send data\n";
    return result;
  }

  return 0;
}

int 
FiberLoad::recvSelf(int commitTag, Channel &theChannel,  FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static Vector vectData(3);

  int result = theChannel.recvVector(dbTag, commitTag, vectData);
  if (result < 0) {
    opserr << "FiberLoad::recvSelf - failed to recv data\n";
    return result;
  }

  this->setTag(vectData(2));
  eps0 = vectData(0);;
  fiber = vectData(1);

  return 0;
}

void 
FiberLoad::Print(OPS_Stream &s, int flag)
{
  s << "FiberLoad - tag " << this->getTag() << endln;
  s << "  Fiber acted on: " << fiber << endln;
  s << "  eps0: " << eps0 << endln;
}

int
FiberLoad::setParameter(const char** argv, int argc, Parameter& param)
{
    if (argc < 1)
        return -1;

    if (strcmp(argv[0], "fiber") == 0) {
        param.setValue(fiber);
            return param.addObject(1, this);
    }
    if (strcmp(argv[0], "eps0") == 0) {
        param.setValue(eps0);
            return param.addObject(2, this);
    }

    return -1;
}

int
FiberLoad::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
  case 1:
    fiber = info.theInt;
    return 0;
  case 2:
    eps0 = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
FiberLoad::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

const Vector&
FiberLoad::getSensitivityData(int gradNumber)
{
  data.Zero();

  switch(parameterID) {
  case 1:
    data(0) = 1.0;
    break;
  case 2:
    data(1) = 1.0;
    break;
  default:
    break;
  }

  return data;
}
