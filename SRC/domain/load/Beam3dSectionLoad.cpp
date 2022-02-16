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
// $Date: 2007-10-17 22:11:35 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/Beam3dSectionLoad.cpp,v $
                                                                        

// Written: fmk 
//
// Purpose: This file contains the class implementation of Beam3dSectionLoad.

#include <Beam3dSectionLoad.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

Vector Beam3dSectionLoad::data(3);

Beam3dSectionLoad::Beam3dSectionLoad(int tag, int _eleTag, int _secTag, int _fibTag, double _eps0)
  :ElementalLoad(tag, LOAD_TAG_Beam3dSectionLoad, _eleTag),
   secTag(_secTag), fibTag(_fibTag), eps0(_eps0)
{

}

Beam3dSectionLoad::Beam3dSectionLoad()
  :ElementalLoad(LOAD_TAG_Beam3dSectionLoad),
    secTag(0), fibTag(0), eps0(0.0)
{

}

Beam3dSectionLoad::~Beam3dSectionLoad()
{

}

const Vector &
Beam3dSectionLoad::getData(int &type, double loadFactor)
{
  type = LOAD_TAG_Beam3dSectionLoad;
  data(0) = secTag;
  data(1) = fibTag;
  data(2) = eps0;
  return data;
}


int 
Beam3dSectionLoad::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  static Vector vectData(5);
  vectData(0) = eleTag;
  vectData(1) = secTag;
  vectData(2) = fibTag;
  vectData(3) = eps0;
  vectData(4) = this->getTag();

  int result = theChannel.sendVector(dbTag, commitTag, vectData);
  if (result < 0) {
    opserr << "Beam3dSectionLoad::sendSelf - failed to send data\n";
    return result;
  }

  return 0;
}

int 
Beam3dSectionLoad::recvSelf(int commitTag, Channel &theChannel,
			    FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static Vector vectData(5);

  int result = theChannel.recvVector(dbTag, commitTag, vectData);
  if (result < 0) {
    opserr << "Beam3dSectionLoad::recvSelf - failed to recv data\n";
    return result;
  }

  eleTag = (int)vectData(0);
  secTag = (int)vectData(1);
  fibTag = (int)vectData(2);
  eps0   = vectData(3);
  this->setTag(vectData(4));

  return 0;
}

void 
Beam3dSectionLoad::Print(OPS_Stream &s, int flag)
{
  s << "Beam3dSectionLoad - Reference load: " << this->getTag() << endln;
  s << "  eps0: "       << eps0     << endln;
  s << "  Element  : "  << eleTag   << endln;
  s << "  Section  : "  << secTag   << endln;
  s << "  Fiber    : "  << fibTag   << endln;
}
