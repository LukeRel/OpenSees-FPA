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

Vector Beam3dSectionLoad::data(4);

Beam3dSectionLoad::Beam3dSectionLoad(int tag, int _eleTag, int _iSec, int _jSec, int _iFib, int _jFib,
    double _eps0, double _fcm, double _RH, double _h, double _t0)
  :ElementalLoad(tag, LOAD_TAG_Beam3dSectionLoad, _eleTag),
   iSec(_iSec), jSec(_jSec), iFib(_iFib), jFib(_jFib), eps0(_eps0), fcm(_fcm), RH(_RH), h(_h), t0(_t0)
{

}

Beam3dSectionLoad::Beam3dSectionLoad()
  :ElementalLoad(LOAD_TAG_Beam3dSectionLoad),
    iSec(0), jSec(0), iFib(0), jFib(0), eps0(0.0), fcm(0.0), RH(0.0), h(0.0), t0(0.0)
{

}

Beam3dSectionLoad::~Beam3dSectionLoad()
{

}

const Vector &
Beam3dSectionLoad::getData(int &type, double loadFactor)
{
  type = LOAD_TAG_Beam3dSectionLoad;
  data(0) = iSec;
  data(1) = jSec;
  data(2) = iFib;
  data(3) = jFib;
  data(4) = eps0;
  data(5) = fcm;
  data(6) = RH;
  data(7) = h;
  data(8) = t0;

  //opserr << "s = "<< secTag <<" f = " << fibTag << " creep = " << creep << endln;

  return data;
}


int 
Beam3dSectionLoad::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  static Vector vectData(5);
  vectData(0) = eleTag;
  vectData(1) = iSec;
  vectData(2) = jSec;
  vectData(3) = iFib;
  vectData(4) = jFib;
  vectData(5) = eps0;
  vectData(6) = fcm;
  vectData(7) = RH;
  vectData(8) = h;
  vectData(9) = t0;
  vectData(10) = this->getTag();

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
  iSec = (int)vectData(1);
  jSec = (int)vectData(2);
  iFib = (int)vectData(3);
  jFib = (int)vectData(4);
  eps0   = vectData(5);
  fcm = vectData(6);
  RH = vectData(7);
  h = vectData(8);
  t0 = vectData(9);
  this->setTag(vectData(10));

  return 0;
}

void 
Beam3dSectionLoad::Print(OPS_Stream &s, int flag)
{
  s << "Beam3dSectionLoad - Reference load: " << this->getTag() << endln;
  s << "  Element  : " << eleTag << endln;
  s << "  Section i: " << iSec << endln;
  s << "  Section j: " << jSec << endln;
  s << "  Fiber i  : " << iFib << endln;
  s << "  Fiber j  : " << jFib << endln;
  s << "  eps0: "       << eps0     << endln;
  s << "  fcm: " << fcm << endln;
  s << "  RH: " << RH << endln;
  s << "  h: " << h << endln;
  s << "  t0: " << t0 << endln;
}
