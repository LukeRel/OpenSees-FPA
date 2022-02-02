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
                                                                        
// $Revision: 1.6 $
// $Date: 2008-08-26 16:52:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/SectionalLoad.cpp,v $
                                                                        
                                                                        
// Written: fmk 11/95
//          modified 11/01 for new design

// Purpose: This file contains the methods for class SectionalLoad.

#include <SectionalLoad.h>
#include <Element.h>
#include <Domain.h>
#include <SectionForceDeformation.h>

SectionalLoad::SectionalLoad(int tag, int cTag, int theSecTag)
  :Load(tag, cTag), secTag(theSecTag), theSection(0)
{

}

SectionalLoad::SectionalLoad(int tag, int cTag)
  :Load(tag, cTag), secTag(0), theSection(0)
{

}


// provided for the FEM_Object broker; the tag and elementTag need
// to be supplied in recvSelf();
SectionalLoad::SectionalLoad(int cTag)
:Load(0, cTag), secTag(0), theSection(0)
{

}


SectionalLoad::~SectionalLoad()
{

}


void
SectionalLoad::setDomain(Domain* theDomain)
{
    /*
    this->DomainComponent::setDomain(theDomain);

    if (theDomain == 0) {
        theSection = 0;
        return;
    }
    */
    // Retrieve section from the model builder	
    SectionForceDeformation* theSection = OPS_getSectionForceDeformation(secTag);
    if (theSection == 0) {
        opserr << "WARNING - SectionalLoad::setDomain - no ele with tag ";
        opserr << secTag << " exists in the domain\n";
    }

}

void 
SectionalLoad::applyLoad(double loadFactor) 
{
  if (theSection != 0)
      theSection->addLoad(this, loadFactor);
}

void 
SectionalLoad::applyLoad(const Vector &loadFactors) 
{
  if (theSection != 0)
      theSection->addLoad(this, loadFactors);
}


const Vector&
SectionalLoad::getSensitivityData(int gradIndex)
{
  static Vector trash(10);

  return trash;
}

int
SectionalLoad::getSectionTag(void) 
{
   return secTag;
}



