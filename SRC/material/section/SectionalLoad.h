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
// $Date: 2010-02-04 19:16:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/SectionalLoad.h,v $
                                                                        
                                                                        
#ifndef SectionalLoad_h
#define SectionalLoad_h

// Written: L.Parente 
//
// Purpose: This file contains the class definition for SectionalLoad.
// SectionalLoad is an abstract class.

#include <Load.h>
#include <Vector.h>

class SectionForceDeformation;

class SectionalLoad : public Load
{
  public:
    SectionalLoad(int tag, int classTag, int _secTag);
    SectionalLoad(int tag, int classTag);
    SectionalLoad(int classTag);    
    ~SectionalLoad();

    virtual void setDomain(Domain *theDomain);
    virtual void applyLoad(double loadfactor);
    virtual void applyLoad(const Vector &loadfactors);
    virtual const Vector &getData(int &type, double loadFactor) = 0;
    virtual const Vector &getSensitivityData(int gradIndex);

    virtual int getSectionTag(void);

  protected:
    int secTag;
    SectionForceDeformation *theSection;  // pointer to associated section

  private:

};

#endif

