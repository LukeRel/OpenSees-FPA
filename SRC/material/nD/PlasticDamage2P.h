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
// $Date: 2006-08-04 18:18:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticIsotropicPlaneStress2D.h,v $

#ifndef PlasticDamage2P_h
#define PlasticDamage2P_h

// Written: Thanh Do
// Created: 07/16
//
// Description: 
//
// What: "@(#) ElasticIsotropicThreeDimesnional.h, revA"

#include <ElasticIsotropicMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class PlasticDamage2P : public NDMaterial
{
  public:

    // Full constructor
    PlasticDamage2P(int tag, double E, double nu, //(CALCOLA BULK E SHEAR)
          double ft, double fc, //(CALCOLA s_y e r)
          double Hk, double Hi, //(CALCOLA H e t)
          double r_bar, double Kinfinity, double Kinit, double d1, double d2, double mDen, //(PARAMETRI CHE SI POTRANNO IGNORARE)
          double Yt0, double bt, double at, double Yc0, double bc, double ac, double kappa); //(IN GATTA)

    // Null constructor
    PlasticDamage2P();

    // Destructor
    ~PlasticDamage2P();

    // Class type return
    const char *getClassType(void) const {return "PlasticDamage2P";};

    // Main methods
    int setTrialStrain (const Vector &v);
    int setTrialStrainIncr (const Vector &v);
    int setTrialStrainIncr (const Vector &v, const Vector &r);

    // Constitutive matrix, stress and strain determination
    const Matrix &getTangent (void);
    const Matrix &getInitialTangent (void);
    const Vector &getStress (void);
    const Vector &getStrain (void);
    
    // Commit updating
    int commitState (void);
    int revertToLastCommit (void);
    int revertToStart (void);

    // Additional functions
    NDMaterial*getCopy(const char *type);
    NDMaterial *getCopy (void);
    const char *getType (void) const;
    int getOrder (void) const;

    // Optional parallel workflow
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    

    // Output flag method
    void Print(OPS_Stream &s, int flag =0);   

    /*

    Response* setResponse(const char** argv, int argc, OPS_Stream& output);
    int getResponse(int responseID, Information& matInformation);

    */


    // Plasticity integration routine
    void plastic_integrator(void);

    // Initialization method
    void initialize(void);

    // Internal functions
    double Kiso(double alpha1);		// isotropic hardening function
    double Kisoprime(double alpha1);	//
    double T(double alpha2);
    double deltaH(double dGamma);

    // Fill vector of state variables for output
    Vector getState();

  protected:

    // Material parameters
    double mKref;			// reference Bulk Modulus 
    double mGref;			// reference Shear Modulus
    
    double mK;			    // bulk modulus 
    double mG;		     	// shear modulus
    double msigma_y;		// yield strength 
    double mrho;			// volumetric term
    double mrho_bar;		// nonassociative flow term
    double mKinf;			// nonlinear isotropic hardening term
    double mKo;			    // nonlinear isotropic hardening term
    double mdelta1; 		// exponential hardening term for drucker prager surface
    double mdelta2;         // exponential hardening term for tension cutoff surface
    double mHard;			// hardening constant
    double mtheta;		    // hardening constant
    double mTo;             // initial tension cutoff strength

    double mKref0;          // record initial Bulk Modulus 
    double mGref0;          // record initial Modulus
    double msigma_y0;       // record initial yield strength
    
    double massDen;         // mass density

    // Internal parameters
    double E;     // elastic modulus
    double nu;    // Poisson ratio 
    double ft;    // tensile yield strength
    double fc;    // compressive yield strength
    double Hk;
    double Hi;
    double r_bar;
    double Kinfinity;
    double Kinit;
    double d1;
    double d2;
    double mDen;
    double Yt0; 
    double bt;
    double at;
    double Yc0;
    double bc;
    double ac;
    double kappa;

    //internal variables
    Vector mEpsilon;
    Vector mEpsilon_n_p;	// plastic strain vector at step n, trail e_p
    Vector mEpsilon_n1_p;	// plastic strain vector at step n+1 
    Vector mSigma;
     
    Vector mBeta_n;		    // backstress at step n, beta_np1_trial = beta_n
    Vector mBeta_n1;		// backstress at step n+1

    double mHprime;	    	// derivative of linear kinematic hardening term 

    double mAlpha1_n;		// alpha1_n
    double mAlpha1_n1;	    // alpha1_n+1
    double mAlpha2_n;		// alpha2_n
    double mAlpha2_n1;	    // alpha2_n+1

    int mElastFlag;    // Flag to determine elastic behavior
    int mFlag;

    Matrix mCe;			// elastic tangent stiffness matrix
    Matrix mCep;			// elastoplastic tangent stiffness matrix
    Vector mI1;			// 2nd Order Identity Tensor	
    Matrix mIIvol;		// IIvol = I1 tensor I1  
    Matrix mIIdev;		// 4th Order Deviatoric Tensor

    Vector mState;		// state vector for output

    // Constants
    static const double one3;
    static const double two3;
    static const double root23;

  private:

    // Variables
    Vector eps;   // strain
    Vector sig;   // stress
    Vector sige;  // effective stress
    Vector eps_p; // plastic strain
    Vector sigeP; // effective stress

    // History variables
    Vector epsCommit;
    Vector sigCommit;
    Vector sigeCommit;
    Vector eps_pCommit;
    Vector sigePCommit;

    // Tangent matrices
    Matrix Ce; 
    Matrix C; 
    Matrix Ccommit; 

};

#endif
