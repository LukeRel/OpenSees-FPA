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

#ifndef DPDamage_h
#define DPDamage_h

// Written in C++: Daniela Fusco, Luca Parente
// Created: 11/21
// 
// Plasticity and damage material based on Di Re et al [2018].
// The plasticity formulation is based on Drucker Prager,
// while the damage formulation is based on Addessi
// 
// The constitutive law computes on three steps:
//  - Elastic prediction
//  - Plastic correction
//  - Damage correction

#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>

#include "DPDamage.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <NDMaterial.h>
#include <ID.h>

class DPDamage : public NDMaterial
{
public:
	// Full Constructor
	DPDamage(int tag, double _E, double _nu, // Parameters
			 double _sig_c, double _sig_t, double _Hk, double _Hi, // Plasticity
			 double _Yt0, double _bt, double _at, double _Yc0, double _bc, double _ac, double _beta); // Damage

	//Null Constructor
	DPDamage();

	//Destructor
	~DPDamage();

	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	NDMaterial* getCopy(void);
	NDMaterial* getCopy(const char* type);

	const char* getType(void) const;
	int getOrder(void) const;

	Response* setResponse(const char** argv, int argc, OPS_Stream& output);
	int getResponse(int responseID, Information& matInformation);

	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

	void Print(OPS_Stream& s, int flag = 0);
	int setParameter(const char** argv, int argc, Parameter& param);
	int updateParameter(int responseID, Information& eleInformation);

	double getRho(void) { return massDen; };

	int setTrialStrain(const Vector& strain_from_element);

	// Unused trialStrain functions
	int setTrialStrain(const Vector& v, const Vector& r);

	//send back the strain
	const Vector& getStrain();

	//send back the stress 
	const Vector& getStress();

	//send back the tangent 
	const Matrix& getTangent();
	const Matrix& getInitialTangent();

protected:

	//material parameters
	double mKref;			// reference Bulk Modulus 
	double mGref;			// reference Shear Modulus
	double mPatm;		    // reference stress first invariant (pressure)
	double mK;			// bulk modulus 
	double mG;			// shear modulus
	double msigma_y;		// yield strength 
	double mrho;			// volumetric term
	double mrho_bar;		// nonassociative flow term
	double mKinf;			// nonlinear isotropic hardening term
	double mKo;			// nonlinear isotropic hardening term
	double mdelta1; 		// exponential hardening term for drucker prager surface
	double mdelta2;       // exponential hardening term for tension cutoff surface
	double mHard;			// hardening constant
	double mtheta;		// hardening constant
	double mTo;           // initial tension cutoff strength
	double E;       // Elastic modulus
	double nu;      // Poisson ratio 
	double K;		// Bulk modulus
	double G;		// Shear modulus
	double sig_c;   // Compressive strength
	double sig_t;   // Tensile strength
	double sig_y;   // Compressive strength
	double mu;		// Friction
	double Hk;      // Kinematic hardening coefficient
	double Hi;      // Isotropic hardening coefficient
	double H;		// Total hardening coefficient
	double theta;	// Relative hardening proportion; full isotropic = 0 < theta < 1 = full kinematic

	double massDen;

	//internal variables
	Vector mEpsilon_n;		// total strain vector at step n
	Vector mEpsilon;		// total strain vector at step n+1
	Vector mEpsilon_n_p;	// plastic strain vector at step n, trail e_p
	Vector mEpsilon_n1_p;	// plastic strain vector at step n+1 
	Vector mEpsilon_e;		// elastic strain vector
	Vector mSigma_n;		// stress at step n
	Vector mSigma;			// stress at step n+1

	Vector mBeta_n;		// backstress at step n, beta_np1_trial = beta_n
	Vector mBeta_n1;		// backstress at step n+1

	double mHprime;		// derivative of linear kinematic hardening term 

	double mAlpha1_n;		// alpha1_n
	double mAlpha1_n1;	// alpha1_n+1
	double mAlpha2_n;		// alpha2_n
	double mAlpha2_n1;	// alpha2_n+1

	int mElastFlag;    // Flag to determine elastic behavior
	int mFlag;

	Matrix mCe;			// elastic tangent stiffness matrix
	Matrix mCep;		// elastoplastic tangent stiffness matrix
	Matrix mCt;			// damage tangent stiffness matrix
	Vector mI1;			// 2nd Order Identity Tensor	
	Matrix mIIvol;		// IIvol = I1 tensor I1  
	Matrix mIIdev;		// 4th Order Deviatoric Tensor

	Vector mState;		// state vector for output

	//functions
	void initialize();	// initializes variables
	int  updateElasticParam(void); //updated Elastic Parameters based on mean stress 

	//plasticity integration routine
	void plastic_integrator(void);

	double Kiso(double alpha1);		// isotropic hardening function
	double Kisoprime(double alpha1);	//
	double T(double alpha2);
	double deltaH(double dGamma);

	Vector getState();		// fills vector of state variables for output

	//vector to tensor notation
	Matrix m;
	const Matrix& tens(const Vector& v);

	//parameters
	static const double one3;
	static const double two3;
	static const double root23;

private:

	// Damage stuff -----------------------------------------------------------------------------------
	void damage();

	// Input parameters
	double Yt0;					// Initial generalized tensile damage variable
	double bt;					// Tensile shape factor
	double at;					// Tensile shape factor
	double Yc0;					// Initial generalized compressive damage variable
	double bc;					// Compression shape factor
	double ac;					// Compression shape factor
	double beta;				// Surface shape factor

	// Damage values
	double Dt;					// Tensile damage at step (k+1)
	double Dc;					// Compressive damage at step (k+1)
	double D;					// Total damage (k+1)
	double Dm1sq;				// [1-D]^2 
	double Dt_n;				// Tensile damage at step (k)
	double Dc_n;				// Compressive damage at step (k)
	double D_n;					// Total damage at previous step (k)

}; //end of DPDamage declarations

#endif
