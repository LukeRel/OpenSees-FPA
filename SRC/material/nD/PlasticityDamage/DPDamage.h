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
		double _Yt0, double _bt, double _at, double _Yc0, double _bc, double _ac, double _beta, // Damage
		double _De); // Degradation

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

	double getRho(void) { return 0; };

	double getDamage(void);

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

	//internal variables
	Vector strain_k;		// total strain vector at step n
	Vector strain;		// total strain vector at step n+1
	Vector strain_k_p;	// plastic strain vector at step n, trail e_p
	Vector strain_p;	// plastic strain vector at step n+1 
	Vector strain_e;		// elastic strain vector
	Vector stress_k;		// stress at step n
	Vector stress;			// stress at step n+1

	Vector zeta_k;		// backstress at step n, beta_np1_trial = beta_n
	Vector zeta;		// backstress at step n+1

	double alpha_k;		// alpha1_n
	double alpha;	// alpha1_n+1

	int mElastFlag;    // Flag to determine elastic behavior
	int mFlag;

	Matrix Ce;			// elastic tangent stiffness matrix
	Matrix Cep;		// elastoplastic tangent stiffness matrix
	Matrix Ct;			// damage tangent stiffness matrix
	Vector I1;			// 2nd Order Identity Tensor	
	Matrix II1;		// 4th Order Identity Tensor
	Matrix IIvol;		// IIvol = I1 tensor I1  
	Matrix IIdev;		// 4th Order Deviatoric Tensor
	Matrix II1T;  // 1*1'

	Vector mState;		// state vector for output

	//functions
	void initialize();	// initializes variables

	//plasticity integration routine
	void plastic_integrator(void);

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
	double Dt_k;				// Tensile damage at step (k)
	double Dc_k;				// Compressive damage at step (k)
	double D_k;					// Total damage at previous step (k)

	// Degradation stuff ------------------------------------------------------------------------------
	double De;

};

#endif
