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

// Written in C++: Luca Parente
// Created: 12/21
// 
// Plasticity and damage material based on Gatta et al [2018].
// The plasticity formulation is based on Von Mises (J2),
// while the damage formulation is based on Addessi
// 
// The constitutive law computes on three steps:
//  - Elastic prediction
//  - Plastic correction
//  - Damage correction


#ifndef J2DamageFiber_h
#define J2DamageFiber_h

#include <NDMaterial.h>

#include <T2Vector.h>
#include <Matrix.h>
#include <Vector.h>


class J2DamageFiber : public NDMaterial
{
public:

	J2DamageFiber(int tag, double _E, double _nu, // Parameters
		double sig_y, double _Hk, double _Hi, // Plasticity
		double _Yt0, double _bt, double _at, double _Yc0, double _bc, double _ac, double _beta, // Damage
		double _De); // Degradation

	J2DamageFiber(const J2DamageFiber&);
	virtual ~J2DamageFiber();

	const char* getClassType(void) const { return "J2DamageFiber"; };
	const char* getType(void) const { return "BeamFiber"; };
	int getOrder() const { return 3; }
	int setTrialStrain(const Vector& fiberStrain);
	int setTrialStrain(const Vector& v, const Vector& r);
	int setTrialStrainIncr(const Vector& v);
	int setTrialStrainIncr(const Vector& v, const Vector& r);

	// Calculates current tangent stiffness.
	const Matrix& getTangent(void);
	const Matrix& getInitialTangent(void);

	// Calculates the corresponding stress increment (rate), for a given strain increment. 
	const Vector& getStress(void);
	const Vector& getStrain(void);
	const Vector& getCommittedStress(void);
	const Vector& getCommittedStrain(void);

	const Vector& getDamage(void);

	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	NDMaterial* getCopy(void);
	NDMaterial* getCopy(const char* code);

	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel,
		FEM_ObjectBroker& theBroker);
	/*
	Response* setResponse(const char** argv, int argc, OPS_Stream& s);
	int getResponse(int responseID, Information& matInformation);
	*/
	void Print(OPS_Stream& s, int flag = 0);

	int setParameter(const char** argv, int argc, Parameter& param);
	int updateParameter(int responseID, Information& eleInformation);

protected:

private:

	// Custom plasticity routine
	int plasticity();

	// Internal parameters
	int ndm;
	double E;       // Elastic modulus
	double nu;      // Poisson ratio
	double K;		// Bulk modulus
	double G;		// Shear modulus
	double sig_y;   // Yield strength
	double Hk;      // Kinematic hardening coefficient
	double Hi;      // Isotropic hardening coefficient

	// Trial variables --------------------------------------------------------------------------------
	Vector stress;				// Stresses vector (k+1)
	Vector strain;				// Total strains vector (k+1)
	Vector strain_e;			// Elastic strains vector (k+1)
	Vector strain_s;			// Vector containing 3 principal total strains
	Vector strain_e_s;			// Vector containint 3 principal elastic strains
	Vector strain_p_dev;		// Deviatoric part of plastic strains (k+1)
	Vector strain_p;			// Total plastic strains (k+1)
	Vector backStress;			// Backstress vector (k+1)
	// Double cumPlastStrainDev;
	
	// ---
	double lambda;				// Positive deviatoric plastic strain increment
	Matrix tangent;				// Material stiffness matrix
	Matrix tangent_e;			// Elastic stiffness matrix
	Matrix tangent_ep;			// Elastoplastic stiffness matrix

	// Define classwide variables
	static Vector tmpVector;
	static Matrix tmpMatrix;

	// Committed variables ----------------------------------------------------------------------------
	Vector stress_k;			// Stresses vector at previous step (k)
	Vector strain_k;			// Total strains vector at previous step (k)
	Vector strain_p_k;			// Plastic strains at previous step (k)
	Vector strain_p_dev_k;		// Deviatoric part of plastic strains at previous step (k)
	Vector backstress_k;		// Total plastic strains at previous step (k)
	double sig_y_k;				// Backstress vector at previous step (k)
	//double CcumPlastStrainDev;
	
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

	// History operators
	double Dt_k;				// Tensile damage at step (k)
	double Dc_k;				// Compressive damage at step (k)
	double D_k;					// Total damage at previous step (k)
	Vector dam;

	// Degradation stuff ------------------------------------------------------------------------------
	double De;

	// Condensation operators -------------------------------------------------------------------------
	int initialize(void); // Condensation operators initialization
	void condensate_consistent(void); // Static condensation
	void condensate_iterative(void); // Static condensation
	void mat3DState(void); // 3D material state determination

	// Final fiber operators
	Vector fiberstress;
	Matrix fibertangent;

	// Set of necessary operators
	Vector strain_m;	// Mantained strains (3)
	Vector strain_c;	// Condensed out strains (3)
	Vector stress_m;	// Mantained stresses (3)
	Vector stress_c;	// Condensed out stresses (3)
	Matrix C;			// 6x6 3D constitutive matrix
	Matrix C_mm;		// 3x3
	Matrix C_mc;		// 3x3
	Matrix C_cm;		// 3x3
	Matrix C_cc;		// 3x3
	Matrix C_mm_e;		// 3x3

	// History operators
	Vector strain_m_k;
	Vector strain_c_k;
	Vector stress_m_k;	// Mantained stresses (3)
	Vector stress_c_k;	// Condensed out stresses (3)
	Matrix C_k;
	Matrix C_mm_k;
	Matrix C_cc_k;

};

#endif
