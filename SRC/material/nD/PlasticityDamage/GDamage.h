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


#ifndef GDamage_h
#define GDamage_h

#include <NDMaterial.h>

#include <T2Vector.h>
#include <Matrix.h>
#include <Vector.h>


class GDamage : public NDMaterial
{
public:

	GDamage(int tag, double _E, double _nu, // Parameters
		double _Yt0, double _bt, double _at, double _Yc0, double _bc, double _ac, double _beta); // Damage

	GDamage(const GDamage&);
	virtual ~GDamage();

	const char* getClassType(void) const { return "GDamage"; };
	const char* getType(void) const { return "ThreeDimensional"; };
	int setTrialStrain(const Vector& strain);
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

	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	NDMaterial* getCopy(void);
	NDMaterial* getCopy(const char* code);

	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel,
		FEM_ObjectBroker& theBroker);

	Response* setResponse(const char** argv, int argc, OPS_Stream& s);
	int getResponse(int responseID, Information& matInformation);
	void Print(OPS_Stream& s, int flag = 0);

	int setParameter(const char** argv, int argc, Parameter& param);
	int updateParameter(int responseID, Information& eleInformation);

protected:
	//vector to tensor notation
	Matrix m;
	const Matrix& tens(const Vector& v);

private:
	// Internal parameters
	int ndm;
	double E;       // Elastic modulus
	double nu;      // Poisson ratio
	double K;		// Bulk modulus
	double G;		// Shear modulus

	// Trial variables --------------------------------------------------------------------------------
	Vector stress;				// Stresses vector (k+1)
	Vector strain;				// Total strains vector (k+1)
	Vector strain_e;			// Elastic strains vector (k+1)
	Vector strain_m;			// Vector containing 3 principal total strains
	Vector strain_e_m;			// Vector containint 3 principal elastic strains
	// Double cumPlastStrainDev;
	
	// ---
	double lambda;				// Positive deviatoric plastic strain increment
	Matrix tangent;				// Material stiffness matrix
	Matrix tangent_e;			// Elastic stiffness matrix

	// Define classwide variables
	static Vector tmpVector;
	static Matrix tmpMatrix;

	// Committed variables ----------------------------------------------------------------------------
	Vector stress_k;			// Stresses vector at previous step (k)
	Vector strain_k;			// Total strains vector at previous step (k)
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
	double Dm1sq;				// [1-D]^2 
	double Dt_k;				// Tensile damage at step (k)
	double Dc_k;				// Compressive damage at step (k)
	double D_k;					// Total damage at previous step (k)
};

#endif
