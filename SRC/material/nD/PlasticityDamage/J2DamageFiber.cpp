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

static int dFlag1 = 0;	// Turn this on for debug
static int dFlag2 = 0;	// Turn this on for debug on damage subroutine

#include <math.h>
#include <stdlib.h>
#include <J2DamageFiber.h>
#include <Information.h>
#include <ID.h>
#include <MaterialResponse.h>
#include <Parameter.h>
#include <elementAPI.h>
#include <algorithm>
#include <Matrix.h>

Matrix J2DamageFiber::tmpMatrix(6, 6);
Vector J2DamageFiber::tmpVector(6);

// --- element: eps(1,1),eps(2,2),eps(3,3),2*eps(1,2),2*eps(2,3),2*eps(1,3) ----
// --- material strain: eps(1,1),eps(2,2),eps(3,3),eps(1,2),eps(2,3),eps(1,3) , same sign ----

void* OPS_J2DamageFiber(void) {

	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 13) {
		opserr << "Want: nDMaterial J2DamageFiber $tag $E $nu $sig_y $Hk $Hi\n";
		opserr << "$Yt0, $bt, $at, $Yc0, $bc, $ac, $beta, <$De>)\n";
		return 0;
	}

	int tag;
	double dData[13]{};
	dData[12] = 0;

	int numData = 1;
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "WARNING invalid integer tag: nDMaterial J2DamageFiber \n";
		return 0;
	}

	numData = numArgs - 1;;
	if (OPS_GetDouble(&numData, dData) != 0) {
		opserr << "WARNING invalid data: nDMaterial J2DamageFiber : " << tag << "\n";
		return 0;
	}

	NDMaterial* theMaterial = new J2DamageFiber(tag,
		dData[0], dData[1], // E and nu
		dData[2], dData[3], dData[4], // Druker Prager plasticity
		dData[5], dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], // Addessi damage
		dData[12]); // Parente degradation

	//opserr<<"J2DamageFiber memory is allocated!"<<endln;
	return theMaterial;
}

// Full constructor
J2DamageFiber::J2DamageFiber(int tag, double _E, double _nu, // Parameters
	double _sig_y, double _Hk, double _Hi, // Plasticity
	double _Yt0, double _bt, double _at, double _Yc0, double _bc, double _ac, double _beta, // Damage
	double _De) // Degradation
	: NDMaterial(tag, ND_TAG_J2DamageFiber),
	E(_E), nu(_nu), sig_y(_sig_y), Hk(_Hk), Hi(_Hi),
	Yt0(_Yt0), bt(_bt), at(_at), Yc0(_Yc0), bc(_bc), ac(_ac), beta(_beta), De(_De),
	tangent(6, 6), tangent_e(6, 6), stress(6), strain(6),
	strain_e(6), strain_p(6), strain_p_dev(6), strain_s(6), strain_e_s(6), backStress(6),
	stress_k(6), strain_k(6), strain_p_k(6), strain_p_dev_k(6), backstress_k(6),
	strain_m(3), strain_c(3), stress_m(3), stress_c(3),
	C(6,6), C_mm(3,3), C_mc(3, 3), C_cm(3, 3), C_cc(3, 3), C_mm_e(3, 3),
	strain_m_k(3), strain_c_k(3), stress_m_k(3), stress_c_k(3), C_k(6,6), C_mm_k(3,3), C_cc_k(3, 3),
	fiberstress(3), fibertangent(3, 3)
{
	// Bulk and shear modulus
	K = E / (3 * (1 - 2 * nu));
	G = E / (2 * (1 + nu));

	this->ndm = 3;

	stress.Zero();
	strain.Zero();
	strain_e.Zero();
	strain_p.Zero();
	strain_e_s.Zero();
	stress_k.Zero();
	strain_k.Zero();
	sig_y_k = sig_y;

	lambda = 0.;

	this->init_condens();

	// Damage variables initialization
	Dt_k = 0.0;
	Dc_k = 0.0;
	D_k = 0.0;
	Dt = Dt_k;
	Dc = Dc_k;
	D = D_k;
	Dm1sq = 1.0;

	//opserr << "Finished constructor." << endln;
}

// Null constructor
J2DamageFiber::J2DamageFiber(const J2DamageFiber& a) : NDMaterial(a.getTag(), ND_TAG_J2DamageFiber),
E(0.0), nu(0.0), sig_y(1e10), Hk(0.0), Hi(0.0),
Yt0(0.0), bt(0.0), at(0.0), Yc0(0.0), bc(0.0), ac(0.0), beta(0.0), De(0.0),
tangent(6, 6), tangent_e(6, 6), stress(6), strain(6),
strain_e(6), strain_p(6), strain_p_dev(6), strain_s(6), strain_e_s(6), backStress(6),
stress_k(6), strain_k(6), strain_p_k(6), strain_p_dev_k(6), backstress_k(6),
strain_m(3), strain_c(3), stress_m(3), stress_c(3),
C(6, 6), C_mm(3, 3), C_mc(3, 3), C_cm(3, 3), C_cc(3, 3), C_mm_e(3, 3),
strain_m_k(3), strain_c_k(3), stress_m_k(3), stress_c_k(3), C_k(6, 6), C_mm_k(3, 3), C_cc_k(3, 3),
fiberstress(3), fibertangent(3,3)
{
	this->ndm = a.ndm;
	this->G = a.G;
	this->K = a.K;
	this->sig_y = a.sig_y;
	this->Hk = a.Hk;
	this->Hi = a.Hi;

	stress.Zero();
	strain.Zero();
	strain_e.Zero();
	strain_s.Zero();
	strain_e_s.Zero();
	sig_y = a.sig_y;
	lambda = 0.;
	stress_k.Zero();
	strain_k.Zero();
	sig_y_k = a.sig_y;

	this->init_condens();

	// Damage variables initialization
	Dt_k = 0.0;
	Dc_k = 0.0;
	D_k = 0.0;
	Dt = Dt_k;
	Dc = Dc_k;
	D = D_k;
	Dm1sq = 1.0;
}

J2DamageFiber::~J2DamageFiber() {
	return;
};

int J2DamageFiber::init_condens(void) {
	// Full strains vector = (eps_11 eps_22 eps_33 gamma_12 gamma_23 gamma_31)
	// Mantained     = 1 4 6 (eps_11 gamma_12 gamma_13) = received
	// Condensed out = 2 3 5 (eps_22 eps_33 eps_23) = must be determined and eventually condensed
	int m[3] = { 0, 3, 5 };
	int c[3] = { 1, 2, 4 };
	int perm[6] = { 0, 3, 5, 1, 2, 4 };

	// Vectors at n+1
	strain_m.Zero();
	strain_c.Zero();
	stress_m.Zero();
	stress_c.Zero();

	// Matrixes at n+1
	C.Zero();
	for (int i = 0; i < 3; i++)	for (int j = 0; j < 3; j++) C(i, j) = K - 2.0 / 3.0 * G;
	for (int i = 0; i < 6; i++) C(i, i) += 2.0 * G;
	for (int i = 0; i < 6; i++) for (int j = 3; j < 6; j++) C(i, j) /= 2.0;
	tangent = C;
	tangent_e = C;
	tangent_ep = C;

	// Other matrixes initialization
	C_mm.Zero(); C_mc.Zero(); C_cm.Zero(); C_cc.Zero();
	int k; int l;
	for (int i = 0;i < 3;i++) for (int j = 0;j < 3;j++) {
		k = m[i]; l = m[j]; C_mm(i, j) = C(k, l);
		k = m[i]; l = c[j]; C_mc(i, j) = C(k, l);
		k = c[i]; l = m[j]; C_cm(i, j) = C(k, l);
		k = c[i]; l = c[j]; C_cc(i, j) = C(k, l);
	}

	// Debug
	if (0 == 1) {
		opserr << "C is: \n"; for (int i = 0;i < 6;i++) { for (int j = 0;j < 6;j++) opserr << C(i, j) << " ";	opserr << endln; }
		opserr << "C_mm is: \n"; for (int i = 0;i < 3;i++) { for (int j = 0;j < 3;j++) opserr << C_mm(i, j) << " ";	opserr << endln; }
		opserr << "C_mc is: \n"; for (int i = 0;i < 3;i++) { for (int j = 0;j < 3;j++) opserr << C_mc(i, j) << " ";	opserr << endln; }
		opserr << "C_cm is: \n"; for (int i = 0;i < 3;i++) { for (int j = 0;j < 3;j++) opserr << C_cm(i, j) << " ";	opserr << endln; }
		opserr << "C_cc is: \n"; for (int i = 0;i < 3;i++) { for (int j = 0;j < 3;j++) opserr << C_cc(i, j) << " ";	opserr << endln; }
	}

	// Saving elastic (3x3) constitutive fiber matrix
	C_mm_e = C_mm;

	// Vectors at n
	strain_m_k.Zero();
	strain_c_k.Zero();
	stress_m_k.Zero();
	stress_c_k.Zero();

	// Matrixes at n
	C_k = C;
	C_mm_k = C_mm;
	C_cc_k = C_cc;

	// Fiber stress and tangent to send back
	fiberstress.Zero();
	fibertangent.Zero();

	return 0;
}

NDMaterial*
J2DamageFiber::getCopy(void)
{
	J2DamageFiber* theCopy =
		new J2DamageFiber(this->getTag(), E, nu, sig_y, Hk, Hi,
			Yt0, bt, at, Yc0, bc, ac, beta, De);

	return theCopy;
}

NDMaterial*
J2DamageFiber::getCopy(const char* type)
{
	if (strcmp(type, this->getType()) == 0)
		return this->getCopy();

	return 0;
}

// Receive fiberstrains(3) -> send to condenssate and make them strain(6) ->
// Send to mat3DState and make them stress(6) and tangent(6,6) ->
// Back to condensate and make them fiberstress(3) and fibertangent(3,3)
// Send back to section
int J2DamageFiber::setTrialStrain(const Vector& fiberStrain)
{
	// Mantained fiber strains
	strain_m = fiberStrain;

	if (dFlag1 == 1) opserr << "New setTrialStrain cycle ----------------------------------------------\n" << endln;

	// Consistent condensation routine [3 > 6 > 3] --------------------------------------
	//
	this->condensate_iterative();
	//this->condensate_consistent();
	//
	// ----------------------------------------------------------------------------------

	// Debug
	if (dFlag1 == 1) {
		opserr << "Condensation complete\n" << endln;
		opserr << "C is: \n"; for (int i = 0;i < 6;i++) { for (int j = 0;j < 6;j++) opserr << C(i, j) << " ";	opserr << endln; }
		opserr << "C_mm is: \n"; for (int i = 0;i < 3;i++) { for (int j = 0;j < 3;j++) opserr << C_mm(i, j) << " ";	opserr << endln; }
		opserr << "C_mc is: \n"; for (int i = 0;i < 3;i++) { for (int j = 0;j < 3;j++) opserr << C_mc(i, j) << " ";	opserr << endln; }
		opserr << "C_cm is: \n"; for (int i = 0;i < 3;i++) { for (int j = 0;j < 3;j++) opserr << C_cm(i, j) << " ";	opserr << endln; }
		opserr << "C_cc is: \n"; for (int i = 0;i < 3;i++) { for (int j = 0;j < 3;j++) opserr << C_cc(i, j) << " ";	opserr << endln; }

		opserr << "\nStress:" << endln;
		opserr << "stress_m =\n"; for (int i = 0;i < 3;i++) opserr << stress_m(i) << " "; opserr << endln;
		opserr << "stress_c =\n"; for (int i = 0;i < 3;i++) opserr << stress_c(i) << " "; opserr << endln;

		opserr << "\nStrain:" << endln;
		opserr << "strain_m =\n"; for (int i = 0;i < 3;i++) opserr << strain_m(i) << " "; opserr << endln;
		opserr << "strain_c =\n"; for (int i = 0;i < 3;i++) opserr << strain_c(i) << " "; opserr << endln;
	}

	// Damage and degradation final correction
	// Performed on condensed in vectors, which size is (3) and (3,3).

	// Optional degradation correction
	D = fmax(D, De);

	// Damage correction in order to avoid singularity
	D = fmin(D,0.95);
	Dm1sq = pow(1.0 - D, 2.0);	// [1-D]^2

	// Elastic fiber strains (3)
	Vector fiberstrain_e(3);
	int m[3] = { 0, 3, 5 };
	for (int i = 0;i < 3;i++) fiberstrain_e(i) = strain(m[i]) - strain_p(m[i]);

	// Nota: convergono ->    tangent = [1-D]^2*Cep,   stress = stress_k + Ce*dstrain
	fibertangent = Dm1sq * C_mm_e;
	fiberstress = Dm1sq * C_mm_e * fiberstrain_e;

	// Debug 3
	if (dFlag1 == 1) {
		opserr << "\nD = " << D << "\n\n";
		opserr << "Outputs after both plasticity and damage:\n";
		opserr << "strain     = [ "; for (int i = 0;i < 6;i++) opserr << strain(i) << " "; opserr << "]\n";
		opserr << "strain_k   = [ "; for (int i = 0;i < 6;i++) opserr << strain_k(i) << " "; opserr << "]\n";
		opserr << "strain_e   = [ "; for (int i = 0;i < 6;i++) opserr << strain_e(i) << " "; opserr << "]\n";
		opserr << "strain_p   = [ "; for (int i = 0;i < 6;i++) opserr << strain_p(i) << " "; opserr << "]\n";
		opserr << "strain_p_k = [ "; for (int i = 0;i < 6;i++) opserr << strain_p_k(i) << " "; opserr << "]\n";
		opserr << "stress     = [ "; for (int i = 0;i < 6;i++) opserr << stress(i) << " "; opserr << "]\n";
		opserr << "stress_k   = [ "; for (int i = 0;i < 6;i++) opserr << stress_k(i) << " "; opserr << "]\n";
		opserr << "tangent_ep:\n[ "; for (int i = 0;i < 6;i++) { for (int j = 0;j < 6;j++) opserr << tangent_ep(j, i) << " "; opserr << "]\n"; }
		opserr << "tangent:\n[ "; for (int i = 0;i < 6;i++) { for (int j = 0;j < 6;j++) opserr << tangent(j, i) << " "; opserr << "]\n"; }
	}

	this->commitState();

	return 0;
};

// Static consistent condensation
void J2DamageFiber::condensate_consistent(void)
{
	// Picking up all the others from the previous steps
	int m[3] = { 0, 3, 5 };
	int c[3] = { 1, 2, 4 };
	strain_c = strain_c_k;
	C_cc = C_cc_k;
	C = C_k;

	// Mantained strain increment
	Vector dstrain_m(3);
	dstrain_m = strain_m - strain_m_k;

	// Condensed out strains
	Vector dstrain_c(3);
	stress_c = C_cc * dstrain_m;
	C_cc.Solve(stress_c, dstrain_c);
	strain_c = -1* dstrain_c;

	// Forming necessary 3D strain vector
	for (int i = 0;i < 3;i++) {
		strain(m[i]) = strain_m(i);
		strain(c[i]) = strain_c(i);
	}

	if (dFlag1 == 1) {
		opserr << "\nCondensate internal debug --------------------------------------------\n" << endln;
		opserr << "Started condensate. Got:" << endln;
		opserr << "strain_m =\n"; for (int i = 0;i < 3;i++) opserr << strain_m(i) << " "; opserr << endln;
		opserr << "strain_c =\n"; for (int i = 0;i < 3;i++) opserr << strain_c(i) << " "; opserr << endln;
		opserr << "strain =\n"; for (int i = 0;i < 6;i++) opserr << strain(i) << " "; opserr << endln;
		opserr << "Moving to mat3DState\n" << endln;
	}

	// Starting mat3DState with "strain(6)" vector ---------------------------------------------
	// 
	// Note that damage final correction must be performed outside!
	// Same goes for degradation. This method sends back:
	//   - stress(6)
	//   - tangent(6,6)
	//   - strain(6)
	this->mat3DState();
	//
	// -----------------------------------------------------------------------------------------

	// Condensation process to revert to fiber size starts here.

	// Condensation matrixes
	C = tangent;
	for (int i = 0;i < 3;i++) {
		// Condensation vectors
		strain_c(i) = strain(c[i]);
		stress_c(i) = stress(c[i]);
		stress_m(i) = stress(m[i]);

		// Condensation matrixes
		for (int j = 0;j < 3;j++) {
			C_mm(i, j) = C(m[i], m[j]);
			C_mc(i, j) = C(m[i], c[j]);
			C_cm(i, j) = C(c[i], m[j]);
			C_cc(i, j) = C(c[i], c[j]);
		}
	}

	// Condensed out strains
	Vector dstress_c(3);
	dstress_c = stress_c - stress_c_k;
	C_cc.Solve(dstress_c, dstrain_c);
	strain_c = -1*dstrain_c;

	// Mantained stresses and tangent
	Matrix C_temp(3, 3);
	C_cc.Solve(C_cm, C_temp);
	stress_m += C_mc * dstrain_c;
	C_mm = C_mm - C_mc * C_temp;

	// Debug
	if (dFlag1 == 1) {
		opserr << "\nDone mat3DState and condensed stuff. Got:\n" << endln;
		opserr << "C is: \n"; for (int i = 0;i < 6;i++) { for (int j = 0;j < 6;j++) opserr << C(i, j) << " ";	opserr << endln; }
		opserr << "C_mm is: \n"; for (int i = 0;i < 3;i++) { for (int j = 0;j < 3;j++) opserr << C_mm(i, j) << " ";	opserr << endln; }
		opserr << "C_mc is: \n"; for (int i = 0;i < 3;i++) { for (int j = 0;j < 3;j++) opserr << C_mc(i, j) << " ";	opserr << endln; }
		opserr << "C_cm is: \n"; for (int i = 0;i < 3;i++) { for (int j = 0;j < 3;j++) opserr << C_cm(i, j) << " ";	opserr << endln; }
		opserr << "C_cc is: \n"; for (int i = 0;i < 3;i++) { for (int j = 0;j < 3;j++) opserr << C_cc(i, j) << " ";	opserr << endln; }

		opserr << "\nStress:" << endln;
		opserr << "stress_m =\n"; for (int i = 0;i < 3;i++) opserr << stress_m(i) << " "; opserr << endln;
		opserr << "stress_c =\n"; for (int i = 0;i < 3;i++) opserr << stress_c(i) << " "; opserr << endln;

		opserr << "\nStrain:" << endln;
		opserr << "strain_m =\n"; for (int i = 0;i < 3;i++) opserr << strain_m(i) << " "; opserr << endln;
		opserr << "strain_c =\n"; for (int i = 0;i < 3;i++) opserr << strain_c(i) << " "; opserr << endln;
		opserr << "\nEnd ------------------------------------------------------------------\n" << endln;
	}

}

// Static consistent condensation
void J2DamageFiber::condensate_iterative(void)
{
	// Picking up all the others from the previous steps
	strain_c = strain_c_k;
	stress_c = stress_c_k;
	C = C_k;

	// Initialize
	Vector dstress_c(3);
	Vector dstrain_c(3);
	int count = 0;
	int maxcount = 50;
	double tol = sig_y*1e-10;
	double norm_sig = 1.0;
	double norm0;

	while ((norm_sig > tol) && (count < maxcount)) {

		count += 1;

		// set norm
		norm_sig = stress_c.Norm();
		if (count == 0)	norm0 = norm_sig;

		// Condensed out strains
		C_cc.Solve(stress_c, dstrain_c);
		strain_c -= dstrain_c;

		// Forming necessary 3D strain vector
		int m[3] = { 0, 3, 5 };
		int c[3] = { 1, 2, 4 };
		for (int i = 0;i < 3;i++) {
			strain(m[i]) = strain_m(i);
			strain(c[i]) = strain_c(i);
		}

		if (dFlag1 == 1) {
			opserr << "\nCondensate internal debug --------------------------------------------\n" << endln;
			opserr << "Started condensate. Got:" << endln;
			opserr << "strain_m =\n"; for (int i = 0;i < 3;i++) opserr << strain_m(i) << " "; opserr << endln;
			opserr << "strain_c =\n"; for (int i = 0;i < 3;i++) opserr << strain_c(i) << " "; opserr << endln;
			opserr << "strain =\n"; for (int i = 0;i < 6;i++) opserr << strain(i) << " "; opserr << endln;
			opserr << "Moving to mat3DState\n" << endln;
		}

		// Starting mat3DState with "strain(6)" vector ---------------------------------------------
		// 
		// Note that damage final correction must be performed outside!
		// Same goes for degradation. This method sends back:
		//   - stress(6)
		//   - tangent(6,6)
		//   - strain(6)
		this->mat3DState();
		//
		// -----------------------------------------------------------------------------------------

		//NDmaterial strain order        = 11, 22, 33, 12, 23, 31  
		//BeamFiberMaterial strain order = 11, 12, 31, 22, 33, 23

		C = tangent;
		for (int i = 0;i < 3;i++) {
			// Condensation vectors
			strain_c(i) = strain(c[i]);
			stress_c(i) = stress(c[i]);
			stress_m(i) = stress(m[i]);

			// Condensation matrixes
			for (int j = 0;j < 3;j++) {
				C_mm(i, j) = C(m[i], m[j]);
				C_mc(i, j) = C(m[i], c[j]);
				C_cm(i, j) = C(c[i], m[j]);
				C_cc(i, j) = C(c[i], c[j]);
			}
		}


	}

	// Debug
	if (1 == 0) {
		opserr << "Done cycles." << endln;
		opserr << "\nDone mat3DState and condensed stuff. Got:\n" << endln;
		opserr << "C is: \n"; for (int i = 0;i < 6;i++) { for (int j = 0;j < 6;j++) opserr << C(i, j) << " ";	opserr << endln; }
		opserr << "C_mm is: \n"; for (int i = 0;i < 3;i++) { for (int j = 0;j < 3;j++) opserr << C_mm(i, j) << " ";	opserr << endln; }
		opserr << "C_mc is: \n"; for (int i = 0;i < 3;i++) { for (int j = 0;j < 3;j++) opserr << C_mc(i, j) << " ";	opserr << endln; }
		opserr << "C_cm is: \n"; for (int i = 0;i < 3;i++) { for (int j = 0;j < 3;j++) opserr << C_cm(i, j) << " ";	opserr << endln; }
		opserr << "C_cc is: \n"; for (int i = 0;i < 3;i++) { for (int j = 0;j < 3;j++) opserr << C_cc(i, j) << " ";	opserr << endln; }

		opserr << "\nStress:" << endln;
		opserr << "stress_m =\n"; for (int i = 0;i < 3;i++) opserr << stress_m(i) << " "; opserr << endln;
		opserr << "stress_c =\n"; for (int i = 0;i < 3;i++) opserr << stress_c(i) << " "; opserr << endln;

		opserr << "\nStrain:" << endln;
		opserr << "strain_m =\n"; for (int i = 0;i < 3;i++) opserr << strain_m(i) << " "; opserr << endln;
		opserr << "strain_c =\n"; for (int i = 0;i < 3;i++) opserr << strain_c(i) << " "; opserr << endln;
		opserr << "\nEnd ------------------------------------------------------------------\n" << endln;
	}

}

// Plastic integration routine
int J2DamageFiber::plasticity() {

	//-----------debug ---------
    /*static int  count =0;

    count ++;

	opserr<<"count: "<<count<<endln;

	if (count == 381)
	opserr<<"count: "<<count<<endln;*/
	// ---------- debug--------

	//double Ctrace = strain_k(0)+strain_k(1)+strain_k(2);
	double trace = strain(0) + strain(1) + strain(2);

	Vector strain_dev(6);
	Vector I2(6);                // unit vector order 2
	I2.Zero();
	for (int i = 0; i < 3; i++)
		I2(i) = 1.0;

	strain_dev = strain;
	strain_dev.addVector(1.0, I2, -trace / 3.0);

	Vector Tstress_dev(6);
	Tstress_dev.addVector(0.0, strain_dev, 2. * G);
	Tstress_dev.addVector(1.0, strain_p_dev_k, -2. * G);

	Vector Teta(6);
	Teta = Tstress_dev;
	Teta.addVector(1.0, backstress_k, -1.0);


	// --- check elastic or plastic--
	double yieldFunction = pow((Teta && Teta), 0.5) - pow(2. / 3, 0.5) * sig_y_k;     // to replace Yn=(2/3)^.5*sig_yn

   // opserr<<"yield function is:"<<yieldFunction<<endln;

	if (yieldFunction > 0) {    // plastic corrector

		lambda = yieldFunction / (2. * G + 2. / 3. * (Hi + Hk));
		//opserr<<"lambda is:"<<lambda<<endln;

		if (lambda < 0) {
			opserr << "Fatal:   J2DamageFiber::lambda is less than zero!" << endln;
			exit(-1);
		}

		sig_y = sig_y_k + pow(2. / 3., 0.5) * Hi * lambda;            //  Note:to replace Y_n+1 = Yn=(2/3)^.5*sig_yn     

		Vector n(6);
		n.addVector(0, Teta, 1. / pow((Teta && Teta), 0.5));

		//Vector eta(6);
		//eta.addVector(0, n, pow( (Teta &&  Teta),0.5)-(2.*G+2./3.*Hk)*lambda);

		backStress.addVector(0.0, backstress_k, 1.0);
		backStress.addVector(1.0, n, 2. / 3. * Hk * lambda);

		strain_p_dev.addVector(0.0, strain_p_dev_k, 1.0);
		strain_p_dev.addVector(1.0, n, lambda);

		//cumPlastStrainDev = CcumPlastStrainDev + pow(2./3.,0.5)*lambda;
		// sig_y = sig_y_k + Hi*pow(2./3., 0.5) * lambda;

		// Stress vector -----------------------------------------------------------
		stress.addVector(0.0, Tstress_dev, 1.0);
		stress.addVector(1.0, n, -2. * G * lambda);
		stress.addVector(1.0, I2, K * trace);

		//double qiu;
		//	qiu=lambda+strain_p_dev_k(0);
		 //opserr<<"qiu is:"<<qiu <<endln;

		// Consistent tangent modulus ----------------------------------------------

		double A = 2. * G / (2. * G + 2. / 3. * Hk + 2. / 3. * Hi);
		double C = 2. * G * lambda / pow(Teta && Teta, 0.5);

		tangent.Zero();

		Matrix I_dev(6, 6);
		I_dev.Zero();

		for (int i = 0; i < 6; i++)  	I_dev(i, i) = 1.0;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				I_dev(i, j) -= 1.0 / 3.0;

		Vector I2(6);   // unit vector order 2
		I2.Zero();

		for (int i = 0; i < 3; i++)
			I2(i) = 1.0;

		tmpMatrix.Zero();

		for (int i = 0; i < 3; i++)       // I2@I2
			for (int j = 0; j < 3; j++)
				tmpMatrix(i, j) = 1.0;

		tangent.addMatrix(0.0, tmpMatrix, K);
		tangent.addMatrix(1.0, I_dev, 2. * G * (1 - C));
		tmpMatrix.Zero();  // n@n

		for (int i = 0; i < 6; i++) {

			for (int j = 0; j < 3; j++)
				tmpMatrix(i, j) = n(i) * n(j);
			for (int j = 3; j < 6; j++)
				tmpMatrix(i, j) = n(i) * n(j) * 2.0;     // To be consistent with the transformation between 4th order tensor and matrix

		}
		tangent.addMatrix(1.0, tmpMatrix, 2. * G * (C - A));

	}
	else {  //elastic case

		// Elastic stress
		sig_y = sig_y_k;
		backStress.addVector(0.0, backstress_k, 1.0);
		strain_p_dev.addVector(0.0, strain_p_dev_k, 1.0);

		// Normal
		Vector n(6);
		n.addVector(0, Teta, 1. / pow((Teta && Teta), 0.5));

		// Stress
		stress.addVector(0.0, Tstress_dev, 1.0);
		stress.addVector(1.0, I2, K * trace);

		// Elastic stiffness
		tangent.Zero();
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				tangent(i, j) = K - 2.0 / 3.0 * G;
		for (int i = 0; i < 6; i++)
			tangent(i, i) += 2.0 * G;

	}


	for (int i = 0; i < 6; i++)
		for (int j = 3; j < 6; j++)
			tangent(i, j) /= 2.0;

	/*	for (int i=0; i<6;i++)
			for(int j=0; j<6; j++)
				 opserr<<"tangent("<<i<<","<<j<<")  is:"<< tangent(i,j)<<endln;
	*/

	return 0;

};

// Damage integration routine
void J2DamageFiber::damage()
{
	//////// 0. Principal strains calculation /////////////////////////////////////////////////////////////////
	// Strain tensors

	// Tolerance
	double tol = 1e-6;

	// Tensor notation
	Matrix epsilon(3, 3);
	Matrix epsilon_e(3, 3);
	for (int i = 0; i < 3;i++) { epsilon(i, i) = strain(i); epsilon_e(i, i) = strain_e(i); };
	epsilon(0, 1) = strain(3);epsilon_e(0, 1) = strain_e(3);
	epsilon(1, 0) = strain(3);epsilon_e(1, 0) = strain_e(3);
	epsilon(1, 2) = strain(4);epsilon_e(1, 2) = strain_e(4);
	epsilon(2, 1) = strain(4);epsilon_e(2, 1) = strain_e(4);
	epsilon(2, 0) = strain(5);epsilon_e(2, 0) = strain_e(5);
	epsilon(0, 2) = strain(5);epsilon_e(0, 2) = strain_e(5);

	// Principal strains
	Vector strain_s(3);
	Vector strain_e_s(3);
	Matrix eVect(3, 3);		// Eigenvectors
	eVect.Eigen3(epsilon); for (int i = 0; i < 3;i++) strain_s(i) = eVect(i, i);
	eVect.Eigen3(epsilon_e); for (int i = 0; i < 3;i++) strain_e_s(i) = eVect(i, i);

	//////// 1. Equivalent strains calculation ////////////////////////////////////////////////////////////////
	// Equivalent total and elastic strains 
	Vector e_t(3);
	Vector e_e(3);
	for (int i = 0;i < 3;i++) e_t[i] = (1.0 - 2.0 * nu) * strain_s[i] + nu * (strain_s[0] + strain_s[1] + strain_s[2]);
	for (int i = 0;i < 3;i++) e_e[i] = (1.0 - 2.0 * nu) * strain_e_s[i] + nu * (strain_e_s[0] + strain_e_s[1] + strain_e_s[2]);

	// Positive and negative equivalent total and elastic strains
	Vector e_tp(3);	Vector e_tn(3);
	Vector e_ep(3);	Vector e_en(3);
	for (int i = 0; i < 3; i++) {
		if (e_t[i] > 0.0) { e_tp[i] = e_t[i];	e_tn[i] = 0.0; }
		else			  { e_tp[i] = 0.0;		e_tn[i] = e_t[i]; }
		if (e_e[i] > 0.0) { e_ep[i] = e_e[i];	e_en[i] = 0.0; }
		else			  { e_ep[i] = 0.0;		e_en[i] = e_e[i]; }
	}

	// Total equivalent tensile and compressive strains
	double Yt = sqrt(e_tp ^ e_tp);
	double Yc = sqrt(e_tn ^ e_tn - beta * ((e_tn[0]) * (e_tn[1]) + (e_tn[1]) * (e_tn[2]) + (e_tn[2]) * (e_tn[0])));

	// Elastic equivalent tensile and compressive strains
	double Yt_e = sqrt(e_ep ^ e_ep);
	double Yc_e = sqrt(e_en ^ e_en - beta * ((e_en[0]) * (e_en[1]) + (e_en[1]) * (e_en[2]) + (e_en[2]) * (e_en[0])));

	////////// 2. Damage coefficients evaluation part ///////////////////////////////////////////////////////////////////////
	// Damage functions at step n
	double Ft = (Yt - Yt0) - Dt_k * (at * Yt + bt);
	double Fc = (Yc - Yc0) - Dc_k * (ac * Yc + bc);

	// Tensile damage parameter at step n+1
	if (Ft < 1e-6) Dt = Dt_k;
	else Dt = (Yt - Yt0) / (at * Yt + bt);
	//Dt = fmax(Dt, Dt_k);
	if (Dt < 0.0)	Dt = 0.0;
	if (Dt > 1.0)	Dt = 1.0;

	// Compressive damage parameter at step n+1
	if (Fc < 0) Dc = Dc_k;
	else Dc = (Yc - Yc0) / (ac * Yc + bc);
	//Dc = fmax(Dc, Dc_k);
	if (Dc < 0.0)	Dc = 0.0;
	if (Dc > 1.0)	Dc = 1.0;

	// Condition for tensile damage Dt > Dc ---> Dt = max(Dc,Dt)
	Dt = fmax(Dc, Dt);

	// Damage 2 parameters coefficients
	int algorithm = 0; // 0 = GATTA, 1 = DI RE
	double alpha = 0.0;

	// Weighting coefficients GATTA
	if (algorithm == 0) {
		if (fabs(Yt_e) / Yt0 + fabs(Yc_e) / Yc0 < tol) D = D_k;
		else {
			alpha = pow(Yt_e / Yt0, 1) / (pow(Yt_e / Yt0, 1) + pow(Yc_e / Yc0, 1));
			if (alpha < 0.0) alpha = 0.0;
			if (alpha > 1.0) alpha = 1.0;

			// Cumulative damage variable
			D = alpha * Dt + (1.0 - alpha) * Dc;
		}
	}
	// Weighting coefficients DI RE
	else {
		double eta_t = Yt_e / (Yt0 + D_k * (at * Yt_e + bt));
		double eta_c = Yc_e / (Yc0 + D_k * (ac * Yc_e + bc));
		if (eta_t < tol) { alpha = 0.0; D = D_k; }
		else {
			alpha = pow(eta_t, 2) / (pow(eta_t, 2) + pow(eta_c, 2));
			if (alpha < 0.0) alpha = 0.0;
			if (alpha > 1.0) alpha = 1.0;

			// Cumulative damage variable
			D = alpha * Dt + (1.0 - alpha) * Dc;
		}
	}

	// unilateral effect trigger
	//D = fmax(D, D_k);

	// Final check
	//D = fmax(D_k, D);

	if (dFlag2 == 1) {
		opserr << "\n--- Damage routine internal debug ----------\n";
		opserr << "\nEquivalent strains:\n";
		opserr << "eps_m  = [ "; for (int i = 0;i < 3;i++) opserr << strain_s(i) << " "; opserr << "]\n";
		opserr << "eps_em = [ "; for (int i = 0;i < 3;i++) opserr << strain_e_s(i) << " "; opserr << "]\n";
		opserr << "e_t    = [ "; for (int i = 0;i < 3;i++) opserr << e_t(i) << " "; opserr << "]\n";
		opserr << "e_e    = [ "; for (int i = 0;i < 3;i++) opserr << e_e(i) << " "; opserr << "]\n";
		opserr << "e_tp   = [ "; for (int i = 0;i < 3;i++) opserr << e_tp(i) << " "; opserr << "]\n";
		opserr << "e_tn   = [ "; for (int i = 0;i < 3;i++) opserr << e_tn(i) << " "; opserr << "]\n";
		opserr << "e_ep   = [ "; for (int i = 0;i < 3;i++) opserr << e_ep(i) << " "; opserr << "]\n";
		opserr << "e_en   = [ "; for (int i = 0;i < 3;i++) opserr << e_en(i) << " "; opserr << "]\n";
		opserr << "Yt   = " << Yt << "\n";
		opserr << "Yc   = " << Yc << "\n";
		opserr << "Yt_e = " << Yt_e << "\n";
		opserr << "Yc_e = " << Yc_e << "\n";

		opserr << "\nDamage variables:\n";
		opserr << "D_c     = " << Dc << "\n";
		opserr << "D_t     = " << Dt << "\n";
		opserr << "alpha_c = " << 1.0-alpha << "\n";
		opserr << "alpha_t = " << alpha << "\n";
		//opserr << "D = " << D << "\n";    -> Moved to setTrialStrain
		opserr << "\n--- End ------------------------------------\n";
	}
}

// 3D Material state determination
void J2DamageFiber::mat3DState(void)
{
	// ----- Plasticity and damage part ---------------------------------------------------------------- //
	// ENG -> SCI strains (x0.5)         gamma_ij --> eps_ij
	for (int i = 3; i < 6; i++) { strain[i] /= 2.0; }

	// Debug 1
	if (dFlag1 == 1) {
		opserr << "\n----------------------------------------------------------------------------------------------------------------------------\n";
		opserr << "\nStarted new setTrialStrain.\n\n";
		opserr << "Inputs before plasticity and damage (initialized at current step):\n";
		opserr << "strain     = [ "; for (int i = 0;i < 6;i++) opserr << strain(i) << " "; opserr << "]\n";
		opserr << "strain_k   = [ "; for (int i = 0;i < 6;i++) opserr << strain_k(i) << " "; opserr << "]\n";
	}

	// Plasticity routine
	this->plasticity();

	// Tangent
	tangent_ep = tangent;

	// Total plastic strains
	strain_p = strain_p_dev; // + strain_vol; no volumetric plastic strains exist in J2.

	// Elastoplastic stiffness matrix
	strain_e = strain - strain_p;

	// Debug 2
	if (dFlag1 == 1) {
		opserr << "\nPlasticity executed!\n\n";
		opserr << "Outputs after plasticity only:\n";
		opserr << "strain     = [ "; for (int i = 0;i < 6;i++) opserr << strain(i) << " "; opserr << "]\n";
		opserr << "strain_e   = [ "; for (int i = 0;i < 6;i++) opserr << strain_e(i) << " "; opserr << "]\n";
		opserr << "strain_p   = [ "; for (int i = 0;i < 6;i++) opserr << strain_p(i) << " "; opserr << "]\n";
		opserr << "strain_p_k = [ "; for (int i = 0;i < 6;i++) opserr << strain_p_k(i) << " "; opserr << "]\n";
		opserr << "stress     = [ "; for (int i = 0;i < 6;i++) opserr << stress(i) << " "; opserr << "]\n";
		opserr << "stress_k   = [ "; for (int i = 0;i < 6;i++) opserr << stress_k(i) << " "; opserr << "]\n";
		opserr << "tangent:\n[ "; for (int i = 0;i < 6;i++) { for (int j = 0;j < 6;j++) opserr << tangent(j, i) << " "; opserr << "]\n"; }
	}

	// Damage routine
	this->damage();

	// SCI -> ENG strains (x2.0)         eps_ij --> gamma_ij
	for (int i = 3; i < 6; i++) { strain[i] *= 2.0;strain_p[i] *= 2.0;strain_e[i] *= 2.0; }
	// ------------------------------------------------------------------------------------------------- //

}

//unused
int J2DamageFiber::setTrialStrain(const Vector& v, const Vector& r) {

	return this->setTrialStrain(v);
};

//unused
int J2DamageFiber::setTrialStrainIncr(const Vector& v) {

	// ----- change to real strain instead of eng. strain
   // ---- since all strain in material is the true strain, not eng.strain. 
	/*
	for (int i = 0; i < 3;i++) {
		tmpVector(i) = v(i);
		tmpVector(i + 3) = v(i + 3) / 2.0;
	}

	if (ndm == 3 && v.Size() == 6)
		strain = strain_k + v;

	else if (ndm == 2 && v.Size() == 3) {
		strain[0] = strain_k[0] + v[0];
		strain[1] = strain_k[1] + v[1];
		strain[2] = 0.0;
		strain[3] = strain_k[2] + v[2];
		strain[4] = 0.0;
		strain[5] = 0.0;
	}
	else {
		opserr << "Fatal:J2DamageFiber:: Material dimension is: " << ndm << endln;
		opserr << "But strain vector size is: " << v.Size() << endln;
		exit(-1);
	}


	this->plasticity();
	*/
	return 0;
};

//unused
int J2DamageFiber::setTrialStrainIncr(const Vector& v, const Vector& r) {

	return this->setTrialStrainIncr(v);
};

// Calculates current tangent stiffness.
const Matrix& J2DamageFiber::getTangent(void) {

	return fibertangent;

};

const Matrix& J2DamageFiber::getInitialTangent(void) {

	return C_mm_e;
};

const Vector& J2DamageFiber::getStress(void) {

	return fiberstress;
};

const Vector& J2DamageFiber::getStrain(void) {

	return strain_m;
};

const Vector& J2DamageFiber::getCommittedStress(void) {

	return stress_m_k;
};

const Vector& J2DamageFiber::getCommittedStrain(void) {

	return strain_m_k;
};

double J2DamageFiber::getDamage(void) {
	
	return D;
}

int J2DamageFiber::commitState(void) {

	// 3D operators
	stress_k = stress;
	strain_k = strain;
	strain_p_k = strain_p;
	strain_p_dev_k = strain_p_dev;
	backstress_k = backStress;
	sig_y_k = sig_y;
	Dc_k = Dc;
	Dt_k = Dt;
	D_k = D;

	// Fiber operators
	strain_m_k = strain_m;
	strain_c_k = strain_c;
	stress_m_k = stress_m;	// Mantained stresses (3)
	stress_c_k = stress_c;	// Condensed out stresses (3)
	C_k = C;
	C_mm_k = C_mm;
	C_cc_k = C_cc;

	return 0;
};

int J2DamageFiber::revertToLastCommit(void) {

	// 3D operators
	stress = stress_k;
	strain = strain_k;
	strain_p = strain_p_k;
	strain_p_dev = strain_p_dev_k;
	backStress = backstress_k;
	sig_y = sig_y_k;
	Dc = Dc_k;
	Dt = Dt_k;
	D = D_k;

	// Fiber operators
	strain_m = strain_m_k;
	strain_c = strain_c_k;
	stress_m = stress_m_k;	// Mantained stresses (3)
	stress_c = stress_c_k;	// Condensed out stresses (3)
	C = C_k;
	C_mm = C_mm_k;
	C_cc = C_cc_k;
	return 0;
};

int J2DamageFiber::revertToStart(void) {
	// -- to be implemented.
	return 0;
}

int J2DamageFiber::sendSelf(int commitTag, Channel& theChannel) {
	// -- to be implemented.

  /*static ID data(1);
  data(0) = tangent;
  return theChannel.sendID(0, commitTag, data);
  */
	return 0;
};

int J2DamageFiber::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker) {
	// -- to be implemented.

  /*static ID data(1);
  int res = theChannel.recvID(0, commitTag, data);
  tangent = data(0);
  return res;
  */
	return 0;
};
/*
Response* J2DamageFiber::setResponse(const char** argv, int argc, OPS_Stream& s) {


	if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0)
		return new MaterialResponse(this, 1, this->getStress());

	else if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0)
		return new MaterialResponse(this, 2, this->getStrain());

	else if (strcmp(argv[0], "tangent") == 0 || strcmp(argv[0], "Tangent") == 0)
		return new MaterialResponse(this, 3, this->getTangent());

	else if (strcmp(argv[0], "damage") == 0 || strcmp(argv[0], "Damage") == 0) {
		return new MaterialResponse(this, 5, this->getDamage());
	}
	else
		return 0;
}



int J2DamageFiber::getResponse(int responseID, Information& matInfo) {

	switch (responseID) {
	case -1:
		return -1;
	case 1:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = this->getStress();
		return 0;

	case 2:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = this->getStrain();
		return 0;

	case 3:
		if (matInfo.theMatrix != 0)
			*(matInfo.theMatrix) = this->getTangent();
		return 0;

	case 5:
		matInfo.setDouble(this->getDamage());
		return 0;

	}
	return 0;
}
*/
void J2DamageFiber::Print(OPS_Stream& s, int flag) {
	// -- to be implemented.
	return;
};


int J2DamageFiber::setParameter(const char** argv, int argc, Parameter& param) {
	// -- to be implemented.


	return 0;
};

int J2DamageFiber::updateParameter(int responseID, Information& eleInformation) {
	// -- to be implemented.
	/*switch (passedParameterID) {
	case -1:
		return -1;

	case 1:
		this->Nd= info.theDouble; //
		break;

	case 2:
		this->G  = info.theDouble;
		break;

	case 3:
		this->K = info.theDouble; //
		break;

	case 4:
		this->SigmaY0  = info.theDouble;
		break;


	case 5:
		this->Hk= info.theDouble; //
		break;

	case 6:
		this->Hi= info.theDouble; //
		break;

	}
	*/


	return 0;
};
