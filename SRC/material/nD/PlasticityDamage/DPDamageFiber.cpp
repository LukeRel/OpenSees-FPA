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

#include <DPDamageFiber.h>

#include <string.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <Information.h>
#include <MaterialResponse.h>
#include <Parameter.h>

static int dFlag1 = 0;	// Turn this on for debug
static int dFlag2 = 0;	// Turn this on for debug on damage subroutine

const double DPDamageFiber::one3 = 1.0 / 3.0;
const double DPDamageFiber::two3 = 2.0 / 3.0;
const double DPDamageFiber::root23 = sqrt(2.0 / 3.0);

#include <elementAPI.h>

static int numDPDamageFiberMaterials = 0;

void* OPS_DPDamageFiber(void)
{
	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 13) {
		opserr << "Want: nDMaterial DPDamageFiber $tag $E $nu $sig_c $sig_t $Hk $Hi\n";
		opserr << "$Yt0, $bt, $at, $Yc0, $bc, $ac, $beta, <$De>)\n";
		return 0;
	}

	int tag;
	double dData[14];
	dData[13] = 0;

	int numData = 1;
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "WARNING invalid integer tag: nDMaterial DPDamageFiber \n";
		return 0;
	}

	numData = numArgs - 1;
	if (OPS_GetDouble(&numData, dData) != 0) {
		opserr << "WARNING invalid data: nDMaterial DPDamageFiber : " << tag << "\n";
		return 0;
	}

	NDMaterial* theMaterial = new DPDamageFiber(tag,
		dData[0], dData[1], // E and nu
		dData[2], dData[3], dData[4], dData[5], // Druker Prager plasticity
		dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12], // Addessi damage
		dData[13]); // Parente degradation

	opserr << "DPDamageFiber memory is allocated!" << endln;
	//for (int i = 0;i < 14;i++) opserr << "dData[" << i << "] = " << dData[i] << endln;
	return theMaterial;
}

//full constructor
DPDamageFiber::DPDamageFiber(int tag, double _E, double _nu, // Parameters
	double _sig_c, double _sig_t, double _Hk, double _Hi, // Plasticity
	double _Yt0, double _bt, double _at, double _Yc0, double _bc, double _ac, double _beta, // Damage
	double _De) // Degradation
	: NDMaterial(tag, ND_TAG_DPDamageFiber),
	E(_E), nu(_nu), sig_c(_sig_c), sig_t(_sig_t), Hk(_Hk), Hi(_Hi),
	Yt0(_Yt0), bt(_bt), at(_at), Yc0(_Yc0), bc(_bc), ac(_ac), beta(_beta), De(_De),
	strain(6), strain_k(6), strain_p_k(6), strain_p(6), strain_e(6),
	stress(6), stress_k(6), zeta_k(6), zeta(6), tangent_e(6, 6), tangent_ep(6, 6),
	I1(6), II1(6, 6), IIvol(6, 6), IIdev(6, 6), II1T(6, 6),
	mState(5), m(3, 3), dam(3),
	strain_m(3), strain_c(3), stress_m(3), stress_c(3),
	C(6, 6), C_mm(3, 3), C_mc(3, 3), C_cm(3, 3), C_cc(3, 3), C_mm_e(3, 3),
	strain_m_k(3), strain_c_k(3), stress_m_k(3), stress_c_k(3), C_k(6, 6), C_mm_k(3, 3), C_cc_k(3, 3),
	fiberstress(3), fibertangent(3, 3)
{
	// Bulk and shear modulus
	K = E / (3.0 * (1.0 - 2.0 * nu));
	G = E / (2.0 * (1.0 + nu));

	// Yielding and friction coefficient
	sig_y = 2.0 * sig_c * sig_t / (sig_c + sig_t);
	mu = root23 * (sig_c - sig_t) / (sig_c + sig_t);

    this->initialize();
}

//null constructor
DPDamageFiber::DPDamageFiber()
    : NDMaterial(),
	E(0.0), nu(0.0), sig_c(0.0), sig_t(0.0), Hk(0.0), Hi(0.0),
	Yt0(0.0), bt(0.0), at(0.0), Yc0(0.0), bc(0.0), ac(0.0), beta(0.0), De(0.0),
    strain(6), strain_k(6), strain_p_k(6), strain_p(6), strain_e(6),
    stress(6), stress_k(6), zeta_k(6), zeta(6), tangent_e(6, 6), tangent_ep(6, 6), tangent(6, 6),
    I1(6), II1(6,6), IIvol(6, 6), IIdev(6, 6), II1T(6, 6),
    mState(5), m(3,3), dam(3),
	strain_m(3), strain_c(3), stress_m(3), stress_c(3),
	C(6, 6), C_mm(3, 3), C_mc(3, 3), C_cm(3, 3), C_cc(3, 3), C_mm_e(3, 3),
	strain_m_k(3), strain_c_k(3), stress_m_k(3), stress_c_k(3), C_k(6, 6), C_mm_k(3, 3), C_cc_k(3, 3),
	fiberstress(3), fibertangent(3, 3)
{
	// Bulk and shear modulus
	K = 0.0;
	G = 0.0;

	// Yielding and friction coefficient
	sig_y = 0.0;
	mu = 0.0;

    this->initialize();
}

//destructor
DPDamageFiber::~DPDamageFiber()
{}

//zero internal variables
void DPDamageFiber::initialize()
{
	strain.Zero();
	strain_k.Zero();
	strain_p_k.Zero();
	strain_p.Zero();
	strain_e.Zero();

	stress.Zero();
	stress_k.Zero();
	zeta_k.Zero();
	zeta.Zero();

	alpha_k = 0.0;
	alpha = 0.0;
	mFlag = 1;

	// Damage variables initialization
	Dt_k = 0.0;
	Dc_k = 0.0;
	D_k = 0.0;
	Dt = Dt_k;
	Dc = Dc_k;
	D = D_k;
	Dm1sq = 1.0;
	dam.Zero();

	// 2nd order Identity Tensor
	I1.Zero();
	I1(0) = 1;
	I1(1) = 1;
	I1(2) = 1;

	// 4th order Identity Tensor
	II1.Zero();
	for (int i = 0;i++;i < 6) II1(i, i) = 1;

	// 4th order Volumetric Tensor
	// IIvol = I1 tensor I1
	IIvol.Zero();
	IIvol(0, 0) = 1;
	IIvol(0, 1) = 1;
	IIvol(0, 2) = 1;
	IIvol(1, 0) = 1;
	IIvol(1, 1) = 1;
	IIvol(1, 2) = 1;
	IIvol(2, 0) = 1;
	IIvol(2, 1) = 1;
	IIvol(2, 2) = 1;

	// 4th order Deviatoric Tensor
	// Note:  this is the contravariant form!
	//        useable for s^a = 2G * IIdev^ab * epsilon_b
	// (Need a different form for s^a = IIdev ^a_b * sigma^a)
	IIdev.Zero();
	IIdev(0, 0) = two3;
	IIdev(0, 1) = -one3;
	IIdev(0, 2) = -one3;
	IIdev(1, 0) = -one3;
	IIdev(1, 1) = two3;
	IIdev(1, 2) = -one3;
	IIdev(2, 0) = -one3;
	IIdev(2, 1) = -one3;
	IIdev(2, 2) = two3;
	IIdev(3, 3) = 0.5;
	IIdev(4, 4) = 0.5;
	IIdev(5, 5) = 0.5;

	C = K * IIvol + 2 * G * IIdev;
	mState.Zero();

	II1T = I1 % I1;

	m.Zero();

	// Full strains vector = (eps_11 eps_22 eps_33 gamma_12 gamma_23 gamma_31)
	// Mantained     = 1 4 6 (eps_11 gamma_12 gamma_13) = received
	// Condensed out = 2 3 5 (eps_22 eps_33 gamma_23) = must be determined and eventually condensed
	int m[3] = { 0, 3, 5 };
	int c[3] = { 1, 2, 4 };
	int perm[6] = { 0, 3, 5, 1, 2, 4 };

	// Vectors at n+1
	strain_m.Zero();
	strain_c.Zero();
	stress_m.Zero();
	stress_c.Zero();

	// Matrixes at n+1
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

}

//make a clone of this material
NDMaterial* DPDamageFiber::getCopy()
{
	DPDamageFiber* clone;
	clone = new DPDamageFiber(this->getTag(), E, nu, sig_c, sig_t, Hk, Hi,
		Yt0, bt, at, Yc0, bc, ac, beta, De);
	return clone;
}

NDMaterial* DPDamageFiber::getCopy(const char* type)
{
	if (strcmp(type, this->getType()) == 0)
		return this->getCopy();

	return 0;
}

//send back type of material
const char* DPDamageFiber::getType() const
{
	return "BeamFiber";
}

//send back order of strain in vector form
int DPDamageFiber::getOrder() const
{
	return 3;
}

//get the strain and integrate plasticity equations
int DPDamageFiber::setTrialStrain(const Vector& strain_from_element)
{
	// Mantained fiber strains
	strain_m = strain_from_element;

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
	D = fmin(D, 0.95);
	Dm1sq = pow(1.0 - D, 2.0);	// [1-D]^2

	// Elastic fiber strains (3)
	Vector fiberstrain_e(3);
	int m[3] = { 0, 3, 5 };
	for (int i = 0;i < 3;i++) fiberstrain_e(i) = strain(m[i]) - strain_p(m[i]);

	// Nota: convergono ->    tangent = [1-D]^2*tangent_ep,   stress = stress_k + Ce*dstrain
	fibertangent = Dm1sq * C_mm;
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
		opserr << "tangent:\n[ "; for (int i = 0;i < 6;i++) { for (int j = 0;j < 6;j++) opserr << C(j, i) << " "; opserr << "]\n"; }
	}

	this->commitState();

	return 0;
}

// Static consistent condensation
void DPDamageFiber::condensate_consistent(void)
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
	strain_c = -1 * dstrain_c;

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
	strain_c = -1 * dstrain_c;

	// Mantained stresses and tangent
	Matrix C_temp(3, 3);
	C_cc.Solve(C_cm, C_temp);
	stress_m += C_mc * strain_c;
	C_mm = C_mm - C_mc * C_temp;

	// Commit fiber operators
	strain_c_k = strain_c;
	stress_m_k = stress_m;	// Mantained stresses (3)
	stress_c_k = stress_c;	// Condensed out stresses (3)
	C_k = C;
	C_mm_k = C_mm;
	C_cc_k = C_cc;

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
void DPDamageFiber::condensate_iterative(void)
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
	double tol = sig_y * 1e-10;
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

//unused trial strain functions
int DPDamageFiber::setTrialStrain(const Vector& v, const Vector& r)
{
	return this->setTrialStrain(v);
}

//plasticity integration routine
void DPDamageFiber::plasticity()
{
	double f1;
	double norm_eta;
	double Invariant_1;
	Vector s(6);
	Vector eta(6);
	Vector dev_ep(6);
	Vector Jact(2);
	Vector n(6);			// normal to the yield surface in strain space

	// Accuracy control
	double tol = 1.0e-6;

	// Elastic strains
	strain_e = strain - strain_p_k;

	// Trial stress
	stress = tangent_e * strain_e;

	// Deviatoric stress
	Invariant_1 = (stress(0) + stress(1) + stress(2));
	s = stress - (Invariant_1 / 3.0) * I1;

	// Vector eta = s-zeta
	eta = s - zeta_k;

	// Norm of eta -> |eta| = |s - zeta|
	norm_eta = sqrt(eta(0) * eta(0) + eta(1) * eta(1) + eta(2) * eta(2)
	+ 2 * (eta(3) * eta(3) + eta(4) * eta(4) + eta(5) * eta(5)));

	// Plastic function
	f1 = norm_eta - root23 * (sig_y + Hi * alpha) + mu * Invariant_1;

	// Check trial state
	if (f1 <= tol) { // Trial state = elastic state - don't need to do any updates.

		// set trial state:
		strain_p = strain_p_k;
		alpha = alpha_k;
		zeta = zeta_k;

		// Elastoplastic matrix
		tangent_ep = tangent_e;

		return;
	}
	else { // Plastic correction required

		// Normal to surface n = eta/|eta|;  (contravaraint)
		if (norm_eta < 1.0e-13) {
			n.Zero();
		}
		else {
			n = eta / norm_eta;
		}

		// Plastic multiplier
		double lambda;
		lambda = f1 / (2 * G + two3 * (Hk + Hi));

		// Plastic strains update
		strain_p = strain_p_k + lambda * n;

		// Back stress variables update
		alpha = alpha_k + root23 * lambda;
		zeta = two3 * Hk * strain_p;

		// Elastic strains
		strain_e = strain_p - strain;

		// Update stress
		stress = tangent_e * (strain_e);

		// Additional terms
		double G2 = pow(G, 2);

		// Additional matrixes
		Matrix nnT(6, 6); // n*n'
		Matrix n1T(6, 6); // n*1'
		nnT = n % n;
		n1T = n % I1;

		// Update tangent
		tangent_ep = tangent_e - lambda * 4 * G / norm_eta * (II1 - one3 * II1T - nnT)
			- (4 * G2 * nnT + 6 * G * K * mu * n1T) / (2 * G + two3 * (Hi + Hk));

	}

	return;
}

// Damage integration routine
void DPDamageFiber::damage()
{
	//////// 1. Principal strains calculation /////////////////////////////////////////////////////////////////
	// Strain tensors from strain vectors (sci)
	Matrix epsilon(3, 3);	epsilon = tens(strain);
	Matrix epsilon_e(3, 3);	epsilon_e = tens(strain_e);

	// Principal strains
	Vector strain_m(3);
	Vector strain_e_m(3);
	Matrix eVect(3, 3);		// Eigenvectors
	eVect.Eigen3(epsilon); for (int i = 0; i < 3;i++) strain_m(i) = eVect(i, i);
	eVect.Eigen3(epsilon_e); for (int i = 0; i < 3;i++) strain_e_m(i) = eVect(i, i);

	//////// 1. Equivalent strains calculation ////////////////////////////////////////////////////////////////
	// Equivalent total and elastic strains 
	Vector e_t(3);
	Vector e_e(3);
	for (int i = 0;i < 3;i++) e_t[i] = (1.0 - 2.0 * nu) * strain_m[i] + nu * (strain_m[0] + strain_m[1] + strain_m[2]);
	for (int i = 0;i < 3;i++) e_e[i] = (1.0 - 2.0 * nu) * strain_e_m[i] + nu * (strain_e_m[0] + strain_e_m[1] + strain_e_m[2]);

	// Positive and negative equivalent total and elastic strains
	Vector e_tp(3);	Vector e_tn(3);
	Vector e_ep(3);	Vector e_en(3);
	for (int i = 0; i < 3; i++) {
		if (e_t[i] > 0.0) { e_tp[i] = e_t[i];	e_tn[i] = 0.0; }
		else { e_tp[i] = 0.0; e_tn[i] = e_t[i]; }
		if (e_e[i] > 0.0) { e_ep[i] = e_e[i]; e_en[i] = 0.0; }
		else { e_ep[i] = 0.0;	e_en[i] = e_e[i]; }
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
	if (Ft < 1e-8) Dt = Dt_k;
	else Dt = (Yt - Yt0) / (at * Yt + bt);
	//Dt = fmax(Dt, Dt_k);
	if (Dt < 0.0)	Dt = 0.0;
	if (Dt > 1.0)	Dt = 1.0;

	// Compressive damage parameter at step n+1
	if (Fc < 1e-8) Dc = Dc_k;
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
		if (fabs(Yt_e) / Yt0 + fabs(Yc_e) / Yc0 < 1.0e-10)  D = D_k;
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
		double eta_t = Yt_e / (Yt0 + Dt_k * (at * Yt_e + bt));
		double eta_c = Yc_e / (Yc0 + Dc_k * (ac * Yc_e + bc));
		if (eta_t == 0) D = D_k;
		else {
			alpha = pow(eta_t, 2) / (pow(eta_t, 2) + pow(eta_c, 2));
			if (alpha < 0.0) alpha = 0.0;
			if (alpha > 1.0) alpha = 1.0;

			// Cumulative damage variable
			D = alpha * Dt + (1.0 - alpha) * Dc;
		}
	}


	// Final check
	//D = fmax(D_k, D);

	if (dFlag2 == 1) {
		opserr << "\n--- Damage routine internal debug ----------\n";
		opserr << "\nEquivalent strains:\n";
		opserr << "eps_m  = [ "; for (int i = 0;i < 3;i++) opserr << strain_m(i) << " "; opserr << "]\n";
		opserr << "eps_em = [ "; for (int i = 0;i < 3;i++) opserr << strain_e_m(i) << " "; opserr << "]\n";
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
		opserr << "alpha_c = " << 1.0 - alpha << "\n";
		opserr << "alpha_t = " << alpha << "\n";
		//opserr << "D = " << D << "\n";    -> Moved to setTrialStrain
		opserr << "\n--- End ------------------------------------\n";
	}
}

// 3D Material state determination
void DPDamageFiber::mat3DState(void)
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

//vector to tensor
const Matrix& DPDamageFiber::tens(const Vector& v)
{
	for (int i = 0; i < 3;i++) m(i, i) = v(i);
	m(0, 1) = v(3);
	m(1, 0) = v(3);
	m(1, 2) = v(4);
	m(2, 1) = v(4);
	m(2, 0) = v(5);
	m(0, 2) = v(5);

	return m;
}

//send back the strain
const Vector& DPDamageFiber::getStrain()
{
	return strain_m;
}

//send back the stress 
const Vector& DPDamageFiber::getStress()
{
	return fiberstress;
}

//send back the tangent 
const Matrix& DPDamageFiber::getTangent()
{
	return fibertangent;
}

//send back the tangent 
const Matrix& DPDamageFiber::getInitialTangent()
{
	return C_mm_e;
}

const Vector& DPDamageFiber::getDamage(void) {
	dam(0) = Dt;
	dam(1) = Dc;
	dam(2) = D;
	return dam;
}

int DPDamageFiber::commitState(void)
{

	// 3D operators
	stress_k = stress;
	strain_k = strain;
	strain_p_k = strain_p;
	alpha_k = alpha;
	zeta_k = zeta;
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
}

int DPDamageFiber::revertToLastCommit(void)
{
	return 0;
}

int DPDamageFiber::revertToStart(void)
{
	if (ops_InitialStateAnalysis) {
		// do nothing, keep state variables from last step
	}
	else {
		// normal call for revertToStart (not initialStateAnalysis)
		this->initialize();
	}

	return 0;
}

Vector DPDamageFiber::getState()
{
	return mState;
}
/*
Response*
DPDamageFiber::setResponse(const char** argv, int argc, OPS_Stream& output)
{
	Response* theResponse = 0;
	const char* matType = this->getType();

	output.tag("NdMaterialOutput");
	output.attr("matType", this->getClassType());
	output.attr("matTag", this->getTag());

	if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0)
		return new MaterialResponse(this, 1, this->getStress());
	else if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0)
		return new MaterialResponse(this, 2, this->getStrain());
	else if (strcmp(argv[0], "state") == 0)
		return new MaterialResponse(this, 3, this->getState());
	else
		return 0;
}

int DPDamageFiber::getResponse(int responseID, Information& matInfo)
{
	switch (responseID) {
	case -1:
		return -1;
	case 1:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getStress();
		return 0;
	case 2:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getStrain();
		return 0;
	case 3:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getState();
		return 0;
	default:
		return -1;
	}
}
*/
int DPDamageFiber::setParameter(const char** argv, int argc, Parameter& param)
{
	if (strcmp(argv[0], "materialState") == 0) {
		// switch elastic/plastic state
		return param.addObject(5, this);
	}
	else if (strcmp(argv[0], "frictionalStrength") == 0) {
		// update rho parameter
		return param.addObject(7, this);
	}
	else if (strcmp(argv[0], "nonassociativeTerm") == 0) {
		// update nonassociative rho_bar parameter
		return param.addObject(8, this);
	}
	else if (strcmp(argv[0], "cohesiveIntercept") == 0) {
		// update zero confinement yield strength
		return param.addObject(9, this);
	}
	else if (strcmp(argv[0], "shearModulus") == 0) {
		// update shear modulus
		return param.addObject(10, this);
	}
	else if (strcmp(argv[0], "bulkModulus") == 0) {
		// update bulk modulus
		return param.addObject(11, this);
	}
	else if (strcmp(argv[0], "updateMaterialStage") == 0) {
		return -1;
	}
	else {
		// invalid parameter type
		opserr << "WARNING: invalid parameter command for DPDamageFiber nDMaterial with tag: " << this->getTag() << endln;
		return -1;
	}

	return -1;
}

int DPDamageFiber::updateParameter(int responseID, Information& info)
{
	if (responseID == 5) {
		// materialState called - update mElasticFlag
		mElastFlag = (int)info.theDouble;
	}
	return 0;
}

int DPDamageFiber::sendSelf(int commitTag, Channel& theChannel)
{
	int res = 0;

	// place data in a vector
	static Vector data(45);
	data(0) = this->getTag();

	data(18) = alpha_k;
	data(20) = mElastFlag;
	data(21) = mFlag;

	data(22) = strain(0);
	data(23) = strain(1);
	data(24) = strain(2);
	data(25) = strain(3);
	data(26) = strain(4);
	data(27) = strain(5);

	data(28) = strain_p_k(0);
	data(29) = strain_p_k(1);
	data(30) = strain_p_k(2);
	data(31) = strain_p_k(3);
	data(32) = strain_p_k(4);
	data(33) = strain_p_k(5);

	data(34) = zeta_k(0);
	data(35) = zeta_k(1);
	data(36) = zeta_k(2);
	data(37) = zeta_k(3);
	data(38) = zeta_k(4);
	data(39) = zeta_k(5);

	data(40) = mState(0);
	data(41) = mState(1);
	data(42) = mState(2);
	data(43) = mState(3);
	data(44) = mState(4);

	res = theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "WARNING: DPDamageFiber::sendSelf - failed to send vector to channel" << endln;
		return -1;
	}

	return 0;
}

int DPDamageFiber::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	int res = 0;

	// receive data
	static Vector data(45);
	res = theChannel.recvVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "WARNING: DPDamageFiber::recvSelf - failed to receive vector from channel" << endln;
		return -1;
	}

	// set member variables
	this->setTag((int)data(0));

	alpha_k = data(18);

	strain(0) = data(22);
	strain(1) = data(23);
	strain(2) = data(24);
	strain(3) = data(25);
	strain(4) = data(26);
	strain(5) = data(27);

	strain_p_k(0) = data(28);
	strain_p_k(1) = data(29);
	strain_p_k(2) = data(30);
	strain_p_k(3) = data(31);
	strain_p_k(4) = data(32);
	strain_p_k(5) = data(33);

	zeta_k(0) = data(34);
	zeta_k(1) = data(35);
	zeta_k(2) = data(36);
	zeta_k(3) = data(37);
	zeta_k(4) = data(38);
	zeta_k(5) = data(39);

	mState(0) = data(40);
	mState(1) = data(41);
	mState(2) = data(42);
	mState(3) = data(43);
	mState(4) = data(44);

	tangent_e = K * IIvol + 2 * G * IIdev;
	tangent_ep = tangent_e;
	tangent = tangent_e;

	return 0;
}

void DPDamageFiber::Print(OPS_Stream& s, int flag)
{
	s << "DPDamageFiber" << endln;
}

