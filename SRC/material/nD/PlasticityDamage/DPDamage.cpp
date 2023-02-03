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

static int dFlag2 = 0;	// Turn this on for debug
static int dFlag3 = 0;	// Turn this on for debug on damage subroutine
static int d_out = 0; // Set to 1 to produce a out_DPDamage.txt output file
static int step = 0;

#include <DPDamage.h>
#include <math.h>
#include <stdlib.h>
#include <Information.h>
#include <ID.h>
#include <MaterialResponse.h>
#include <Parameter.h>
#include <elementAPI.h>
#include <algorithm>

Matrix DPDamage::tmpMatrix(6, 6);
Vector DPDamage::tmpVector(6);

// --- element: eps(1,1),eps(2,2),eps(3,3),2*eps(1,2),2*eps(2,3),2*eps(1,3) ----
// --- material strain: eps(1,1),eps(2,2),eps(3,3),eps(1,2),eps(2,3),eps(1,3) , same sign ----

void* OPS_DPDamage(void) {

	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 14) {
		opserr << "Want: nDMaterial DPDamage $tag $E $nu $sig_t $sig_c $Hk $Hi\n";
		opserr << "$Yt0, $bt, $at, $Yc0, $bc, $ac, $beta, <$De>)\n";
		return 0;
	}

	int tag;
	double dData[14];
	dData[13] = 0;

	int numData = 1;
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "WARNING invalid integer tag: nDMaterial DPDamage \n";
		return 0;
	}

	numData = numArgs - 1;;
	if (OPS_GetDouble(&numData, dData) != 0) {
		opserr << "WARNING invalid data: nDMaterial DPDamage : " << tag << "\n";
		return 0;
	}

	NDMaterial* theMaterial = new DPDamage(tag,
		dData[0], dData[1], // E and nu
		dData[2], dData[3], dData[4], dData[5], // Druker Prager plasticity
		dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12], // Addessi damage
		dData[13]); // Parente degradation

	//opserr<<"DPDamage memory is allocated!"<<endln;
	return theMaterial;
}

// Full constructor
DPDamage::DPDamage(int tag, double _E, double _nu, // Parameters
	double _sig_t, double _sig_c, double _Hk, double _Hi, // Plasticity
	double _Yt0, double _bt, double _at, double _Yc0, double _bc, double _ac, double _beta, // Damage
	double _De) // Degradation
	: NDMaterial(tag, ND_TAG_DPDamage),
	E(_E), nu(_nu), sig_t(_sig_t), sig_c(_sig_c), Hk(_Hk), Hi(_Hi),
	Yt0(_Yt0), bt(_bt), at(_at), Yc0(_Yc0), bc(_bc), ac(_ac), beta(_beta), De(_De),
	strain(6), strain_e(6), strain_p(6), strain_p_dev(6), strain_m(6), strain_e_m(6),
	strain_k(6), strain_p_k(6), strain_p_dev_k(6),
	tangent(6, 6), tangent_e(6, 6),
	stress(6), stress_k(6), backStress(6), backstress_k(6),
	I2(6), Idev(6, 6), Ivol(6, 6), dam(3)
{
	// Bulk and shear modulus
	K = E / (3. * (1 - 2 * nu));
	G = E / (2. * (1 + nu));

	// Yield and friction
	sig_y = 2. * sig_c * sig_t / (sig_c + sig_t);
	mu = sqrt(2. / 3) * (sig_c - sig_t) / (sig_c + sig_t);

	this->ndm = 3;

	strain.Zero();
	strain_p.Zero();
	strain_e.Zero();
	strain_m.Zero();
	strain_e_m.Zero();
	strain_k.Zero();
	strain_p_k.Zero();
	stress.Zero();
	stress_k.Zero();
	sig_y_k = sig_y;

	lambda = 0.;

	// Elastic tangent
	for (int i = 0; i < 3; i++)	for (int j = 0; j < 3; j++) tangent(i, j) = K - 2.0 / 3.0 * G;
	for (int i = 0; i < 6; i++) tangent(i, i) += 2.0 * G;
	for (int i = 0; i < 6; i++) for (int j = 3; j < 6; j++) tangent(i, j) /= 2.0;
	tangent_e = tangent;
	tangent_ep = tangent;

	// Unit vector order 2
	I2.Zero();
	for (int i = 0; i < 3; i++)	I2(i) = 1.0;

	// Volumetric tensor
	Ivol.Zero();
	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++)	Ivol(i, j) = 1.0;

	// Deviatoric tensor
	Idev.Zero();
	for (int i = 0; i < 6; i++) Idev(i, i) = 1.0;
	for (int i = 0; i < 3; i++)	for (int j = 0; j < 3; j++)	Idev(i, j) -= 1.0 / 3.0;

	// Damage variables initialization
	Dt_k = 0.0;
	Dc_k = 0.0;
	D_k = 0.0;
	Dt = Dt_k;
	Dc = Dc_k;
	D = D_k;
	Dm1sq = 1.0;
	dam.Zero();

	// Initializing damage output file
	if (d_out == 1) {
		using namespace std;
		ofstream outdata;
		outdata.open("out_DPDamage.txt");
		outdata << "Step eps11 eps22 eps33 gam12 gam13 gam23 sig11 sig22 sig33 tau12 tau13 tau33 Dt Dc D" << endln;
		outdata.close();
	}
}

// Null constructor
DPDamage::DPDamage(const DPDamage& a) : NDMaterial(a.getTag(), ND_TAG_DPDamage),
E(0.0), nu(0.0), sig_t(1e10), sig_c(1e10), Hk(0.0), Hi(0.0),
Yt0(0.0), bt(0.0), at(0.0), Yc0(0.0), bc(0.0), ac(0.0), beta(0.0), De(0.0),
tangent(6, 6), tangent_e(6, 6), stress(6), strain(6), I2(6), Idev(6,6), Ivol(6,6),
strain_e(6), strain_p(6), strain_p_dev(6), strain_m(6), strain_e_m(6), backStress(6),
stress_k(6), strain_k(6), strain_p_k(6), strain_p_dev_k(6), backstress_k(6), dam(3)
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
	strain_m.Zero();
	strain_e_m.Zero();
	sig_y = a.sig_y;
	lambda = 0.;
	stress_k.Zero();
	strain_k.Zero();
	sig_y_k = a.sig_y;

	// Elastic tangent
	for (int i = 0; i < 3; i++)	for (int j = 0; j < 3; j++)	tangent(i, j) = K - 2.0 / 3.0 * G;
	for (int i = 0; i < 6; i++)	tangent(i, i) += 2.0 * G;
	for (int i = 0; i < 6; i++)	for (int j = 3; j < 6; j++)	tangent(i, j) /= 2.0;
	tangent_e = tangent;
	tangent_ep = tangent;

	// Unit vector order 2
	I2.Zero();
	for (int i = 0; i < 3; i++)	I2(i) = 1.0;

	// Volumetric tensor
	Ivol.Zero();
	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++)	Ivol(i, j) = 1.0;

	// Deviatoric tensor
	Idev.Zero();
	for (int i = 0; i < 6; i++) Idev(i, i) = 1.0;
	for (int i = 0; i < 3; i++)	for (int j = 0; j < 3; j++)	Idev(i, j) -= 1.0 / 3.0;

	// Damage variables initialization
	Dt_k = 0.0;
	Dc_k = 0.0;
	D_k = 0.0;
	Dt = Dt_k;
	Dc = Dc_k;
	D = D_k;
	Dm1sq = 1.0;
	dam.Zero();
}

NDMaterial*
DPDamage::getCopy(void)
{
	DPDamage* theCopy =
		new DPDamage(this->getTag(), E, nu, sig_t, sig_c, Hk, Hi,
			Yt0, bt, at, Yc0, bc, ac, beta, De);

	return theCopy;
}

NDMaterial*
DPDamage::getCopy(const char* type)
{
	if (strcmp(type, this->getType()) == 0)
		return this->getCopy();

	return 0;
}

int DPDamage::setTrialStrain(const Vector& pStrain)
{
	// Debug
	if (D_k < 0.1) { dFlag2 = 0; dFlag3 = 0; }
	else { dFlag2 = 0; dFlag3 = 0; }

	// Strains from the analysis
	strain = pStrain;

	// ----- Plasticity and damage part ---------------------------------------------------------------- //
	// ENG -> SCI strains (x0.5)         gamma_ij --> eps_ij
	for (int i = 3; i < 6; i++) { strain[i] /= 2.0; }

	// Debug 1
	if (dFlag2 == 1) {
		opserr << "\n----------------------------------------------------------------------------------------------------------------------------\n";
		opserr << "\nStarted new setTrialStrain.\n\n";
		opserr << "Inputs before plasticity and damage (initialized at current step):\n";
		opserr << "strain     = [ "; for (int i = 0; i < 6; i++) opserr << strain(i) << " "; opserr << "]\n";
		opserr << "strain_k   = [ "; for (int i = 0; i < 6; i++) opserr << strain_k(i) << " "; opserr << "]\n";
	}

	// Plasticity routine
	this->plasticity();

	// Tangent
	tangent_ep = tangent;

	// Total plastic strains
	strain_p = strain_p_dev; // + strain_vol; no volumetric plastic strains exist in J2.

	for (int i = 3; i < 6; i++) { strain_p[i] /= 2.0; }

	// Elastoplastic stiffness matrix
	strain_e = strain - strain_p;

	// Debug 2
	if (dFlag2 == 1) {
		opserr << "\nPlasticity executed!\n\n";
		opserr << "Outputs after plasticity only:\n";
		opserr << "strain     = [ "; for (int i = 0; i < 6; i++) opserr << strain(i) << " "; opserr << "]\n";
		opserr << "strain_e   = [ "; for (int i = 0; i < 6; i++) opserr << strain_e(i) << " "; opserr << "]\n";
		opserr << "strain_p   = [ "; for (int i = 0; i < 6; i++) opserr << strain_p(i) << " "; opserr << "]\n";
		opserr << "strain_p_k = [ "; for (int i = 0; i < 6; i++) opserr << strain_p_k(i) << " "; opserr << "]\n";
		opserr << "stress     = [ "; for (int i = 0; i < 6; i++) opserr << stress(i) << " "; opserr << "]\n";
		opserr << "stress_k   = [ "; for (int i = 0; i < 6; i++) opserr << stress_k(i) << " "; opserr << "]\n";
		opserr << "tangent:\n[ "; for (int i = 0; i < 6; i++) { for (int j = 0; j < 6; j++) opserr << tangent(j, i) << " "; opserr << "]\n"; }
	}

	// Damage routine
	this->damage();

	// SCI -> ENG strains (x2.0)         eps_ij --> gamma_ij
	for (int i = 3; i < 6; i++) {
		strain[i] *= 2.0;
		strain_p[i] *= 2.0;
		strain_e[i] *= 2.0;
		}
	// ------------------------------------------------------------------------------------------------- //

	// Incremental quantities
	//D = 0;	// Damage trigger
	//Vector dstress = stress - stress_k;
	//Vector dstrain = strain - strain_k;
	//Vector dstrain_p = strain_p - strain_p_k;
	//Vector dstrain_e = dstrain - dstrain_p;
	//double dD = D - D_k;
	//Matrix d_dstrain_e(6,6);
	//for (int i = 0;i < 6;i++) d_dstrain_e(i,i) = dstrain_e(i) / strain_e(i);

	// Optional degradation correction
	//D = fmax(D, De);

	// Damage correction in order to avoid singularity
	//D = fmin(D, 0.999);

	// Constitutive matrix and stress
	Dm1sq = pow(1.0 - D, 2.0);	// [1-D]^2
	tangent = fmax(1.0e-3,Dm1sq) * tangent_ep;
	stress.addMatrixVector(0, tangent_e, strain_e, Dm1sq);

	// Debug 3
	if (dFlag2 == 1) {
		opserr << "\nD = " << D << "\n\n";
		opserr << "Outputs after both plasticity and damage:\n";
		opserr << "strain     = [ "; for (int i = 0; i < 6; i++) opserr << strain(i) << " "; opserr << "]\n";
		opserr << "strain_k   = [ "; for (int i = 0; i < 6; i++) opserr << strain_k(i) << " "; opserr << "]\n";
		opserr << "strain_e   = [ "; for (int i = 0; i < 6; i++) opserr << strain_e(i) << " "; opserr << "]\n";
		opserr << "strain_p   = [ "; for (int i = 0; i < 6; i++) opserr << strain_p(i) << " "; opserr << "]\n";
		opserr << "strain_p_k = [ "; for (int i = 0; i < 6; i++) opserr << strain_p_k(i) << " "; opserr << "]\n";
		opserr << "stress     = [ "; for (int i = 0; i < 6; i++) opserr << stress(i) << " "; opserr << "]\n";
		opserr << "stress_k   = [ "; for (int i = 0; i < 6; i++) opserr << stress_k(i) << " "; opserr << "]\n";
		opserr << "tangent_ep:\n[ "; for (int i = 0; i < 6; i++) { for (int j = 0; j < 6; j++) opserr << tangent_ep(j, i) << " "; opserr << "]\n"; }
		opserr << "tangent:\n[ "; for (int i = 0; i < 6; i++) { for (int j = 0; j < 6; j++) opserr << tangent(j, i) << " "; opserr << "]\n"; }
	}

	// Internal commits

	// Damage output
	if (d_out == 1) {
		step++;

		// Damage output
		using namespace std;
		ofstream outdata;
		outdata.open("out_DPDamage.txt", ios::app);
		outdata << step << " ";
		outdata << strain(0) << " " << strain(1) << " " << strain(2) << " " << strain(3) << " " << strain(4) << " " << strain(5) << " ";
		outdata << stress(0) << " " << stress(1) << " " << stress(2) << " " << stress(3) << " " << stress(4) << " " << stress(5) << " ";
		outdata << Dt << " " << Dc << " " << D << endln;
		outdata.close();
	}

	return 0;
};

// Plastic integration routine
void DPDamage::plasticity() {

	// Trace = eps_11 + eps_22 + eps_33
	double trace = strain(0) + strain(1) + strain(2);

	// Deviatoric strains = strain - 1/3*trace*I2
	Vector strain_dev(6);
	strain_dev = strain;
	strain_dev.addVector(1.0, I2, -trace / 3.0);

	// Deviatoric stress = 2*G*(strain_dev - strain_p_dev_k)
	Vector Tstress_dev(6);
	Tstress_dev.addVector(0.0, strain_dev, 2. * G);
	Tstress_dev.addVector(1.0, strain_p_dev_k, -2. * G);

	// First invariant = sigt_11 + sigt_22 + sigt_33
	double I1 = 0;
	if (sig_c != sig_t) for (int i = 0; i < 3; i++)	for (int j = 0; j < 3; j++) I1 += tangent_e(i, j) * strain(j);

	// Eta vector = stress_dev_trial - backstress_k
	Vector Teta(6);
	Teta = Tstress_dev;
	Teta.addVector(1.0, backstress_k, -1.0);
	double eta_m = sqrt(Teta && Teta);

	// Yield function = |eta| - sqrt(2/3)*sig_y_k + mu*I1
	double F = eta_m - sqrt(2. / 3) * sig_y_k + mu*I1;

	// Plastic ----------------------------------------------------------------------------------------------------------
	if (F > 0) {

		// Plastic multiplier = F/(2*G + 2/3*(Hi + Hk))
		lambda = F / (2. * G + 2. / 3. * (Hi + Hk));

		if (lambda < 0) {
			opserr << "Fatal:   DPDamage::lambda is less than zero!" << endln;
			exit(-1);
		}

		// New yielding threshold
		sig_y = sig_y_k + pow(2. / 3., 0.5) * Hi * lambda;

		// Normal vector n = eta/|eta|
		Vector n(6);
		n.addVector(0, Teta, 1. / eta_m);

		// Backstress = backstress_k + 2/3*Hk*lambda*n
		backStress.addVector(0.0, backstress_k, 1.0);
		backStress.addVector(1.0, n, 2. / 3. * Hk * lambda);

		// Deviatoric plastic strains = strain_p_dev_k + lambda*n
		strain_p_dev.addVector(0.0, strain_p_dev_k, 1.0);
		strain_p_dev.addVector(1.0, n, lambda);

		// Stress vector = stress_dev_trial - 2*G*lambda*n + K*trace*I2 ----------------------
		stress.addVector(0.0, Tstress_dev, 1.0);
		stress.addVector(1.0, n, -2. * G * lambda);
		stress.addVector(1.0, I2, K * trace);

		// Consistent tangent modulus --------------------------------------------------------
		// C^ep = K*Ivol + 2*G*(1-B)*Idev + 2*G*(B-A)*nn - 2*G*C*n*I2'

		double A = 2. * G / (2 * G + 2. / 3 * (Hk + Hi));
		double B = 2. * G * lambda / eta_m;

		tangent.Zero();
		tangent.addMatrix(0.0, Ivol, K);
		tangent.addMatrix(1.0, Idev, 2. * G * (1 - B));

		tmpMatrix.Zero();  // n@n
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 3; j++)	tmpMatrix(i, j) = n(i) * n(j);
			for (int j = 3; j < 6; j++) tmpMatrix(i, j) = n(i) * n(j) * 2.0;
		}

		tangent.addMatrix(1.0, tmpMatrix, 2. * G * (B - A));

		// If Drucker Prager plasticity
		if (sig_c != sig_t) {

			double C = 3. * K * mu / (2 * G + 2. / 3 * (Hk + Hi));

			tmpMatrix.Zero();  // n@I2
			for (int i = 0; i < 6; i++) {
				for (int j = 0; j < 3; j++)	tmpMatrix(i, j) = n(i) * I2(j);
				for (int j = 3; j < 6; j++) tmpMatrix(i, j) = n(i) * I2(j) * 2.0;
			}

			tangent.addMatrix(1.0, tmpMatrix, -2. * G * C);
		}
	}
	// Elastic ----------------------------------------------------------------------------------------------------------
	else {

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
		for (int i = 0; i < 3; i++)	for (int j = 0; j < 3; j++)	tangent(i, j) = K - 2.0 / 3.0 * G;
		for (int i = 0; i < 6; i++)	tangent(i, i) += 2.0 * G;

	}

	for (int i = 0; i < 6; i++)	for (int j = 3; j < 6; j++)	tangent(i, j) /= 2.0;

};

// Damage integration routine
void DPDamage::damage()
{
	//////// 0. Principal strains calculation /////////////////////////////////////////////////////////////////
	// Strain tensors
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
	double Ft = (Yt - Yt0) - D_k * (at * Yt + bt);
	double Fc = (Yc - Yc0) - D_k * (ac * Yc + bc);

	if (Ft < 0 && Fc < 0) { Dt = Dt_k; Dc = Dc_k; }
	else {

		// Compressive damage parameter at step n+1
		Dc = (Yc - Yc0) / (ac * Yc + bc);
		Dc = fmin(fmax(Dc, Dc_k), 1.0);

		// Tensile damage parameter at step n+1
		// Condition for tensile damage Dt > Dc ---> Dt = max(Dc,Dt)
		Dt = (Yt - Yt0) / (at * Yt + bt);
		Dt = fmax(fmin(fmax(Dt, Dt_k), 1.0), Dc);
	}

	// Damage 2 parameters coefficients
	int algorithm = 1; // 0 = GATTA, 1 = DI RE
	double alpha = 0.0;

	// Weighting coefficients GATTA
	if (algorithm == 0) {
		if (fabs(Yt_e) / Yt0 + fabs(Yc_e) / Yc0 < 1.0e-10) { alpha = 0.0; D = D_k; }
		else {
			alpha = pow(Yt_e / Yt0, 1) / (pow(Yt_e / Yt0, 1) + pow(Yc_e / Yc0, 1));
			if (alpha < 0.0) alpha = 0.0;
			if (alpha > 1.0) alpha = 1.0;
		}
	}
	// Weighting coefficients DI RE
	else {
		double eta_t = Yt_e / (Yt0 + D_k * (at * Yt_e + bt));
		double eta_c = Yc_e / (Yc0 + D_k * (ac * Yc_e + bc));
		if (eta_t == 0 && eta_c == 0) { alpha = 0.5; }
		else {
			alpha = pow(eta_t, 2) / (pow(eta_t, 2) + pow(eta_c, 2));
			alpha = fmax(0.0,alpha);
			alpha = fmin(1.0,alpha);
		}
	}

	// Cumulative damage variable
	D = alpha * Dt + (1.0 - alpha) * Dc;

	if (dFlag3 == 1) {
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
		opserr << "alpha_c = " << 1.0-alpha << "\n";
		opserr << "alpha_t = " << alpha << "\n";
		//opserr << "D = " << D << "\n";    -> Moved to setTrialStrain
		opserr << "\n--- End ------------------------------------\n";
	}
}

// Returns current tangent stiffness.
const Matrix& DPDamage::getTangent(void) {

	if (ndm == 3)
		return tangent;
	else {
		static Matrix workM(3, 3);
		workM(0, 0) = tangent(0, 0);
		workM(0, 1) = tangent(0, 1);
		workM(0, 2) = tangent(0, 3);
		workM(1, 0) = tangent(1, 0);
		workM(1, 1) = tangent(1, 1);
		workM(1, 2) = tangent(1, 3);
		workM(2, 0) = tangent(3, 0);
		workM(2, 1) = tangent(3, 1);
		workM(2, 2) = tangent(3, 3);
		return workM;

	}

};

const Matrix& DPDamage::getInitialTangent(void) {

	return tangent_e;
};

const Vector& DPDamage::getStress(void) {

	return stress;
};

const Vector& DPDamage::getStrain(void) {

	return strain;
};

const Vector& DPDamage::getDamage(void) {
	dam[0] = Dt;
	dam[1] = Dc;
	dam[2] = D;
	return dam;
}

int DPDamage::commitState(void) {

	stress_k = stress;
	strain_k = strain;
	strain_p_k = strain_p;
	strain_p_dev_k = strain_p_dev;
	backstress_k = backStress;
	sig_y_k = sig_y;
	D_k = D;
	Dc_k = Dc;
	Dt_k = Dt;
	//Dc_commit = Dc;
	//Dt_commit = Dt;
	//D_commit = D;

	return 0;
};

int DPDamage::revertToLastCommit(void) {
	//opserr << "Using DPDamage::revertToLastCommit" << endln;
	stress = stress_k;
	strain = strain_k;
	strain_p = strain_p_k;
	strain_p_dev = strain_p_dev_k;
	backStress = backstress_k;
	sig_y = sig_y_k;
	Dc = Dc_commit;
	Dt = Dt_commit;
	D = D_commit;
	return 0;
};

int DPDamage::revertToStart(void) {

	stress.Zero();
	strain.Zero();
	strain_p.Zero();
	strain_p_dev.Zero();
	backStress.Zero();
	sig_y = 2. * sig_c * sig_t / (sig_c + sig_t);
	Dc = 0;
	Dt = 0;
	D = 0;

	stress_k.Zero();
	strain_k.Zero();
	strain_p_k.Zero();
	strain_p_dev_k.Zero();
	backstress_k.Zero();
	sig_y_k = sig_y;
	D_k = 0;
	Dc_k = 0;
	Dt_k = 0;
	Dc_commit = 0;
	Dt_commit = 0;
	D_commit = 0;

	return 0;
}

int DPDamage::sendSelf(int commitTag, Channel& theChannel) {
	// -- to be implemented.

  /*static ID data(1);
  data(0) = tangent;
  return theChannel.sendID(0, commitTag, data);
  */
	return 0;
};

int DPDamage::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker) {
	// -- to be implemented.

  /*static ID data(1);
  int res = theChannel.recvID(0, commitTag, data);
  tangent = data(0);
  return res;
  */
	return 0;
};

Response* DPDamage::setResponse(const char** argv, int argc, OPS_Stream& s) {


	if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0)
		return new MaterialResponse(this, 1, stress);

	else if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0)
		return new MaterialResponse(this, 2, strain);

	else if (strcmp(argv[0], "tangent") == 0 || strcmp(argv[0], "Tangent") == 0)
		return new MaterialResponse(this, 3, tangent);

	else if (strcmp(argv[0], "plasticStrain") == 0 || strcmp(argv[0], "plasticStrain") == 0)
		return new MaterialResponse(this, 4, strain_p);

	else if (strcmp(argv[0], "damage") == 0 || strcmp(argv[0], "Damage") == 0)
		return new MaterialResponse(this, 5, this->getDamage());

	else
		return 0;
}

int DPDamage::getResponse(int responseID, Information& matInfo) {



	switch (responseID) {
	case -1:
		return -1;
	case 1:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = stress;
		return 0;

	case 2:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = strain;
		return 0;

	case 3:
		if (matInfo.theMatrix != 0)
			*(matInfo.theMatrix) = tangent;
		return 0;

	case 4:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = strain_p;
		return 0;

	case 5:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = this->getDamage();
		return 0;

	}



	return 0;
}

void DPDamage::Print(OPS_Stream& s, int flag) {
	// -- to be implemented.
	return;
};

int DPDamage::setParameter(const char** argv, int argc, Parameter& param) {
	// -- to be implemented.


	return 0;
};

int DPDamage::updateParameter(int responseID, Information& eleInformation) {
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
