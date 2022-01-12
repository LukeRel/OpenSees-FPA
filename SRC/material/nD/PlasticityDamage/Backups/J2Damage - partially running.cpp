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

static int of = 0;	// Turn this on for debug
static int df = 0;	// Turn this off for debug in plastic state only

#include <math.h>
#include <stdlib.h>
#include <J2Damage.h>
#include <Information.h>
#include <ID.h>
#include <MaterialResponse.h>
#include <Parameter.h>

Matrix J2Damage::tmpMatrix(6, 6);
Vector J2Damage::tmpVector(6);

// --- element: eps(1,1),eps(2,2),eps(3,3),2*eps(1,2),2*eps(2,3),2*eps(1,3) ----
// --- material strain: eps(1,1),eps(2,2),eps(3,3),eps(1,2),eps(2,3),eps(1,3) , same sign ----


#include <elementAPI.h>

void* OPS_J2Damage(void) {

	//opserr<<"OPS_J2Damage is called!"<<endln;

	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 13) {
		opserr << "Want: nDMaterial J2Damage $tag $E $nu $sig_y $Hk $Hi\n";
		opserr << "$Yt0, $bt, $at, $Yc0, $bc, $ac, $beta)\n";
		return 0;
	}

	int tag;
	double dData[12];

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

	NDMaterial* theMaterial = new J2Damage(tag,
		dData[0], dData[1], // E and nu
		dData[2], dData[3], dData[4], // Druker Prager plasticity
		dData[5], dData[6], dData[7], dData[8], dData[9], dData[10], dData[11]); // Addessi damage

	opserr<<"J2Damage memory is allocated!"<<endln;
	return theMaterial;
}

// Full constructor
J2Damage::J2Damage(int tag, double _E, double _nu, // Parameters
	double _sig_y, double _Hk, double _Hi, // Plasticity
	double _Yt0, double _bt, double _at, double _Yc0, double _bc, double _ac, double _beta) // Damage
	: NDMaterial(tag, ND_TAG_J2Damage),
	E(_E), nu(_nu), sig_y(_sig_y), Hk(_Hk), Hi(_Hi),
	Yt0(_Yt0), bt(_bt), at(_at), Yc0(_Yc0), bc(_bc), ac(_ac), beta(_beta),
	tangent(6, 6), tangent_e(6, 6), stress(6), strain(6),
	strain_e(6), strain_p(6), strain_p_dev(6), strain_m(6), strain_e_m(6), backStress(6),
	stress_k(6), strain_k(6), strain_p_k(6), strain_p_dev_k(6), backstress_k(6)
{
	// Bulk and shear modulus
	K = E / (3 * (1 - 2 * nu));
	G = E / (2 * (1 + nu));

	this->ndm = 3;

	stress.Zero();
	strain.Zero();
	strain_e.Zero();
	strain_m.Zero();
	strain_e_m.Zero();
	stress_k.Zero();
	strain_k.Zero();
	sig_y_k = sig_y;

	lambda = 0.;

	// Elastic tangent
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			tangent(i, j) = K - 2.0 / 3.0 * G;
	for (int i = 0; i < 6; i++)
		tangent(i, i) += 2.0 * G;
	for (int i = 0; i < 6; i++)
		for (int j = 3; j < 6; j++)
			tangent(i, j) /= 2.0;
	tangent_e = tangent;
	tangent_ep = tangent;

	// Damage variables initialization
	Dt_k = 0.0;
	Dc_k = 0.0;
	D_k = 0.0;
	Dt = Dt_k;
	Dc = Dc_k;
	D = D_k;
	Dm1sq = 0.0;

}

// Null constructor
J2Damage::J2Damage(const J2Damage& a) : NDMaterial(a.getTag(), ND_TAG_J2Damage),
E(0.0), nu(0.0), sig_y(1e10), Hk(0.0), Hi(0.0),
Yt0(0.0), bt(0.0), at(0.0), Yc0(0.0), bc(0.0), ac(0.0), beta(0.0),
tangent(6, 6), tangent_e(6, 6), stress(6), strain(6),
strain_e(6), strain_p(6), strain_p_dev(6), strain_m(6), strain_e_m(6), backStress(6),
stress_k(6), strain_k(6), strain_p_k(6), strain_p_dev_k(6), backstress_k(6)
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
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			tangent(i, j) = K - 2.0 / 3.0 * G;
	for (int i = 0; i < 6; i++)
		tangent(i, i) += 2.0 * G;
	for (int i = 0; i < 6; i++)
		for (int j = 3; j < 6; j++)
			tangent(i, j) /= 2.0;
	tangent_e = tangent;
	tangent_ep = tangent;

	// Damage variables initialization
	Dt_k = 0.0;
	Dc_k = 0.0;
	D_k = 0.0;
	Dt = Dt_k;
	Dc = Dc_k;
	D = D_k;
	Dm1sq = 0.0;
}


J2Damage::~J2Damage() {
	return;
};

// Plastic integration routine
int J2Damage::plasticity() {

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
			opserr << "Fatal:   J2Damage::lambda is less than zero!" << endln;
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
void J2Damage::damage()
{
	//////// 1. Equivalent strains part ///////////////////////////////////////////////////////////////////////
	// Equivalent total and elastic strains 
	Vector e_tot(3);
	Vector e_el(3);
	for (int i = 0;i < 3;i++) e_tot[i] = (1.0 - 2.0 * nu) * strain_m[i] + nu * (strain_m[0] + strain_m[1] + strain_m[2]);
	for (int i = 0;i < 3;i++) e_el[i] = (1.0 - 2.0 * nu) * strain_e_m[i] + nu * (strain_e_m[0] + strain_e_m[1] + strain_e_m[2]);

	// Positive and negative equivalent total and elastic strains
	Vector e_tp(3);	Vector e_tn(3);
	Vector e_ep(3);	Vector e_en(3);
	for (int i = 0; i < 3; i++) {
		if (e_tot[i] > 0.0) { e_tp[i] = e_tot[i];	e_tn[i] = 0.0; }
		else { e_tp[i] = 0.0; e_tn[i] = e_tot[i]; }
		if (e_el[i] > 0.0) { e_ep[i] = e_el[i]; e_en[i] = 0.0; }
		else { e_ep[i] = 0.0;	e_en[i] = e_el[i]; }
	}

	// Total equivalent tensile and compressive strains
	double Yt = sqrt((e_tp[0]) * (e_tp[0]) + (e_tp[1]) * (e_tp[1]) + (e_tp[2]) * (e_tp[2]));
	double Yc = sqrt((e_tn[0]) * (e_tn[0]) + (e_tn[1]) * (e_tn[1]) + (e_tn[2]) * (e_tn[2]) +
			 beta * ((e_tn[0]) * (e_tn[1]) + (e_tn[1]) * (e_tn[2]) + (e_tn[2]) * (e_tn[0])));

	// Elastic equivalent tensile and compressive strains
	double Yt_el = sqrt((e_ep[0]) * (e_ep[0]) + (e_ep[1]) * (e_ep[1]) + (e_ep[2]) * (e_ep[2]));
	double Yc_el = sqrt((e_en[0]) * (e_en[0]) + (e_en[1]) * (e_en[1]) + (e_en[2]) * (e_en[2]) +
				beta * ((e_en[0]) * (e_en[1]) + (e_en[1]) * (e_en[2]) + (e_en[2]) * (e_en[0])));

	////////// 2. Damage coefficients evaluation part ///////////////////////////////////////////////////////////////////////
	// Damage functions at step n
	double Ft = (Yt - Yt0) - Dt_k * (at * Yt + bt);
	double Fc = (Yc - Yc0) - Dc_k * (ac * Yc + bc);

	// Tensile damage parameter at step n+1
	if (Ft <= 0.0) Dt = Dt_k; // Damage not increasing
	else Dt = (Yt - Yt0) / (at * Yt + bt); // Tensile damage

	// Compressive damage parameter at step n+1
	if (Fc <= 0.0) Dc = Dc_k; // Damage not increasing
	else Dc = (Yc - Yc0) / (ac * Yc + bc); // Compressive damage

	// Condition for tensile damage Dt > Dc ---> Dt = max(Dc,Dt)
	if (Dt <= Dc) Dt = Dc;
	if (Dt > 1.0) Dt = 1.0;
	if (Dc > 1.0)	Dc = 1.0;

	// Weighting coefficients
	double alpt; double alpc;
	if (Yt_el == 0.0 && Yc_el == 0.0) { alpt = 0.0; alpc = 0.0; }
	else { alpt = (fabs(Yt_el)/Yt0) / (fabs(Yt_el)/Yt0 + fabs(Yc_el)/Yc0); alpc = 1.0 - alpt; }

	// Cumulative damage variable Dm1sq = [1-D]^2
	D = alpt * Dt + alpc * Dc;
	Dm1sq = pow(((1.0 - Dt) * alpt + (1.0 - Dc) * alpc), 2.0);

	if (df == 1) {
		opserr << "\n---------- Damage internal terms --------------------------------\n";
		opserr << "\nEquivalent strains:\n";
		opserr << "e_tot = [ "; for (int i = 0;i < 3;i++) opserr << e_tot(i) << " "; opserr << "]\n";
		opserr << "e_el = [ "; for (int i = 0;i < 3;i++) opserr << e_el(i) << " "; opserr << "]\n";
		opserr << "e_tp = [ "; for (int i = 0;i < 3;i++) opserr << e_tp(i) << " "; opserr << "]\n";
		opserr << "e_tn = [ "; for (int i = 0;i < 3;i++) opserr << e_tn(i) << " "; opserr << "]\n";
		opserr << "e_ep = [ "; for (int i = 0;i < 3;i++) opserr << e_ep(i) << " "; opserr << "]\n";
		opserr << "e_en = [ "; for (int i = 0;i < 3;i++) opserr << e_en(i) << " "; opserr << "]\n";
		opserr << "Yt = " << Yt << "\n";
		opserr << "Yc = " << Yc << "\n";
		opserr << "Yt_el = " << Yt_el << "\n";
		opserr << "Yc_el = " << Yc_el << "\n";

		opserr << "\nDamage variables:\n";
		opserr << "Dc = " << Dc << "\n";
		opserr << "Dt = " << Dt << "\n";
		opserr << "alpha_c = " << alpc << "\n";
		opserr << "alpha_t = " << alpt << "\n";
		opserr << "D = " << D << "\n";
		opserr << "\n---------- End --------------------------------------------------\n";
	}

}

int J2Damage::setTrialStrain(const Vector& pStrain)
{
	if (D_k > 0) df = 0; // Debug for D > 0 at previous step - Set df = 1 for debug.

	if (ndm == 3 && pStrain.Size() == 6)
		strain = pStrain;
	else if (ndm == 2 && pStrain.Size() == 3) {
		strain[0] = pStrain[0];
		strain[1] = pStrain[1];
		strain[2] = 0.0;
		strain[3] = pStrain[2];
		strain[4] = 0.0;
		strain[5] = 0.0;
	}
	else {
		opserr << "Fatal:J2Damage:: Material dimension is: " << ndm << endln;
		opserr << "But strain vector size is: " << pStrain.Size() << endln;
		exit(-1);
	}

	// Change to real strain instead of eng. strain
	for (int i = 3; i < 6; i++) {
		strain[i] /= 2.0;
	}

	// Volumetric part of strains
	double trace = strain(0) + strain(1) + strain(2);
	Vector I2(6);                // unit vector order 2
	I2.Zero();
	for (int i = 0; i < 3; i++)
		I2(i) = 1.0;
	Vector strain_vol = 1.0 / 3.0 * trace * I2;

	// Debug 1
	if (df == 1) {
		opserr << "\n----------------------------------------------------------------------------------------------------------------------------\n";
		opserr << "\nStarted new setTrialStrain.\n\n";
		opserr << "Inputs before plasticity and damage (initialized at current step):\n";
		opserr << "strain     = [ "; for (int i = 0;i < 6;i++) opserr << strain(i) << " "; opserr << "]\n";
		opserr << "strain_k   = [ "; for (int i = 0;i < 6;i++) opserr << strain_k(i) << " "; opserr << "]\n";
	}

	// Plasticity routine
	this->plasticity();
	// ------------------------------------------------------------------------------------------------- //
	
	// Total plastic strains
	strain_p = strain_p_dev;// +strain_vol;
	if (df == 1) {
		opserr << "Deviatoric strains:\n";
		opserr << "strain_dev   = [ "; for (int i = 0;i < 6;i++) opserr << strain(i) - strain_vol(i) << " "; opserr << "]\n";
		opserr << "Volumetric strains:\n";
		opserr << "strain_vol   = [ "; for (int i = 0;i < 6;i++) opserr << strain_vol(i) << " "; opserr << "]\n";
		opserr << "Deviatoric plastic strains:\n";
		opserr << "strain_p_dev = [ "; for (int i = 0;i < 6;i++) opserr << strain_p(i) << " "; opserr << "]\n";
		opserr << "Resulting elastic strains:\n";
		opserr << "strain_e     = [ "; for (int i = 0;i < 6;i++) opserr << strain(i) - strain_p(i) << " "; opserr << "]\n";
		opserr << "Total plastic strains:\n";
		opserr << "strain_p     = [ "; for (int i = 0;i < 6;i++) opserr << strain_p(i) + strain_vol(i) << " "; opserr << "]\n";
		opserr << "Resulting elastic strains:\n";
		opserr << "strain_e     = [ "; for (int i = 0;i < 6;i++) opserr << strain(i) + strain_vol(i) - strain_p(i) << " "; opserr << "]\n\n";
	}

	// Elastoplastic stiffness matrix
	strain_e = strain - strain_p;
	Vector dstress = stress - stress_k;
	Vector dstrain = strain - strain_k;
	/*
	for (int i = 0;i < 6;i++)	for (int j = 0;j < 6;j++) {
		if (dstrain(j) == 0) tangent_ep(i, j) = tangent(i,j);
		else tangent_ep(i, j) = dstress(i) / dstrain(j);
	}
	*/

	// Elastic strains
	/*
	Matrix tangent_e_i(6, 6);
	tangent_e_i.Invert(tangent_e);
	strain_e = tangent_e_i * stress;	// eps_e = Ce^-1*sig
	*/
	//opserr << "tangent_e_i = \n["; for (int i = 0;i < 6;i++) { for (int j = 0;j < 6;j++) opserr << tangent_e_i(j, i) << " "; opserr << "]\n"; }

	// Strain Tensors
	Matrix Eps(3, 3);		// Total strain tensor
	for (int i = 0; i < 3;i++) Eps(i, i) = strain(i);
	Eps(0, 1) = strain(3); Eps(1, 2) = strain(4); Eps(2, 0) = strain(5);
	Eps(1, 0) = strain(3); Eps(2, 1) = strain(4); Eps(0, 2) = strain(5);
	Matrix Eps_e(3, 3);		// Elastic strain tensor
	for (int i = 0; i < 3;i++) Eps_e(i, i) = strain_e(i);
	Eps_e(0, 1) = strain_e(3); Eps_e(1, 2) = strain_e(4); Eps_e(2, 0) = strain_e(5);
	Eps_e(1, 0) = strain_e(3); Eps_e(2, 1) = strain_e(4); Eps_e(0, 2) = strain_e(5);

	// Principal strains
	Matrix eVect(3, 3);		// Eigenvectors
	eVect.Eigen3(Eps); for (int i = 0; i < 3;i++) strain_m(i) = eVect(i, i);
	eVect.Eigen3(Eps_e); for (int i = 0; i < 3;i++) strain_e_m(i) = eVect(i, i);

	// Debug 2
	if (df == 1) {
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
	// ------------------------------------------------------------------------------------------------- //

	// Change back to eng. strains
	for (int i = 3; i < 6; i++) {
		strain[i] *= 2.0;
		strain_p[i] *= 2.0;
		strain_e[i] *= 2.0;
	}

	// Stress damage correction
	//D = 0;	// Damage trigger
	dstrain = strain - strain_k;
	Vector dstrain_p = strain_p - strain_p_k;
	Vector dstrain_e = dstrain - dstrain_p;
	double dD = D - D_k;

	// Constitutive matrix and stress
	//tangent.addMatrix(0,tangent,Dam); // tangent = [1-D]^2*tangent
	//tangent = Dam *tangent_e;
	//stress = tangent*(strain_e);

	//Nota: stress = tangent_e * strain_e restituisce le tensioni giuste post plasticità

	/* Procedura Addessi
	tangent = pow(1.0 - D, 2.0) * tangent;	// Secant damage
	stress = stress_k + pow(1.0 - D, 2.0) * tangent_e * dstrain_e - 2.0 * (1.0 - dD) * tangent_e * dstrain;	// Exact stress

	// Procedura Di Re
	Matrix Ct(6, 6);
	Vector dD_dstrain(6);
	double dD_e = 0.0;
	for (int i = 0;i < 6;i++) {
		if (dstrain(i) == 0) dD_dstrain(i) = 0.0;
		else dD_dstrain(i) = dD/dstrain(i);
		dD_e += dD_dstrain(i) * strain_e(i);
	}
	tangent = pow(1 - D, 2.0) * tangent -2 * (1 - D) * tangent_e * dD_e;
	stress = stress_k + tangent * dstrain;
	*/

	// Nota: convergono ->    tangent = [1-D]^2*Cep,   stress = stress_k + Ce*dstrain
	tangent = pow(1-D,2)*tangent;
	stress = stress_k + tangent * dstrain;

	// Commits
	stress_k = stress;
	strain_k = strain;
	/*
	strain_p_k = strain_p;
	strain_p_dev_k = strain_p_dev;
	backstress_k = backStress;
	sig_y_k = sig_y;
	Dc_k = Dc;
	Dt_k = Dt;
	D_k = D;
	*/

	// Debug 3
	if (df == 1) {
		opserr << "\nDamage executed!\n\n";
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

	return 0;
};

int J2Damage::setTrialStrain(const Vector& v, const Vector& r) {

	return this->setTrialStrain(v);
};

int J2Damage::setTrialStrainIncr(const Vector& v) {

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
		opserr << "Fatal:J2Damage:: Material dimension is: " << ndm << endln;
		opserr << "But strain vector size is: " << v.Size() << endln;
		exit(-1);
	}


	this->plasticity();
	*/
	return 0;
};

int J2Damage::setTrialStrainIncr(const Vector& v, const Vector& r) {



	return this->setTrialStrainIncr(v);


};

// Calculates current tangent stiffness.
const Matrix& J2Damage::getTangent(void) {

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

const Matrix& J2Damage::getInitialTangent(void) {

	return tangent_e;
};

const Vector& J2Damage::getStress(void) {

	return stress;
};

const Vector& J2Damage::getStrain(void) {

	return strain;
};

const Vector& J2Damage::getCommittedStress(void) {

	return stress_k;
};

const Vector& J2Damage::getCommittedStrain(void) {

	return strain_k;
};


int J2Damage::commitState(void) {

	stress_k = stress;
	strain_k = strain;
	strain_p_k = strain_p;
	strain_p_dev_k = strain_p_dev;
	backstress_k = backStress;
	sig_y_k = sig_y;
	Dc_k = Dc;
	Dt_k = Dt;
	D_k = D;
	//CcumPlastStrainDev = cumPlastStrainDev;

	return 0;
};

int J2Damage::revertToLastCommit(void) {

	// -- to be implemented.
	return 0;
};

int J2Damage::revertToStart(void) {
	// -- to be implemented.
	return 0;
}



NDMaterial*
J2Damage::getCopy(void)
{
	J2Damage* theCopy =
		new J2Damage(this->getTag(), E, nu, sig_y, Hk, Hi,
			Yt0, bt, at, Yc0, bc, ac, beta);

	return theCopy;
}

NDMaterial*
J2Damage::getCopy(const char* type)
{
	if (strcmp(type, this->getType()) == 0)
		return this->getCopy();

	return 0;
}


int J2Damage::sendSelf(int commitTag, Channel& theChannel) {
	// -- to be implemented.

  /*static ID data(1);
  data(0) = tangent;
  return theChannel.sendID(0, commitTag, data);
  */
	return 0;
};

int J2Damage::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker) {
	// -- to be implemented.

  /*static ID data(1);
  int res = theChannel.recvID(0, commitTag, data);
  tangent = data(0);
  return res;
  */
	return 0;
};


Response* J2Damage::setResponse(const char** argv, int argc, OPS_Stream& s) {


	if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0)
		return new MaterialResponse(this, 1, stress);

	else if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0)
		return new MaterialResponse(this, 2, strain);

	else if (strcmp(argv[0], "tangent") == 0 || strcmp(argv[0], "Tangent") == 0)
		return new MaterialResponse(this, 3, tangent);


	else if (strcmp(argv[0], "plasticStrainDev") == 0 || strcmp(argv[0], "plasticStrainDevs") == 0)
		return new MaterialResponse(this, 4, strain_p_dev);

	else if (strcmp(argv[0], "damage") == 0 || strcmp(argv[0], "Damage") == 0)
		return new MaterialResponse(this, 5, D);

	else
		return 0;

}



int J2Damage::getResponse(int responseID, Information& matInfo) {



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
			*(matInfo.theVector) = strain_p_dev;
		return 0;

	case 5:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = D;
		return 0;

	}



	return 0;
};

void J2Damage::Print(OPS_Stream& s, int flag) {
	// -- to be implemented.
	return;
};


int J2Damage::setParameter(const char** argv, int argc, Parameter& param) {
	// -- to be implemented.


	return 0;
};

int J2Damage::updateParameter(int responseID, Information& eleInformation) {
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
