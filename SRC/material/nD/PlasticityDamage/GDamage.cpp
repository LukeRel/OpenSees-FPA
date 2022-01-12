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

static int dFlag1 = 1;	// Turn this on for debug
static int dFlag2 = 1;	// Turn this on for debug on damage subroutine

#include <math.h>
#include <stdlib.h>
#include <GDamage.h>
#include <Information.h>
#include <ID.h>
#include <MaterialResponse.h>
#include <Parameter.h>
#include <elementAPI.h>
#include <algorithm>

Matrix GDamage::tmpMatrix(6, 6);
Vector GDamage::tmpVector(6);

// --- element: eps(1,1),eps(2,2),eps(3,3),2*eps(1,2),2*eps(2,3),2*eps(1,3) ----
// --- material strain: eps(1,1),eps(2,2),eps(3,3),eps(1,2),eps(2,3),eps(1,3) , same sign ----

void* OPS_GDamage(void) {

	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 10) {
		opserr << "Want: nDMaterial GDamage $tag $E $nu\n";
		opserr << "$Yt0, $bt, $at, $Yc0, $bc, $ac, $beta)\n";
		return 0;
	}

	int tag;
	double dData[9];

	int numData = 1;
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "WARNING invalid integer tag: nDMaterial GDamage \n";
		return 0;
	}

	numData = numArgs - 1;;
	if (OPS_GetDouble(&numData, dData) != 0) {
		opserr << "WARNING invalid data: nDMaterial GDamage : " << tag << "\n";
		return 0;
	}

	NDMaterial* theMaterial = new GDamage(tag,
		dData[0], dData[1], // E and nu
		dData[2], dData[3], dData[4], dData[5], dData[6], dData[7], dData[8]); // Addessi damage

	opserr<<"GDamage memory is allocated!"<<endln;
	return theMaterial;
}

// Full constructor
GDamage::GDamage(int tag, double _E, double _nu, // Parameters
	double _Yt0, double _bt, double _at, double _Yc0, double _bc, double _ac, double _beta) // Damage
	: NDMaterial(tag, ND_TAG_GDamage),
	E(_E), nu(_nu),
	Yt0(_Yt0), bt(_bt), at(_at), Yc0(_Yc0), bc(_bc), ac(_ac), beta(_beta),
	tangent(6, 6), tangent_e(6, 6), stress(6), strain(6),
	strain_e(6), strain_m(6), strain_e_m(6),
	stress_k(6), strain_k(6),
	m(3, 3)
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

	m.Zero();

	lambda = 0.;

	// Elastic tangent
	for (int i = 0; i < 3; i++)	for (int j = 0; j < 3; j++) tangent(i, j) = K - 2.0 / 3.0 * G;
	for (int i = 0; i < 6; i++) tangent(i, i) += 2.0 * G;
	for (int i = 0; i < 6; i++) for (int j = 3; j < 6; j++) tangent(i, j) /= 2.0;
	tangent_e = tangent;

	// Damage variables initialization
	Dt_k = 0.0;
	Dc_k = 0.0;
	D_k = 0.0;
	Dt = Dt_k;
	Dc = Dc_k;
	D = D_k;
	Dm1sq = 1.0;

}

// Null constructor
GDamage::GDamage(const GDamage& a) : NDMaterial(a.getTag(), ND_TAG_GDamage),
E(0.0), nu(0.0),
Yt0(0.0), bt(0.0), at(0.0), Yc0(0.0), bc(0.0), ac(0.0), beta(0.0),
tangent(6, 6), tangent_e(6, 6), stress(6), strain(6),
strain_e(6), strain_m(6), strain_e_m(6),
stress_k(6), strain_k(6),
m(3, 3)
{
	// Bulk and shear modulus
	K = E / (3 * (1 - 2 * nu));
	G = E / (2 * (1 + nu));

	this->ndm = a.ndm;

	stress.Zero();
	strain.Zero();
	strain_e.Zero();
	strain_m.Zero();
	strain_e_m.Zero();
	lambda = 0.;
	stress_k.Zero();
	strain_k.Zero();

	m.Zero();

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

	// Damage variables initialization
	Dt_k = 0.0;
	Dc_k = 0.0;
	D_k = 0.0;
	Dt = Dt_k;
	Dc = Dc_k;
	D = D_k;
	Dm1sq = 1.0;
}

GDamage::~GDamage() {
	return;
};

NDMaterial*
GDamage::getCopy(void)
{
	GDamage* theCopy =
		new GDamage(this->getTag(), E, nu, Yt0, bt, at, Yc0, bc, ac, beta);

	return theCopy;
}

NDMaterial*
GDamage::getCopy(const char* type)
{
	if (strcmp(type, this->getType()) == 0)
		return this->getCopy();

	return 0;
}

// Damage integration routine
void GDamage::damage()
{
	//Tolerance
	double tol = 1.0e-10;
	
	//////// 0. Principal strains calculation /////////////////////////////////////////////////////////////////
	// Strain tensors
	// Strain tensors from strain vectors (sci)
	Matrix epsilon(3, 3);	 epsilon = tens(strain);
	Matrix epsilon_e(3, 3);	 epsilon_e = tens(strain_e);

	// Principal strains
	Vector strain_m(3);
	Vector strain_e_m(3);
	Matrix eVect(3, 3);		 // Eigenvectors
	eVect.Eigen3(epsilon);   for (int i = 0; i < 3;i++) strain_m(i) = eVect(i, i);
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
	//double Ft = (Yt - Yt0) - Dt_k * (at * Yt + bt);
	//double Fc = (Yc - Yc0) - Dc_k * (ac * Yc + bc); -> not really necessary for algorithm

	// Tensile damage parameter at step n+1
	Dt = (Yt - Yt0) / (at * Yt + bt);
	Dt = fmax(Dt, Dt_k);
	if (Dt < 0.0)	Dt = 0.0;
	if (Dt > 1.0)	Dt = 1.0;

	// Compressive damage parameter at step n+1
	Dc = (Yc - Yc0) / (ac * Yc + bc);
	Dc = fmax(Dc,Dc_k);
	if (Dc < 0.0)	Dc = 0.0;
	if (Dc > 1.0)	Dc = 1.0;

	// Condition for tensile damage Dt > Dc ---> Dt = max(Dc,Dt)
	Dt = fmax(Dc, Dt);

	// Weighting coefficients
	double alpha;
	if (Yt_e == 0.0 && Yc_e == 0.0) alpha = 0.0;
	else alpha = (fabs(Yt_e)/Yt0) / (fabs(Yt_e)/Yt0 + fabs(Yc_e)/Yc0);
	if (alpha < 0.0) alpha = 0.0;
	if (alpha > 1.0) alpha = 1.0;

	// Cumulative damage variable Dm1sq = [1-D]^2
	D = alpha * Dt + (1.0-alpha) * Dc;

	if (fabs(Yt_e) / Yt0 + fabs(Yc_e) / Yc0 < tol) D = D_k;

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
		opserr << "alpha_c = " << 1.0-alpha << "\n";
		opserr << "alpha_t = " << alpha << "\n";
		//opserr << "D = " << D << "\n";    -> Moved to setTrialStrain
		opserr << "\n--- End ------------------------------------\n";
	}
}

//vector to tensor
const Matrix& GDamage::tens(const Vector& v)
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

int GDamage::setTrialStrain(const Vector& pStrain)
{
	// Debug
	if (D_k < 0.1) {dFlag1 = 0; dFlag2 = 0;}
	else { dFlag1 = 0; dFlag2 = 0; }

	// Strains from the analysis
	strain = pStrain;
	
	// ----- Plasticity and damage part ---------------------------------------------------------------- //
	// ENG -> SCI strains (x0.5)         gamma_ij --> eps_ij
	for (int i = 3; i < 6; i++) {strain[i] /= 2.0;}

	// Debug 1
	if (dFlag1 == 1) {
		opserr << "\n----------------------------------------------------------------------------------------------------------------------------\n";
		opserr << "\nStarted new setTrialStrain.\n\n";
		opserr << "Inputs before plasticity and damage (initialized at current step):\n";
		opserr << "strain     = [ "; for (int i = 0;i < 6;i++) opserr << strain(i) << " "; opserr << "]\n";
		opserr << "strain_k   = [ "; for (int i = 0;i < 6;i++) opserr << strain_k(i) << " "; opserr << "]\n";
	}

	// All strains are elastic before damage occurs
	strain_e = strain;
	
	// Damage routine
	this->damage();
	
	// SCI -> ENG strains (x2.0)         eps_ij --> gamma_ij
	for (int i = 3; i < 6; i++) { strain[i] *= 2.0; strain_e[i] *= 2.0; }
	// ------------------------------------------------------------------------------------------------- //

	// Incremental quantities
	//D = 0;	// Damage trigger
	//Vector dstress = stress - stress_k;
	Vector dstrain = strain - strain_k;
	Vector dstrain_e = dstrain;
	double dD = D - D_k;

	// Constitutive matrix and stress ---------------------------------------------------------------- //
	// Damage correction in order to avoid singularity
	Dm1sq = pow(1.0 - D, 2.0);	// [1-D]^2

	// Nota: convergono ->    tangent = [1-D]^2*Cep,   stress = stress_k + Ce*dstrain
	tangent = Dm1sq * tangent_e;
	//stress = stress_k + Dm1sq * tangent_e * dstrain_e - 2.0 * (1.0 - D) * dD * tangent_e * strain_e;
	stress = tangent*strain_e;

	// Debug 2
	if (dFlag1 == 1) {
		opserr << "\nD = " << D << "\n\n";
		opserr << "Outputs after both plasticity and damage:\n";
		opserr << "strain     = [ "; for (int i = 0;i < 6;i++) opserr << strain(i) << " "; opserr << "]\n";
		opserr << "strain_k   = [ "; for (int i = 0;i < 6;i++) opserr << strain_k(i) << " "; opserr << "]\n";
		opserr << "strain_e   = [ "; for (int i = 0;i < 6;i++) opserr << strain_e(i) << " "; opserr << "]\n";
		opserr << "stress     = [ "; for (int i = 0;i < 6;i++) opserr << stress(i) << " "; opserr << "]\n";
		opserr << "stress_k   = [ "; for (int i = 0;i < 6;i++) opserr << stress_k(i) << " "; opserr << "]\n";
		opserr << "tangent:\n[ "; for (int i = 0;i < 6;i++) { for (int j = 0;j < 6;j++) opserr << tangent(j, i) << " "; opserr << "]\n"; }
	}

	return 0;
};

//unused
int GDamage::setTrialStrain(const Vector& v, const Vector& r) {

	return this->setTrialStrain(v);
};

//unused
int GDamage::setTrialStrainIncr(const Vector& v) {

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
		opserr << "Fatal:GDamage:: Material dimension is: " << ndm << endln;
		opserr << "But strain vector size is: " << v.Size() << endln;
		exit(-1);
	}


	this->plasticity();
	*/
	return 0;
};

//unused
int GDamage::setTrialStrainIncr(const Vector& v, const Vector& r) {

	return this->setTrialStrainIncr(v);
};

// Calculates current tangent stiffness.
const Matrix& GDamage::getTangent(void) {

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

const Matrix& GDamage::getInitialTangent(void) {

	return tangent_e;
};

const Vector& GDamage::getStress(void) {

	return stress;
};

const Vector& GDamage::getStrain(void) {

	return strain;
};

const Vector& GDamage::getCommittedStress(void) {

	return stress_k;
};

const Vector& GDamage::getCommittedStrain(void) {

	return strain_k;
};

int GDamage::commitState(void) {

	stress_k = stress;
	strain_k = strain;
	Dc_k = Dc;
	Dt_k = Dt;
	D_k = D;
	//CcumPlastStrainDev = cumPlastStrainDev;

	return 0;
};

int GDamage::revertToLastCommit(void) {

	// -- to be implemented.
	return 0;
};

int GDamage::revertToStart(void) {
	// -- to be implemented.
	return 0;
}

int GDamage::sendSelf(int commitTag, Channel& theChannel) {
	// -- to be implemented.

  /*static ID data(1);
  data(0) = tangent;
  return theChannel.sendID(0, commitTag, data);
  */
	return 0;
};

int GDamage::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker) {
	// -- to be implemented.

  /*static ID data(1);
  int res = theChannel.recvID(0, commitTag, data);
  tangent = data(0);
  return res;
  */
	return 0;
};


Response* GDamage::setResponse(const char** argv, int argc, OPS_Stream& s) {


	if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0)
		return new MaterialResponse(this, 1, stress);

	else if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0)
		return new MaterialResponse(this, 2, strain);

	else if (strcmp(argv[0], "tangent") == 0 || strcmp(argv[0], "Tangent") == 0)
		return new MaterialResponse(this, 3, tangent);

	else if (strcmp(argv[0], "damage") == 0 || strcmp(argv[0], "Damage") == 0)
		return new MaterialResponse(this, 4, D);

	else
		return 0;

}



int GDamage::getResponse(int responseID, Information& matInfo) {



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
			*(matInfo.theVector) = D;
		return 0;

	}



	return 0;
};

void GDamage::Print(OPS_Stream& s, int flag) {
	// -- to be implemented.
	return;
};


int GDamage::setParameter(const char** argv, int argc, Parameter& param) {
	// -- to be implemented.


	return 0;
};

int GDamage::updateParameter(int responseID, Information& eleInformation) {
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
