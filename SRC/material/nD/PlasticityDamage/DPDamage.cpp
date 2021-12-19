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

#include <DPDamage.h>

#include <string.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <Information.h>
#include <MaterialResponse.h>
#include <Parameter.h>

static int dFlag1 = 0;	// Turn this on for debug
static int dFlag2 = 0;	// Turn this on for debug on damage subroutine

const double DPDamage::one3 = 1.0 / 3.0;
const double DPDamage::two3 = 2.0 / 3.0;
const double DPDamage::root23 = sqrt(2.0 / 3.0);

#include <elementAPI.h>

static int numDPDamageMaterials = 0;

void* OPS_DPDamage(void)
{
	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 13) {
		opserr << "Want: nDMaterial DPDamage $tag $E $nu $sig_c $sig_t $Hk $Hi\n";
		opserr << "$Yt0, $bt, $at, $Yc0, $bc, $ac, $beta)\n";
		return 0;
	}

	int tag;
	double dData[13];

	int numData = 1;
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "WARNING invalid integer tag: nDMaterial DPDamage \n";
		return 0;
	}

	numData = numArgs - 1;
	if (OPS_GetDouble(&numData, dData) != 0) {
		opserr << "WARNING invalid data: nDMaterial DPDamage : " << tag << "\n";
		return 0;
	}

	NDMaterial* theMaterial = new DPDamage(tag,
		dData[0], dData[1], // E and nu
		dData[2], dData[3], dData[4], dData[5], // Druker Prager plasticity
		dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12]); // Addessi damage

	opserr << "DPDamage memory is allocated!" << endln;
	for (int i = 0;i < 13;i++) opserr << "dData[" << i << "] = " << dData[i] << endln;
	return theMaterial;
}

//full constructor
DPDamage::DPDamage(int tag, double _E, double _nu, // Parameters
	double _sig_c, double _sig_t, double _Hk, double _Hi, // Plasticity
	double _Yt0, double _bt, double _at, double _Yc0, double _bc, double _ac, double _beta) // Damage
    : NDMaterial(tag, ND_TAG_DPDamage),
	E(_E), nu(_nu), sig_c(_sig_c), sig_t(_sig_t), Hk(_Hk), Hi(_Hi),
	Yt0(_Yt0), bt(_bt), at(_at), Yc0(_Yc0), bc(_bc), ac(_ac), beta(_beta),
    mEpsilon(6),
	mEpsilon_n(6),
    mEpsilon_n_p(6),
    mEpsilon_n1_p(6),
	mEpsilon_e(6),
    mSigma(6),
	mSigma_n(6),
    mBeta_n(6),
    mBeta_n1(6),
    mCe(6, 6),
    mCep(6, 6),
	mCt(6, 6),
    mI1(6),
    mIIvol(6, 6),
    mIIdev(6, 6),
    mState(5),
	m(3, 3)
{
	// Bulk and shear modulus
	K = E / (3.0 * (1.0 - 2.0 * nu));
	G = E / (2.0 * (1.0 + nu));

	// Yielding and friction coefficient
	sig_y = 2.0 * sig_c * sig_t / (sig_c + sig_t);
	mu = root23 * (sig_c - sig_t) / (sig_c + sig_t);

	// Total hardening and proportion
	H = Hi + Hk;
	theta = Hk / (Hi + Hk);

	massDen = 0.0;
    mKref = K;
    mGref = G;
    mPatm = 101.0;
    mK = K;
    mG = G;
    msigma_y = sig_y;
    mrho = mu;
    mrho_bar = mu;
    mKinf = K;
    mKo = K;
    mdelta1 = 0.0;
    mdelta2 = 0.0;
    mHard = H;
    mtheta = theta;

    if (mrho == 0.0) {
        mTo = 1e10;
    }
    else {
        mTo = root23 * msigma_y / mrho;
    }
    // set the elastic flag
    //  0 = elastic+no param update; 1 = elastic+param update; 2 = elastoplastic+no param update (default)
    mElastFlag = 2;

    // Use these values to deactivate yield surface 1 - Create Pure Tension Cutoff
    //msigma_y = 1e10;
    //mTo      = 100;

    this->initialize();
}

//null constructor
DPDamage::DPDamage()
    : NDMaterial(),
	E(0.0), nu(0.0), sig_c(0.0), sig_t(0.0), Hk(0.0), Hi(0.0),
	Yt0(0.0), bt(0.0), at(0.0), Yc0(0.0), bc(0.0), ac(0.0), beta(0.0),
    mEpsilon(6),
	mEpsilon_n(6),
    mEpsilon_n_p(6),
    mEpsilon_n1_p(6),
	mEpsilon_e(6),
    mSigma(6),
	mSigma_n(6),
    mBeta_n(6),
    mBeta_n1(6),
    mCe(6, 6),
    mCep(6, 6),
	mCt(6, 6),
    mI1(6),
    mIIvol(6, 6),
    mIIdev(6, 6),
    mState(5),
	m(3,3)
{
	// Bulk and shear modulus
	K = 0.0;
	G = 0.0;

	// Yielding and friction coefficient
	sig_y = 0.0;
	mu = 0.0;

	// Total hardening and proportion
	H = Hi + Hk;
	theta = Hk / (Hi+Hk);

    massDen = 0.0;
    mKref = 0.0;
    mGref = 0.0;
    mPatm = 101.0;
    mK = 0.0;
    mG = 0.0;
    msigma_y = 1e+10;
    mrho = 0.0;
    mrho_bar = 0.0;
    mKinf = 0.0;
    mKo = 0.0;
    mdelta1 = 0.0;
    mdelta2 = 0.0;
    mHard = 0.0;
    mtheta = 0.0;
    mTo = 0.0;

    mElastFlag = 2;

    this->initialize();
}

//destructor
DPDamage::~DPDamage()
{}

//zero internal variables
void DPDamage::initialize()
{
    mEpsilon.Zero();
	mEpsilon_n.Zero();
    mEpsilon_n_p.Zero();
    mEpsilon_n1_p.Zero();
	mEpsilon_e.Zero();

    mSigma.Zero();
	mSigma_n.Zero();
    mBeta_n.Zero();
    mBeta_n1.Zero();

    mAlpha1_n = 0.0;
    mAlpha1_n1 = 0.0;
    mAlpha2_n = 0.0;
    mAlpha2_n1 = 0.0;
    mFlag = 1;

    mHprime = (1 - mtheta) * mHard;

	// Damage variables initialization
	Dt_n = 0.0;
	Dc_n = 0.0;
	D_n = 0.0;
	Dt = Dt_n;
	Dc = Dc_n;
	D = D_n;
	Dm1sq = 1.0;

    // 2nd order Identity Tensor
    mI1.Zero();
    mI1(0) = 1;
    mI1(1) = 1;
    mI1(2) = 1;

    // 4th order Volumetric Tensor
    // IIvol = I1 tensor I1
    mIIvol.Zero();
    mIIvol(0, 0) = 1;
    mIIvol(0, 1) = 1;
    mIIvol(0, 2) = 1;
    mIIvol(1, 0) = 1;
    mIIvol(1, 1) = 1;
    mIIvol(1, 2) = 1;
    mIIvol(2, 0) = 1;
    mIIvol(2, 1) = 1;
    mIIvol(2, 2) = 1;

    // 4th order Deviatoric Tensor
    // Note:  this is the contravariant form!
    //        useable for s^a = 2G * IIdev^ab * epsilon_b
    // (Need a different form for s^a = IIdev ^a_b * sigma^a)
    mIIdev.Zero();
    mIIdev(0, 0) = two3;
    mIIdev(0, 1) = -one3;
    mIIdev(0, 2) = -one3;
    mIIdev(1, 0) = -one3;
    mIIdev(1, 1) = two3;
    mIIdev(1, 2) = -one3;
    mIIdev(2, 0) = -one3;
    mIIdev(2, 1) = -one3;
    mIIdev(2, 2) = two3;
    mIIdev(3, 3) = 0.5;
    mIIdev(4, 4) = 0.5;
    mIIdev(5, 5) = 0.5;

    mCe = mK * mIIvol + 2 * mG * mIIdev;
    mCep = mCe;
	mCt = mCe;
    mState.Zero();

	m.Zero();
}

//make a clone of this material
NDMaterial* DPDamage::getCopy()
{
	DPDamage* clone;
	clone = new DPDamage(this->getTag(), E, nu, sig_c, sig_t, Hk, Hi,
		Yt0, bt, at, Yc0, bc, ac, beta);
	return clone;
}

NDMaterial* DPDamage::getCopy(const char* type)
{
	if (strcmp(type, this->getType()) == 0)
		return this->getCopy();

	return 0;
}

//send back type of material
const char* DPDamage::getType() const
{
	return "ThreeDimensional";
}

//send back order of strain in vector form
int DPDamage::getOrder() const
{
	return 6;
}

//get the strain and integrate plasticity equations
int DPDamage::setTrialStrain(const Vector& strain_from_element)
{
	// Debug

	mEpsilon = strain_from_element;

	// Debug 1
	if (dFlag1 == 1) {
		opserr << "\n----------------------------------------------------------------------------------------------------------------------------\n";
		opserr << "\nStarted new setTrialStrain.\n\n";
		opserr << "Inputs before plasticity and damage (initialized at current step):\n";
		opserr << "strain     = [ "; for (int i = 0;i < 6;i++) opserr << mEpsilon(i) << " "; opserr << "]\n";
		opserr << "strain_k   = [ "; for (int i = 0;i < 6;i++) opserr << mEpsilon_n(i) << " "; opserr << "]\n";
	}

	// Plasticity
	this->plastic_integrator();

	// Updated elastic strains
	mEpsilon_e = mEpsilon - mEpsilon_n1_p;

	// Debug 2
	if (dFlag1 == 1) {
		opserr << "\nPlasticity executed!\n\n";
		opserr << "Outputs after plasticity only:\n";
		opserr << "strain     = [ "; for (int i = 0;i < 6;i++) opserr << mEpsilon(i) << " "; opserr << "]\n";
		opserr << "strain_e   = [ "; for (int i = 0;i < 6;i++) opserr << mEpsilon_e(i) << " "; opserr << "]\n";
		opserr << "strain_p   = [ "; for (int i = 0;i < 6;i++) opserr << mEpsilon_n1_p(i) << " "; opserr << "]\n";
		opserr << "strain_p_k = [ "; for (int i = 0;i < 6;i++) opserr << mEpsilon_n_p(i) << " "; opserr << "]\n";
		opserr << "stress     = [ "; for (int i = 0;i < 6;i++) opserr << mSigma(i) << " "; opserr << "]\n";
		opserr << "stress_k   = [ "; for (int i = 0;i < 6;i++) opserr << mSigma_n(i) << " "; opserr << "]\n";
		opserr << "elastic tangent:\n[ "; for (int i = 0;i < 6;i++) { for (int j = 0;j < 6;j++) opserr << mCep(j, i) << " "; opserr << "]\n"; }
		opserr << "elastoplastic tangent:\n[ "; for (int i = 0;i < 6;i++) { for (int j = 0;j < 6;j++) opserr << mCe(j, i) << " "; opserr << "]\n"; }
	}

	// Damage routine
	for (int i = 3; i < 6; i++) { mEpsilon[i] /= 2.0;mEpsilon_n1_p[i] /= 2.0;mEpsilon_e[i] /= 2.0; }
	this->damage();
	for (int i = 3; i < 6; i++) { mEpsilon[i] *= 2.0;mEpsilon_n1_p[i] *= 2.0;mEpsilon_e[i] *= 2.0; }
	// ------------------------------------------------------------------------------------------------- //

	// Incremental quantities
	//D = 0;	// Damage trigger
	//Vector dstress = stress - stress_k;
	Vector dstrain = mEpsilon - mEpsilon_n;
	Vector dstrain_p = mEpsilon_n1_p - mEpsilon_n_p;
	Vector dstrain_e = dstrain - dstrain_p;
	double dD = D - D_n;

	// Constitutive matrix and stress ---------------------------------------------------------------- //
	mCt = Dm1sq * mCep;
	mSigma = mSigma_n + mCep * dstrain_e - 2 * (1 - D) * dD * mCe * mEpsilon_e;

	// Debug 3
	if (D > 1) {
		opserr << "\nD = " << D << "\n\n";
		opserr << "Outputs after both plasticity and damage:\n";
		opserr << "mEpsilon     = [ "; for (int i = 0;i < 6;i++) opserr << mEpsilon(i) << " "; opserr << "]\n";
		opserr << "mEpsilon_n   = [ "; for (int i = 0;i < 6;i++) opserr << mEpsilon_n(i) << " "; opserr << "]\n";
		opserr << "mEpsilon_e   = [ "; for (int i = 0;i < 6;i++) opserr << mEpsilon_e(i) << " "; opserr << "]\n";
		opserr << "mEpsilon_p   = [ "; for (int i = 0;i < 6;i++) opserr << mEpsilon_n1_p(i) << " "; opserr << "]\n";
		opserr << "mEpsilon_p_n = [ "; for (int i = 0;i < 6;i++) opserr << mEpsilon_n_p(i) << " "; opserr << "]\n";
		opserr << "mSigma       = [ "; for (int i = 0;i < 6;i++) opserr << mSigma(i) << " "; opserr << "]\n";
		opserr << "mSigma_n     = [ "; for (int i = 0;i < 6;i++) opserr << mSigma_n(i) << " "; opserr << "]\n";
		opserr << "tangent_ep:\n[ "; for (int i = 0;i < 6;i++) { for (int j = 0;j < 6;j++) opserr << mCep(j, i) << " "; opserr << "]\n"; }
		opserr << "tangent:\n[ "; for (int i = 0;i < 6;i++) { for (int j = 0;j < 6;j++) opserr << mCt(j, i) << " "; opserr << "]\n"; }
	}

	/* Nota: Ho ottenuto un ottimo risultato con i seguenti (commit interno attivo)
	tangent = pow(1 - D, 2) * tangent; <- qui tangent è elastoplastica
	stress = stress_k + tangent * dstrain_e - 2 * (1 - D) * dD * tangent * strain_e;
	*/

	// Commits
	this->commitState();

	// Debug 3
	if (dFlag1 == 1) {
		opserr << "\nDamage executed!\n\n";
		opserr << "Outputs after both plasticity and damage:\n";
		opserr << "strain     = [ "; for (int i = 0;i < 6;i++) opserr << mEpsilon(i) << " "; opserr << "]\n";
		opserr << "strain_k   = [ "; for (int i = 0;i < 6;i++) opserr << mEpsilon_n(i) << " "; opserr << "]\n";
		opserr << "strain_e   = [ "; for (int i = 0;i < 6;i++) opserr << mEpsilon_e(i) << " "; opserr << "]\n";
		opserr << "strain_p   = [ "; for (int i = 0;i < 6;i++) opserr << mEpsilon_n_p(i) << " "; opserr << "]\n";
		opserr << "strain_p_k = [ "; for (int i = 0;i < 6;i++) opserr << mEpsilon_n1_p(i) << " "; opserr << "]\n";
		opserr << "stress     = [ "; for (int i = 0;i < 6;i++) opserr << mSigma(i) << " "; opserr << "]\n";
		opserr << "stress_k   = [ "; for (int i = 0;i < 6;i++) opserr << mSigma_n(i) << " "; opserr << "]\n";
		opserr << "tangent_ep:\n[ "; for (int i = 0;i < 6;i++) { for (int j = 0;j < 6;j++) opserr << mCep(j, i) << " "; opserr << "]\n"; }
		opserr << "tangent:\n[ "; for (int i = 0;i < 6;i++) { for (int j = 0;j < 6;j++) opserr << mCe(j, i) << " "; opserr << "]\n"; }
	}

	return 0;
}

//unused trial strain functions
int DPDamage::setTrialStrain(const Vector& v, const Vector& r)
{
	return this->setTrialStrain(v);
}

//plasticity integration routine
void DPDamage::plastic_integrator()
{
	bool okay;		// boolean variable to ensure satisfaction of multisurface kuhn tucker conditions
	double f1;
	double f2;
	double norm_eta;
	double Invariant_1;
	double Invariant_ep;
	double norm_ep;
	double norm_dev_ep;
	//Vector epsilon_e(6); --> Moved to memory
	Vector s(6);
	Vector eta(6);
	Vector dev_ep(6);
	Vector Jact(2);

	double NormCep;

	double alpha1;			// hardening parameter for DP surface
	double alpha2;			// hardening parameter for tension cut-off
	Vector n(6);			// normal to the yield surface in strain space
	Vector R(2);			// residual vector
	Vector gamma(2);		// vector of consistency parameters
	Vector dgamma(2);		// incremental vector of consistency parameters
	Matrix g(2, 2);			// jacobian of the corner region (return map)
	Matrix g_contra(2, 2);	// inverse of jacobian of the corner region

	// Accuracy control
	double tol = 1.0e-6;
	int mMax = 2;
	int countMax = 100;
	double fTOL;
	double gTOL;
	fTOL = 0.0;
	gTOL = -1.0e-6;

	// set trial state:

	// epsilon_n1_p_trial = ..._n1_p  = ..._n_p
	mEpsilon_n1_p = mEpsilon_n_p;

	// alpha1_n+1_trial
	mAlpha1_n1 = mAlpha1_n;
	// alpha2_n+1_trial
	mAlpha2_n1 = mAlpha2_n;

	// beta_n+1_trial
	mBeta_n1 = mBeta_n;

	// epsilon_elastic = epsilon_n+1 - epsilon_n_p
	mEpsilon_e = mEpsilon - mEpsilon_n1_p;

	// trial stress
	mSigma = mCe * mEpsilon_e;

	// deviator stress tensor: s = 2G * IIdev * epsilon_e
	//I1_trial
	Invariant_1 = (mSigma(0) + mSigma(1) + mSigma(2));

	// s_n+1_trial
	s = mSigma - (Invariant_1 / 3.0) * mI1;

	//eta_trial = s_n+1_trial - beta_n;
	eta = s - mBeta_n;

	// compute yield function value (contravariant norm)
	norm_eta = sqrt(eta(0) * eta(0) + eta(1) * eta(1) + eta(2) * eta(2) + 2 * (eta(3) * eta(3) + eta(4) * eta(4) + eta(5) * eta(5)));

	// f1_n+1_trial
	f1 = norm_eta + mrho * Invariant_1 - root23 * Kiso(mAlpha1_n1);

	// f2_n+1_trial
	f2 = Invariant_1 - T(mAlpha2_n1);

	// update elastic bulk and shear moduli 
	this->updateElasticParam();

	// check trial state
	int count = 1;
	if ((f1 <= fTOL) && (f2 <= fTOL) || mElastFlag < 2) {

		okay = true;
		// trial state = elastic state - don't need to do any updates.
		mCep = mCe;
		count = 0;

		// set state variables for recorders
		Invariant_ep = mEpsilon_n1_p(0) + mEpsilon_n1_p(1) + mEpsilon_n1_p(2);

		norm_ep = sqrt(mEpsilon_n1_p(0) * mEpsilon_n1_p(0) + mEpsilon_n1_p(1) * mEpsilon_n1_p(1) + mEpsilon_n1_p(2) * mEpsilon_n1_p(2)
			+ 0.5 * (mEpsilon_n1_p(3) * mEpsilon_n1_p(3) + mEpsilon_n1_p(4) * mEpsilon_n1_p(4) + mEpsilon_n1_p(5) * mEpsilon_n1_p(5)));

		dev_ep = mEpsilon_n1_p - one3 * Invariant_ep * mI1;

		norm_dev_ep = sqrt(dev_ep(0) * dev_ep(0) + dev_ep(1) * dev_ep(1) + dev_ep(2) * dev_ep(2)
			+ 0.5 * (dev_ep(3) * dev_ep(3) + dev_ep(4) * dev_ep(4) + dev_ep(5) * dev_ep(5)));

		mState(0) = Invariant_1;
		mState(1) = norm_eta;
		mState(2) = Invariant_ep;
		mState(3) = norm_dev_ep;
		mState(4) = norm_ep;
		return;
	}
	else {
		// plastic correction required
		okay = false;

		// determine number of active surfaces.  size & fill Jact
		if ((f1 > fTOL) && (f2 <= fTOL)) {
			// f1 surface only
			Jact(0) = 1;
			Jact(1) = 0;
		}
		else if ((f1 <= fTOL) && (f2 > fTOL)) {
			// f2 surface only
			Jact(0) = 0;
			Jact(1) = 1;
		}
		else if ((f1 > fTOL) && (f2 > fTOL)) {
			// both surfaces active
			Jact(0) = 1;
			Jact(1) = 1;
		}
	}

	//-----------------MultiSurface Placity Return Map--------------------------------------
	while (!okay && count < countMax+1) {

		alpha1 = mAlpha1_n;
		alpha2 = mAlpha2_n;

		//  n = eta / norm_eta;  (contravaraint)
		if (norm_eta < 1.0e-13) {
			n.Zero();
		}
		else {
			n = eta / norm_eta;
		}

		// initialize R, gamma1, gamma2, dgamma1, dgamma2 = 0
		R.Zero();
		gamma.Zero();
		dgamma.Zero();
		// initialize g such that det(g) = 1
		g(0, 0) = 1;
		g(1, 1) = 1;
		g(1, 0) = 0;
		g(0, 1) = 0;

		// Newton procedure to compute nonlinear gamma1 and gamma2
		//initialize terms
		if (Jact(0) == 1) {
			R(0) = norm_eta - (2 * mG + two3 * mHprime) * gamma(0) + mrho * Invariant_1
				- 9 * mK * mrho * mrho_bar * gamma(0) - 9 * mK * mrho * gamma(1) - root23 * Kiso(alpha1);
			g(0, 0) = -2 * mG - two3 * (mHprime + Kisoprime(alpha1)) - 9 * mK * mrho * mrho_bar;
		}
		else if (Jact(1) == 1) {
			R(1) = Invariant_1 - 9 * mK * mrho_bar * gamma(0) - 9 * mK * gamma(1) - T(alpha2);
			g(1, 1) = -9 * mK + mdelta2 * T(alpha2);
		}
		if (Jact(0) == 1 && Jact(1) == 1) {
			g(0, 1) = -9 * mK * mrho;
			g(1, 0) = mrho_bar * (-9 * mK + mdelta2 * T(alpha2));
		}
		g.Invert(g_contra);

		// iteration counter
		int m = 0;

		//iterate
		while ((fabs(R.Norm()) > tol) && (m < mMax)) {

			dgamma = -1 * g_contra * R;
			gamma += dgamma;

			//update alpha1 and alpha2
			alpha1 = mAlpha1_n + root23 * gamma(0);
			alpha2 = mAlpha2_n + mrho_bar * gamma(0) + gamma(1);

			// reset g & R matrices
			g(0, 0) = 1;
			g(1, 1) = 1;
			g(1, 0) = 0;
			g(0, 1) = 0;
			R.Zero();
			if (Jact(0) == 1) {
				R(0) = norm_eta - (2 * mG + two3 * mHprime) * gamma(0) + mrho * Invariant_1
					- 9 * mK * mrho * mrho_bar * gamma(0) - 9 * mK * mrho * gamma(1) - root23 * Kiso(alpha1);
				g(0, 0) = -2 * mG - two3 * (mHprime + Kisoprime(alpha1)) - 9 * mK * mrho * mrho_bar;
			}
			else if (Jact(1) == 1) {
				R(1) = Invariant_1 - 9 * mK * mrho_bar * gamma(0) - 9 * mK * gamma(1) - T(alpha2);
				g(1, 1) = -9 * mK + mdelta2 * T(alpha2);
			}
			if (Jact(0) == 1 && Jact(1) == 1) {
				g(0, 1) = -9 * mK * mrho;
				g(1, 0) = mrho_bar * (-9 * mK + mdelta2 * T(alpha2));
			}
			g.Invert(g_contra);

			m++;
		}

		// check maintain Kuhn-Tucker conditions
		f1 = norm_eta - (2 * mG + two3 * mHprime) * gamma(0) + mrho * Invariant_1
			- 9 * mK * mrho * mrho_bar * gamma(0) - 9 * mK * mrho * gamma(1) - root23 * Kiso(alpha1);
		f2 = Invariant_1 - 9 * mK * mrho_bar * gamma(0) - 9 * mK * gamma(1) - T(alpha2);
		
		if (count > countMax) {
			okay = true;
			break;
		}
		
		// check active surfaces
		if ((Jact(0) == 1) && (Jact(1) == 0)) {
			// f2 may be > or < f2_tr because of softening of f2 related to alpha1
			if (f2 >= fTOL) {
				// okay = false;
				Jact(0) = 1;
				Jact(1) = 1;
				count += 1;

			}
			else {
				okay = true;

			}
		}
		else if ((Jact(0) == 0) && (Jact(1) == 1)) {
			// f1 will always be less than f1_tr
			okay = true;
		}
		else if ((Jact(0) == 1) && (Jact(1) == 1)) {
			if ((gamma(0) <= gTOL) && (gamma(1) > gTOL)) {
				// okay = false;
				Jact(0) = 0;
				Jact(1) = 1;
				count += 1;
			}
			else if ((gamma(0) > gTOL) && (gamma(1) <= gTOL)) {
				// okay = false;
				Jact(0) = 1;
				Jact(1) = 0;
				count += 1;
			}
			else if ((gamma(0) > gTOL) && (gamma(1) > gTOL)) {
				okay = true;
			}
		}

		if ((count > 3) && (!okay)) {
			Jact(0) = 1;
			Jact(1) = 1;
			count += 100;
		}

		if (count > 3) {
			opserr << "Jact = " << Jact;
			opserr << "count = " << count << endln;
		}

	} // end of while(!okay) loop


	//update everything and exit!

	Vector b1(6);
	Vector b2(6);
	Vector n_covar(6);
	Vector temp1(6);
	Vector temp2(6);

	// update alpha1 and alpha2
	mAlpha1_n1 = alpha1;
	mAlpha2_n1 = alpha2;

	//update epsilon_n1_p
	//first calculate n_covar
	// n_a = G_ab * n^b = covariant
	n_covar(0) = n(0);
	n_covar(1) = n(1);
	n_covar(2) = n(2);
	n_covar(3) = 2 * n(3);
	n_covar(4) = 2 * n(4);
	n_covar(5) = 2 * n(5);
	mEpsilon_n1_p = mEpsilon_n_p + (mrho_bar * gamma(0) + gamma(1)) * mI1 + gamma(0) * n_covar;


	Invariant_ep = mEpsilon_n1_p(0) + mEpsilon_n1_p(1) + mEpsilon_n1_p(2);

	norm_ep = sqrt(mEpsilon_n1_p(0) * mEpsilon_n1_p(0) + mEpsilon_n1_p(1) * mEpsilon_n1_p(1) + mEpsilon_n1_p(2) * mEpsilon_n1_p(2)
		+ 0.5 * (mEpsilon_n1_p(3) * mEpsilon_n1_p(3) + mEpsilon_n1_p(4) * mEpsilon_n1_p(4) + mEpsilon_n1_p(5) * mEpsilon_n1_p(5)));

	dev_ep = mEpsilon_n1_p - one3 * Invariant_ep * mI1;

	norm_dev_ep = sqrt(dev_ep(0) * dev_ep(0) + dev_ep(1) * dev_ep(1) + dev_ep(2) * dev_ep(2)
		+ 0.5 * (dev_ep(3) * dev_ep(3) + dev_ep(4) * dev_ep(4) + dev_ep(5) * dev_ep(5)));

	// update sigma
	mSigma -= (3 * mK * mrho_bar * gamma(0) + 3 * mK * gamma(1)) * mI1 + 2 * mG * gamma(0) * n;

	s -= 2 * mG * gamma(0) * n;
	Invariant_1 -= 9 * mK * mrho_bar * gamma(0) + 9 * mK * gamma(1);
	//mSigma        = s + Invariant_1/3.0 * mI1;

	//update beta_n1
	mBeta_n1 = mBeta_n - (two3 * mHprime * gamma(0)) * n;

	//eta_n+1 = s_n+1 - beta_n+1;
	eta = s - mBeta_n1;
	norm_eta = sqrt(eta(0) * eta(0) + eta(1) * eta(1) + eta(2) * eta(2) + 2 * (eta(3) * eta(3) + eta(4) * eta(4) + eta(5) * eta(5)));

	// update Cep
	// note: Cep is contravariant
	if ((Jact(0) == 1) && (Jact(1) == 0)) {
		b1 = 2 * mG * n + 3 * mK * mrho * mI1;
		b2.Zero();
	}
	else if ((Jact(0) == 0) && (Jact(1) == 1)) {
		b1.Zero();
		b2 = 3 * mK * mI1;
	}
	else if ((Jact(0) == 1) && (Jact(1) == 1)) {
		b1 = 2 * mG * n + 3 * mK * mrho * mI1;
		b2 = 3 * mK * mI1;
	}

	temp1 = g_contra(0, 0) * b1 + g_contra(0, 1) * b2;
	temp2 = mrho_bar * temp1 + g_contra(1, 0) * b1 + g_contra(1, 1) * b2;

	NormCep = 0.0;
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			mCep(i, j) = mCe(i, j)
				+ 3 * mK * mI1(i) * temp2(j)
				+ 2 * mG * n(i) * temp1(j)
				- 4 * mG * mG / norm_eta * gamma(0) * (mIIdev(i, j) - n(i) * n(j));
			NormCep += mCep(i, j) * mCep(i, j);
		}
	}

	if (NormCep < tol) {
		mCep = 1.0e-3 * mCe;
		opserr << "NormCep = " << NormCep << endln;
	}

	mState(0) = Invariant_1;
	mState(1) = norm_eta;
	mState(2) = Invariant_ep;
	mState(3) = norm_dev_ep;
	mState(4) = norm_ep;

	return;
}

// Damage integration routine
void DPDamage::damage()
{
	//////// 1. Principal strains calculation /////////////////////////////////////////////////////////////////
	// Strain tensors from strain vectors (sci)
	Matrix epsilon(3, 3);	epsilon = tens(mEpsilon);
	Matrix epsilon_e(3, 3);	epsilon_e = tens(mEpsilon_e);

	// Principal strains
	Vector strain_m(3);
	Vector strain_e_m(3);
	Matrix eVect(3, 3);		// Eigenvectors
	eVect.Eigen3(epsilon); for (int i = 0; i < 3;i++) strain_m(i) = eVect(i, i);
	eVect.Eigen3(epsilon_e); for (int i = 0; i < 3;i++) strain_e_m(i) = eVect(i, i);

	//////// 2. Equivalent strains calculation ////////////////////////////////////////////////////////////////
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

	////////// 3. Damage coefficients evaluation part ///////////////////////////////////////////////////////////////////////
	// Damage functions at step n
	double Ft = (Yt - Yt0) - Dt_n * (at * Yt + bt);
	double Fc = (Yc - Yc0) - Dc_n * (ac * Yc + bc);

	// Tensile damage parameter at step n+1
	Dt = (Yt - Yt0) / (at * Yt + bt);
	if (Ft > 0.0)	Dt = fmax(Dt, Dt_n);
	else			Dt = Dt_n;
	if (Dt < 0.0)	Dt = 0.0;
	if (Dt > 1.0)	Dt = 1.0;

	// Compressive damage parameter at step n+1
	Dc = (Yc - Yc0) / (ac * Yc + bc);
	if (Fc > 0.0)	Dc = fmax(Dc, Dc_n);
	else			Dc = Dc_n;
	if (Dc < 0.0)	Dc = 0.0;
	if (Dc > 1.0)	Dc = 1.0;

	// Condition for tensile damage Dt > Dc ---> Dt = max(Dc,Dt)
	Dt = fmax(Dc, Dt);

	// Weighting coefficients
	double alpha;
	if (Yt_e == 0.0 && Yc_e == 0.0) alpha = 0.0;
	else alpha = (fabs(Yt_e) / Yt0) / (fabs(Yt_e) / Yt0 + fabs(Yc_e) / Yc0);
	if (alpha < 0.0) alpha = 0.0;
	if (alpha > 1.0) alpha = 1.0;

	if (fabs(Yt_e)/Yt0 + fabs(Yc_e)/Yc0 < 1.0e-10) D = D_n;

	// Cumulative damage variable Dm1sq = [1-D]^2
	D = alpha * Dt + (1.0 - alpha) * Dc;
	Dm1sq = pow(1.0 - D, 2.0);

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

int DPDamage::updateElasticParam()
{
	double Sigma_mean = 0.0;
	if (mElastFlag == 1 && mFlag == 1) {
		Sigma_mean = -one3 * (mSigma(0) + mSigma(1) + mSigma(2));
		if (Sigma_mean < 0.0) Sigma_mean = 0.0;  // prevents modulus update for cases where tension exists 
		mK = mKref * pow(1 + (Sigma_mean / mPatm), 0.5);
		mG = mGref * pow(1 + (Sigma_mean / mPatm), 0.5);
		mCe = mK * mIIvol + 2 * mG * mIIdev;
		mFlag = 0;
		//opserr << "Plastic Integrator -->" << "K = " << mK  << "  G =" << mG << endln;
	}
	else if (mElastFlag != 1) {
		mFlag = 1;
	}

	return 0;
}

double DPDamage::Kiso(double alpha1)
{
	return msigma_y + mtheta * mHard * alpha1;// +(mKinf - mKo) * (1 - exp(-mdelta1 * alpha1));
}

double DPDamage::Kisoprime(double alpha1)
{
	return mtheta * mHard;// +(mKinf - mKo) * mdelta1 * exp(-mdelta1 * alpha1);
}

double DPDamage::T(double alpha2)
{
	return mTo * exp(-mdelta2 * alpha2);
}

double DPDamage::deltaH(double dGamma)
{
	return mHprime * root23 * dGamma;
}

//vector to tensor
const Matrix& DPDamage::tens(const Vector& v)
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
const Vector& DPDamage::getStrain()
{
	return mEpsilon;
}

//send back the stress 
const Vector& DPDamage::getStress()
{
	return mSigma;
}

//send back the tangent 
const Matrix& DPDamage::getTangent()
{
	return mCt;
}

//send back the tangent 
const Matrix& DPDamage::getInitialTangent()
{
	return mCe;
}

int DPDamage::commitState(void)
{
	mEpsilon_n_p = mEpsilon_n1_p;
	mAlpha1_n = mAlpha1_n1;
	mAlpha2_n = mAlpha2_n1;
	mBeta_n = mBeta_n1;

	mSigma_n = mSigma;
	mEpsilon_n = mEpsilon;
	Dc_n = Dc;
	Dt_n = Dt;
	D_n = D;

	return 0;
}

int DPDamage::revertToLastCommit(void)
{
	return 0;
}

int DPDamage::revertToStart(void)
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

Vector DPDamage::getState()
{
	return mState;
}

Response*
DPDamage::setResponse(const char** argv, int argc, OPS_Stream& output)
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

int DPDamage::getResponse(int responseID, Information& matInfo)
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

int DPDamage::setParameter(const char** argv, int argc, Parameter& param)
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
		opserr << "WARNING: invalid parameter command for DPDamage nDMaterial with tag: " << this->getTag() << endln;
		return -1;
	}

	return -1;
}

int DPDamage::updateParameter(int responseID, Information& info)
{
	if (responseID == 5) {
		// materialState called - update mElasticFlag
		mElastFlag = (int)info.theDouble;
	}
	else if (responseID == 7) {
		// frictionalStrength called - update rho and tension cutoff
		mrho = info.theDouble;
		if (mrho == 0.0) {
			mTo = 1e10;
		}
		else {
			mTo = root23 * msigma_y / mrho;
		}
	}
	else if (responseID == 8) {
		// nonassociativeTerm called - update rho_bar
		mrho_bar = info.theDouble;
	}
	else if (responseID == 9) {
		// cohesiveIntercept called - update sigma_y and tension cutoff
		msigma_y = info.theDouble;
		if (mrho == 0.0) {
			mTo = 1e10;
		}
		else {
			mTo = root23 * msigma_y / mrho;
		}
	}
	else if (responseID == 10) {
		// shearModulus called - update G and Ce
		mG = info.theDouble;
		mCe = mK * mIIvol + 2 * mG * mIIdev;
	}
	else if (responseID == 11) {
		// bulkModulus called - update K and Ce
		mK = info.theDouble;
		mCe = mK * mIIvol + 2 * mG * mIIdev;
	}

	return 0;
}

int DPDamage::sendSelf(int commitTag, Channel& theChannel)
{
	int res = 0;

	// place data in a vector
	static Vector data(45);
	data(0) = this->getTag();
	data(1) = mKref;
	data(2) = mGref;
	data(3) = mK;
	data(4) = mG;
	data(5) = msigma_y;
	data(6) = mrho;
	data(7) = mrho_bar;
	data(8) = mKinf;
	data(9) = mKo;
	data(10) = mdelta1;
	data(11) = mdelta2;
	data(12) = mHard;
	data(13) = mtheta;
	data(14) = massDen;
	data(15) = mPatm;
	data(16) = mTo;
	data(17) = mHprime;
	data(18) = mAlpha1_n;
	data(19) = mAlpha2_n;
	data(20) = mElastFlag;
	data(21) = mFlag;

	data(22) = mEpsilon(0);
	data(23) = mEpsilon(1);
	data(24) = mEpsilon(2);
	data(25) = mEpsilon(3);
	data(26) = mEpsilon(4);
	data(27) = mEpsilon(5);

	data(28) = mEpsilon_n_p(0);
	data(29) = mEpsilon_n_p(1);
	data(30) = mEpsilon_n_p(2);
	data(31) = mEpsilon_n_p(3);
	data(32) = mEpsilon_n_p(4);
	data(33) = mEpsilon_n_p(5);

	data(34) = mBeta_n(0);
	data(35) = mBeta_n(1);
	data(36) = mBeta_n(2);
	data(37) = mBeta_n(3);
	data(38) = mBeta_n(4);
	data(39) = mBeta_n(5);

	data(40) = mState(0);
	data(41) = mState(1);
	data(42) = mState(2);
	data(43) = mState(3);
	data(44) = mState(4);

	res = theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "WARNING: DPDamage::sendSelf - failed to send vector to channel" << endln;
		return -1;
	}

	return 0;
}

int DPDamage::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	int res = 0;

	// receive data
	static Vector data(45);
	res = theChannel.recvVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "WARNING: DPDamage::recvSelf - failed to receive vector from channel" << endln;
		return -1;
	}

	// set member variables
	this->setTag((int)data(0));
	mKref = data(1);
	mGref = data(2);
	mK = data(3);
	mG = data(4);
	msigma_y = data(5);
	mrho = data(6);
	mrho_bar = data(7);
	mKinf = data(8);
	mKo = data(9);
	mdelta1 = data(10);
	mdelta2 = data(11);
	mHard = data(12);
	mtheta = data(13);
	massDen = data(14);
	mPatm = data(15);
	mTo = data(16);
	mHprime = data(17);
	mAlpha1_n = data(18);
	mAlpha2_n = data(19);
	mElastFlag = (int)data(20);
	mFlag = (int)data(21);

	mEpsilon(0) = data(22);
	mEpsilon(1) = data(23);
	mEpsilon(2) = data(24);
	mEpsilon(3) = data(25);
	mEpsilon(4) = data(26);
	mEpsilon(5) = data(27);

	mEpsilon_n_p(0) = data(28);
	mEpsilon_n_p(1) = data(29);
	mEpsilon_n_p(2) = data(30);
	mEpsilon_n_p(3) = data(31);
	mEpsilon_n_p(4) = data(32);
	mEpsilon_n_p(5) = data(33);

	mBeta_n(0) = data(34);
	mBeta_n(1) = data(35);
	mBeta_n(2) = data(36);
	mBeta_n(3) = data(37);
	mBeta_n(4) = data(38);
	mBeta_n(5) = data(39);

	mState(0) = data(40);
	mState(1) = data(41);
	mState(2) = data(42);
	mState(3) = data(43);
	mState(4) = data(44);

	mCe = mK * mIIvol + 2 * mG * mIIdev;
	mCep = mCe;
	mCt = mCe;

	return 0;
}

void DPDamage::Print(OPS_Stream& s, int flag)
{
	s << "DPDamage" << endln;
}

