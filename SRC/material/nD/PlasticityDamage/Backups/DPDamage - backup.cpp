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
#include <Channel.h>
#include <cmath>
#include <elementAPI.h>
#include <Matrix.h>
#include <Matrix.cpp>

const double DPDamage::one3 = 1.0 / 3.0;
const double DPDamage::two3 = 2.0 / 3.0;
const double DPDamage::root23 = sqrt(2.0 / 3.0);

// Starting kudos message
static int numPD2P = 0;

// Output debug flag
int df = 1;

void *
OPS_DPDamage(void)
{
	if (numPD2P == 0) { numPD2P++;
		opserr << "Two parameters damage and plasticity 3D constitutive law based on Addessi et al.\nBy LP and DF.\n";
	}

	NDMaterial* theMaterial = 0;

	int numArgs = OPS_GetNumRemainingInputArgs();

	if (numArgs < 14) {
		opserr << "Want: nDMaterial DPDamage $tag $E $nu $ft $fc $Hk $Hi\n";
		opserr << "$Yt0, $bt, $at, $Yc0, $bc, $ac, $beta)\n";
		return 0;
	}

	int tag;
	double dData[13];
	// Consiglia dati  così:  dData[4] = 0.6;

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

	theMaterial = new DPDamage(tag,
		dData[0], dData[1], // E and nu
		dData[2], dData[3], dData[4], dData[5], // Druger Prager plasticity
		dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12]); // Addessi damage

	opserr << "Tcl script read and material created succesfully.\n";
	return theMaterial;
}

// Full constructor - Initial material terms determination
DPDamage::DPDamage(int tag, double _E, double _nu, // Parameters
	double _ft, double _fc, double _Hk, double _Hi, // Plasticity
	double _Yt0, double _bt, double _at, double _Yc0, double _bc, double _ac, double _beta) // Damage
	:NDMaterial(tag, ND_TAG_DPDamage),
	E(_E), nu(_nu),	ft(_ft), fc(_fc), Hk(_Hk), Hi(_Hi),
	Yt0(_Yt0), bt(_bt), at(_at), Yc0(_Yc0), bc(_bc), ac(_ac), beta(_beta),
	Eps(6),	Eps_n_p(6),	Eps_n1_p(6), Sig(6),
	Zeta_n(6), Zeta_n1(6),
	Ce(6, 6), Cep(6, 6),
	I1(6), IIvol(6, 6),	IIdev(6, 6),
	State(5)
{
	// Bulk and shear modulus
	K = E / (3 * (1 - 2 * nu));
	G = E / (2 * (1 + nu));

	// Mohr Coulomb plasticity parameters
	sigma_y = 2 * fc * ft / (fc + ft);
	mu = root23 * (fc - ft) / (fc + ft);

	// Condensed hardening parameters
	H = Hi + Hk;
	theta = Hk / (Hi + Hk);

	// Set the elastic flag
	// 0 = elastic+no param update; 1 = elastic+param update; 2 = elastoplastic+no param update (default)
	ElastFlag = 2;
	
	this->initialize();

	opserr << "Full constructor executed!\n";
}


// Null constructor - Empty terms determination
DPDamage::DPDamage()
  :NDMaterial (0, ND_TAG_DPDamage),
	E(0.0), nu(0.0), ft(0.0), fc(0.0), Hk(0.0), Hi(0.0),
	Yt0(0.0), bt(0.0), at(0.0), Yc0(0.0), bc(0.0), ac(0.0), beta(0.0),
	Eps(6),	Eps_n_p(6),	Eps_n1_p(6),
	Sig(6),	Zeta_n(6),	Zeta_n1(6),
	Ce(6, 6),Cep(6, 6),
	I1(6),IIvol(6, 6),IIdev(6, 6),
	State(5)
{
	E = 0.0;
	nu = 0.0;
	K = 0.0;
	G = 0.0;
	sigma_y = 1e+10;
	mu = 0.0;

	ElastFlag = 2;

	this->initialize();
}

// Destructor
DPDamage::~DPDamage ()
{}

// Zero internal variables
void DPDamage::initialize()
{
	// Strain vectors
	Eps.Zero();
	Eps_n_p.Zero();
	Eps_n1_p.Zero();

	// Stress and backstress vectors
	Sig.Zero();
	Zeta_n.Zero();
	Zeta_n1.Zero();

	// Backstress variables
	Alpha1_n = 0.0;
	Alpha1_n1 = 0.0;
	Alpha2_n = 0.0;
	Alpha2_n1 = 0.0;
	Flag = 1;

	// 2nd order Identity Tensor
	I1.Zero();
	I1(0) = 1;
	I1(1) = 1;
	I1(2) = 1;

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

	// Elastic and elastoplastic constitutive matrix
	Ce = K * IIvol + 2 * G * IIdev;
	Cep = Ce;
	State.Zero();

	// Damage variables initialization
	Dt_n = 0;
	Dc_n = 0;
}

int
DPDamage::setTrialStrain(Vector const& v1, Vector const& v2) {
	return this->setTrialStrain(v1);
}

int
DPDamage::setTrialStrain(const Vector& strain)
{
	// Strains vector from global analysis
	Eps = strain;

	if (df == 1) {
		opserr << "--------------------------------------------------------------------------\n";
		opserr << "\nTotal strains:     e = [ ";
		for (int i = 0;i < 6;i++) opserr << Eps(i) << " ";
		opserr << "]\n";
	}

	/* PLASTICITY -------------------------------------- */
	this->plastic_integrator();

	// Elastic strains
	Vector Eps_e(6);
	Eps_e = Eps - Eps_n1_p;

	if (df == 1) {
		opserr << "\nPlastic routine executed succesfully:";
		opserr << "\nTotal strains: e = ";
		for (int i = 0;i < 6;i++) opserr << Eps(i) << " ";
		opserr << "\nElastic strains: ee = ";
		for (int i = 0;i < 6;i++) opserr << Eps_e(i) << " ";
		opserr << "\nPlastic strains: ep = ";
		for (int i = 0;i < 6;i++) opserr << Eps_n1_p(i) << " ";
		opserr << "\nStresses: sig = ";
		for (int i = 0;i < 6;i++) opserr << Sig(i) << " ";
		opserr << "\n";
	}

	/* DAMAGE -------------------------------------------*/
	this->damage_integrator();

	
	return 0;
}

// Plasticity integration routine
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
	Vector Eps_e(6);
	Vector s(6);
	Vector eta(6);
	Vector dev_ep(6);
	Vector Jact(2);

	double fTOL;
	double gTOL;
	fTOL = 0.0;
	gTOL = -1.0e-10;

	double NorCep;

	double alpha1;			// hardening parameter for DP surface
	double alpha2;			// hardening parameter for tension cut-off
	Vector n(6);			// normal to the yield surface in strain space
	Vector R(2);			// residual vector
	Vector gamma(2);		// vector of consistency parameters
	Vector dgamma(2);		// incremental vector of consistency parameters
	Matrix g(2, 2);			// jacobian of the corner region (return map)
	Matrix g_contra(2, 2);	// inverse of jacobian of the corner region

	// Calcolo variabili di storia
	Eps_n1_p = Eps_n_p;	// eps^p_n+1_trial
	Alpha1_n1 = Alpha1_n;			// alpha1_n+1_trial
	Alpha2_n1 = Alpha2_n;			// alpha2_n+1_trial
	Zeta_n1 = Zeta_n;				// beta_n+1_trial

	// Deformazioni elastiche
	Eps_e = Eps - Eps_n1_p;

	// Tensioni elastiche trial
	Sig = Ce * Eps_e;

	// deviator stress tensor: s = 2G * IIdev * Eps_e
	//I1_trial
	Invariant_1 = (Sig(0) + Sig(1) + Sig(2));

	// s_n+1_trial
	s = Sig - (Invariant_1 / 3.0) * I1;

	//eta_trial = s_n+1_trial - beta_n;
	eta = s - Zeta_n;

	// compute yield function value (contravariant norm)
	norm_eta = sqrt(eta(0) * eta(0) + eta(1) * eta(1) + eta(2) * eta(2) + 2 * (eta(3) * eta(3) + eta(4) * eta(4) + eta(5) * eta(5)));

	// f1_n+1_trial
	f1 = norm_eta + mu * Invariant_1 - root23 * qiso(Alpha1_n1);

	// f2_n+1_trial
	f2 = Invariant_1;

	
	// check trial state
	int count = 1;
	if ((f1 <= fTOL) && (f2 <= fTOL) || ElastFlag < 2) {

		okay = true;
		// trial state = elastic state - don't need to do any updates.
		Cep = Ce;
		count = 0;

		// set state variables for recorders
		Invariant_ep = Eps_n1_p(0) + Eps_n1_p(1) + Eps_n1_p(2);

		norm_ep = sqrt(Eps_n1_p(0) * Eps_n1_p(0) + Eps_n1_p(1) * Eps_n1_p(1) + Eps_n1_p(2) * Eps_n1_p(2)
			+ 0.5 * (Eps_n1_p(3) * Eps_n1_p(3) + Eps_n1_p(4) * Eps_n1_p(4) + Eps_n1_p(5) * Eps_n1_p(5)));

		dev_ep = Eps_n1_p - one3 * Invariant_ep * I1;

		norm_dev_ep = sqrt(dev_ep(0) * dev_ep(0) + dev_ep(1) * dev_ep(1) + dev_ep(2) * dev_ep(2)
			+ 0.5 * (dev_ep(3) * dev_ep(3) + dev_ep(4) * dev_ep(4) + dev_ep(5) * dev_ep(5)));

		State(0) = Invariant_1;
		State(1) = norm_eta;
		State(2) = Invariant_ep;
		State(3) = norm_dev_ep;
		State(4) = norm_ep;
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
	while (!okay) {

		alpha1 = Alpha1_n;
		alpha2 = Alpha2_n;

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
			R(0) = norm_eta - (2 * G + two3 * Hk) * gamma(0) + mu * Invariant_1
				- 9 * K * mu * mu * gamma(0) - 9 * K * mu * gamma(1) - root23 * qiso(alpha1);
			g(0, 0) = -2 * G - two3 * (Hk + Hi) - 9 * K * mu * mu;
		}
		if (Jact(0) == 1 && Jact(1) == 1) {
			g(0, 1) = -9 * K * mu;
			g(1, 0) = mu * (-9 * K);
		}
		g.Invert(g_contra);

		// iteration counter
		int m = 0;

		//iterate
		while ((fabs(R.Norm()) > 1e-10) && (m < 10)) {

			dgamma = -1 * g_contra * R;
			gamma += dgamma;

			//update alpha1 and alpha2
			alpha1 = Alpha1_n + root23 * gamma(0);
			alpha2 = Alpha2_n + mu * gamma(0) + gamma(1);

			// reset g & R matrices
			g(0, 0) = 1;
			g(1, 1) = 1;
			g(1, 0) = 0;
			g(0, 1) = 0;
			R.Zero();
			if (Jact(0) == 1) {
				R(0) = norm_eta - (2 * G + two3 * Hk) * gamma(0) + mu * Invariant_1
					- 9 * K * mu * mu * gamma(0) - 9 * K * mu * gamma(1) - root23 * qiso(alpha1);
				g(0, 0) = -2 * G - two3 * (Hk + Hi) - 9 * K * mu * mu;
			}
			if (Jact(0) == 1 && Jact(1) == 1) {
				g(0, 1) = -9 * K * mu;
				g(1, 0) = mu * (-9 * K);
			}
			g.Invert(g_contra);

			m++;
		}

		// check maintain Kuhn-Tucker conditions
		f1 = norm_eta - (2 * G + two3 * Hk) * gamma(0) + mu * Invariant_1
			- 9 * K * mu * mu * gamma(0) - 9 * K * mu * gamma(1) - root23 * qiso(alpha1);
		f2 = Invariant_1 - 9 * K * mu * gamma(0) - 9 * K * gamma(1);

		if (count > 100) {
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
	Alpha1_n1 = alpha1;
	Alpha2_n1 = alpha2;

	//update epsilon_n1_p
	//first calculate n_covar
	// n_a = G_ab * n^b = covariant
	n_covar(0) = n(0);
	n_covar(1) = n(1);
	n_covar(2) = n(2);
	n_covar(3) = 2 * n(3);
	n_covar(4) = 2 * n(4);
	n_covar(5) = 2 * n(5);
	Eps_n1_p = Eps_n_p + (mu * gamma(0) + gamma(1)) * I1 + gamma(0) * n_covar;


	Invariant_ep = Eps_n1_p(0) + Eps_n1_p(1) + Eps_n1_p(2);

	norm_ep = sqrt(Eps_n1_p(0) * Eps_n1_p(0) + Eps_n1_p(1) * Eps_n1_p(1) + Eps_n1_p(2) * Eps_n1_p(2)
		+ 0.5 * (Eps_n1_p(3) * Eps_n1_p(3) + Eps_n1_p(4) * Eps_n1_p(4) + Eps_n1_p(5) * Eps_n1_p(5)));

	dev_ep = Eps_n1_p - one3 * Invariant_ep * I1;

	norm_dev_ep = sqrt(dev_ep(0) * dev_ep(0) + dev_ep(1) * dev_ep(1) + dev_ep(2) * dev_ep(2)
		+ 0.5 * (dev_ep(3) * dev_ep(3) + dev_ep(4) * dev_ep(4) + dev_ep(5) * dev_ep(5)));

	// update sigma
	Sig -= (3 * K * mu * gamma(0) + 3 * K * gamma(1)) * I1 + 2 * G * gamma(0) * n;

	s -= 2 * G * gamma(0) * n;
	Invariant_1 -= 9 * K * mu * gamma(0) + 9 * K * gamma(1);
	//Sig        = s + Invariant_1/3.0 * I1;

	//update beta_n1
	Zeta_n1 = Zeta_n - (two3 * Hk * gamma(0)) * n;

	//eta_n+1 = s_n+1 - beta_n+1;
	eta = s - Zeta_n1;
	norm_eta = sqrt(eta(0) * eta(0) + eta(1) * eta(1) + eta(2) * eta(2) + 2 * (eta(3) * eta(3) + eta(4) * eta(4) + eta(5) * eta(5)));

	// update Cep
	// note: Cep is contravariant
	if ((Jact(0) == 1) && (Jact(1) == 0)) {
		b1 = 2 * G * n + 3 * K * mu * I1;
		b2.Zero();
	}
	else if ((Jact(0) == 0) && (Jact(1) == 1)) {
		b1.Zero();
		b2 = 3 * K * I1;
	}
	else if ((Jact(0) == 1) && (Jact(1) == 1)) {
		b1 = 2 * G * n + 3 * K * mu * I1;
		b2 = 3 * K * I1;
	}

	temp1 = g_contra(0, 0) * b1 + g_contra(0, 1) * b2;
	temp2 = mu * temp1 + g_contra(1, 0) * b1 + g_contra(1, 1) * b2;

	NorCep = 0.0;
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			Cep(i, j) = Ce(i, j)
				+ 3 * K * I1(i) * temp2(j)
				+ 2 * G * n(i) * temp1(j)
				- 4 * G * G / norm_eta * gamma(0) * (IIdev(i, j) - n(i) * n(j));
			NorCep += Cep(i, j) * Cep(i, j);
		}
	}

	if (NorCep < 1e-10) {
		Cep = 1.0e-3 * Ce;
		opserr << "NorCep = " << NorCep << endln;
	}

	State(0) = Invariant_1;
	State(1) = norm_eta;
	State(2) = Invariant_ep;
	State(3) = norm_dev_ep;
	State(4) = norm_ep;

	return;
}

// Damage integration routine
void DPDamage::damage_integrator()
{
	if (df == 1) opserr << "Started damage routine.\n";
	//////// 1. Principal total and elastic strains //////////////////////////////////////////////////////////

	// Vector and matrices to be used within the method
	static Vector Eps_e(6);			// Elastic strains vector - eps^e

	// Elastic strain
	Eps_e = Eps - Eps_n1_p;

	// Principal and elastic strains from total strains
	Vector Eps_m = this->Eigen3(Eps);
	Vector Eps_e_m = this->Eigen3(Eps_e);

	if (df == 1) {
		opserr << "\nGot through the eigen3 method:\n";
		opserr << "Principal total strains: epr = ";
		for (int i = 0;i < 3;i++) opserr << Eps_m(i) << " ";
		opserr << "\nPrincipal elastic strains: epr_e = ";
		for (int i = 0;i < 3;i++) opserr << Eps_e_m(i) << " ";
		opserr << "\n";
	}

	//////// 2. Equivalent strains part ///////////////////////////////////////////////////////////////////////

	// Equivalent total strains 
	Vector e_tot(3);
	for (int i = 0;i < 3;i++)
		e_tot[i] = (1 - 2 * nu) * Eps_m[i] + nu * (Eps_m[0] + Eps_m[1] + Eps_m[2]);

	// Equivalent elastic strains
	Vector e_el(3);
	for (int i = 0;i < 3;i++)
		e_el[i] = (1 - 2 * nu) * Eps_e_m[i] + nu * (Eps_e_m[0] + Eps_e_m[1] + Eps_e_m[2]);

	// Positive and negative equivalent total strains
	Vector e_tp(3);
	Vector e_tn(3);
	for (int i = 0; i < 3; i++) {
		if (e_tot[i] > 0) {
			e_tp[i] = e_tot[i];
			e_tn[i] = 0;
		}
		else {
			e_tp[i] = 0;
			e_tn[i] = e_tot[i];
		}
	}

	// Total equivalent tensile strains
	double Yt = sqrt((e_tp[0]) * (e_tp[0]) + (e_tp[1]) * (e_tp[1]) + (e_tp[2]) * (e_tp[2]));

	// Total equivalent compressive strains
	double Yc = sqrt((e_tn[0]) * (e_tn[0]) + (e_tn[1]) * (e_tn[1]) + (e_tn[2]) * (e_tn[2]) +
		beta * ((e_tn[0]) * (e_tn[1]) + (e_tn[1]) * (e_tn[2]) + (e_tn[2]) * (e_tn[0])));

	// Positive and negative equivalent elastic strains
	Vector e_ep(3);
	Vector e_en(3);
	for (int i = 0; i < 3; i++) {
		if (e_el[i] > 0) {
			e_ep[i] = e_el[i];
			e_en[i] = 0;
		}
		else {
			e_ep[i] = 0;
			e_en[i] = e_el[i];
		}
	}
	// Elastic equivalent tensile strains
	double Yt_el = sqrt((e_ep[0]) * (e_ep[0]) + (e_ep[1]) * (e_ep[1]) + (e_ep[2]) * (e_ep[2]));

	// Elastic equivalent compressive strains
	double Yc_el = sqrt((e_en[0]) * (e_en[0]) + (e_en[1]) * (e_en[1]) + (e_en[2]) * (e_en[2]) +
		beta * ((e_en[0]) * (e_en[1]) + (e_en[1]) * (e_en[2]) + (e_en[2]) * (e_en[0])));

	if (df == 1) {
		opserr << "\nEvaluated all equivalent strains.\n";
		opserr << "e_tp = [ ";
		for (int i = 0;i < 3;i++) opserr << e_tp(i) << " ";
		opserr << "]\ne_tn = [ ";
		for (int i = 0;i < 3;i++) opserr << e_tn(i) << " ";
		opserr << "]\ne_ep = [ ";
		for (int i = 0;i < 3;i++) opserr << e_ep(i) << " ";
		opserr << "]\ne_en = [ ";
		for (int i = 0;i < 3;i++) opserr << e_en(i) << " ";
		opserr << "]\nYt = " << Yt;
		opserr << "\nYc = " << Yc;
		opserr << "\nYt_el = " << Yt_el;
		opserr << "\nYc_el = " << Yc_el;
		opserr << "\n";
	}

	////////// 3. Damage coefficients evaluation part ///////////////////////////////////////////////////////////////////////

	// Damage functions at step n
	double Ft = (Yt - Yt0) - Dt_n * (at * Yt + bt);
	double Fc = (Yc - Yc0) - Dc_n * (ac * Yc + bc);

	// Tensile damage parameter at step n+1
	if (Ft <= 0) {                // No damage
		Dt_n1 = Dt_n;
	}
	else {                       // Tensile damage
		Dt_n1 = (Yt - Yt0) / (at * Yt + bt);
	}

	// Compressive damage parameter at step n+1
	if (Fc <= 0) { // No damage
		Dc_n1 = Dc_n;
	}
	else { // Compressive damage
		Dc_n1 = (Yc - Yc0) / (ac * Yc + bc);
	}

	// Condition for traction damage Dt>Dc ---> Dt=max(Dc,Dt)
	if (Dt_n1 <= Dc_n1) {
		Dt_n1 = Dc_n1;
	}
	else {
		Dt_n1 = Dt_n1;
	}

	// Weighting coefficients
	double alpt;
	double alpc;
	if (Yt_el == 0 && Yc_el == 0) {
		alpt = 0;
		alpc = 1;
	} else {
		alpt = (Yt_el / Yt0) / ((Yt_el / Yt0) + (Yc_el / Yc0));
		alpc = 1 - alpt;
	}

	// Cumulative damage variable Dam = [1-D]^2
	double Dam = pow(((1 - Dt_n1) * alpt + (1 - Dc_n1) * alpc), 2);

	if (df == 1) {
		opserr << "\nEvaluated all damage variables.\n";
		opserr << "Dc = " << Dc_n1 << "\n";
		opserr << "Dt = " << Dt_n1 << "\n";
		opserr << "alpha_c = " << alpc << "\n";
		opserr << "alpha_t = " << alpt << "\n";
		opserr << "[1-D]^2 = " << Dam << "\n";
	}

	// Stress at n+1
	Sig = Dam * Cep * (Eps - Eps_n1_p);
	// Di Re 2018: Ct = [(1-D)^2*Cep - 2(1-D)*C*Eps_e*dD/dEps]

	if (df == 1) {
		opserr << "\nEvaluated stresses.\n";
		opserr << "Result: Sig = [ ";
		for (int i = 0;i < 6;i++) opserr << Sig(i) << " ";
		opserr << "]\n";
	}

	// Secant damage and tangent plastic constitutive matrix
	C = Dam * Cep;

	if (df == 1) {
		opserr << "\nEvaluated tangent matrix. All clear! Congrats!\n";
		opserr << "\nC matrix:\n";
		for (int i = 0;i < 6;i++) {
			for (int j = 0;j < 6;j++) {
				opserr << Cep(j, i) << " ";
			}
			opserr << "\n";
		}
		opserr << "\n";
	}
}

// Functions for plasticity
double DPDamage::qiso(double alpha1)
{
	return sigma_y + Hi * alpha1;
}

// Functions for damage
Vector DPDamage::principal_values(Vector e)
{
	// Invariants
	double I1 = e[0]+e[1]+e[2];
	double I2 = 0.5*(e[0]*e[1]+e[1]*e[2]+e[2]*e[0]-(e[3]*e[3]+e[4]*e[4]+e[5]*e[5]));
	double I3 = e[0]*e[1]*e[2]+2.0*e[3]*e[4]*e[5] -(e[3]*e[3]*e[2]+e[4]*e[4]*e[0]+e[5]*e[5]*e[1]);

	// Check if 1-2-3-4-5-6 position correspond to 11-22-33-12-23-31.
	// Evaluate if 12-23-31 correspond to eps=gamma/2 or gamma.

	// Intermediate quantities
	double Q = (3*I2-I1*I1)/9;
	double R = (2*pow(I1,3)-9*I1*I2+27*I3)/54;
	double arg = R/sqrt(pow(-Q,3));
	double theta = acos(arg);

	opserr << "Q=" << Q << " R=" << R << " arg=" << arg << " theta=" << theta << "\n";

	// Principal values
	Vector em(3);
	if (Q == 0 && R == 0) {
		em.Zero();
	}
	else {
		em[0] = 2.0 * sqrt(-Q) * cos(theta / 3) + I1 / 3;
		em[1] = 2.0 * sqrt(-Q) * cos((theta + 2 * 3.14157) / 3) + I1 / 3;
		em[2] = 2.0 * sqrt(-Q) * cos((theta + 4 * 3.14157) / 3) + I1 / 3;
	}
	return em;
}

// Function for principal strains
Vector DPDamage::Eigen3(const Vector e)
{
	//.... compute eigenvalues and vectors for a 3 x 3 symmetric matrix
	//
	//.... INPUTS:
	//        M(3,3) - matrix with initial values (only upper half used)
	//
	//.... OUTPUTS
	//        v(3,3) - matrix of eigenvectors (by column)
	//        d(3)   - eigenvalues associated with columns of v
	//        rot    - number of rotations to diagonalize
	//
	//---------------------------------------------------------------eig3==

	//.... Storage done as follows:
	//
	//       | v(1,1) v(1,2) v(1,3) |     |  d(1)  a(1)  a(3)  |
	//       | v(2,1) v(2,2) v(2,3) |  =  |  a(1)  d(2)  a(2)  |
	//       | v(3,1) v(3,2) v(3,3) |     |  a(3)  a(2)  d(3)  |
	//
	//        Transformations performed on d(i) and a(i) and v(i,j) become
	//        the eigenvectors.  
	//
	//---------------------------------------------------------------eig3==

	int     rot, its, i, j, k;
	double  g, h, aij, sm, thresh, t, c, s, tau;

	static Matrix  v(3, 3);
	static Vector  d(3);
	static Vector  a(3);
	static Vector  b(3);
	static Vector  z(3);

	static const double tol = 1.0e-08;

	//.... move array into one-d arrays
	a(0) = e(3);
	a(1) = e(4);
	a(2) = e(5);


	for (i = 0; i < 3; i++) {
		d(i) = e(i);
		b(i) = e(i);
		z(i) = 0.0;

		for (j = 0; j < 3; j++)
			v(i, j) = 0.0;

		v(i, i) = 1.0;

	} //end for i

	rot = 0;
	its = 0;

	sm = fabs(a(0)) + fabs(a(1)) + fabs(a(2));

	while (sm > tol) {
		//.... set convergence test and threshold
		if (its < 3)
			thresh = 0.011 * sm;
		else
			thresh = 0.0;

		//.... perform sweeps for rotations
		for (i = 0; i < 3; i++) {

			j = (i + 1) % 3;
			k = (j + 1) % 3;

			aij = a(i);

			g = 100.0 * fabs(aij);

			if (fabs(d(i)) + g != fabs(d(i)) ||
				fabs(d(j)) + g != fabs(d(j))) {

				if (fabs(aij) > thresh) {

					a(i) = 0.0;
					h = d(j) - d(i);

					if (fabs(h) + g == fabs(h))
						t = aij / h;
					else {
						//t = 2.0 * sign(h/aij) / ( fabs(h/aij) + sqrt(4.0+(h*h/aij/aij)));
						double hDIVaij = h / aij;
						if (hDIVaij > 0.0)
							t = 2.0 / (hDIVaij + sqrt(4.0 + (hDIVaij * hDIVaij)));
						else
							t = -2.0 / (-hDIVaij + sqrt(4.0 + (hDIVaij * hDIVaij)));
					}

					//.... set rotation parameters

					c = 1.0 / sqrt(1.0 + t * t);
					s = t * c;
					tau = s / (1.0 + c);

					//.... rotate diagonal terms

					h = t * aij;
					z(i) = z(i) - h;
					z(j) = z(j) + h;
					d(i) = d(i) - h;
					d(j) = d(j) + h;

					//.... rotate off-diagonal terms

					h = a(j);
					g = a[k];
					a(j) = h + s * (g - h * tau);
					a(k) = g - s * (h + g * tau);

					//.... rotate eigenvectors

					for (k = 0; k < 3; k++) {
						g = v(k, i);
						h = v(k, j);
						v(k, i) = g - s * (h + g * tau);
						v(k, j) = h + s * (g - h * tau);
					} // end for k

					rot = rot + 1;

				} // end if fabs > thresh 
			} //else
			else
				a(i) = 0.0;

		}  // end for i

		//.... update the diagonal terms
		for (i = 0; i < 3; i++) {
			b(i) = b(i) + z(i);
			d(i) = b(i);
			z(i) = 0.0;
		} // end for i

		its += 1;

		sm = fabs(a(0)) + fabs(a(1)) + fabs(a(2));

	} //end while sm
	static Vector  dd(3);
	if (d(0) > d(1))
	{
		if (d(0) > d(2))
		{
			dd(0) = d(0);
			if (d(1) > d(2))
			{
				dd(1) = d(1);
				dd(2) = d(2);
			}
			else
			{
				dd(1) = d(2);
				dd(2) = d(1);
			}
		}
		else
		{
			dd(0) = d(2);
			dd(1) = d(0);
			dd(2) = d(1);
		}
	}
	else
	{
		if (d(1) > d(2))
		{
			dd(0) = d(1);
			if (d(0) > d(2))
			{
				dd(1) = d(0);
				dd(2) = d(2);
			}
			else
			{
				dd(1) = d(2);
				dd(2) = d(0);
			}
		}
		else
		{
			dd(0) = d(2);
			dd(1) = d(1);
			dd(2) = d(0);
		}
	}

	return d;
}


Vector DPDamage::getState()
{
	return State;
}

// Other methods:
//	  - setTrialStrainIncr: used for cases where stress depends by
//		the rate of change of the strains.
//	  - getTangent: returns constitutive C matrix.
//	  - getInitialTangent: returns initial constitutive C matrix.
//	  - getStress: returns the stress vector (sigma)
//	  - getStrain: returns the strain vector (epsilon)

const Matrix&
DPDamage::getTangent (void)
{
  return C;
}

const Matrix&
DPDamage::getInitialTangent (void)
{
  return Ce;
}

const Vector&
DPDamage::getStress (void)
{
  return Sig;
}

const Vector&
DPDamage::getStrain (void)
{
  return Eps;
}

// History methods:
//	  - commitState: invoked at the end of every n-th cycle.
//		Resets n+1 values to n vaues.
//	  - revertToLastCommit: reverts n+1 values to n values.
//	  - revertToStart: reverts histry variables at n-th cycle to
//		initial value at n=1.
//	  - getCopy: creates a copy for every instance of the material.
//

int
DPDamage::commitState(void)
{
	// Recursive vectors
	Eps_n_p = Eps_n1_p;
	Alpha1_n = Alpha1_n1;
	Alpha2_n = Alpha2_n1;
	Zeta_n = Zeta_n1;

	// Damage parameters
	Dt_n = Dt_n1;
	Dc_n = Dc_n1;

	return 0;
}

int DPDamage::revertToLastCommit(void)
{
	return 0;
}

int
DPDamage::revertToStart (void)
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

NDMaterial*
DPDamage::getCopy(void)
{
	DPDamage* theCopy =
		new DPDamage(this->getTag(), E, nu, ft, fc, Hk, Hi,
			Yt0, bt, at, Yc0, bc, ac, beta);

	return theCopy;
}

NDMaterial*
DPDamage::getCopy(const char* type)
{
	if (strcmp(type, this->getType()) == 0)
		return this->getCopy();

	return 0;
}

const char*
DPDamage::getType (void) const
{
  return "ThreeDimensional";
}

int
DPDamage::getOrder (void) const
{
  return 6;
}


int 
DPDamage::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(10);

  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "DPDamage::sendSelf -- could not send Vector\n";
    return res;
  }
  
  return res;
}

int 
DPDamage::recvSelf(int commitTag, Channel &theChannel, 
					FEM_ObjectBroker &theBroker)
{
  static Vector data(10);
  
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "DPDamage::sendSelf -- could not send Vector\n";
    return res;
  }

  return res;
}

void 
DPDamage::Print(OPS_Stream &s, int flag) {
  opserr << "DPDamage: " << this->getTag();
  opserr << "strain: " << Eps;
  opserr << "strain: " << Sig;
  opserr << "tangent: " << C;
}