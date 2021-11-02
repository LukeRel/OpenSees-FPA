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

// Written in C++: Daniela Fusco and Luca Parente
// Created: 11/21
                                                                        
#include <PlasticDamage2P.h>           
#include <Channel.h>
#include <cmath>
#include <elementAPI.h>

const double PlasticDamage2P::one3 = 1.0 / 3.0;
const double PlasticDamage2P::two3 = 2.0 / 3.0;
const double PlasticDamage2P::root23 = sqrt(2.0 / 3.0);

static Vector Iv6(6); 
static Matrix Ivp(6,6); 
static Matrix Idp(6,6); 
static Matrix I(6,6);
static Matrix Id(6,6); 
static Matrix P(6, 6);   //Matrice deviatorica
static Matrix l(6, 6);   //Matrice volumetrica

void *
OPS_NewPlasticDamage2P(void)
{
  NDMaterial *theMaterial = 0;
  
  int numArgs = OPS_GetNumRemainingInputArgs();
  
  if (numArgs < 20 || numArgs > 30) {
    opserr << "Want: nDMaterial PlasticDamage2P $tag $E $nu $ft $fc $Hk $Hi\n";
	opserr << "$r_bar $Kinfinity $Kinit, $d1, $d2, $mDen\n";
	opserr << "$Yt0, $bt, $at, $Yc0, $bc, $ac, $kappa)\n";
    return 0;	
  }
  
  int tag;
  double dData[19];
  //Consiglia dati  così:  dData[4] = 0.6;
  
  
  int numData = 1;
  if (OPS_GetInt(&numData, &tag) != 0) {
    opserr << "WARNING invalid integer tag: nDMaterial EasticIsotropic \n";
    return 0;
  }
  
  numData = numArgs - 1;;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid data: nDMaterial EasticIsotropic : " << tag <<"\n";
    return 0;
  }  
  
  theMaterial = new PlasticDamage2P(tag, 
					    dData[0], dData[1], dData[2],dData[3], 
					    dData[4], dData[5], dData[6],dData[7],
	                    dData[8], dData[9], dData[10], dData[11],
	                    dData[12], dData[13], dData[14], dData[15],
	                    dData[16], dData[17], dData[18]);

  return theMaterial;
}

// Full constructor _ Input
PlasticDamage2P::PlasticDamage2P(int tag, double E, double nu, //(CALCOLA BULK E SHEAR)
	double ft, double fc, //(CALCOLA s_y e r)
	double Hk, double Hi, //(CALCOLA H e t)
	double r_bar, double Kinfinity, double Kinit, double d1, double d2, double mDen, //(PARAMETRI CHE SI POTRANNO IGNORARE)
	double Yt0, double bt, double at, double Yc0, double bc, double ac, double kappa) //(IN GATTA)
	:NDMaterial(tag, ND_TAG_PlasticDamage2P),
	eps(6), sig(6), sige(6), eps_p(6), sigeP(6),
	epsCommit(6), sigCommit(6), sigeCommit(6), eps_pCommit(6), sigePCommit(6),
	Ce(6, 6), C(6, 6), Ccommit(6, 6),
	mEpsilon(6),
	mEpsilon_n_p(6),
	mEpsilon_n1_p(6),
	mSigma(6),
	mBeta_n(6),
	mBeta_n1(6),
	mCe(6, 6),
	mCep(6, 6),
	mI1(6),
	mIIvol(6, 6),
	mIIdev(6, 6),
	mState(5)

	
{
	//Calcola Bulk e Shear
	mK = E / (3 * (1 - 2 * nu));
	mG = E / (2 * (1 + nu));

	//Calcola parametri plasticità
	msigma_y = 2 * fc * ft / (fc + ft);
	mrho = root23 * (fc - ft) / (fc + ft);

	//Parametri di harding
	mHard = Hi + Hk;
	mtheta = Hi / (Hi + Hk);

	//Parametri aggiuntivi plasticità
	mrho_bar = r_bar;
	mKinf = Kinfinity;
	mKo = Kinit;
	mdelta1 = d1;
	mdelta2 = d2;
	massDen = mDen;

	//FORSE DA CANCELLARE
	mKref0 = mK;
	mGref0 = mG;
	msigma_y0 = msigma_y;


	if (mrho == 0.0) {
		mTo = 1e10;
	}
	else {
		mTo = root23 * msigma_y / mrho;
	}
	// set the elastic flag
	// 0 = elastic+no param update; 1 = elastic+param update; 2 = elastoplastic+no param update (default)
	mElastFlag = 2;

	// Use these values to deactivate yield surface 1 - Create Pure Tension Cutoff
	// msigma_y = 1e10;
	// mTo      = 100;

	//Inizializzazione
	eps.Zero();
	sig.Zero();
	sige.Zero();
	eps_p.Zero();
	sigeP.Zero();
	Ce.Zero();

	Iv6.Zero(); Iv6(0) = 1.;Iv6(1) = 1.;Iv6(2) = 1.;

	// Tensore Volumetrico (1, in Di Re)
	Ivp.Zero();
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			Ivp(i, j) = 1.;

	// Tensore Deviatorico (P, in Di Re)
	Idp.Zero();
	I.Zero();
	Id.Zero();
	for (int i = 0; i < 6; i++) {
		Idp(i, i) = 1.;
		if (i < 3) {
			I(i, i) = 1.0;
			Id(i, i) = 1.0;
		}
		else {
			I(i, i) = 0.5;
			Id(i, i) = 0.5;
		}
	}
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++) {
			Id(i, j) = Idp(i, j) - 1 / 3.;
			Idp(i, j) = Id(i, j);
		}

	Ce.addMatrix(0.0, Ivp, mK);
	Ce.addMatrix(1.0, Id, 2. * mG);

	C = Ce;
	P = Id;              // Matrice deviatorica
	l = Ivp;             // Matrice volumetrica



	
	// Inizializza variabili di danno
	Dt_n = 0;
	Dc_n = 0;
	

	this->initialize();
}


// Null constructor _ Input
PlasticDamage2P::PlasticDamage2P()
  :NDMaterial (0, ND_TAG_PlasticDamage2P),
	eps(6), sig(6), sige(6), eps_p(6), sigeP(6),
	epsCommit(6), sigCommit(6), sigeCommit(6), eps_pCommit(6), sigePCommit(6),
	Ce(6, 6), C(6, 6), Ccommit(6, 6),
	mEpsilon(6),
	mEpsilon_n_p(6),
	mEpsilon_n1_p(6),
	mSigma(6),
	mBeta_n(6),
	mBeta_n1(6),
	mCe(6, 6),
	mCep(6, 6),
	mI1(6),
	mIIvol(6, 6),
	mIIdev(6, 6),
	mState(5)
{

}

PlasticDamage2P::~PlasticDamage2P ()
{

}




//zero internal variables
void PlasticDamage2P::initialize()
{
	mEpsilon.Zero();
	mEpsilon_n_p.Zero();
	mEpsilon_n1_p.Zero();

	mSigma.Zero();
	mBeta_n.Zero();
	mBeta_n1.Zero();

	mAlpha1_n = 0.0;
	mAlpha1_n1 = 0.0;
	mAlpha2_n = 0.0;
	mAlpha2_n1 = 0.0;
	mFlag = 1;

	mHprime = (1 - mtheta) * mHard;

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
	mState.Zero();
}





// Plastic - Damage routine_________________________________
int
PlasticDamage2P::setTrialStrain(const Vector& strain)
{
	// Vector and matrices to be used within the method
	static Vector epsilon_e(6);			// Elastic strains vector - eps^e
	

	// Strains vector from global analysis
	eps = strain;

	// Variable renaming
	mEpsilon_n_p = eps_p;
	mEpsilon = eps;



	/* PLASTICITY -------------------------------------- */
	this->plastic_integrator();



	/* DAMAGE -------------------------------------------*/
	
	// Elastic strain
	epsilon_e = eps - mEpsilon_n1_p;

	// Principal strains from total strains
	Vector eps_princ = principal_values(eps);

	// Principal strains from elastic strains
	Vector epse_princ = principal_values(epsilon_e);

	// Equivalent total strains 
	Vector e_tot(3);
	e_tot[0] = (1 - 2 * nu) * eps_princ[0] + nu * (eps_princ[0] + eps_princ[1] + eps_princ[2]);
	e_tot[1] = (1 - 2 * nu) * eps_princ[1] + nu * (eps_princ[0] + eps_princ[1] + eps_princ[2]);
	e_tot[2] = (1 - 2 * nu) * eps_princ[2] + nu * (eps_princ[0] + eps_princ[1] + eps_princ[2]);

	// Equivalent elastic strains
	Vector e_el(3);
	e_el[0] = (1 - 2 * nu) * epse_princ[0] + nu * (epse_princ[0] + epse_princ[1] + epse_princ[2]);
	e_el[1] = (1 - 2 * nu) * epse_princ[1] + nu * (epse_princ[0] + epse_princ[1] + epse_princ[2]);
	e_el[2] = (1 - 2 * nu) * epse_princ[2] + nu * (epse_princ[0] + epse_princ[1] + epse_princ[2]);

	// Equivalent total strains positive and negative
	Vector e_tot_pos(3);
	Vector e_tot_neg(3);
	for(int i = 0; i < 2; i++) {
		if (e_tot[i] > 0) {
			e_tot_pos[i]= e_tot[i];
			e_tot_neg[i] = 0;
		}
		else {
			e_tot_pos[i] = 0;
			e_tot_neg[i] = e_tot[i];
		}
	}
	//Equivalent Deformation Total Traction
	double Yt = sqrt((e_tot_pos[0]) * (e_tot_pos[0]) + (e_tot_pos[1]) * (e_tot_pos[1]) + (e_tot_pos[2]) * (e_tot_pos[2]));
	//Equivalent Deformation Total Compression
	double Yc = sqrt((e_tot_neg[0]) * (e_tot_neg[0]) + (e_tot_neg[1]) * (e_tot_neg[1]) + (e_tot_neg[2]) * (e_tot_neg[2])) +
		        (kappa/2)*((e_tot_neg[0]) * (e_tot_neg[0]) + (e_tot_neg[0]) * (e_tot_neg[1]) + (e_tot_neg[0]) * (e_tot_neg[2])+
					       (e_tot_neg[1]) * (e_tot_neg[0]) + (e_tot_neg[1]) * (e_tot_neg[1]) + (e_tot_neg[1]) * (e_tot_neg[2])+
				           (e_tot_neg[2]) * (e_tot_neg[0]) + (e_tot_neg[2]) * (e_tot_neg[1]) + (e_tot_neg[2]) * (e_tot_neg[2]));




	// Equivalent elastic strains positive and negative
	Vector e_el_pos(3);
	Vector e_el_neg(3);
	for (int i = 0; i < 2; i++) {
		if (e_el[i] > 0) {
			e_el_pos[i] = e_el[i];
			e_el_neg[i] = 0;
		}
		else {
			e_el_pos[i] = 0;
			e_el_neg[i] = e_el[i];
		}
	}
	//Equivalent Deformation Elastic Traction
	double Yt_el = sqrt((e_el_pos[0]) * (e_el_pos[0]) + (e_el_pos[1]) * (e_el_pos[1]) + (e_el_pos[2]) * (e_el_pos[2]));
	//Equivalent Deformation Elastic Compression
	double Yc_el = sqrt((e_el_neg[0]) * (e_el_neg[0]) + (e_el_neg[1]) * (e_el_neg[1]) + (e_el_neg[2]) * (e_el_neg[2])) +
		          (kappa / 2) * ((e_el_neg[0]) * (e_el_neg[0]) + (e_el_neg[0]) * (e_el_neg[1]) + (e_el_neg[0]) * (e_el_neg[2]) +
			                     (e_el_neg[1]) * (e_el_neg[0]) + (e_el_neg[1]) * (e_el_neg[1]) + (e_el_neg[1]) * (e_el_neg[2]) +
			                     (e_el_neg[2]) * (e_el_neg[0]) + (e_el_neg[2]) * (e_el_neg[1]) + (e_el_neg[2]) * (e_el_neg[2]));



	// Historical damage parameters
	// Dt_n Variabile a trazione di danno allo step precedente 
	// Dc_n Variabile a compressione di danno allo step precedente
	// Dt_n1 Variabile a trazione di danno allo step corrente
	// Dc_n1 Variabile a compressione di danno allo step corrente


	// Yield Functions at step n
	double Ft = (Yt - Yt0) - Dt_n * (at * Yt + bt);
	double Fc = (Yc - Yc0) - Dc_n * (ac * Yc + bc);



	// Traction Damage parameter at step n+1
	if (Ft < 0) {                // No damage
		Dt_n1= Dt_n;
	}
	else {                       // Damage Traction
		Dt_n1 = (Yt - Yt0) / (at * Yt + bt);
	}


	// Compression Damage parameter at step n+1
	if (Fc < 0) {                // No damage
		Dc_n1 = Dc_n;
	}
	else {                       // Damage Compression
		Dc_n1 = (Yc - Yc0) / (ac * Yc + bc);
	}


	// Condition for traction damage Dt>Dc ---> Dt=max(Dc,Dt)
	if (Dt_n1 < Dc_n1) {                
		Dt_n1 = Dc_n1;
	}
	else {                      
		Dt_n1 = Dt_n1;
	}

	// Weighting coefficients
	double alpt = (Yt_el / Yt0) / ((Yt_el / Yt0) + (Yc_el / Yc0));
	double alpc = 1 - alpt;


	//Stress at n+1 !!!Attenzione lui in plastic integrator non lo calcola come cep(eps-epsp) quindi vedi se posso scrivere cosi' oppure lui ha solo complicato le cose
	mSigma = pow(((1 - Dt_n1) * alpt + (1 - Dc_n1) * alpc), 2) * mCep * (eps - mEpsilon_n1_p);
	//devi inserire nel codice un get stress e un get tangent!!

	// Matrice secante di danno e tangente in plasticità
	C = pow(((1 - Dt_n1) * alpt + (1 - Dc_n1) * alpc), 2) * mCep;
	//devi inserire nel codice un get stress e un get tangent!!


	return 0;
}


// Plasticity integration routine___________________________
void PlasticDamage2P::plastic_integrator()
{
	bool okay;		// boolean variable to ensure satisfaction of multisurface kuhn tucker conditions
	double f1;
	double f2;
	double norm_eta;
	double Invariant_1;
	double Invariant_ep;
	double norm_ep;
	double norm_dev_ep;
	Vector epsilon_e(6);
	Vector s(6);
	Vector eta(6);
	Vector dev_ep(6);
	Vector Jact(2);

	double fTOL;
	double gTOL;
	fTOL = 0.0;
	gTOL = -1.0e-10;

	double NormCep;

	double alpha1;			// hardening parameter for DP surface
	double alpha2;			// hardening parameter for tension cut-off
	Vector n(6);			// normal to the yield surface in strain space
	Vector R(2);			// residual vector
	Vector gamma(2);		// vector of consistency parameters
	Vector dgamma(2);		// incremental vector of consistency parameters
	Matrix g(2, 2);			// jacobian of the corner region (return map)
	Matrix g_contra(2, 2);	// inverse of jacobian of the corner region

	// Calcolo variabili di storia
	mEpsilon_n1_p = mEpsilon_n_p;	// eps^p_n+1_trial
	mAlpha1_n1 = mAlpha1_n;			// alpha1_n+1_trial
	mAlpha2_n1 = mAlpha2_n;			// alpha2_n+1_trial
	mBeta_n1 = mBeta_n;				// beta_n+1_trial

	// Deformazioni elastiche
	epsilon_e = mEpsilon - mEpsilon_n1_p;

	// Tensioni elastiche trial
	mSigma = mCe * epsilon_e;

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
	while (!okay) {

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
		for (int i = 0; i < 2; i++) {
			if (Jact(i) == 1) {
				R(0) = norm_eta - (2 * mG + two3 * mHprime) * gamma(0) + mrho * Invariant_1
					- 9 * mK * mrho * mrho_bar * gamma(0) - 9 * mK * mrho * gamma(1) - root23 * Kiso(alpha1);
				g(0, 0) = -2 * mG - two3 * (mHprime + Kisoprime(alpha1)) - 9 * mK * mrho * mrho_bar;
			}
			else if (Jact(i) == 2) {
				R(1) = Invariant_1 - 9 * mK * mrho_bar * gamma(0) - 9 * mK * gamma(1) - T(alpha2);
				g(1, 1) = -9 * mK + mdelta2 * T(alpha2);
			}
		}
		if (Jact(0) == 1 && Jact(1) == 1) {
			g(0, 1) = -9 * mK * mrho;
			g(1, 0) = mrho_bar * (-9 * mK + mdelta2 * T(alpha2));
		}
		g.Invert(g_contra);

		// iteration counter
		int m = 0;

		//iterate
		while ((fabs(R.Norm()) > 1e-10) && (m < 10)) {

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
			for (int i = 0; i < 2; i++) {
				if (Jact(i) == 1) {
					R(0) = norm_eta - (2 * mG + two3 * mHprime) * gamma(0) + mrho * Invariant_1
						- 9 * mK * mrho * mrho_bar * gamma(0) - 9 * mK * mrho * gamma(1) - root23 * Kiso(alpha1);
					g(0, 0) = -2 * mG - two3 * (mHprime + Kisoprime(alpha1)) - 9 * mK * mrho * mrho_bar;
				}
				else if (Jact(i) == 2) {
					R(1) = Invariant_1 - 9 * mK * mrho_bar * gamma(0) - 9 * mK * gamma(1) - T(alpha2);
					g(1, 1) = -9 * mK + mdelta2 * T(alpha2);
				}
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

	if (NormCep < 1e-10) {
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





// Functions for plasticity___________________________________
double PlasticDamage2P::Kiso(double alpha1)
{
	return msigma_y + mtheta * mHard * alpha1 + (mKinf - mKo) * (1 - exp(-mdelta1 * alpha1));
}

double PlasticDamage2P::Kisoprime(double alpha1)
{
	return mtheta * mHard + (mKinf - mKo) * mdelta1 * exp(-mdelta1 * alpha1);
}

double PlasticDamage2P::T(double alpha2)
{
	return mTo * exp(-mdelta2 * alpha2);
}

double PlasticDamage2P::deltaH(double dGamma)
{
	return mHprime * root23 * dGamma;
}




// Functions for damage________________________________________
Vector PlasticDamage2P::principal_values(Vector e)
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

	// Principal values
	Vector em(6);
	em[0] = 2.0*sqrt(-Q)*cos(theta/3) + I1/3;
	em[1] = 2.0*sqrt(-Q)*cos((theta+2*3.14157)/3) + I1/3;
	em[2] = 2.0*sqrt(-Q)*cos((theta+4*3.14157)/3) + I1/3;
	return em;
}

Vector PlasticDamage2P::getState()
{
	return mState;
}

// Other methods:
//	  - setTrialStrainIncr: used for cases where stress depends by
//		the rate of change of the strains.
//	  - getTangent: returns constitutive C matrix.
//	  - getInitialTangent: returns initial constitutive C matrix.
//	  - getStress: returns the stress vector (sigma)
//	  - getStrain: returns the strain vector (epsilon)

int
PlasticDamage2P::setTrialStrainIncr (const Vector &strain)
{
  eps += strain;
  this->setTrialStrain(eps);
  return 0;
}

int
PlasticDamage2P::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  eps += strain;
  this->setTrialStrain(eps);
  return 0;
}

const Matrix&
PlasticDamage2P::getTangent (void)
{
  return C;
}

const Matrix&
PlasticDamage2P::getInitialTangent (void)
{
  return Ce;
}

const Vector&
PlasticDamage2P::getStress (void)
{
  return sig;
}

const Vector&
PlasticDamage2P::getStrain (void)
{
  return eps;
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
PlasticDamage2P::commitState(void)
{
	mEpsilon_n_p = mEpsilon_n1_p;
	mAlpha1_n = mAlpha1_n1;
	mAlpha2_n = mAlpha2_n1;
	mBeta_n = mBeta_n1;


	Dt_n = Dt_n1;
	Dc_n = Dc_n1;


	sigCommit = sig;
	sigeCommit = sige;
	sigePCommit = sige;

	C = Ccommit;// cancella


	return 0;
}

int PlasticDamage2P::revertToLastCommit(void)
{
	return 0;
}

int
PlasticDamage2P::revertToStart (void)
{
	if (ops_InitialStateAnalysis) {
		// do nothing, keep state variables from last step
	}
	else {
		// normal call for revertToStart (not initialStateAnalysis)
		this->initialize();
	}

  eps.Zero();
  sig.Zero();
  sige.Zero();
  eps_p.Zero();
  sigeP.Zero();
  Ce.Zero();

  return 0;
}

NDMaterial*
PlasticDamage2P::getCopy(const char *type)
{
  if (strcmp(type,"ThreeDimensional") == 0 || strcmp(type,"3D") == 0) {
	  PlasticDamage2P* theCopy =
		  new PlasticDamage2P(this->getTag(), E, nu, ft, fc, 
			  Hk, Hi, r_bar,
			  Kinfinity, Kinit, d1, d2, mDen,
			  Yt0, bt, at, Yc0, bc, ac, kappa);
  
    return theCopy;  
  } else {
    return 0;
  }
}

NDMaterial*
PlasticDamage2P::getCopy (void)
{
  PlasticDamage2P *theCopy =
    new PlasticDamage2P(this->getTag(), E, nu, ft, fc,
		Hk, Hi, r_bar,
		Kinfinity, Kinit, d1, d2, mDen,
		Yt0, bt, at, Yc0, bc, ac, kappa);
  return theCopy;
}

const char*
PlasticDamage2P::getType (void) const
{
  return "ThreeDimensional";
}

int
PlasticDamage2P::getOrder (void) const
{
  return 6;
}


int 
PlasticDamage2P::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(10);

  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "PlasticDamage2P::sendSelf -- could not send Vector\n";
    return res;
  }
  
  return res;
}

int 
PlasticDamage2P::recvSelf(int commitTag, Channel &theChannel, 
					FEM_ObjectBroker &theBroker)
{
  static Vector data(10);
  
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "PlasticDamage2P::sendSelf -- could not send Vector\n";
    return res;
  }

  return res;
}

void 
PlasticDamage2P::Print(OPS_Stream &s, int flag) {
  opserr << "PlasticDamage2P: " << this->getTag();
  opserr << "strain: " << eps;
  opserr << "strain: " << sig;
  opserr << "tangent: " << C;
}