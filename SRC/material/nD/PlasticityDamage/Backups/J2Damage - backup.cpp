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
                                                                        
#include <J2Damage.h>           
#include <Channel.h>
#include <cmath>
#include <elementAPI.h>
#include <Matrix.h>
#include <Matrix.cpp>

double J2Damage::IIdev[3][3][3][3]; //rank 4 deviatoric 
double J2Damage::IIvol[3][3][3][3]; //rank 4 I bun I

// Static vectors and matrices
Vector J2Damage::eps(6);
Vector J2Damage::sig(6);
Matrix J2Damage::Cep(6,6);

const double J2Damage::one3 = 1.0 / 3.0;
const double J2Damage::two3 = 2.0 / 3.0;
const double J2Damage::root23 = sqrt(2.0 / 3.0);

// Starting kudos message
static int numPD2P = 0;

// Output debug flag
static int df = 0;

void *
OPS_NewJ2Damage(void)
{
	if (numPD2P == 0) { numPD2P++;
		opserr << "Two parameters damage and plasticity 3D constitutive law based on Addessi et al.\nBy LP and DF.\n";
	}

	NDMaterial* theMaterial = 0;

	int numArgs = OPS_GetNumRemainingInputArgs();

	if (numArgs < 14) {
		opserr << "Want: nDMaterial J2Damage $tag $E $nu $ft $fc $Hk $Hi\n";
		opserr << "$Yt0, $bt, $at, $Yc0, $bc, $ac, $beta)\n";
		return 0;
	}

	int tag;
	double dData[13];
	// Consiglia dati  così:  dData[4] = 0.6;

	int numData = 1;
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "WARNING invalid integer tag: nDMaterial J2Damage \n";
		return 0;
	}

	numData = numArgs - 1;;
	if (OPS_GetDouble(&numData, dData) != 0) {
		opserr << "WARNING invalid data: nDMaterial J2Damage: " << tag << "\n";
		return 0;
	}

	theMaterial = new J2Damage(tag,
		dData[0], dData[1], // E and nu
		dData[2], dData[3], dData[4], dData[5], // Druger Prager plasticity
		dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12]); // Addessi damage

	opserr << "Tcl script read and material created succesfully.\n";
	return theMaterial;
}

// Full constructor - Initial material terms determination
J2Damage::J2Damage(int tag, double _E, double _nu, // Parameters
	double _ft, double _fc, double _Hk, double _Hi, // Plasticity
	double _Yt0, double _bt, double _at, double _Yc0, double _bc, double _ac, double _beta) // Damage
	:NDMaterial(tag, ND_TAG_J2Damage),
	E(_E), nu(_nu),	ft(_ft), fc(_fc), Hk(_Hk), Hi(_Hi),
	Yt0(_Yt0), bt(_bt), at(_at), Yc0(_Yc0), bc(_bc), ac(_ac), beta(_beta),
	Eps(3,3), Eps_n_p(3,3),	Eps_n1_p(3,3), Sig(3,3),
	State(5)
{
	// Bulk and G modulus
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

	opserr << "eps: ";
	for (int i = 0;i < 6;i++) opserr << eps(i) << " ";
	opserr << "\nCt:\n";
	for (int i = 0;i < 6;i++) {
		for (int j = 0;j < 6;j++)
			opserr << Ct(j, i) << " ";
		opserr << "\n";
	}
	opserr << "sig: ";
	for (int i = 0;i < 6;i++) opserr << sig(i) << " ";
	opserr << "\nFull constructor executed. --------------------------------------------------\n";
	opserr << "\n";
}


// Null constructor - Empty terms determination
J2Damage::J2Damage()
  :NDMaterial (0, ND_TAG_J2Damage),
	E(0.0), nu(0.0), ft(0.0), fc(0.0), Hk(0.0), Hi(0.0),
	Yt0(0.0), bt(0.0), at(0.0), Yc0(0.0), bc(0.0), ac(0.0), beta(0.0),
	Eps(3, 3), Eps_n_p(3, 3), Eps_n1_p(3, 3), Sig(3, 3),
	State(5)
{
	K = 0.0;
	G = 0.0;
	sigma_y = 1e+10;
	mu = 0.0;

	H = 0.0;
	theta = 0.0;

	ElastFlag = 2;

	this->initialize();

	opserr << "Null constructor executed.\n";
}

// Destructor
J2Damage::~J2Damage ()
{}

// Zero internal variables
void J2Damage::initialize()
{
	// Strain vectors
	Eps.Zero();
	Eps_n_p.Zero();
	Eps_n1_p.Zero();

	// Stress and backSig tensors
	Sig.Zero();
	Flag = 1;

	alpha_n = 0.0;
	alpha_n1 = 0.0;

	//zero rank4 IIdev and IIvol
	int i, j, k, l;
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			for (k = 0; k < 3; k++) {
				for (l = 0; l < 3; l++) {

					IIvol[i][j][k][l] = 0.0;

					IIdev[i][j][k][l] = 0.0;

				} // end for l
			} // end for k
		} // end for j
	} // end for i

	//form rank4 IIvol 
	IIvol[0][0][0][0] = 1.0;
	IIvol[0][0][1][1] = 1.0;
	IIvol[0][0][2][2] = 1.0;
	IIvol[1][1][0][0] = 1.0;
	IIvol[1][1][1][1] = 1.0;
	IIvol[1][1][2][2] = 1.0;
	IIvol[2][2][0][0] = 1.0;
	IIvol[2][2][1][1] = 1.0;
	IIvol[2][2][2][2] = 1.0;

	//form rank4 IIdev

	IIdev[0][0][0][0] = two3; // 0.666667 
	IIdev[0][0][1][1] = -one3; //-0.333333 
	IIdev[0][0][2][2] = -one3; //-0.333333 
	IIdev[0][1][0][1] = 0.5;
	IIdev[0][1][1][0] = 0.5;
	IIdev[0][2][0][2] = 0.5;
	IIdev[0][2][2][0] = 0.5;
	IIdev[1][0][0][1] = 0.5;
	IIdev[1][0][1][0] = 0.5;
	IIdev[1][1][0][0] = -one3; //-0.333333 
	IIdev[1][1][1][1] = two3; // 0.666667 
	IIdev[1][1][2][2] = -one3; //-0.333333 
	IIdev[1][2][1][2] = 0.5;
	IIdev[1][2][2][1] = 0.5;
	IIdev[2][0][0][2] = 0.5;
	IIdev[2][0][2][0] = 0.5;
	IIdev[2][1][1][2] = 0.5;
	IIdev[2][1][2][1] = 0.5;
	IIdev[2][2][0][0] = -one3; //-0.333333 
	IIdev[2][2][1][1] = -one3; //-0.333333 
	IIdev[2][2][2][2] = two3; // 0.666667 

	plastic_integrator();

	// Elastic and elastoplastic constitutive matrix
	Ce = Cep;
	Ct = Ce;
	State.Zero();

	// Damage variables initialization
	Dt_n = 0;
	Dc_n = 0;

}

int
J2Damage::setTrialStrain(Vector const& v1, Vector const& v2) {
	return this->setTrialStrain(v1);
}

int
J2Damage::setTrialStrain(const Vector& strain)
{
	// Strain tensor initialization
	Eps.Zero();
	Eps = strainTensorMap(strain);

	if (df == 1) {
		opserr << "--------------------------------------------------------------------------\n";
		opserr << "\nTotal strains:     e = [ ";
		for (int i = 0;i < 6;i++) for (int j = 0;j < 6;j++) opserr << Eps(j,i) << " ";
		opserr << "]\n";
	}

	/* PLASTICITY -------------------------------------- */
	this->plastic_integrator();

	// Elastic strains
	Eps_e = Eps - Eps_n1_p;

	if (df == 1) {
		opserr << "\nPlastic routine executed succesfully:";
		opserr << "\nTotal strains: e = ";
		for (int i = 0;i < 6;i++) for (int j = 0;j < 6;j++) opserr << Eps(j,i) << " ";
		opserr << "\nElastic strains: ee = ";
		for (int i = 0;i < 6;i++) for (int j = 0;j < 6;j++) opserr << Eps_e(j, i) << " ";
		opserr << "\nPlastic strains: ep = ";
		for (int i = 0;i < 6;i++) for (int j = 0;j < 6;j++) opserr << Eps_n1_p(j, i) << " ";
		opserr << "\nStresses: sig = ";
		for (int i = 0;i < 6;i++) for (int j = 0;j < 6;j++) opserr << Sig(j, i) << " ";
		opserr << "\n";
	}

	/* DAMAGE -------------------------------------------*/
	this->damage_integrator();
	
	return 0;
}

// Plasticity integration routine
void J2Damage::plastic_integrator()
{
	const double tolerance = (1.0e-8) * sigma_y;
	const double dt = ops_Dt; //time step

	static Matrix dev_strain(3, 3); //deviatoric strain
	static Matrix dev_Sig(3, 3); //deviatoric Str
	static Matrix normal(3, 3);     //normal to yield surface

	double NbunN; //normal bun normal 

	double norm_tau = 0.0;   //norm of deviatoric Sig 
	double inv_norm_tau = 0.0;
	double phi = 0.0; //trial value of yield function
	double trace = 0.0; //trace of strain
	double gamma = 0.0; //consistency parameter
	double resid = 1.0;
	double tang = 0.0;
	double theta = 0.0;
	double theta_inv = 0.0;

	double c1 = 0.0;
	double c2 = 0.0;
	double c3 = 0.0;

	int i, j, k, l;
	int ii, jj;

	int iteration_counter;
	const int max_iterations = 25;

	//compute the deviatoric strains

	trace = Eps(0, 0) + Eps(1, 1) + Eps(2, 2);

	dev_strain = Eps;
	for (i = 0; i < 3; i++)
		dev_strain(i, i) -= (one3 * trace);

	//compute the trial deviatoric stress
	//   dev_Sig = (2.0*G) * ( dev_strain - Eps_n_p ) ;
	dev_Sig = 2.0*G*(dev_strain-Eps_n_p);

	//compute norm of deviatoric stress
	norm_tau = 0.0;
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++)
			norm_tau += dev_Sig(i, j) * dev_Sig(i, j);
	} //end for i 

	norm_tau = sqrt(norm_tau);

	if (norm_tau > tolerance) {
		inv_norm_tau = 1.0 / norm_tau;
		normal = inv_norm_tau * dev_Sig;
	}
	else {
		normal.Zero();
		inv_norm_tau = 0.0;
	} //end if 

	//compute trial value of yield function

	phi = norm_tau - root23 * q(alpha_n);

	// check if phi > 0 

	if (phi > 0.0) { //plastic

	   //solve for gamma 
		gamma = 0.0;
		resid = 1.0;
		iteration_counter = 0;
		while (fabs(resid) > tolerance) {

			resid = norm_tau - (2.0 * G) * gamma - root23 * q(alpha_n + root23 * gamma);

			tang = -(2.0 * G) - two3 * qprime(alpha_n + root23 * gamma);

			gamma -= (resid / tang);

			iteration_counter++;

			if (iteration_counter > max_iterations) {
				opserr << "More than " << max_iterations;
				opserr << " iterations in constituive subroutine J2-plasticity \n";
				break;
			} //end if 

		} //end while resid

		gamma *= (1.0 - 1e-08);

		//update plastic internal variables

		Eps_n1_p = Eps_n_p + gamma * normal;

		alpha_n1 = alpha_n + root23 * gamma;

		//recompute deviatoric Stress

		dev_Sig = (2.0 * G) * (dev_strain - Eps_n1_p);

		//compute the terms for plastic part of tangent

		theta = (2.0 * G) + two3 * qprime(alpha_n1);

		theta_inv = 1.0 / theta;

	}
	else { //elastic 

	  //update history variables -- they remain unchanged

		Eps_n1_p = Eps_n_p;
		alpha_n1 = alpha_n;

		//no extra tangent terms to compute 

		gamma = 0.0;
		theta = 0.0;
		theta_inv = 0.0;

	} //end if phi > 0

	//add on K part of stress
	Sig = dev_Sig;
	for (i = 0; i < 3; i++)
		Sig(i, i) += K * trace;

	//compute the tangent
	c1 = -4.0 * G * G;
	c2 = c1 * theta_inv;
	c3 = c1 * gamma * inv_norm_tau;

	for (ii = 0; ii < 6; ii++) {
		for (jj = 0; jj < 6; jj++) {

			index_map(ii, i, j);
			index_map(jj, k, l);

			NbunN = normal(i, j) * normal(k, l);

			//elastic terms
			tangent[i][j][k][l] = K * IIvol[i][j][k][l];

			tangent[i][j][k][l] += (2.0 * G) * IIdev[i][j][k][l];

			//plastic terms 
			tangent[i][j][k][l] += c2 * NbunN;

			tangent[i][j][k][l] += c3 * (IIdev[i][j][k][l] - NbunN);

			//minor symmetries 
			tangent[j][i][k][l] = tangent[i][j][k][l];
			tangent[i][j][l][k] = tangent[i][j][k][l];
			tangent[j][i][l][k] = tangent[i][j][k][l];

		} // end for jj
	} // end for ii

	// Form plasticity matrix "Cep"
	for (ii = 0; ii < 6; ii++) {
		for (jj = 0; jj < 6; jj++) {

			index_map(ii, i, j);
			index_map(jj, k, l);

			Cep(ii, jj) = tangent[i][j][k][l];

		} //end for j
	} //end for i

	// Form stress vector "sig"
	sig(0) = Sig(0, 0);
	sig(1) = Sig(1, 1);
	sig(2) = Sig(2, 2);
	sig(3) = Sig(0, 1);
	sig(4) = Sig(1, 2);
	sig(5) = Sig(2, 0);

	return;
}

// Damage integration routine
void J2Damage::damage_integrator()
{
	if (df == 1) opserr << "Started damage routine.\n";
	//////// 1. Principal total and elastic strains //////////////////////////////////////////////////////////

	// Strains vector
	eps = strainVectorMap(Eps);
	eps_e = strainVectorMap(Eps_e);

	// Principal and elastic strains from total strains
	Vector eps_m = eigen3(eps);
	Vector eps_e_m = eigen3(eps_e);

	if (df == 1) {
		opserr << "\nGot through the eigen3 method:\n";
		opserr << "Principal total strains: epr = ";
		for (int i = 0;i < 3;i++) opserr << eps_m(i) << " ";
		opserr << "\nPrincipal elastic strains: epr_e = ";
		for (int i = 0;i < 3;i++) opserr << eps_e_m(i) << " ";
		opserr << "\n";
	}

	//////// 2. Equivalent strains part ///////////////////////////////////////////////////////////////////////

	// Equivalent total strains 
	Vector e_tot(3);
	for (int i = 0;i < 3;i++)
		e_tot[i] = (1 - 2 * nu) * eps_m[i] + nu * (eps_m[0] + eps_m[1] + eps_m[2]);

	// Equivalent elastic strains
	Vector e_el(3);
	for (int i = 0;i < 3;i++)
		e_el[i] = (1 - 2 * nu) * eps_e_m[i] + nu * (eps_e_m[0] + eps_e_m[1] + eps_e_m[2]);

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
	}
	else {
		alpt = (Yt_el / Yt0) / ((Yt_el / Yt0) + (Yc_el / Yc0));
		alpc = 1 - alpt;
	}

	// Cumulative damage variable Dam = [1-D]^2
	double Dam = pow(((1 - Dt_n1) * alpt + (1 - Dc_n1) * alpc), 2);
	Dam = 1; // Remove damage

	if (df == 1) {
		opserr << "\nEvaluated all damage variables.\n";
		opserr << "Dc = " << Dc_n1 << "\n";
		opserr << "Dt = " << Dt_n1 << "\n";
		opserr << "alpha_c = " << alpc << "\n";
		opserr << "alpha_t = " << alpt << "\n";
		opserr << "[1-D]^2 = " << Dam << "\n";
	}

	// Stress at n+1
	sig = Dam * Cep * (eps - eps_n1_p);
	// Di Re 2018: Ct = [(1-D)^2*Cep - 2(1-D)*C*eps_e*dD/deps]

	if (df == 1) {
		opserr << "\nEvaluated Stres.\n";
		opserr << "Result: Sig = [ ";
		for (int i = 0;i < 6;i++) opserr << sig(i) << " ";
		opserr << "]\n";
	}

	// Secant damage and tangent plastic constitutive matrix
	Ct = Dam * Cep;

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

// matrix_index ---> tensor indices i,j
void J2Damage::index_map(int matrix_index, int& i, int& j) {
	switch (matrix_index + 1) { //add 1 for standard tensor indices

	case 1:
		i = 1;
		j = 1;
		break;

	case 2:
		i = 2;
		j = 2;
		break;

	case 3:
		i = 3;
		j = 3;
		break;

	case 4:
		i = 1;
		j = 2;
		break;

	case 5:
		i = 2;
		j = 3;
		break;

	case 6:
		i = 3;
		j = 1;
		break;


	default:
		i = 1;
		j = 1;
		break;

	} //end switch

	i--; //subtract 1 for C-indexing
	j--;

	return;
}

// Strain tensor to strain vector and viceversa
Vector J2Damage::strainVectorMap(Matrix E) {
	Vector e(6);
	e(0) = E(0, 0);
	e(1) = E(1, 1);
	e(2) = E(2, 2);
	e(3) = 2.0 * E(0, 1);
	e(4) = 2.0 * E(1, 2);
	e(5) = 2.0 * E(2, 0);
	return e;
}
Matrix J2Damage::strainTensorMap(Vector e) {
	Matrix E(3,3);
	E(0,0) = e(0);
	E(1,1) = e(1);
	E(2,2) = e(2);
	E(0,1) = e(3)/2;
	E(1,2) = e(4)/2;
	E(2,0) = e(5)/2;
	return E;
}

//hardening function
double J2Damage::q(double alpha1)
{
	return    sigma_y + H * alpha1;
}

//hardening function derivative
double J2Damage::qprime(double alpha1)
{
	return  H;
}

// Functions for plasticity
double J2Damage::qiso(double alpha1)
{
	return sigma_y + Hi * alpha1;
}

// Function for principal strains
Vector J2Damage::eigen3(const Vector e)
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

Vector J2Damage::getState()
{
	return State;
}

// Other methods:
//	  - setTrialStrainIncr: used for cases where Sig depends by
//		the rate of change of the strains.
//	  - getTangent: returns constitutive C matrix.
//	  - getInitialTangent: returns initial constitutive C matrix.
//	  - getStress: returns the Sig vector (sigma)
//	  - getStrain: returns the strain vector (epsilon)

const Matrix&
J2Damage::getTangent (void)
{
  return Ct;
}

const Matrix&
J2Damage::getInitialTangent (void)
{
  return Ce;
}

const Vector&
J2Damage::getStress (void)
{
  return sig;
}

const Vector&
J2Damage::getStrain (void)
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
J2Damage::commitState(void)
{
	// Recursive vectors
	Eps_n_p = Eps_n1_p;
	alpha_n = alpha_n1;

	// Damage parameters
	Dt_n = Dt_n1;
	Dc_n = Dc_n1;

	return 0;
}

int J2Damage::revertToLastCommit(void)
{
	return 0;
}

int
J2Damage::revertToStart (void)
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
J2Damage::getCopy(void)
{
	J2Damage* theCopy =
		new J2Damage(this->getTag(), E, nu, ft, fc, Hk, Hi,
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

const char*
J2Damage::getType (void) const
{
  return "ThreeDimensional";
}

int
J2Damage::getOrder (void) const
{
  return 6;
}


int 
J2Damage::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(10);

  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "J2Damage::sendSelf -- could not send Vector\n";
    return res;
  }
  
  return res;
}

int 
J2Damage::recvSelf(int commitTag, Channel &theChannel, 
					FEM_ObjectBroker &theBroker)
{
  static Vector data(10);
  
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "J2Damage::sendSelf -- could not send Vector\n";
    return res;
  }

  return res;
}

void 
J2Damage::Print(OPS_Stream &s, int flag) {
  opserr << "J2Damage: " << this->getTag();
  opserr << "strain: " << Eps;
  opserr << "strain: " << Sig;
  opserr << "tangent: " << Ct;
}