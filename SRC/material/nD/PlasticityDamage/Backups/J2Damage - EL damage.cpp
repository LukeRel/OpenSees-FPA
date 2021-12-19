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
** ****************************************************************** */

// $Revision: 1.7 $
// $Date: 2008-10-20 22:23:03 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/J2Damage.cpp,v $

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
#include <FEM_ObjectBroker.h>

#include <Information.h>
#include <Parameter.h>

#include <string.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>

//parameters
const double J2Damage::one3 = 1.0 / 3.0;
const double J2Damage::two3 = 2.0 / 3.0;
const double J2Damage::four3 = 4.0 / 3.0;
const double J2Damage::root23 = sqrt(2.0 / 3.0);

double J2Damage::initialTangent[3][3][3][3];   //material tangent
double J2Damage::IIdev[3][3][3][3]; //rank 4 deviatoric 
double J2Damage::IbunI[3][3][3][3]; //rank 4 I bun I

//static vectors and matrices
Vector J2Damage::v(6);
Vector J2Damage::s(6);
Matrix J2Damage::m(6,6);
Vector J2Damage::strain_vec(6);
Vector J2Damage::stress_vec(6);
Vector J2Damage::strain_k(6);
Vector J2Damage::stress_k(6);
Matrix J2Damage::tangent_e(6, 6);
Matrix J2Damage::tangent_ep(6, 6);
Matrix J2Damage::tangent_matrix(6, 6);

//debug flags
int dFlag1 = 0; //everything except damage
int dFlag2 = 0; //damage

void* OPS_J2Damage()
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 13) {
        opserr << "Want: nDMaterial J2Damage $tag $E $nu $sig_y $Hk $Hi\n";
        opserr << "$Yt0, $bt, $at, $Yc0, $bc, $ac, $beta)\n";
        return 0;
    }

    int tag;
    numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
        opserr << "WARNING invalid J2Damage tag\n";
        return 0;
    }

    double data[12] = { 0,0,0,0,0,0,0,0,0,0,0,0 };
    numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 12) {
        numdata = 12;
    }
    if (OPS_GetDoubleInput(&numdata, data)) {
        opserr << "WARNING invalid J2Damage double inputs\n";
        return 0;
    }

    NDMaterial* mat = new J2Damage(tag,
        data[0], data[1], // E and nu
        data[2], data[3], data[4], // Druker Prager plasticity
        data[5], data[6], data[7], data[8], data[9], data[10], data[11]); // Addessi damage
    if (mat == 0) {
        opserr << "WARNING: failed to create J2Damage material\n";
        return 0;
    }

    if (1 == 1) {
        opserr << "------------------------------------------------------------------------------------" << endln;
        opserr << "Created J2Damage succesfully. Based on Von Mises (J2) plasticity and Addessi damage." << endln;
        opserr << "Inputs:" << endln;
        opserr << "E = " << data[0] << "\nnu = " << data[1] << "\nsig_y = " << data[2] << "\nHki = " << data[3] << endln;
        opserr << "Hi = " << data[4] << "\nYt0 = " << data[5] << "\nbt = " << data[6] << "\nat = " << data[7] << endln;
        opserr << "Yc0 = " << data[8] << "\nbc = " << data[9] << "\nac = " << data[10] << "\nbeta = " << data[11] << endln;
        opserr << "------------------------------------------------------------------------------------\n" << endln;
    }

    return mat;
}

//zero internal variables
void J2Damage::zero()
{
    xi_n = 0.0;
    xi_nplus1 = 0.0;

    epsilon_e.Zero();
    epsilon_p_n.Zero();
    epsilon_p_nplus1.Zero();

    stress.Zero();
    stress_k.Zero();

    strain_k.Zero();
    strain.Zero();

    D = 0.0;
    D_k = 0.0;

    Dt = 0.0;
    Dt_k = 0.0;

    Dc = 0.0;
    Dc_k = 0.0;
}

//null constructor
J2Damage::J2Damage() :
    NDMaterial(0, ND_TAG_J2Damage),
    epsilon_e(3,3),
    epsilon_p_n(3, 3),
    epsilon_p_nplus1(3, 3),
    stress(3, 3),
    strain(3, 3),
    parameterID(0)
{
    K = 0.0;
    G = 0.0;
    sig_y = 0.0;
    sigma_infty = 0.0;
    delta = 0.0;
    H = 0.0;
    eta = 0.0;
    rho = 0.0;

    this->zero();     // or (*this).zero( ) 

    int i, j, k, l;

    //zero rank4 IIdev and IbunI 
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < 3; k++) {
                for (l = 0; l < 3; l++) {

                    IbunI[i][j][k][l] = 0.0;

                    IIdev[i][j][k][l] = 0.0;

                } // end for l
            } // end for k
        } // end for j
    } // end for i


    //form rank4 IbunI 

    IbunI[0][0][0][0] = 1.0;
    IbunI[0][0][1][1] = 1.0;
    IbunI[0][0][2][2] = 1.0;
    IbunI[1][1][0][0] = 1.0;
    IbunI[1][1][1][1] = 1.0;
    IbunI[1][1][2][2] = 1.0;
    IbunI[2][2][0][0] = 1.0;
    IbunI[2][2][1][1] = 1.0;
    IbunI[2][2][2][2] = 1.0;

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
}

//full constructor
J2Damage::J2Damage(int tag, double _E, double _nu, // Parameters
    double _sig_y, double _Hk, double _Hi, // Plasticity
    double _Yt0, double _bt, double _at, double _Yc0, double _bc, double _ac, double _beta) // Damage
    : NDMaterial(tag, ND_TAG_J2Damage),
    E(_E), nu(_nu), sig_y(_sig_y), Hk(_Hk), Hi(_Hi),
    Yt0(_Yt0), bt(_bt), at(_at), Yc0(_Yc0), bc(_bc), ac(_ac), beta(_beta),
    epsilon_e(3,3),
    epsilon_p_n(3, 3),
    epsilon_p_nplus1(3, 3),
    stress(3, 3),
    strain(3, 3),
    parameterID(0)
{
    // Bulk and shear modulus
    K = E / (3 * (1 - 2 * nu));
    G = E / (2 * (1 + nu));

    K = K;
    G = G;
    sigma_infty = sig_y;
    delta = 0.0;
    H = Hk + Hi;
    eta = 0.0;
    rho = 0.0;
    
    if (0 == 1) {
        opserr << "Loaded full constructor. Material inputs: " << endln;
        opserr << "K = " << K << endln;
        opserr << "G = " << G << endln;
        opserr << "sig_y = " << sig_y << endln;
        opserr << "sig_inf = " << sig_y << endln;
        opserr << "delta = " << delta << endln;
        opserr << "H = " << H << endln;
        opserr << "eta = " << eta << endln;
        opserr << "rho = " << rho << endln;
    }
    
    this->zero();

    int i, j, k, l;

    //zero rank4 IIdev and IbunI 
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < 3; k++) {
                for (l = 0; l < 3; l++) {

                    IbunI[i][j][k][l] = 0.0;

                    IIdev[i][j][k][l] = 0.0;

                } // end for l
            } // end for k
        } // end for j
    } // end for i


    //form rank4 IbunI 

    IbunI[0][0][0][0] = 1.0;
    IbunI[0][0][1][1] = 1.0;
    IbunI[0][0][2][2] = 1.0;
    IbunI[1][1][0][0] = 1.0;
    IbunI[1][1][1][1] = 1.0;
    IbunI[1][1][2][2] = 1.0;
    IbunI[2][2][0][0] = 1.0;
    IbunI[2][2][1][1] = 1.0;
    IbunI[2][2][2][2] = 1.0;

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
}

//elastic constructor
J2Damage ::
J2Damage(int    tag,
    double _E,
    double _nu) :
    NDMaterial(tag, ND_TAG_J2Damage),
    epsilon_e(3,3),
    epsilon_p_n(3, 3),
    epsilon_p_nplus1(3, 3),
    stress(3, 3),
    strain(3, 3),
    parameterID(0)
{
    // Bulk and shear modulus
    K = E / (3 * (1 - 2 * nu));
    G = E / (2 * (1 + nu));

    sig_y = 1.0e16 * G;
    sigma_infty = sig_y;
    delta = 0.0;
    H = 0.0;
    eta = 0.0;

    this->zero();

    int i, j, k, l;

    //zero rank4 IIdev and IbunI 
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < 3; k++) {
                for (l = 0; l < 3; l++) {

                    IbunI[i][j][k][l] = 0.0;

                    IIdev[i][j][k][l] = 0.0;

                } // end for l
            } // end for k
        } // end for j
    } // end for i


    //form rank4 IbunI 

    IbunI[0][0][0][0] = 1.0;
    IbunI[0][0][1][1] = 1.0;
    IbunI[0][0][2][2] = 1.0;
    IbunI[1][1][0][0] = 1.0;
    IbunI[1][1][1][1] = 1.0;
    IbunI[1][1][2][2] = 1.0;
    IbunI[2][2][0][0] = 1.0;
    IbunI[2][2][1][1] = 1.0;
    IbunI[2][2][2][2] = 1.0;

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
}

//destructor
J2Damage :: ~J2Damage()
{  }

//get the strain and integrate plasticity equations
int J2Damage::setTrialStrain(const Vector& strain_from_element)
{
    // Strains in science notation
    strain_vec = sci(strain_from_element);

    strain.Zero();

    for(int i = 0;i<3;i++) strain(i, i) = strain_vec(i);
    strain(0, 1) = strain_vec(3);
    strain(1, 0) = strain_vec(3);
    strain(1, 2) = strain_vec(4);
    strain(2, 1) = strain_vec(4);
    strain(2, 0) = strain_vec(5);
    strain(0, 2) = strain_vec(5);

    // Elastic tangent
    this->getETangent();

    // Debug 1
    if (dFlag1 == 1) {
        opserr << "\n------------------------------------------------------------------------------------------------------------------------\n";
        opserr << "\nStarted new setTrialStrain.\n\n";
        opserr << "Inputs before plasticity and damage (initialized at current step):\n";
        opserr << "strain (eng) = [ "; for (int i = 0;i < 6;i++) opserr << strain_from_element(i) << " "; opserr << "]\n";
        opserr << "strain (sci) = [ "; for (int i = 0;i < 6;i++) opserr << strain_vec(i) << " "; opserr << "]\n";
        opserr << "tangent_e:\n[ "; for (int i = 0;i < 6;i++) { for (int j = 0;j < 6;j++) opserr << tangent_e(j, i) << " "; opserr << "]\n"; }
    }

    // Plasticity routine
    this->plastic_integrator();

    // Plasticity stress and tangent
    this->getEpTangent();
    stress_vec = vec(stress);
    tangent_matrix = tangent_ep;

    // Vectorize
    Vector strain_e(6);
    Vector strain_p(6);
    Vector strain_p_k(6);
    Vector stress_t(6);
    Vector Deps(6);
    Vector Deps_p(6);
    Vector Deps_e(6);
    strain_vec = eng(strain_vec);
    strain_e = eng(vec(epsilon_e));
    strain_p = eng(vec(epsilon_p_nplus1));
    strain_p_k = eng(vec(epsilon_p_n));
    stress_t = tangent_e * strain_e;
    Deps = strain_vec - strain_k;           // Deps  = (eps^n+1 - eps^n)
    Deps_p = strain_p - strain_p_k;         // Deps_p = (eps_p^n+1 - eps_p^n)
    Deps_e = Deps - Deps_p;                 // Deps_e = Deps - Deps_p

    // Debug 2
    if (dFlag1 == 1) {
        opserr << "\nPlasticity executed!\n\n";
        opserr << "Outputs after plasticity only:\n";
        opserr << "strain     = [ "; for (int i = 0;i < 6;i++) opserr << strain_vec(i) << " "; opserr << "]\n";
        opserr << "strain_e   = [ "; for (int i = 0;i < 6;i++) opserr << strain_e(i) << " "; opserr << "]\n";
        opserr << "strain_p   = [ "; for (int i = 0;i < 6;i++) opserr << strain_p(i) << " "; opserr << "]\n";
        opserr << "strain_p_k = [ "; for (int i = 0;i < 6;i++) opserr << strain_p_k(i) << " "; opserr << "]\n";
        opserr << "stress     = [ "; for (int i = 0;i < 6;i++) opserr << stress_vec(i) << " "; opserr << "]\n";
        opserr << "stress_t   = [ "; for (int i = 0;i < 6;i++) opserr << stress_t(i) << " "; opserr << "]\n";
        opserr << "tangent_ep:\n[ "; for (int i = 0;i < 6;i++) { for (int j = 0;j < 6;j++) opserr << tangent_ep(j, i) << " "; opserr << "]\n"; }
    }

    // Damage routine [local inputs: strain, strain_p]
    this->damage_integrator();
    double dD = D - D_k;

    // Damage stress and tangent
    //tangent_matrix = tangent_ep;
    //stress_vec = stress_k + tangent_ep *Deps;
    /*
    Vector CeDeps(6);
    Vector CeDeps_e(6);
    Vector CeEps_e(6);
    CeDeps = tangent_e * Deps;
    CeDeps_e = tangent_e * Deps_e;
    CeEps_e = tangent_e * strain_e;
    stress_vec.Zero();
    tangent_matrix = pow((1.0 - D),2.0) * tangent_e;
    for (int i = 0;i < 6;i++)
        stress_vec(i) = stress_k(i) + pow((1.0 - D), 2.0) * CeDeps_e(i);// -2.0 * (1.0 - D) * dD * CeEps_e(i);
    */
    /*
    opserr << "Cep*Deeng  = [ "; for (int i = 0;i < 6;i++) opserr << CepDepsEng(i) << " "; opserr << "]\n";
    opserr << "stress_k   = [ "; for (int i = 0;i < 6;i++) opserr << stress_k(i) << " "; opserr << "]\n";
    opserr << "stress     = [ "; for (int i = 0;i < 6;i++) opserr << stress_vec(i) << " "; opserr << "]\n";
    opserr << "sk + eng   = [ "; for (int i = 0;i < 6;i++) opserr << stress_k(i) + CepDepsEng(i) << " "; opserr << "]\n";
    */
    tangent_matrix = Dm1sq * tangent_e;
    stress_vec = Dm1sq * stress_vec;
    this->commitState();
    // Commits
    /*
    strain_k = strain_vec;
    stress_k = stress_vec;
    epsilon_p_n = epsilon_p_nplus1;
    xi_n = xi_nplus1;
    Dc_k = Dc;
    Dt_k = Dt;
    D_k = D;
    */

    // Debug 3
    if (dFlag1 == 1) {
        opserr << "\nDamage executed!\n\n";
        opserr << "Outputs after both plasticity and damage:\n";
        opserr << "D       = " << D << endln;
        opserr << "[1-D]^2 = " << Dm1sq << endln;
        opserr << "strain     = [ "; for (int i = 0;i < 6;i++) opserr << strain_from_element(i) << " "; opserr << "]\n";
        opserr << "strain_e   = [ "; for (int i = 0;i < 6;i++) opserr << strain_e(i) << " "; opserr << "]\n";
        opserr << "strain_p   = [ "; for (int i = 0;i < 6;i++) opserr << strain_p(i) << " "; opserr << "]\n";
        opserr << "strain_p_k = [ "; for (int i = 0;i < 6;i++) opserr << strain_p_k(i) << " "; opserr << "]\n";
        opserr << "stress     = [ "; for (int i = 0;i < 6;i++) opserr << stress_vec(i) << " "; opserr << "]\n";
        opserr << "tangent:\n[ "; for (int i = 0;i < 6;i++) { for (int j = 0;j < 6;j++) opserr << tangent_ep(j, i) << " "; opserr << "]\n"; }
    }

    return 0;
}

//--------------------Plasticity-------------------------------------

//plasticity integration routine
void J2Damage::plastic_integrator()
{
    const double tolerance = (1.0e-08) * sig_y;
    const double dt = ops_Dt; //time step
    static Matrix dev_strain(3, 3); //deviatoric strain
    static Matrix dev_stress(3, 3); //deviatoric stress
    static Matrix normal(3, 3);     //normal to yield surface
    double NbunN; //normal bun normal
    double norm_tau = 0.0;   //norm of deviatoric stress 
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
    const int max_iterations = 50;

    //compute the deviatoric strains
    trace = strain(0, 0) + strain(1, 1) + strain(2, 2);
    dev_strain = strain;
    for (i = 0; i < 3; i++) dev_strain(i, i) -= (one3 * trace);

    //compute the trial deviatoric stresses
    //   dev_stress = (2.0*G) * ( dev_strain - epsilon_p_n ) ;
    dev_stress = dev_strain;
    dev_stress -= epsilon_p_n;
    dev_stress *= 2.0 * G;

    //compute norm of deviatoric stress
    norm_tau = 0.0;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) norm_tau += dev_stress(i, j) * dev_stress(i, j);
    } //end for i 

    norm_tau = sqrt(norm_tau);

    if (norm_tau > tolerance) {
        inv_norm_tau = 1.0 / norm_tau;
        normal = inv_norm_tau * dev_stress;
    }
    else {
        normal.Zero();
        inv_norm_tau = 0.0;
    } //end if 

    //compute trial value of yield function
    phi = norm_tau - root23 * q(xi_n);

    // check if phi > 0 
    if (phi > 0.0) { //plastic

       //solve for gamma 
        gamma = 0.0;
        resid = 1.0;
        iteration_counter = 0;
        while (fabs(resid) > tolerance) {

            resid = norm_tau - (2.0 * G) * gamma - root23 * q(xi_n + root23 * gamma);
            tang = -(2.0 * G) - two3 * qprime(xi_n + root23 * gamma);

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
        epsilon_p_nplus1 = epsilon_p_n + gamma * normal;
        xi_nplus1 = xi_n + root23 * gamma;

        //recompute deviatoric stresses 
        dev_stress = (2.0 * G) * (dev_strain - epsilon_p_nplus1);

        //compute the terms for plastic part of tangent
        theta = (2.0 * G) + two3 * qprime(xi_nplus1);
        if (eta > 0.0 && dt > 0.0) theta += (eta / dt);
        theta_inv = 1.0 / theta;

    }
    else { //elastic 

      //update history variables -- they remain unchanged

        epsilon_p_nplus1 = epsilon_p_n;

        xi_nplus1 = xi_n;

        //no extra tangent terms to compute 

        gamma = 0.0;
        theta = 0.0;
        theta_inv = 0.0;

    } //end if phi > 0


    //add on K part of stress
    stress = dev_stress;
    for (i = 0; i < 3; i++)
        stress(i, i) += K * trace;

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
            tangent[i][j][k][l] = K * IbunI[i][j][k][l];

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

    return;
}

// set up for initial elastic
void J2Damage::doInitialTangent()
{
    int ii, jj, i, j, k, l;

    //compute the deviatoric strains
    for (ii = 0; ii < 6; ii++) {
        for (jj = 0; jj < 6; jj++) {

            index_map(ii, i, j);
            index_map(jj, k, l);

            //elastic terms
            initialTangent[i][j][k][l] = K * IbunI[i][j][k][l];
            initialTangent[i][j][k][l] += (2.0 * G) * IIdev[i][j][k][l];

            //minor symmetries 
            //minor symmetries 
            initialTangent[j][i][k][l] = initialTangent[i][j][k][l];
            initialTangent[i][j][l][k] = initialTangent[i][j][k][l];
            initialTangent[j][i][l][k] = initialTangent[i][j][k][l];

        } // end for jj
    } // end for ii

    return;
}

//Hardening function
double J2Damage::q(double xi)
{
    //  q(xi) = simga_infty + (sig_y - sigma_infty)*exp(-delta*xi) + H*xi 

    return    sigma_infty
        + (sig_y - sigma_infty) * exp(-delta * xi) + H * xi;
}

//Hardening function derivative
double J2Damage::qprime(double xi)
{
    return  (sig_y - sigma_infty) * (-delta) * exp(-delta * xi) + H;
}

//matrix_index ---> tensor indices i,j
void J2Damage::index_map(int matrix_index, int& i, int& j)
{
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

//--------------------Damage----------------------------------------

//damage integration routine
void J2Damage::damage_integrator()
{
    //////// 0. Principal strains calculation /////////////////////////////////////////////////////////////////
    // Elastic strain vector
    epsilon_e = strain - epsilon_p_nplus1;

    // Principal strains
    Vector strain_m(3);
    Vector strain_e_m(3);
    Matrix eVect(3, 3);		// Eigenvectors
    eVect.Eigen3(strain); for (int i = 0; i < 3;i++) strain_m(i) = eVect(i, i);
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
    double Yc = sqrt(e_tn ^ e_tn + beta * ((e_tn[0]) * (e_tn[1]) + (e_tn[1]) * (e_tn[2]) + (e_tn[2]) * (e_tn[0])));

    // Elastic equivalent tensile and compressive strains
    double Yt_e = sqrt(e_ep ^ e_ep);
    double Yc_e = sqrt(e_en ^ e_en + beta * ((e_en[0]) * (e_en[1]) + (e_en[1]) * (e_en[2]) + (e_en[2]) * (e_en[0])));

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

    // Boundary limits
    if (Dt < 0.0) Dt = 0.0; if (Dt > 1.0) Dt = 1.0;
    if (Dc < 0.0) Dc = 0.0; if (Dc > 1.0) Dc = 1.0;

    // Weighting coefficients
    double alpt; double alpc;
    if (Yt_e == 0.0 && Yc_e == 0.0) { alpt = 0.0; alpc = 0.0; }
    else { alpt = (fabs(Yt_e) / Yt0) / (fabs(Yt_e) / Yt0 + fabs(Yc_e) / Yc0); alpc = 1.0 - alpt; }

    // Cumulative damage variable Dm1sq = [1-D]^2
    D = alpt * Dt + alpc * Dc;
    Dm1sq = pow(((1.0 - Dt) * alpt + (1.0 - Dc) * alpc), 2.0);

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
        opserr << "alpha_c = " << alpc << "\n";
        opserr << "alpha_t = " << alpt << "\n";
        //opserr << "D = " << D << "\n";    -> Moved to setTrialStrain
        opserr << "\n--- End ------------------------------------\n";
    }
}

//--------------------Other methods----------------------------------

//rank2 tensor(3,3) to vector(6) notation
const Vector&
J2Damage::vec(const Matrix& m)
{
    for (int i = 0;i < 3;i++) v(i) = m(i, i);
    v(3) = m(0, 1); v(4) = m(1, 2); v(5) = m(2, 0);
    return v;
}

const Vector&
J2Damage::sci(const Vector& v) //science to engineering notation for strains
{
    for (int i = 0;i < 3;i++) s(i) = v(i);
    for (int i = 3;i < 6;i++) s(i) = v(i)/2.0;
    return s;
}

const Vector&
J2Damage::eng(const Vector& v) //science to engineering notation for strains
{
    for (int i = 0;i < 3;i++) s(i) = v(i);
    for (int i = 3;i < 6;i++) s(i) = v(i)*2.0;
    return s;
}

//unused trial strain functions
int J2Damage::setTrialStrain(const Vector& v, const Vector& r)
{
    return this->setTrialStrain(v);
}

int J2Damage::setTrialStrainIncr(const Vector& v)
{
    static Vector newStrain(6);
    newStrain(0) = strain(0, 0) + v(0);
    newStrain(1) = strain(1, 1) + v(1);
    newStrain(2) = strain(2, 2) + v(2);
    newStrain(3) = 2.0 * strain(0, 1) + v(3);
    newStrain(4) = 2.0 * strain(1, 2) + v(4);
    newStrain(5) = 2.0 * strain(2, 0) + v(5);

    return this->setTrialStrain(newStrain);
}

int J2Damage::setTrialStrainIncr(const Vector& v, const Vector& r)
{
    return this->setTrialStrainIncr(v);
}

//send back the strain
const Vector& J2Damage::getStrain()
{ return strain_vec; }

//send back the stress 
const Vector& J2Damage::getStress()
{ return stress_vec; }

//send back the tangent 
const Matrix& J2Damage::getInitialTangent()
{
    // matrix to tensor mapping
    //  Matrix      Tensor
    // -------     -------
    //   0           0 0
    //   1           1 1
    //   2           2 2   
    //   3           0 1  ( or 1 0 )
    //   4           1 2  ( or 2 1 )
    //   5           2 0  ( or 0 2 ) 

    int ii, jj;
    int i, j, k, l;

    this->doInitialTangent();

    for (ii = 0; ii < 6; ii++) {
        for (jj = 0; jj < 6; jj++) {

            index_map(ii, i, j);
            index_map(jj, k, l);

            tangent_matrix(ii, jj) = initialTangent[i][j][k][l];

        } //end for j
    } //end for i

    return tangent_matrix;
}

//send back the elastic tangent 
const Matrix& J2Damage::getETangent()
{
    // matrix to tensor mapping
    //  Matrix      Tensor
    // -------     -------
    //   0           0 0
    //   1           1 1
    //   2           2 2   
    //   3           0 1  ( or 1 0 )
    //   4           1 2  ( or 2 1 )
    //   5           2 0  ( or 0 2 ) 

    int ii, jj;
    int i, j, k, l;

    this->doInitialTangent();

    for (ii = 0; ii < 6; ii++) {
        for (jj = 0; jj < 6; jj++) {

            index_map(ii, i, j);
            index_map(jj, k, l);

            tangent_e(ii, jj) = initialTangent[i][j][k][l];

        } //end for j
    } //end for i


    return tangent_e;
}

//send back the elastoplastic tangent 
const Matrix& J2Damage::getEpTangent()
{
    // matrix to tensor mapping
    //  Matrix      Tensor
    // -------     -------
    //   0           0 0
    //   1           1 1
    //   2           2 2   
    //   3           0 1  ( or 1 0 )
    //   4           1 2  ( or 2 1 )
    //   5           2 0  ( or 0 2 ) 

    int ii, jj;
    int i, j, k, l;

    for (ii = 0; ii < 6; ii++) {
        for (jj = 0; jj < 6; jj++) {

            index_map(ii, i, j);
            index_map(jj, k, l);

            tangent_ep(ii, jj) = tangent[i][j][k][l];

        } //end for j
    } //end for i


    return tangent_ep;
}

//send back the tangent 
const Matrix& J2Damage::getTangent()
{ return tangent_matrix; }

int
J2Damage::commitState()
{
    // Commits
    strain_k = strain_vec;
    stress_k = stress_vec;

    // Plasticity
    epsilon_p_n = epsilon_p_nplus1;
    xi_n = xi_nplus1;

    // Damage
    Dc_k = Dc;
    Dt_k = Dt;
    D_k = D;

    return 0;
}

int
J2Damage::revertToLastCommit()
{
    // Commits
    strain_vec = strain_k;
    stress_vec = stress_k;

    // Plasticity
    epsilon_p_nplus1 = epsilon_p_n;
    xi_nplus1 = xi_n;

    // Damage
    Dc = Dc_k;
    Dt = Dt_k;
    D = D_k;

    return 0;
}

int
J2Damage::revertToStart() {

    // added: C.McGann, U.Washington for InitialStateAnalysis
    if (ops_InitialStateAnalysis) {
        // do nothing, keep state variables from last step
    }
    else {
        // normal call for revertToStart (not initialStateAnalysis)
        this->zero();
    }

    return 0;
}

int
J2Damage::setParameter(const char** argv, int argc,
    Parameter& param)
{
    if (strcmp(argv[0], "K") == 0)
        return param.addObject(1, this);

    else if (strcmp(argv[0], "G") == 0 || strcmp(argv[0], "mu") == 0)
        return param.addObject(2, this);

    else if (strcmp(argv[0], "rho") == 0)
        return param.addObject(3, this);

    return -1;
}

int
J2Damage::updateParameter(int parameterID, Information& info)
{
    switch (parameterID) {
    case 1:
        K = info.theDouble;
        return 0;
    case 2:
        G = info.theDouble;
        return 0;
    case 3:
        rho = info.theDouble;
        return 0;
    default:
        return -1;
    }
}

int
J2Damage::activateParameter(int paramID)
{
    parameterID = paramID;

    return 0;
}

int
J2Damage::sendSelf(int commitTag, Channel& theChannel)
{
    // we place all the data needed to define material and it's state
    // int a vector object
    static Vector data(10 + 9);
    int cnt = 0;
    data(cnt++) = this->getTag();
    data(cnt++) = K;
    data(cnt++) = G;
    data(cnt++) = sig_y;
    data(cnt++) = sigma_infty;
    data(cnt++) = delta;
    data(cnt++) = H;
    data(cnt++) = eta;
    data(cnt++) = rho;

    data(cnt++) = xi_n;

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            data(cnt++) = epsilon_p_n(i, j);


    // send the vector object to the channel
    if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
        opserr << "J2Damage::sendSelf - failed to send vector to channel\n";
        return -1;
    }

    return 0;
}

int
J2Damage::recvSelf(int commitTag, Channel& theChannel,
    FEM_ObjectBroker& theBroker)
{
    // recv the vector object from the channel which defines material param and state
    static Vector data(10 + 9);
    if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
        opserr << "J2Damage::recvSelf - failed to recv vector from channel\n";
        return -1;
    }

    // set the material parameters and state variables
    int cnt = 0;
    this->setTag(data(cnt++));
    K = data(cnt++);
    G = data(cnt++);
    sig_y = data(cnt++);
    sigma_infty = data(cnt++);
    delta = data(cnt++);
    H = data(cnt++);
    eta = data(cnt++);
    rho = data(cnt++);

    xi_n = data(cnt++);

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            epsilon_p_n(i, j) = data(cnt++);

    epsilon_p_nplus1 = epsilon_p_n;
    xi_nplus1 = xi_n;

    return 0;
}

//make a clone of this material
NDMaterial*
J2Damage::getCopy()
{
    J2Damage * clone;
    clone = new J2Damage(this->getTag(), E, nu, sig_y, Hk, Hi,
        Yt0, bt, at, Yc0, bc, ac, beta);
    return clone;
}

NDMaterial*
J2Damage::getCopy(const char* type)
{
    if (strcmp(type, this->getType()) == 0)
        return this->getCopy();
    return 0;
}

//print out material data
void J2Damage::Print(OPS_Stream& s, int flag)
{
    s << endln;
    s << "J2-Plasticity : ";
    s << this->getType() << endln;
    s << "K Modulus =   " << K << endln;
    s << "G Modulus =  " << G << endln;
    s << "sig_y =        " << sig_y << endln;
    s << "Sigma_infty =    " << sigma_infty << endln;
    s << "Delta =          " << delta << endln;
    s << "H =              " << H << endln;
    s << "Eta =            " << eta << endln;
    s << "Rho =            " << rho << endln;
    s << endln;
}

//send back type of material
const char* J2Damage::getType() const
{
    return "ThreeDimensional";
}

//send back order of strain in vector form
int J2Damage::getOrder() const
{
    return 6;
}