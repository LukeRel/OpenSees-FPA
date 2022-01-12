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

// Written: Ed "C++" Love
// Do not ask Prashant about this code.  He has no clue. 
//
// J2Damage isotropic hardening material class
// 
//  Elastic Model
//  sigma = K*trace(epsilion_elastic) + (2*G)*dev(epsilon_elastic)
//
//  Yield Function
//  phi(sigma,q) = || dev(sigma) ||  - sqrt(2/3)*q(xi) 
//
//  Saturation Isotropic Hardening with linear term
//  q(xi) = simga_0 + (sigma_infty - sigma_0)*exp(-delta*xi) + H*xi 
//
//  Flow Rules
//  \dot{epsilon_p} =  gamma * d_phi/d_sigma
//  \dot{xi}        = -gamma * d_phi/d_q 
//
//  Linear Viscosity 
//  gamma = phi / eta  ( if phi > 0 ) 
//
//  Backward Euler Integration Routine 
//  Yield condition enforced at time n+1 
//
//  Send strains in following format :
// 
//     strain_vec = {   eps_00
//                      eps_11
//	                eps_22 		      
//                    2 eps_01   
//            	      2 eps_12   
//		      2 eps_20    }   <--- note the 2
// 
//  set eta := 0 for rate independent case
//

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
Vector J2Damage::strain_vec(6);
Vector J2Damage::stress_vec(6);
Matrix J2Damage::tangent_matrix(6, 6);

void* OPS_J2Damage()
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 7) {
        opserr << "WARNING: Insufficient arguments\n";
        opserr << "Want: nDMaterial J2Damage tag? K? G? sig0? sigInf? delta? H? <eta?>\n";
        return 0;
    }

    int tag;
    numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
        opserr << "WARNING invalid J2Damage tag\n";
        return 0;
    }

    double data[8] = { 0,0,0,0,0,0,0,0 };
    numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 7) {
        numdata = 7;
    }
    if (OPS_GetDoubleInput(&numdata, data)) {
        opserr << "WARNING invalid J2Damage double inputs\n";
        return 0;
    }

    NDMaterial* mat = new J2Damage(tag, data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7]);
    if (mat == 0) {
        opserr << "WARNING: failed to create J2Damage material\n";
        return 0;
    }

    return mat;
}

//zero internal variables
void J2Damage::zero()
{
    xi_n = 0.0;
    xi_nplus1 = 0.0;

    epsilon_p_n.Zero();
    epsilon_p_nplus1.Zero();

    stress.Zero();
    strain.Zero();
}

//null constructor
J2Damage::J2Damage() :
    NDMaterial(0, ND_TAG_J2Damage),
    epsilon_p_n(3, 3),
    epsilon_p_nplus1(3, 3),
    stress(3, 3),
    strain(3, 3),
    parameterID(0)
{
    bulk = 0.0;
    shear = 0.0;
    sigma_0 = 0.0;
    sigma_infty = 0.0;
    delta = 0.0;
    Hard = 0.0;
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
J2Damage::J2Damage(int    tag,
    double K,
    double G,
    double yield0,
    double yield_infty,
    double d,
    double H,
    double viscosity,
    double r)
    :
    NDMaterial(tag, ND_TAG_J2Damage),
    epsilon_p_n(3, 3),
    epsilon_p_nplus1(3, 3),
    stress(3, 3),
    strain(3, 3),
    parameterID(0)
{
    bulk = K;
    shear = G;
    sigma_0 = yield0;
    sigma_infty = yield_infty;
    delta = d;
    Hard = H;
    eta = viscosity;
    rho = r;
    
    if (0 == 1) {
        opserr << "Loaded full constructor. Material inputs: " << endln;
        opserr << "K = " << K << endln;
        opserr << "G = " << G << endln;
        opserr << "sig_y = " << yield0 << endln;
        opserr << "sig_inf = " << yield_infty << endln;
        opserr << "delta = " << d << endln;
        opserr << "Hard = " << H << endln;
        opserr << "eta = " << viscosity << endln;
        opserr << "rho = " << r << endln;
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
    double K,
    double G) :
    NDMaterial(tag, ND_TAG_J2Damage),
    epsilon_p_n(3, 3),
    epsilon_p_nplus1(3, 3),
    stress(3, 3),
    strain(3, 3),
    parameterID(0)
{
    bulk = K;
    shear = G;
    sigma_0 = 1.0e16 * shear;
    sigma_infty = sigma_0;
    delta = 0.0;
    Hard = 0.0;
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
    strain.Zero();

    strain(0, 0) = strain_from_element(0);
    strain(1, 1) = strain_from_element(1);
    strain(2, 2) = strain_from_element(2);

    strain(0, 1) = 0.50 * strain_from_element(3);
    strain(1, 0) = strain(0, 1);

    strain(1, 2) = 0.50 * strain_from_element(4);
    strain(2, 1) = strain(1, 2);

    strain(2, 0) = 0.50 * strain_from_element(5);
    strain(0, 2) = strain(2, 0);


    this->plastic_integrator();

    this->getStress();
    this->getTangent();

    if (0 == 1) {
        opserr << "strain = [ "; for (int i = 0;i < 6;i++) opserr << strain_from_element(i) << " "; opserr << "]\n";
        opserr << "stress = [ "; for (int i = 0;i < 6;i++) opserr << stress_vec(i) << " "; opserr << "]\n";
        opserr << "tangent:\n[ "; for (int i = 0;i < 6;i++) { for (int j = 0;j < 6;j++) opserr << tangent_matrix(j, i) << " "; opserr << "]\n"; }
    }

    return 0;
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
{
    strain_vec(0) = strain(0, 0);
    strain_vec(1) = strain(1, 1);
    strain_vec(2) = strain(2, 2);

    strain_vec(3) = 2.0 * strain(0, 1);
    strain_vec(4) = 2.0 * strain(1, 2);
    strain_vec(5) = 2.0 * strain(2, 0);

    return strain_vec;
}

//send back the stress 
const Vector& J2Damage::getStress()
{
    stress_vec(0) = stress(0, 0);
    stress_vec(1) = stress(1, 1);
    stress_vec(2) = stress(2, 2);

    stress_vec(3) = stress(0, 1);
    stress_vec(4) = stress(1, 2);
    stress_vec(5) = stress(2, 0);

    return stress_vec;
}

//send back the tangent 
const Matrix& J2Damage::getTangent()
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

            tangent_matrix(ii, jj) = tangent[i][j][k][l];

        } //end for j
    } //end for i


    return tangent_matrix;
}

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

//--------------------Plasticity-------------------------------------

//plasticity integration routine
void J2Damage::plastic_integrator()
{
    const double tolerance = (1.0e-8) * sigma_0;

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
    const int max_iterations = 25;

    //compute the deviatoric strains

    trace = strain(0, 0) + strain(1, 1) + strain(2, 2);

    dev_strain = strain;
    for (i = 0; i < 3; i++)
        dev_strain(i, i) -= (one3 * trace);

    //compute the trial deviatoric stresses

    //   dev_stress = (2.0*shear) * ( dev_strain - epsilon_p_n ) ;
    dev_stress = dev_strain;
    dev_stress -= epsilon_p_n;
    dev_stress *= 2.0 * shear;

    //compute norm of deviatoric stress

    norm_tau = 0.0;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++)
            norm_tau += dev_stress(i, j) * dev_stress(i, j);
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

            resid = norm_tau
                - (2.0 * shear) * gamma
                - root23 * q(xi_n + root23 * gamma);
            if (eta > 0.0 && dt > 0.0)
                resid -= (eta / dt) * gamma;

            tang = -(2.0 * shear)
                - two3 * qprime(xi_n + root23 * gamma);
            if (eta > 0.0 && dt > 0.0)
                tang -= (eta / dt);

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

        dev_stress = (2.0 * shear) * (dev_strain - epsilon_p_nplus1);

        //compute the terms for plastic part of tangent

        theta = (2.0 * shear)
            + two3 * qprime(xi_nplus1);
        if (eta > 0.0 && dt > 0.0)
            theta += (eta / dt);

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


    //add on bulk part of stress

    stress = dev_stress;
    for (i = 0; i < 3; i++)
        stress(i, i) += bulk * trace;

    //compute the tangent

    c1 = -4.0 * shear * shear;
    c2 = c1 * theta_inv;
    c3 = c1 * gamma * inv_norm_tau;

    for (ii = 0; ii < 6; ii++) {
        for (jj = 0; jj < 6; jj++) {

            index_map(ii, i, j);
            index_map(jj, k, l);

            NbunN = normal(i, j) * normal(k, l);

            //elastic terms
            tangent[i][j][k][l] = bulk * IbunI[i][j][k][l];

            tangent[i][j][k][l] += (2.0 * shear) * IIdev[i][j][k][l];

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
            initialTangent[i][j][k][l] = bulk * IbunI[i][j][k][l];
            initialTangent[i][j][k][l] += (2.0 * shear) * IIdev[i][j][k][l];

            //minor symmetries 
            //minor symmetries 
            initialTangent[j][i][k][l] = initialTangent[i][j][k][l];
            initialTangent[i][j][l][k] = initialTangent[i][j][k][l];
            initialTangent[j][i][l][k] = initialTangent[i][j][k][l];

        } // end for jj
    } // end for ii

    return;
}

//hardening function
double J2Damage::q(double xi)
{
    //  q(xi) = simga_infty + (sigma_0 - sigma_infty)*exp(-delta*xi) + H*xi 

    return    sigma_infty
        + (sigma_0 - sigma_infty) * exp(-delta * xi)
        + Hard * xi;
}

//hardening function derivative
double J2Damage::qprime(double xi)
{
    return  (sigma_0 - sigma_infty) * (-delta) * exp(-delta * xi)
        + Hard;
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

int
J2Damage::commitState()
{
    epsilon_p_n = epsilon_p_nplus1;
    xi_n = xi_nplus1;

    return 0;
}

int
J2Damage::revertToLastCommit()
{
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
        bulk = info.theDouble;
        return 0;
    case 2:
        shear = info.theDouble;
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
    data(cnt++) = bulk;
    data(cnt++) = shear;
    data(cnt++) = sigma_0;
    data(cnt++) = sigma_infty;
    data(cnt++) = delta;
    data(cnt++) = Hard;
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
    bulk = data(cnt++);
    shear = data(cnt++);
    sigma_0 = data(cnt++);
    sigma_infty = data(cnt++);
    delta = data(cnt++);
    Hard = data(cnt++);
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
    clone = new J2Damage(this->getTag(), bulk, shear, sigma_0, sigma_infty, delta, Hard, eta, rho);
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
    s << "Bulk Modulus =   " << bulk << endln;
    s << "Shear Modulus =  " << shear << endln;
    s << "Sigma_0 =        " << sigma_0 << endln;
    s << "Sigma_infty =    " << sigma_infty << endln;
    s << "Delta =          " << delta << endln;
    s << "H =              " << Hard << endln;
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