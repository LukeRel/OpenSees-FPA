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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/J2Damage.h,v $

// Written: Ed "C++" Love

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

#ifndef J2Damage_h
#define J2Damage_h

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>
#include <NDMaterial.h>

class J2Damage : public NDMaterial {

    //-------------------Declarations-------------------------------

public:

    //null constructor
    J2Damage();

    //full constructor
    J2Damage(int    tag,
        double K,
        double G,
        double yield0,
        double yield_infty,
        double d,
        double H,
        double viscosity = 0,
        double rho = 0.0);

    //elastic constructor
    J2Damage(int tag, double K, double G);

    //destructor
    virtual ~J2Damage();

    const char* getClassType(void) const { return "J2Damage"; };

    //swap history variables
    int commitState();

    //revert to last saved state
    int revertToLastCommit();

    //revert to start
    int revertToStart();

    //sending and receiving
    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel,
        FEM_ObjectBroker& theBroker);

    //print out material data
    void Print(OPS_Stream& s, int flag = 0);

    NDMaterial* getCopy(const char* type);
    NDMaterial* getCopy(void);
    const char* getType(void) const;
    int getOrder(void) const;

    double getRho(void) { return rho; }

    int setParameter(const char** argv, int argc, Parameter& param);
    int updateParameter(int parameterID, Information& info);
    int activateParameter(int paramID);

    //get the strain and integrate plasticity equations
    int setTrialStrain(const Vector& strain_from_element);

    //unused trial strain functions
    int setTrialStrain(const Vector& v, const Vector& r);
    int setTrialStrainIncr(const Vector& v);
    int setTrialStrainIncr(const Vector& v, const Vector& r);

    //send back the strain
    const Vector& getStrain();

    //send back the stress 
    const Vector& getStress();

    //send back the tangent 
    const Matrix& getTangent();
    const Matrix& getInitialTangent();

protected:

    //material parameters
    double bulk;        //bulk modulus
    double shear;       //shear modulus
    double sigma_0;     //initial yield stress
    double sigma_infty; //final saturation yield stress
    double delta;       //exponential hardening parameter
    double Hard;        //linear hardening parameter
    double eta;         //viscosity

    //internal variables
    Matrix epsilon_p_n;       // plastic strain time n
    Matrix epsilon_p_nplus1;  // plastic strain time n+1
    double xi_n;              // xi time n
    double xi_nplus1;         // xi time n+1

    //material response 
    Matrix stress;                //stress tensor
    double tangent[3][3][3][3];   //material tangent
    static double initialTangent[3][3][3][3];   //material tangent
    static double IIdev[3][3][3][3]; //rank 4 deviatoric 
    static double IbunI[3][3][3][3]; //rank 4 I bun I 

    //material input
    Matrix strain;               //strain tensor

    //parameters
    static const double one3;
    static const double two3;
    static const double four3;
    static const double root23;

    //zero internal variables
    void zero();

    //plasticity integration routine
    void plastic_integrator();

    void doInitialTangent();

    //hardening function
    double q(double xi);

    //hardening function derivative
    double qprime(double xi);

    //matrix index to tensor index mapping
    void index_map(int matrix_index, int& i, int& j);

    double rho;

    int parameterID;

private:

    //static vectors and matrices
    static Vector strain_vec;     //strain in vector notation
    static Vector stress_vec;     //stress in vector notation
    static Matrix tangent_matrix; //material tangent in matrix notation

}; //end of J2Damage declarations


#endif
