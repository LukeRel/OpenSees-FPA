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
    J2Damage(int tag, double _E, double _nu, // Parameters
        double _sig_y, double _Hk, double _Hi, // Plasticity
        double _Yt0, double _bt, double _at, double _Yc0, double _bc, double _ac, double _beta); // Damage

    //elastic constructor
    J2Damage(int tag, double _E, double _nu);

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
    const Matrix& getInitialTangent();
    const Matrix& getETangent();
    const Matrix& getEpTangent();
    const Matrix& getTangent();

protected:

    //internal variables
    Matrix epsilon_e;         // elastic strains time n+1
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

    //Integration routines
    void plastic_integrator();          //plasticity integration routine
    void damage_integrator();           //damage integration routine

    // Plasticity functions
    double q(double xi);                //hardening function
    double qprime(double xi);           //hardening function derivative

    //matrix index to tensor index mapping

    // Utilities
    void index_map(int matrix_index, int& i, int& j);
    void doInitialTangent();
    const Vector& vec(const Matrix& m); //tensor to vector notation for strains and stresses
    const Vector& eng(const Vector& v); //science to engineering notation for strains
    const Vector& sci(const Vector& v); //engineering to science notation for strains

    double rho;
    int parameterID;

private:
    // Internal parameters
    int ndm;
    double E;       // Elastic modulus
    double nu;      // Poisson ratio 
    double K;		// Bulk modulus
    double G;		// Shear modulus
    double sig_y;   // Yield strength
    double Hk;      // Kinematic hardening coefficient
    double Hi;      // Isotropic hardening coefficient
    double Yt0;     // Initial generalized tensile damage variable
    double bt;		// Tensile shape factor
    double at;		// Tensile shape factor
    double Yc0;     // Initial generalized compressive damage variable
    double bc;		// Compression shape factor
    double ac;		// Compression shape factor
    double beta;	// Surface shape factor
    double sigma_infty; //final saturation yield stress
    double delta;       //exponential hardening parameter
    double H;        //linear hardening parameter
    double eta;         //viscosity

    // Static vectors and matrices
    static Vector s;              //general size 6 vector
    static Vector v;              //general size 6 vector
    static Matrix m;              //general size 6 matrix
    static Vector strain_vec;     //strain in vector notation (k+1)
    static Vector stress_vec;     //stress in vector notation (k+1)
    static Matrix tangent_matrix; //material tangent in matrix notation
    static Matrix tangent_e;      //elastic tangent in matrix notation
    static Matrix tangent_ep;     //elastoplastic tangent in matrix notation

    // Committed stuff
    static Vector strain_k;       //strain in vector notation (k)
    static Vector stress_k;       //stress in vector notation (k)

    // Current damage variables
    double Dt;					// Tensile damage at step (k+1)
    double Dc;					// Compressive damage at step (k+1)
    double D;					// Total damage (k+1)
    double Dm1sq;				// [1-D]^2

    // Committed damage variables
    double Dt_k;				// Tensile damage at step (k)
    double Dc_k;				// Compressive damage at step (k)
    double D_k;					// Total damage at previous step (k)

}; //end of J2Damage declarations


#endif
