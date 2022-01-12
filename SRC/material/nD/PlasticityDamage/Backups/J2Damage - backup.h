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
                                                                        
// $Revision: 1.6 $
// $Date: 2006-08-04 18:18:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticIsotropicPlaneStress2D.h,v $

#ifndef J2Damage_h
#define J2Damage_h

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

#include <ElasticIsotropicMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class J2Damage : public NDMaterial
{
  public:

    // Full constructor
    J2Damage(int _tag, double _E, double _nu, // Parameters
          double _ft, double _fc, double _Hk, double _Hi, // Plasticity
          double _Yt0, double _bt, double _at, double _Yc0, double _bc, double _ac, double _beta); // Damage

    // Null constructor
    J2Damage();

    // Destructor
    ~J2Damage();

    // Class type return
    const char *getClassType(void) const {return "J2Damage";};

    // Main methods
    int setTrialStrain (const Vector &v);
    int setTrialStrain(const Vector& v, const Vector& r);

    // Constitutive matrix, stress and strain determination
    const Matrix &getTangent (void);
    const Matrix &getInitialTangent (void);
    const Vector &getStress (void);
    const Vector &getStrain (void);
    
    // Commit updating
    int commitState (void);
    int revertToLastCommit (void);
    int revertToStart (void);

    // Additional functions
    NDMaterial *getCopy(const char *type);
    NDMaterial *getCopy(void);
    const char *getType(void) const;
    int getOrder(void) const;

    // Optional parallel workflow
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    

    // Output flag method
    void Print(OPS_Stream &s, int flag =0);   

    /*

    Response* setResponse(const char** argv, int argc, OPS_Stream& output);
    int getResponse(int responseID, Information& matInformation);

    */

    // Fill vector of state variables for output
    Vector getState();

  protected:

    // Plasticity integration routine
    void plastic_integrator(void);
    void index_map(int matrix_index, int& i, int& j);

    Vector strainVectorMap(Matrix E);
    Matrix strainTensorMap(Vector e);

    // Damage integration routine
    void damage_integrator(void);

    // Initialization method
    void initialize(void);

    // Principal values
    Vector eigen3(const Vector e);

    // Internal plasticity functions
    double qiso(double alpha1);		    // Isotropic hardening function: q = Hi*alpha
    double q(double alpha1);		    // Hardening function: q = H*alpha
    double qprime(double alpha1);		// Hardening function derivate: H

    // Material parameters
    double K;			    // bulk modulus 
    double G;		     	// shear modulus
    double sigma_y;		// yield strength 
    double mu;			// volumetric term
    double H;			// hardening constant
    double theta;		    // hardening constant

    // Internal parameters
    double E;       // Elastic modulus
    double nu;      // Poisson ratio 
    double ft;      // Tensile yield strength
    double fc;      // Compressive yield strength
    double Hk;      // Kinematic hardening coefficient
    double Hi;      // Isotropic hardening coefficient
    double Yt0;     // Initial generalized tensile damage variable
    double bt;  
    double at;
    double Yc0;     // Initial generalized compressive damage variable
    double bc;
    double ac;
    double beta;

    // Damage variables
    double Dt_n;            // Tensile damage at step n
    double Dc_n;            // Compressive damage at step n
    double Dt_n1;           // Tensile damage at step n+1
    double Dc_n1;           // Compressive damage at step n+1

    // Rank 2 tensor notation
    Matrix Eps;             // strain tensor
    Matrix Eps_e;           // strain tensor
    Matrix Sig;             // stress tensor
    Matrix Eps_n_p;	        // Plastic strains vector at step n, trial e_p
    Matrix Eps_n1_p;	    // Plastic strains vector at step n+1

    // Rank 4 tensor notation
    double tangent[3][3][3][3];   //material tangent
    static double IIdev[3][3][3][3]; //rank 4 deviatoric 
    static double IIvol[3][3][3][3]; //rank 4 I bun I

    // Back stress coefficients
    double alpha_n;		    // First hardening isotropic variable at step n
    double alpha_n1;	    // First hardening isotropic variable at step n+1

    // Flags
    int ElastFlag;          // Flag to determine elastic behavior
    int Flag;

    // Constitutive matrixes
    Matrix Ct;              // Tangent elastoplastic matrix
    Matrix Ce;			    // elastic tangent stiffness matrix
    Vector I1;			    // 2nd Order Identity Tensor

    Vector State;		    // state vector for output

    // Constants
    static const double one3;
    static const double two3;
    static const double root23;

  private:
    // Vector notation
    static Vector eps;             // Total strains vector
    static Vector sig;             // Total stress vector
    static Matrix Cep;             // Elastoplastic tangent matrix

    Vector eps_e;           // Elastic strains vector
    Vector eps_n_p;	        // Plastic strains vector at step n, trial e_p
    Vector eps_n1_p;	    // Plastic strains vector at step n+1 

};

#endif
