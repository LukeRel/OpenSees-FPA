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

#ifndef DPDamage_h
#define DPDamage_h

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

#include <ElasticIsotropicMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class DPDamage : public NDMaterial
{
  public:

    // Full constructor
    DPDamage(int tag, double E, double nu, // Parameters
          double ft, double fc, double Hk, double Hi, // Plasticity
          double Yt0, double bt, double at, double Yc0, double bc, double ac, double beta); // Damage

    // Null constructor
    DPDamage();

    // Destructor
    ~DPDamage();

    // Class type return
    const char *getClassType(void) const {return "DPDamage";};

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

    // Plasticity integration routine
    void plastic_integrator(void);

    // Damage integration routine
    void damage_integrator(void);

    // Initialization method
    void initialize(void);

    // Principal values
    Vector principal_values(Vector e);
    Vector Eigen3(const Vector e);

    // Internal functions
    double qiso(double alpha1);		    // Isotropic hardening function: q = Hi*alpha

    // Fill vector of state variables for output
    Vector getState();

  protected:

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

    // Internal variables
    Vector Eps;          // Total strains vector
    Vector Eps_e;        // Elastic strains vector
    Vector Eps_n_p;	     // Plastic strains vector at step n, trial e_p
    Vector Eps_n1_p;	 // Plastic strains vector at step n+1 
    Vector Sig;          // Tension vector
     
    Vector Zeta_n;		    // Backstress at step n, beta_np1_trial = beta_n
    Vector Zeta_n1;		    // Backstress at step n+1

    double Alpha1_n;		// First hardening isotropic variable at step n
    double Alpha1_n1;	    // First hardening isotropic variable at step n+1
    double Alpha2_n;		// Second hardening isotropic variable at step n
    double Alpha2_n1;	    // Second hardening isotropic variable at step n+1

    int ElastFlag;         // Flag to determine elastic behavior
    int Flag;

    Matrix C;               // Tangent elastoplastic matrix
    Matrix Ce;			    // elastic tangent stiffness matrix
    Matrix Cep;		    	// elastoplastic tangent stiffness matrix
    Vector I1;			    // 2nd Order Identity Tensor	
    Matrix IIvol;		    // IIvol = I1 tensor I1  
    Matrix IIdev;		    // 4th Order Deviatoric Tensor

    Vector State;		    // state vector for output

    // Constants
    static const double one3;
    static const double two3;
    static const double root23;

  private:

};

#endif
