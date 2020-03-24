#ifndef CONSTANTS_H
#define CONSTANTS_H

//#define NOISE_ON
//#define RAND_FLAG (100)
//#define BIN_FILE
//#define FROM_INPUT_FILE
//#define DISLOCATION_ACTIVITY_ON
//#define ISOTROPIC_ELASTICITY
//#define DEBUG_ON
//#define INCLUDE_MISFIT_WORK_TERM

#define TOTAL_TIMESTEP (2)
#define SAVE_CONFIG (2)
#define dt (0.177)
//System Size
#define L (128)
#define M (128)
#define N (128)


#define MISFIT_STRAIN (-0.003)

#ifdef DISLOCATION_ACTIVITY_ON
	//No of slip modes
	#define P (8)
#endif

//No of Variants of the precipitate
#define V (4)

//Equilibrium concentration of gamma phase
#define cg_e (0.160)

//Equilibrium concentration of gamma prime phase
#define cgp_e (0.229)

//Equilibrium lattice parameter
//#define alatt (1)

//Grid size
//#define lo (20e-9)
#define dx (20e-9)
#define dy (20e-9)
#define dz (20e-9)

//Atomic mobility
#define Mc (1.544e-16) // mol^2/(Jms)
#define interface_omega (3.8996e6) // J/m^3 

#define Lphi (5.79e-9) //m^2/N.s
//#define Leta (2.46e-6) //m^2/N.s
#define Leta (1.4486e-11) //tuned value

//#define keta (7.8e-4)  //J/m
  #define keta (7.8e-6) //tuned value

#define kphi (9.36e-10) //J/m

// fo in
#define fo (3.2175e4) // in J./mol

//Molar volume
#define Vm (1e-5) // in m^3/mol

#ifdef  ISOTROPIC_ELASTICITY
	#define G 74.8e9  // in Pascal
	#define v_poisson 0.276
#endif
	#define Cijkla (399.19e9)  // C11
	#define Cijklb (64.9e9)  // C12
	#define Cijklc (117e9)  // C44


#define alpha_cross (10)
//alpha_cross is the theta value 

#define s_app_mag (0)


//Basic constants
#define SQRT_SIX (2.449489742783178)


#endif


