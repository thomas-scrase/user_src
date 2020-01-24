/* 
	RICEMYOFILAMENTCONSTANTS_HPP
	Constants variables for Rice et al. 2008 myofilament model
	Some parameters have been tuned to match experimental data;

	Haibo Ni
	qiangzi.ni@gmail.com
*/


#ifndef RICEMYOFILAMENTCONSTANTS_HPP
#define RICEMYOFILAMENTCONSTANTS_HPP


#define p0 1.754251548964904
#define p1 0.905622641626625
#define p2 0.499437793063966
#define p3 0.400000001127317
#define p4 1.000000000000000
#define p5 1.000000000000000
#define p6 1.981229252026256
#define p7 0.511387864324941
#define p8 0.023420000000000
// rolling back to its original value.
static const double SLmax   = 2.4;// belus 2010. fig6, was 2.4;        //   (um) maximum sarcomere length
static const double SLmin   = 1.4;//1.4;        //   (um) minimum sarcomere length
static const double len_thin  = 1.2;//1.2;      //   (um) thin filament length
static const double len_thick = 1.65; // //1.65;     //   (um) thick filament length
static const double len_hbare = 0.1;      //   (um) length of bare portion of thick filament
//   Temperature Dependence
static const double Qkon = 1.5;
static const double Qkoff = 1.3;
// static const double Qkoff = 1.4;
static const double Qkn_p = 1.6;
static const double Qkp_n = 1.6;
static const double Qfapp = 6.25;
// static const double Qgapp = 6.25;
static const double Qgapp = 2.5;
static const double Qhf = 6.25;
static const double Qhb = 6.25;
static const double Qgxb = 6.25;

//   Ca binding to troponin
static const double kon     = p0 * 50/*e-3*/;    //   (1/[ms uM])
static const double koffL   = p1 * p0 * 250e-3; //   (1/ms)
static const double koffH   = p1 * p0 * 25e-3;  //   (1/ms)
static const double perm50  = p6 * 0.5;      //   perm variable that controls n to p transition
static const double nperm   = p7 * 15;       //     in Hill-like fashion
static const double kn_p    = 500e-3;     //   (1/ms)
static const double kp_n    = 50e-3;      //   (1/ms)
static const double koffmod = 1.0;        //   mod to change species

//   Thin filament regulation and crossbridge cycling
static const double fapp    = p2 * 500e-3;   //   (1/ms) XB on rate
static const double gapp    = p3 * 70e-3;    //   (1/ms) XB off rate
static const double gslmod  = 6;          //   controls SL effect on gapp
static const double hf      = p2 * 2000e-3;  //   (1/ms) rate between pre-force and force states
static const double hfmdc   = 5;          //
static const double hb      = p3 * 400e-3;   //   (1/ms) rate between pre-force and force states
static const double hbmdc   = 0;          //
static const double gxb     = p3 * 70e-3;    //   (1/ms) ATP consuming transition rate
static const double sigmap  = 8;          //   distortion dependence of STP using transition gxb
static const double sigman  = 1;          //
static const double xbmodsp = 4.0 / 3.0;      //   mouse specific modification for XB cycling rates

//   Mean strain of strongly-bound states
static const double x_0     = 0.007;      //   (um) strain induced by head rotation
static const double xPsi    = 2;          //   scaling factor balancing SL motion and XB cycling

//   Normalized active a nd passive force

static const double PCon_t  = 0.002;      //   (norm Force) passive force due to titin
static const double PExp_t  = 10;         //     these apply to trabeculae and single cells only
static const double SL_c    = 2.25;       //   (um) resting length for collagen
static const double PCon_c  = 0.02;       //   (norm Force) passive force due to collagen
static const double PExp_c  = 70;         //     these apply to trabeculae and single cells only

//   Calculation of complete muscle response
static const double massf   = 0.000025e6;  //   ([norm Force ms^2]/um) muscle mass
static const double visc    = 0.003e3;    //   ([norm Force ms]/um) muscle viscosity
static const double KSE     = 1;          //   (norm Force/um) series elastic element
static const double kxb     = 120;        //   (mN/mm^2) maximal force
static const double Trop_conc = 70e-3;       //   (uM) troponin concentration
static const double Temp = 310.0/*-17.0*/;

#endif 