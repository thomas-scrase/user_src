#include "TNNP_06_MC_IKr.h"

namespace oomph{
	
TNNP06MCIKr::TNNP06MCIKr(){


	Ko = 5.4;
    Cao = 2.0;
    Nao = 140.0;
    pKNa = 0.03;
    R = 8314.472;
    T = 310.0;
    F = 96485.3415;
    RTONF = R*T/F;
    NFORT = F/(R*T);
    CAPACITANCE = 0.185;

    double BufMg = 5.67;
    double MgiTotal = 5.0;
    double K_MgBuf = 0.174;
    double bMg = BufMg - MgiTotal + K_MgBuf;
    double cMg = K_MgBuf * MgiTotal;
    Mgi = 0.5 * (sqrt(bMg * bMg + 4 * cMg) - bMg);

    TAU_hNaL = 600;







    KpCa = 0.0005;
    GbNa = 0.00029;
    GbCa = 0.000592;
    maxsr = 2.5;
    minsr = 1.0;
    EC = 1.5;
    k1_ = 0.15;
    k2_ = 0.045;
    k4 = 0.005;
    k3 = 0.060;
    Vrel = 0.102;
    Vmaxup = 0.006375;
    Kup = 0.00025;
    Vxfer = 0.0038;
    Bufsr = 10.0;
    Kbufsr = 0.3;
    GCaL = 0.00003980;
    epi_factor   = 2.3; // * 1.5;
    endo_factor  = 2.3;
    mcell_factor = 2.3;
    GK1 = 5.405;
    knaca = 1000.0;
    KmNai = 12.3;
    KmCa = 1.38;
    ksat = 0.27;
    n = 0.35;
    knak = 2.724;
    KmK = 1.0;
    KmNa = 40.0;
    Vleak = 0.00036;
    Bufss = 0.4;
    Kbufss = 0.00025;
    Vc = 16404;
    Vss = 54.68;
    Vsr = 1094;
    inverseVcF2 = 1 / (2 * Vc * F);
    inverseVcF = 1. / (Vc * F);
    inversevssF2 = 1 / (2 * Vss * F);
    Bufc = 0.2;
    Kbufc = 0.001;
    GNa = 14.838;
    s = 1.0648;
    SPM = 0.0014613;
    phi = 0.8838;





	//Assign the names of the variables used by this model
	Names_Of_Cell_Variables =
	{
        "sm",
        "sh",
        "sj",
        "sxs",
        "ss",
        "sr",
        "sd",
        "sf",
        "sf2",
        "sfcass",
        "mNaL",
        "hNaL",
        "f_sd",
        "f_sf",
        "f_sf2",
        "f_sfcass",
        "Ki",
        "Nai",
        "Cai",
        "CaSS",
        "CaSR",
        "sRR",
        // "XR1",
		// "XR2",
		"wt_O",
        "wt_C1",
        "wt_C2",
        "wt_I",
        "wt_C3"
	};
	Names_Of_Other_Parameters =
	{

	};
	Names_Of_Other_Variables =
	{

	};
	Names_Of_Output_Data =
	{

	};

    FinalizeConstruction();
}


// sm = 0.0,
// sh = 0.75,
// sj = 0.75,
// sxs = 0.0,
// ss = 1.0,
// sr = 0.0,
// sd = 0.0,
// sf = 1.0,
// sf2 = 1.0,
// sfcass = 1.0,
// mNaL = 0.0,
// hNaL = 0.0,
// f_sd = 0.0,
// f_sf = 0.0,
// f_sf2 = 0.0,
// f_sfcass = 0.0,
// Ki = 138.3
// Nai = 7.67
// Cai = 0.00007
// CaSS = 0.00007
// CaSR = 1.3
// sRR = 1.0,
// wt_O = 0.0,
// wt_C1 = 0.0,
// wt_C2 = 0.0,
// wt_I = 0.0,
// wt_C3 = 1.0,


double TNNP06MCIKr::return_initial_state_variable(const unsigned &v, const unsigned &cell_type)
{
    switch(v){
        case sm_TT:
            return /*sm =*/ 0.0;
        case sh_TT:
            return /*sh =*/ 0.75;
        case sj_TT:
            return /*sj =*/ 0.75;
        case sxs_TT:
            return /*sxs =*/ 0.0;
        case ss_TT:
            return /*ss =*/ 1.0;
        case sr_TT:
            return /*sr =*/ 0.0;
        case sd_TT:
            return /*sd =*/ 0.0;
        case sf_TT:
            return /*sf =*/ 1.0;
        case sf2_TT:
            return /*sf2 =*/ 1.0;
        case sfcass_TT:
            return /*sfcass =*/ 1.0;
        case mNaL_TT:
            return /*mNaL =*/ 0.0;
        case hNaL_TT:
            return /*hNaL =*/ 0.0;
        case f_sd_TT:
            return /*f_sd =*/ 0.0;
        case f_sf_TT:
            return /*f_sf =*/ 0.0;
        case f_sf2_TT:
            return /*f_sf2 =*/ 0.0;
        case f_sfcass_TT:
            return /*f_sfcass =*/ 0.0;
        case Ki_TT:
            return /*Ki =*/ 138.3;
        case Nai_TT:
            return /*Nai =*/ 7.67;
        case Cai_TT:
            return /*Cai =*/ 0.00007;
        case CaSS_TT:
            return /*CaSS =*/ 0.00007;
        case CaSR_TT:
            return /*CaSR =*/ 1.3;
        case sRR_TT:
            return /*sRR =*/ 1.0;
        case wt_O_TT:
            return /*wt_O =*/ 0.0;
        case wt_C1_TT:
            return /*wt_C1 =*/ 0.0;
        case wt_C2_TT:
            return /*wt_C2 =*/ 0.0;
        case wt_I_TT:
            return /*wt_I =*/ 0.0;
        case wt_C3_TT:
            return /*wt_C3 =*/ 1.0;
    }
}

double TNNP06MCIKr::return_initial_membrane_potential(const unsigned &cell_type)
{
	return -86.2;
}

void TNNP06MCIKr::Calculate_Derivatives(const double &Vm,
                                            const Vector<double> &CellVariables,
                                            const double &t,
                                            const unsigned &cell_type,
                                            const double &Istim,
                                            const Vector<double> &Other_Parameters,
                                            const Vector<double> &Other_Variables,

                                            Vector<double> &Variable_Derivatives,
                                            double &Iion)
{
 //    for(unsigned i=0; i<Num_Cell_Vars; i++){
 //        Variable_Derivatives[i] = -1e9;
 //    }

 //    // std::cout << "BOOM" << std::endl;
	// //Unpack the cell variables
	// const double sm	= CellVariables[sm_TT];
 //    const double sh	= CellVariables[sh_TT];
 //    const double sj	= CellVariables[sj_TT];
 //    const double sxs	= CellVariables[sxs_TT];
 //    const double ss	= CellVariables[ss_TT];
 //    const double sr	= CellVariables[sr_TT];
 //    const double sd	= CellVariables[sd_TT];
 //    const double sf	= CellVariables[sf_TT];
 //    const double sf2	= CellVariables[sf2_TT];
 //    const double sfcass	= CellVariables[sfcass_TT];
 //    const double mNaL	= CellVariables[mNaL_TT];
 //    const double hNaL	= CellVariables[hNaL_TT];
 //    const double f_sd	= CellVariables[f_sd_TT];
 //    const double f_sf	= CellVariables[f_sf_TT];
 //    const double f_sf2	= CellVariables[f_sf2_TT];
 //    const double f_sfcass	= CellVariables[f_sfcass_TT];
 //    const double Ki	= CellVariables[Ki_TT];
 //    const double Nai	= CellVariables[Nai_TT];
 //    const double Cai	= CellVariables[Cai_TT];
 //    const double CaSS	= CellVariables[CaSS_TT];
 //    const double CaSR	= CellVariables[CaSR_TT];
 //    const double sRR	= CellVariables[sRR_TT];
 //    // const double XR1	= CellVariables[XR1_TT];
	// // const double XR2	= CellVariables[XR2_TT];
	// const double wt_O	= CellVariables[wt_O_TT];
 //    const double wt_C1	= CellVariables[wt_C1_TT];
 //    const double wt_C2	= CellVariables[wt_C2_TT];
 //    const double wt_I	= CellVariables[wt_I_TT];
 //    const double wt_C3	= CellVariables[wt_C3_TT];

 //    //Whatever this is
 //    const double test = 0.0;



 //    const double Ek = RTONF * (log((Ko / Ki)));
 //    const double Ena = RTONF * (log((Nao / Nai)));
 //    const double Eks = RTONF * (log((Ko + pKNa * Nao) / (Ki + pKNa * Nai)));
 //    const double Eca = 0.5 * RTONF * (log((Cao / Cai)));
 //    const double Ak1 = 0.1 / (1. + exp(0.06 * (Vm - Ek - 200)));
 //    const double Bk1 = (3. * exp(0.0002 * (Vm - Ek + 100)) + exp(0.1 * (Vm - Ek - 10))) / (1. + exp(-0.5 * (Vm - Ek)));
 //    const double rec_iK1 = Ak1 / (Ak1 + Bk1);
 //    const double rec_iNaK = (1. / (1. + 0.1245 * exp(-0.1 * Vm * F / (R * T)) + 0.0353 * exp(-Vm * F / (R * T))));
 //    const double rec_ipK = 1. / (1. + exp((25 - Vm) / 5.98));


 //    double GNaL;
 //    if (cell_type == 100 || cell_type == 103)        GNaL = 0.0065; // * 10.5;
 //    else if (cell_type == 102 || cell_type == 105)  GNaL = 0.0065; // * 10.5;
 //    else if (cell_type == 101 || cell_type == 104) GNaL = 0.0095; // * 10.5;

 //    const double INaL = GNaL * mNaL * mNaL * mNaL * hNaL * (Vm - Ena);

 //    // // if (ikrFormulation.compare("markov") == 0)
 //    //     if (cell_type == 100 || cell_type == 103)        GNaL = 0.0065; // * 10.5;
 //    //     else if (cell_type == 102 || cell_type == 105)  GNaL = 0.0065; // * 10.5;
 //    //     else if (cell_type == 101 || cell_type == 104) GNaL = 0.0095; // * 10.5;

 //    //     const double INaL = GNaL * mNaL * mNaL * mNaL * hNaL * (Vm - Ena);
 //    // // else
 //    // //     INaL = 0;




 //    //if (ikrFormulation.compare("fink") == 0)
 //    //   CalculateFINKICaL();
 //    //else
    
 //    const double ICaL = GCaL * sd * sf * sf2 * sfcass * 4 * (Vm - 15)*(F * NFORT)*(0.25 * exp(2 * (Vm - 15) * NFORT) * CaSS - Cao) / (exp(2 * (Vm - 15) * NFORT) - 1.);



 //    double ito_right_ventricle_factor;
 //    if (cell_type == 103 || cell_type == 105 || cell_type == 104)
	// 	ito_right_ventricle_factor = 4.0;
	// else 
	// 	ito_right_ventricle_factor = 1.0;
	
	// double Gto;
 //    if (cell_type == 100 || cell_type == 103)		Gto = 0.294 * ito_right_ventricle_factor; 
 //    else if (cell_type == 102 || cell_type == 105)	Gto = 0.073 * ito_right_ventricle_factor;
 //    else if (cell_type == 101 || cell_type == 104) Gto = 0.294 * ito_right_ventricle_factor;

 //    const double Ito = Gto * sr * ss * (Vm - Ek);




 //    // if (_mutantType == WT) {
		
		
		

 //        const double wt_a1 = 2.172;
 //        const double wt_b1 = 1.077;
 //        const double wt_a2 = 0.00655  * exp(0.5 * 0.05547153 * (Vm - 36.));
 //        const double wt_a  = 0.00555  * exp(0.05547153 * (Vm - 12.));
 //        const double wt_b  = 0.002357 * exp(-0.036588 * (Vm));
 //        const double wt_b2 = 0.65     * 0.0029357 * exp(0.69 * -0.02158 * (Vm));
 //        const double wt_ai = 0.11     * 0.439 * exp(1.7 * -0.02352 * (Vm + 25.))*(4.5 / Ko);
 //        const double wt_bi = 0.4      * 0.656 * exp(0.000942 * (Vm))*((pow((4.5 / Ko), 0.3)));
 //        const double wt_mu = (wt_ai * wt_b2) / wt_bi;

	// 	const double dwt_C3 = (wt_b * wt_C2)-(wt_a * wt_C3);
 //        const double dwt_C2 = -((wt_b + wt_a1) * wt_C2)+(wt_a * wt_C3)+(wt_b1 * wt_C1);
 //        const double dwt_C1 = -((wt_b1 + wt_a2 + wt_a2) * wt_C1) + (wt_a1 * wt_C2) + (wt_b2 * wt_O) + (wt_mu * wt_I);
 //        const double dwt_O = -((wt_b2 + wt_bi) * wt_O) + (wt_a2 * wt_C1) + (wt_ai * wt_I);
 //        const double dwt_I = -((wt_mu + wt_ai) * wt_I) + (wt_a2 * wt_C1) + (wt_bi * wt_O);

 //        // wt_O  += dt*dO;
 //        // wt_C1 += dt*dC1;
 //        // wt_C2 += dt*dC2;
 //        // wt_I  += dt*dI;
 //        // const double wt_C3 = 1 - (wt_O + wt_C2 + wt_C1 + wt_I); //+= dt*dC3;
 //        double Gkr;
 //        if (cell_type == 100 || cell_type == 103)
 //        	Gkr = 0.0135 * pow(Ko, 0.59) * epi_factor;   //2.3;
 //        else if (cell_type == 102 || cell_type == 105)
 //        	Gkr = 0.0135 * pow(Ko, 0.59) * endo_factor; //2.3;
 //        else
 //        	Gkr = 0.0135 * pow(Ko, 0.59) * mcell_factor; //2.3;

 //        const double IKr = Gkr * wt_O * (Vm - Ek);
	// // }
	// // else if (_mutantType == HOMOZYGOUS) {
	// // 	//alpha_percentage = 1;
	// // 	epi_factor   = 2.3;//(2.3*3.0)/2.0; //(2.3*2.8)/2.0;	//2.3;
	// // 	endo_factor  = 2.3; //(2.3*1.8)/2.0; //(2.3*1.8)/2.0;	//1.8;
	// // 	mcell_factor = 2.3; //2.3;				//2.0;

	// // 	//CORRECT
	// // 	sqt1_a1 = 2.172;
	// // 	sqt1_b1 = 0.5 *  1.077;
	// // 	sqt1_a2 = 0.3 *  0.00655 * exp(0.05547153 * (Vm - 36. + 15.));
	// // 	sqt1_a  =		0.00555 * exp(0.05547153 * (Vm - 12. + 15.));
	// // 	sqt1_b  =		0.002357 * exp(-0.036588 * (Vm));
	// // 	sqt1_b2 = 0.00077 * 0.0029357 * exp(1.3 * 3.3 * -0.02158 * (Vm)); //0.0029357*exp(2.5*-0.02158*(Vm));
	// // 	sqt1_ai = 0.439 * exp(-0.02352 * (Vm + 25. + 15.))*(4.5 / Ko);
	// // 	sqt1_bi = 0.025 * 0.656 * exp(0.000942 * (Vm + 15.))*((pow((4.5 / Ko), 0.3)));
	// // 	sqt1_mu = (sqt1_ai * sqt1_b2) / sqt1_bi;

	// // 	sqt1_a1 = 			2.172;
	// // 	sqt1_b1 = 0.5*		1.077;
	// // 	sqt1_a2 = 3*			0.00655*exp(0.05547153*(Vm-36.));
	// // 	sqt1_a =  			0.00555*exp(0.05547153*(Vm-12.));
	// // 	sqt1_b =  			0.002357*exp(-0.036588*(Vm));
	// // 	sqt1_b2 = 0.0047*	0.0029357*exp(3.3*-0.02158*(Vm));//0.0029357*exp(2.5*-0.02158*(Vm));
	// // 	sqt1_ai = 			0.439*exp(-0.02352*(Vm+25.))*(4.5/Ko);
	// // 	sqt1_bi = 0.015*		0.656*exp(0.000942*(Vm))*((pow((4.5/Ko),0.3)));
	// // 	sqt1_mu = 			(sqt1_ai*sqt1_b2*sqt1_a2)/(sqt1_a2*sqt1_bi);

	// // 	//dC3 = (sqt1_b * sqt1_C2)-(sqt1_a * sqt1_C3);
 // //        dC2 = -((sqt1_b + sqt1_a1) * sqt1_C2)+(sqt1_a * sqt1_C3)+(sqt1_b1 * sqt1_C1);
 // //        dC1 = -((sqt1_b1 + sqt1_a2 + sqt1_a2) * sqt1_C1) + (sqt1_a1 * sqt1_C2) + (sqt1_b2 * sqt1_O) + (sqt1_mu * sqt1_I);
 // //        dO = -((sqt1_b2 + sqt1_bi) * sqt1_O) + (sqt1_a2 * sqt1_C1) + (sqt1_ai * sqt1_I);
 // //        dI = -((sqt1_mu + sqt1_ai) * sqt1_I) + (sqt1_a2 * sqt1_C1) + (sqt1_bi * sqt1_O);

 // //        sqt1_O += dt*dO;
 // //        sqt1_C1 += dt*dC1;
 // //        sqt1_C2 += dt*dC2;
 // //        sqt1_I += dt*dI;
 // //        sqt1_C3 = 1 - (sqt1_O + sqt1_C2 + sqt1_C1 + sqt1_I); //+= dt*dC3;

 // //        if (cell_type == EPI)
 // //        	Gkr = 0.0135 * pow(Ko, 0.59) * epi_factor;   //2.3;
 // //        else if (cell_type == ENDO)
 // //        	Gkr = 0.0135 * pow(Ko, 0.59) * endo_factor; //2.3;
 // //        else
 // //        	Gkr = 0.0135 * pow(Ko, 0.59) * mcell_factor; //2.3;

 // //        return (SQT1_IKr = Gkr * sqt1_O * (Vm - Ek));
	// // }
	// // else {// HETEROZYGOUS
	// // 	alpha_percentage = 0.5;
	// // 	//Sum of above two formulations
	// // 	IKr = alpha_percentage*Calculate_SQT1_MarkovIKr() + (1-alpha_percentage)*Calculate_WT_MarkovIKr();
	// // }

	// // CalculateNVNKIKr();

 //    double iks_right_ventricle_factor;
 //    if (cell_type == 103 || cell_type == 105 || cell_type == 104)
	// 	iks_right_ventricle_factor = 2.0;
	// else 
	// 	iks_right_ventricle_factor = 1.0;
	
 //    //OURS
 //    double Gks;
 //    if (cell_type == 100 || cell_type == 103)		Gks = 0.392 * iks_right_ventricle_factor; 
 //    else if (cell_type == 102 || cell_type == 105)  Gks = 0.392 * iks_right_ventricle_factor;
 //    else if (cell_type == 101 || cell_type == 104) Gks = 0.098 * iks_right_ventricle_factor; 
	
 //    const double IKs = Gks * sxs * sxs * (Vm - Eks);


    
 //    const double IK1 = GK1 * rec_iK1 * (Vm - Ek);


    
    
    
    
    
 //    const double INaCa = knaca * (1. / (KmNai * KmNai * KmNai + Nao * Nao * Nao))*(1. / (KmCa + Cao)) * (1. / (1 + ksat * exp((n - 1) * Vm * NFORT)))*
 //                    (exp(n * Vm * NFORT) * Nai * Nai * Nai * Cao - exp((n - 1) * Vm * NFORT) * Nao * Nao * Nao * Cai * 2.5);
    

    
    
    
 //    const double INaK = knak * (Ko / (Ko + KmK))*(Nai / (Nai + KmNa)) * rec_iNaK;
    

 //    double GpCa;
 //    if (test == 0.7)
 //        GpCa = 0.0619;
 //    else if (test == 1.4)
 //        GpCa = 0.3714;
 //    else if (test == 1.8)
 //        GpCa = 0.8666;
 //    else
 //        GpCa = 0.1238;

    
 //    const double IpCa = GpCa * Cai / (KpCa + Cai);


 //    double GpK;
 //    if (test == 0.7)
	//     GpK = 0.0730;
	// else if (test == 1.4)
	//     GpCa = 0.0073;
	// else if (test == 1.8)
	//     GpCa = 0.00219;
	// else
	//     GpCa = 0.0146;

	// const double IpK = GpK * rec_ipK * (Vm - Ek);

	
 //    const double IbNa = GbNa * (Vm - Ena);
    
    
 //    const double IbCa = GbCa * (Vm - Eca);

    
    
    
 //    const double kCaSR = maxsr - ((maxsr - minsr) / (1 + (EC / CaSR)*(EC / CaSR)));

    
 //    const double k1 = k1_ / kCaSR;

    
 //    const double k2 = k2_*kCaSR;

    
 //    const double dsRR = k4 * (1 - sRR) - k2 * CaSS*sRR;
 //    // sRR += dt*dRR;
    
 //    const double sOO = k1 * CaSS * CaSS * sRR / (k3 + k1 * CaSS * CaSS);

    
 //    const double Irel = Vrel * sOO * (CaSR - CaSS);

    
 //    const double Ileak = Vleak * (CaSR - Cai);

    
    
 //    const double Iup = Vmaxup / (1. + ((Kup * Kup) / (Cai * Cai)));

    
 //    const double Ixfer = Vxfer * (CaSS - Cai);

    
    
 //    const double CaCSQN = Bufsr * CaSR / (CaSR + Kbufsr);
 //    const double dCaSR = dt * (Iup - Irel - Ileak);
 //    const double bjsr = Bufsr - CaCSQN - dCaSR - CaSR + Kbufsr;
 //    const double cjsr = Kbufsr * (CaCSQN + dCaSR + CaSR);
 //    // CaSR = (sqrt(bjsr * bjsr + 4 * cjsr) - bjsr) / 2;

    
    
    
    
    
    
    
    

 //    const double CaSSBuf = Bufss * CaSS / (CaSS + Kbufss);
 //    const double dCaSS = dt * (-Ixfer * (Vc / Vss) + Irel * (Vsr / Vss)+(-ICaL * inversevssF2 * CAPACITANCE));
 //    const double bcss = Bufss - CaSSBuf - dCaSS - CaSS + Kbufss;
 //    const double ccss = Kbufss * (CaSSBuf + dCaSS + CaSS);
 //    // CaSS = (sqrt(bcss * bcss + 4 * ccss) - bcss) / 2;

    
    
 //    const double CaBuf = Bufc * Cai / (Cai + Kbufc);
 //    const double dCai = dt * ((-(IbCa + IpCa - 2 * INaCa) * inverseVcF2 * CAPACITANCE)-(Iup - Ileak)*(Vsr / Vc) + Ixfer);
 //    const double bc = Bufc - CaBuf - dCai - Cai + Kbufc;
 //    const double cc = Kbufc * (CaBuf + dCai + Cai);
 //    // Cai = (sqrt(bc * bc + 4 * cc) - bc) / 2;

    
 //    const double INa = GNa * sm * sm * sm * sh * sj * (Vm - Ena);

 //    const double dNai = -(INa + IbNa + 3 * INaK + 3 * INaCa) * inverseVcF*CAPACITANCE;
 //    // Nai += dt*dNai;

 //    const double dKi = -(Istim + IK1 + Ito + IKr + IKs - 2 * INaK + IpK) * inverseVcF*CAPACITANCE;
 //    // Ki += dt*dKi;



 //    //compute steady state values and time constants
 //    const double AM = 1. / (1. + exp((-60. - Vm) / 5.));
 //    const double BM = 0.1 / (1. + exp((Vm + 35.) / 5.)) + 0.10 / (1. + exp((Vm - 50.) / 200.));
 //    const double TAU_M = AM*BM;
 //    const double M_INF = 1. / ((1. + exp((-56.86 - Vm) / 9.03))*(1. + exp((-56.86 - Vm) / 9.03)));

 //    double TAU_H;
 //    if (Vm >= -40.0) {
 //        const double AH_1 = 0.;
 //        const double BH_1 = (0.77 / (0.13 * (1. + exp(-(Vm + 10.66) / 11.1))));
 //        TAU_H = 1.0 / (AH_1 + BH_1);
 //    } else {
 //        const double AH_2 = (0.057 * exp(-(Vm + 80.) / 6.8));
 //        const double BH_2 = (2.7 * exp(0.079 * Vm)+(3.1e5) * exp(0.3485 * Vm));
 //        TAU_H = 1.0 / (AH_2 + BH_2);
 //    }

 //    const double H_INF = 1. / ((1. + exp((Vm + 71.55) / 7.43))*(1. + exp((Vm + 71.55) / 7.43)));

 //    double TAU_J;
 //    if (Vm >= -40.) {
 //        const double AJ_1 = 0.;
 //        const double BJ_1 = (0.6 * exp((0.057) * Vm) / (1. + exp(-0.1 * (Vm + 32.))));
 //        TAU_J = 1.0 / (AJ_1 + BJ_1);
 //    } else {
 //        const double AJ_2 = (((-2.5428e4) * exp(0.2444 * Vm)-(6.948e-6) * exp(-0.04391 * Vm))*(Vm + 37.78) / (1. + exp(0.311 * (Vm + 79.23))));
 //        const double BJ_2 = (0.02424 * exp(-0.01052 * Vm) / (1. + exp(-0.1378 * (Vm + 40.14))));
 //        TAU_J = 1.0 / (AJ_2 + BJ_2);
 //    }

 //    const double J_INF = H_INF;

 //    //======================================================================
 //    //    // IKr area (ORIGINAL TEN TUSSCHER)
 //    //    Xr1_INF=1./(1.+exp((-26.-Vm)/7.));
 //    //    axr1=450./(1.+exp((-45.-Vm)/10.));
 //    //    bxr1=6./(1.+exp((Vm-(-30.))/11.5));
 //    //    TAU_Xr1=axr1*bxr1;
 //    //    Xr2_INF=1./(1.+exp((Vm-(-88.))/24.));
 //    //    axr2=3./(1.+exp((-60.-Vm)/20.));
 //    //    bxr2=1.12/(1.+exp((Vm-60.)/20.));
 //    //    TAU_Xr2=axr2*bxr2;



 //    // if (_mutant) {//SQT1
 //    //     axr1 = 450. / (1. + exp((-45. - Vm) / 10.));
 //    //     bxr1 = 6. / (1. + exp((Vm - (-30.)) / 11.5));
 //    //     TAU_Xr1 = axr1*bxr1;
 //    //     Xr1_INF = 2.0 * 1. / (1. + exp((-26. - Vm) / 7.));

 //    //     axr2 = 3. / (1. + exp((-60. - Vm) / 20.));
 //    //     bxr2 = 1.12 / (1. + exp((Vm - 60.) / 20.));
 //    //     TAU_Xr2 = axr2*bxr2;
 //    //     Xr2_INF = 2.0 * 1. / (1. + exp((0.7 * Vm - (-88.)) / 24.)); //* 0.7

 //    // } else {//Wild Type	[+3 to Vm]
 //        const double axr1 = 450. / (1. + exp((-45. - Vm + 3) / 10.));
 //        const double bxr1 = 6. / (1. + exp((Vm - (-30.) + 3) / 11.5));
 //        const double TAU_Xr1 = axr1*bxr1;
 //        const double Xr1_INF = 1. / (1. + exp((-26. - Vm + 3) / 7.));

 //        const double axr2 = 3. / (1. + exp((-60. - Vm + 3) / 20.));
 //        const double bxr2 = 1.12 / (1. + exp((Vm - 60. + 3) / 20.));
 //        const double TAU_Xr2 = axr2*bxr2;
 //        const double Xr2_INF = 1. / (1. + exp((0.7 * Vm - (-88.) + 3) / 24.)); //* 0.7
 //    // }
 //    //======================================================================

 //    const double Xs_INF = 1. / (1. + exp((-5. - Vm) / 14.));
 //    const double Axs = (1400. / (sqrt(1. + exp((5. - Vm) / 6))));
 //    const double Bxs = (1. / (1. + exp((Vm - 35.) / 15.)));
 //    const double TAU_Xs = Axs * Bxs + 80;


 //    double R_INF;
 //    double S_INF;
 //    double TAU_R;
 //    double TAU_S;
 //    if (cell_type == 100 || cell_type == 103) {
 //        R_INF = 1. / (1. + exp((20 - Vm) / 6.));
 //        S_INF = 1. / (1. + exp((Vm + 20) / 5.));
 //        TAU_R = 9.5 * exp(-(Vm + 40.)*(Vm + 40.) / 1800.) + 0.8;
 //        TAU_S = 85. * exp(-(Vm + 45.)*(Vm + 45.) / 320.) + 5. / (1. + exp((Vm - 20.) / 5.)) + 3.;
 //    }else if (cell_type == 102 || cell_type == 105) {
 //        R_INF = 1. / (1. + exp((20 - Vm) / 6.));
 //        S_INF = 1. / (1. + exp((Vm + 28) / 5.));
 //        TAU_R = 9.5 * exp(-(Vm + 40.)*(Vm + 40.) / 1800.) + 0.8;
 //        TAU_S = 1000. * exp(-(Vm + 67)*(Vm + 67) / 1000.) + 8.;
 //    }else if (cell_type == 101 || cell_type == 104) {
 //        R_INF = 1. / (1. + exp((20 - Vm) / 6.));
 //        S_INF = 1. / (1. + exp((Vm + 20) / 5.));
 //        TAU_R = 9.5 * exp(-(Vm + 40.)*(Vm + 40.) / 1800.) + 0.8;
 //        TAU_S = 85. * exp(-(Vm + 45.)*(Vm + 45.) / 320.) + 5. / (1. + exp((Vm - 20.) / 5.)) + 3.;
 //    }

 //    const double D_INF = 1. / (1. + exp((-8 - Vm) / 7.5));
 //    const double Ad = 1.4 / (1. + exp((-35 - Vm) / 13)) + 0.25;
 //    const double Bd = 1.4 / (1. + exp((Vm + 5) / 5));
 //    const double Cd = 1. / (1. + exp((50 - Vm) / 20));
 //    const double TAU_D = Ad * Bd + Cd;
 //    const double F_INF = 1. / (1. + exp((Vm + 20) / 7));
 //    const double Af = 1102.5 * exp(-(Vm + 27)*(Vm + 27) / 225);
 //    const double Bf = 200. / (1 + exp((13 - Vm) / 10.));
 //    const double Cf = (180. / (1 + exp((Vm + 30) / 10))) + 20;

 //    double TAU_F;
 //    if (test == 0.7)
 //        TAU_F = (Af + Bf + Cf)*0.6;
 //    else if (test == 1.4)
 //        TAU_F = (Af + Bf + Cf)*1.5;
 //    else if (test == 1.8)
 //        TAU_F = (Af + Bf + Cf)*2.0;
 //    else
 //        TAU_F = Af + Bf + Cf;

 //    const double F2_INF = 0.67 / (1. + exp((Vm + 35) / 7)) + 0.33;
 //    const double Af2 = 600 * exp(-(Vm + 25)*(Vm + 25) / 170);
 //    const double Bf2 = 31 / (1. + exp((25 - Vm) / 10));
 //    const double Cf2 = 16 / (1. + exp((Vm + 30) / 10));
 //    const double TAU_F2 = Af2 + Bf2 + Cf2;
 //    const double FCaSS_INF = 0.6 / (1 + (CaSS / 0.05)*(CaSS / 0.05)) + 0.4;
 //    const double TAU_FCaSS = 80. / (1 + (CaSS / 0.05)*(CaSS / 0.05)) + 2.;

 //    //INaL Area
 //    const double alpha_mNaL = (0.32 * (Vm + 47.13)) / (1 - exp(-0.1 * (Vm + 47.13)));
 //    const double beta_mNaL = 0.08 * exp(-Vm / 11.0);
 //    const double mNaL_INF = alpha_mNaL / (alpha_mNaL + beta_mNaL);
 //    const double TAU_mNaL = 1.0 / (alpha_mNaL + beta_mNaL); //TAU_hNaL initialised in constructor.
 //    const double hNaL_INF = 1.0 / (1 + exp((Vm + 91.0) / 6.1));

 //    //FINK IK1
    
 //    const double K1Mg = 2.8 * exp(-(Vm - s * Ek) / 180.);
 //    const double KBMg = 0.45 * exp(-(Vm - s * Ek) / 20.);
 //    const double Kd1SPM = 0.0007 * exp(-(Vm - s * Ek + 8 * Mgi) / 4.8);
 //    const double Kd2SPM = 0.04 * exp(-(Vm - s * Ek) / 9.1);

 //    const double psi = 1 + Mgi / KBMg;
    
 //    const double r1_INF = (psi * psi) / ((SPM / Kd1SPM) + (Mgi / K1Mg) + psi * psi * psi);
 //    const double r2_INF = 1. / (1. + (SPM / Kd2SPM));
    
 //    const double GoverGmax = phi * r1_INF + (1 - phi) * r2_INF;

 //    //Fink ICaL
 //    const double fD_INF = 1. / (1. + exp((5 - Vm) / 7.5));
 //    const double fTAU_F = (Af + Bf + Cf) / 4.;
 //    const double fF2_INF = 0.75 / (1. + exp((Vm + 35) / 7)) + 0.25;
 //    const double fTAU_F2 = (Af2 + Bf2 + Cf2) / 2.;
 //    const double fFCaSS_INF = 0.4 / (1 + (CaSS / 0.05)*(CaSS / 0.05)) + 0.6;



 //    //Update gates
 //    // sm = M_INF - (M_INF - sm) * exp(-dt / TAU_M);
 //    // sh = H_INF - (H_INF - sh) * exp(-dt / TAU_H);
 //    // sj = J_INF - (J_INF - sj) * exp(-dt / TAU_J);

 //    // //Markov IKr Gates
 //    // // sxr1 = Xr1_INF - (Xr1_INF - sxr1) * exp(-dt / TAU_Xr1);
 //    // // sxr2 = Xr2_INF - (Xr2_INF - sxr2) * exp(-dt / TAU_Xr2);

 //    // //ten Tusscher IKr Gates
 //    // // tt_sxr1 = tt_Xr1_INF - (tt_Xr1_INF - tt_sxr1) * exp(-dt / tt_TAU_Xr1);
 //    // // tt_sxr2 = tt_Xr2_INF - (tt_Xr2_INF - tt_sxr2) * exp(-dt / tt_TAU_Xr2);

 //    // sxs = Xs_INF - (Xs_INF - sxs) * exp(-dt / TAU_Xs);
 //    // ss = S_INF - (S_INF - ss) * exp(-dt / TAU_S);
 //    // sr = R_INF - (R_INF - sr) * exp(-dt / TAU_R);
 //    // sd = D_INF - (D_INF - sd) * exp(-dt / TAU_D);
 //    // sf = F_INF - (F_INF - sf) * exp(-dt / TAU_F);
 //    // sf2 = F2_INF - (F2_INF - sf2) * exp(-dt / TAU_F2);
 //    // sfcass = FCaSS_INF - (FCaSS_INF - sfcass) * exp(-dt / TAU_FCaSS);

 //    // //INaL gates
 //    // mNaL = mNaL_INF - (mNaL_INF - mNaL) * exp(-dt / TAU_mNaL);
 //    // hNaL = hNaL_INF - (hNaL_INF - hNaL) * exp(-dt / TAU_hNaL);

 //    // //Fink ICaL
 //    // f_sd = fD_INF - (fD_INF - f_sd) * exp(-dt / TAU_D);
 //    // f_sf = F_INF - (F_INF - f_sf) * exp(-dt / TAU_F);
 //    // f_sf2 = fF2_INF - (F2_INF - f_sf2) * exp(-dt / TAU_F2);
 //    // f_sfcass = fFCaSS_INF - (fFCaSS_INF - f_sfcass) * exp(-dt / TAU_FCaSS);

 //    // const double dsm = ( M_INF - (M_INF - sm) * exp(-dt / TAU_M) - sm )/dt;
 //    // const double dsh = ( H_INF - (H_INF - sh) * exp(-dt / TAU_H) - sh )/dt;
 //    // const double dsj = ( J_INF - (J_INF - sj) * exp(-dt / TAU_J) - sj )/dt;
 //    // const double dsxs = ( Xs_INF - (Xs_INF - sxs) * exp(-dt / TAU_Xs) - sxs )/dt;
 //    // const double dss = ( S_INF - (S_INF - ss) * exp(-dt / TAU_S) - ss )/dt;
 //    // const double dsr = ( R_INF - (R_INF - sr) * exp(-dt / TAU_R) - sr )/dt;
 //    // const double dsd = ( D_INF - (D_INF - sd) * exp(-dt / TAU_D) - sd )/dt;
 //    // const double dsf = ( F_INF - (F_INF - sf) * exp(-dt / TAU_F) - sf )/dt;
 //    // const double dsf2 = ( F2_INF - (F2_INF - sf2) * exp(-dt / TAU_F2) - sf2 )/dt;
 //    // const double dsfcass = ( FCaSS_INF - (FCaSS_INF - sfcass) * exp(-dt / TAU_FCaSS) - sfcass )/dt;
 //    // const double dmNaL = ( mNaL_INF - (mNaL_INF - mNaL) * exp(-dt / TAU_mNaL) - mNaL )/dt;
 //    // const double dhNaL = ( hNaL_INF - (hNaL_INF - hNaL) * exp(-dt / TAU_hNaL) - hNaL )/dt;
 //    // const double df_sd = ( fD_INF - (fD_INF - f_sd) * exp(-dt / TAU_D) - f_sd )/dt;
 //    // const double df_sf = ( F_INF - (F_INF - f_sf) * exp(-dt / TAU_F) - f_sf )/dt;
 //    // const double df_sf2 = ( fF2_INF - (F2_INF - f_sf2) * exp(-dt / TAU_F2) - f_sf2 )/dt;
 //    // const double df_sfcass = ( fFCaSS_INF - (fFCaSS_INF - f_sfcass) * exp(-dt / TAU_FCaSS) - f_sfcass )/dt;


 //    std:: cout << "sm_TT " << TAU_M << std::endl;
 //    std:: cout << "sh_TT " << TAU_H << std::endl;
 //    std:: cout << "sj_TT " << TAU_J << std::endl;
 //    std:: cout << "sxs_TT " << TAU_Xs << std::endl;
 //    std:: cout << "ss_TT " << TAU_S << std::endl;
 //    std:: cout << "sr_TT " << TAU_R << std::endl;
 //    std:: cout << "sd_TT " << TAU_D << std::endl;
 //    std:: cout << "sf_TT " << TAU_F << std::endl;
 //    std:: cout << "sf2_TT " << TAU_F2 << std::endl;
 //    std:: cout << "sfcass_TT " << TAU_FCaSS << std::endl;
 //    std:: cout << "mNaL_TT " << TAU_mNaL << std::endl;
 //    std:: cout << "hNaL_TT " << TAU_hNaL << std::endl;
 //    std:: cout << "f_sd_TT " << TAU_D << std::endl;
 //    std:: cout << "f_sf_TT " << TAU_F << std::endl;
 //    std:: cout << "f_sf2_TT " << TAU_F2 << std::endl;
 //    std:: cout << "f_sfcass_TT " << TAU_FCaSS << std::endl;

 //    std::cout << "sm_TT " <<  ( M_INF - (M_INF - sm) * exp(-dt / TAU_M) - sm )/dt << std::endl;
 //    std::cout << "sh_TT " <<  ( H_INF - (H_INF - sh) * exp(-dt / TAU_H) - sh )/dt << std::endl;
 //    std::cout << "sj_TT " <<  ( J_INF - (J_INF - sj) * exp(-dt / TAU_J) - sj )/dt << std::endl;
 //    std::cout << "sxs_TT " <<  ( Xs_INF - (Xs_INF - sxs) * exp(-dt / TAU_Xs) - sxs )/dt << std::endl;
 //    std::cout << "ss_TT " <<  ( S_INF - (S_INF - ss) * exp(-dt / TAU_S) - ss )/dt << std::endl;
 //    std::cout << "sr_TT " <<  ( R_INF - (R_INF - sr) * exp(-dt / TAU_R) - sr )/dt << std::endl;
 //    std::cout << "sd_TT " <<  ( D_INF - (D_INF - sd) * exp(-dt / TAU_D) - sd )/dt << std::endl;
 //    std::cout << "sf_TT " <<  ( F_INF - (F_INF - sf) * exp(-dt / TAU_F) - sf )/dt << std::endl;
 //    std::cout << "sf2_TT " <<  ( F2_INF - (F2_INF - sf2) * exp(-dt / TAU_F2) - sf2 )/dt << std::endl;
 //    std::cout << "sfcass_TT " <<  ( FCaSS_INF - (FCaSS_INF - sfcass) * exp(-dt / TAU_FCaSS) - sfcass )/dt << std::endl;
 //    std::cout << "mNaL_TT " <<  ( mNaL_INF - (mNaL_INF - mNaL) * exp(-dt / TAU_mNaL) - mNaL )/dt << std::endl;
 //    std::cout << "hNaL_TT " <<  ( hNaL_INF - (hNaL_INF - hNaL) * exp(-dt / TAU_hNaL) - hNaL )/dt << std::endl;
 //    std::cout << "f_sd_TT " <<  ( fD_INF - (fD_INF - f_sd) * exp(-dt / TAU_D) - f_sd )/dt << std::endl;
 //    std::cout << "f_sf_TT " <<  ( F_INF - (F_INF - f_sf) * exp(-dt / TAU_F) - f_sf )/dt << std::endl;
 //    std::cout << "f_sf2_TT " <<  ( fF2_INF - (F2_INF - f_sf2) * exp(-dt / TAU_F2) - f_sf2 )/dt << std::endl;
 //    std::cout << "f_sfcass_TT " <<  ( fFCaSS_INF - (fFCaSS_INF - f_sfcass) * exp(-dt / TAU_FCaSS) - f_sfcass )/dt << std::endl;

 //    // if(TAU_M<1e-12){exit(0);}
 //    // if(TAU_H<1e-12){exit(0);}
 //    // if(TAU_J<1e-12){exit(0);}
 //    // if(TAU_Xs<1e-12){exit(0);}
 //    // if(TAU_S<1e-12){exit(0);}
 //    // if(TAU_R<1e-12){exit(0);}
 //    // if(TAU_D<1e-12){exit(0);}
 //    // if(TAU_F<1e-12){exit(0);}
 //    // if(TAU_F2<1e-12){exit(0);}
 //    // if(TAU_FCaSS<1e-12){exit(0);}
 //    // if(TAU_mNaL<1e-12){exit(0);}
 //    // if(TAU_hNaL<1e-12){exit(0);}
 //    // if(TAU_D<1e-12){exit(0);}
 //    // if(TAU_F<1e-12){exit(0);}
 //    // if(TAU_F2<1e-12){exit(0);}
 //    // if(TAU_FCaSS<1e-12){exit(0);}

 //    //Package up the derivatives
 //    Variable_Derivatives[sm_TT] =  ( M_INF - (M_INF - sm) * exp(-(dt*1.0e-6) / TAU_M) - sm )/(dt*1.0e-6);
 //    Variable_Derivatives[sh_TT] =  ( H_INF - (H_INF - sh) * exp(-(dt*1.0e-6) / TAU_H) - sh )/(dt*1.0e-6);
 //    Variable_Derivatives[sj_TT] =  ( J_INF - (J_INF - sj) * exp(-(dt*1.0e-6) / TAU_J) - sj )/(dt*1.0e-6);
 //    Variable_Derivatives[sxs_TT] =  ( Xs_INF - (Xs_INF - sxs) * exp(-(dt*1.0e-6) / TAU_Xs) - sxs )/(dt*1.0e-6);
 //    Variable_Derivatives[ss_TT] =  ( S_INF - (S_INF - ss) * exp(-(dt*1.0e-6) / TAU_S) - ss )/(dt*1.0e-6);
 //    Variable_Derivatives[sr_TT] =  ( R_INF - (R_INF - sr) * exp(-(dt*1.0e-6) / TAU_R) - sr )/(dt*1.0e-6);
 //    Variable_Derivatives[sd_TT] =  ( D_INF - (D_INF - sd) * exp(-(dt*1.0e-6) / TAU_D) - sd )/(dt*1.0e-6);
 //    Variable_Derivatives[sf_TT] =  ( F_INF - (F_INF - sf) * exp(-(dt*1.0e-6) / TAU_F) - sf )/(dt*1.0e-6);
 //    Variable_Derivatives[sf2_TT] =  ( F2_INF - (F2_INF - sf2) * exp(-(dt*1.0e-6) / TAU_F2) - sf2 )/(dt*1.0e-6);
 //    Variable_Derivatives[sfcass_TT] =  ( FCaSS_INF - (FCaSS_INF - sfcass) * exp(-(dt*1.0e-6) / TAU_FCaSS) - sfcass )/(dt*1.0e-6);
 //    Variable_Derivatives[mNaL_TT] =  ( mNaL_INF - (mNaL_INF - mNaL) * exp(-(dt*1.0e-6) / TAU_mNaL) - mNaL )/(dt*1.0e-6);
 //    Variable_Derivatives[hNaL_TT] =  ( hNaL_INF - (hNaL_INF - hNaL) * exp(-(dt*1.0e-6) / TAU_hNaL) - hNaL )/(dt*1.0e-6);
 //    Variable_Derivatives[f_sd_TT] =  ( fD_INF - (fD_INF - f_sd) * exp(-(dt*1.0e-6) / TAU_D) - f_sd )/(dt*1.0e-6);
 //    Variable_Derivatives[f_sf_TT] =  ( F_INF - (F_INF - f_sf) * exp(-(dt*1.0e-6) / TAU_F) - f_sf )/(dt*1.0e-6);
 //    Variable_Derivatives[f_sf2_TT] =  ( fF2_INF - (F2_INF - f_sf2) * exp(-(dt*1.0e-6) / TAU_F2) - f_sf2 )/(dt*1.0e-6);
 //    Variable_Derivatives[f_sfcass_TT] =  ( fFCaSS_INF - (fFCaSS_INF - f_sfcass) * exp(-(dt*1.0e-6) / TAU_FCaSS) - f_sfcass )/(dt*1.0e-6);
 //    Variable_Derivatives[Ki_TT] = dKi;
 //    Variable_Derivatives[Nai_TT] = dNai;
 //    Variable_Derivatives[Cai_TT] = ( (sqrt(bc * bc + 4 * cc) - bc) / 2 - Cai )/dt;
 //    Variable_Derivatives[CaSS_TT] = ( (sqrt(bcss * bcss + 4 * ccss) - bcss) / 2 - CaSS )/dt;
 //    Variable_Derivatives[CaSR_TT] = ( (sqrt(bjsr * bjsr + 4 * cjsr) - bjsr) / 2 - CaSR)/dt;
 //    Variable_Derivatives[sRR_TT] = dsRR;
 //    // Variable_Derivatives["XR1"] = dXR1;
	// // Variable_Derivatives["XR2"] = dXR2;
	// Variable_Derivatives[wt_O_TT] = dwt_O;
 //    Variable_Derivatives[wt_C1_TT] = dwt_C1;
 //    Variable_Derivatives[wt_C2_TT] = dwt_C2;
 //    Variable_Derivatives[wt_I_TT] = dwt_I;
 //    Variable_Derivatives[wt_C3_TT] = dwt_C3;

 //    // std::cout << Variable_Derivatives[sm_TT] << std::endl;
 //    // std::cout << Variable_Derivatives[sh_TT] << std::endl;
 //    // std::cout << Variable_Derivatives[sj_TT] << std::endl;
 //    // std::cout << Variable_Derivatives[sxs_TT] << std::endl;
 //    // std::cout << Variable_Derivatives[ss_TT] << std::endl;
 //    // std::cout << Variable_Derivatives[sr_TT] << std::endl;
 //    // std::cout << Variable_Derivatives[sd_TT] << std::endl;
 //    // std::cout << Variable_Derivatives[sf_TT] << std::endl;
 //    // std::cout << Variable_Derivatives[sf2_TT] << std::endl;
 //    // std::cout << Variable_Derivatives[sfcass_TT] << std::endl;
 //    // std::cout << Variable_Derivatives[mNaL_TT] << std::endl;
 //    // std::cout << Variable_Derivatives[hNaL_TT] << std::endl;
 //    // std::cout << Variable_Derivatives[f_sd_TT] << std::endl;
 //    // std::cout << Variable_Derivatives[f_sf_TT] << std::endl;
 //    // std::cout << Variable_Derivatives[f_sf2_TT] << std::endl;
 //    // std::cout << Variable_Derivatives[f_sfcass_TT] << std::endl;
 //    // std::cout << Variable_Derivatives[Ki_TT] << std::endl;
 //    // std::cout << Variable_Derivatives[Nai_TT] << std::endl;
 //    // std::cout << Variable_Derivatives[Cai_TT] << std::endl;
 //    // std::cout << Variable_Derivatives[CaSS_TT] << std::endl;
 //    // std::cout << Variable_Derivatives[CaSR_TT] << std::endl;
 //    // std::cout << Variable_Derivatives[sRR_TT] << std::endl;
 //    // std::cout << Variable_Derivatives[wt_O_TT] << std::endl;
 //    // std::cout << Variable_Derivatives[wt_C1_TT] << std::endl;
 //    // std::cout << Variable_Derivatives[wt_C2_TT] << std::endl;
 //    // std::cout << Variable_Derivatives[wt_I_TT] << std::endl;
 //    // std::cout << Variable_Derivatives[wt_C3_TT] << std::endl;

 //    for(unsigned i=0; i<Num_Cell_Vars; i++){
 //        std::cout << Names_Of_Cell_Variables[i] << " " << Variable_Derivatives[i] << std::endl;
 //        if(!std::isfinite(Variable_Derivatives[i])){
 //            exit(0);
 //        }
 //    }

 //    Iion = -(IKr + IKs + IK1 + Ito + INa + IbNa + ICaL + IbCa + INaK + INaCa + IpCa + IpK + INaL + Istim);
 //    // std::cout << "IKr " << IKr << std::endl;
 //    // std::cout << "IKs " << IKs << std::endl;
 //    // std::cout << "IK1 " << IK1 << std::endl;
 //    // std::cout << "Ito " << Ito << std::endl;
 //    // std::cout << "INa " << INa << std::endl;
 //    // std::cout << "IbNa " << IbNa << std::endl;
 //    // std::cout << "ICaL " << ICaL << std::endl;
 //    // std::cout << "IbCa " << IbCa << std::endl;
 //    // std::cout << "INaK " << INaK << std::endl;
 //    // std::cout << "INaCa " << INaCa << std::endl;
 //    // std::cout << "IpCa " << IpCa << std::endl;
 //    // std::cout << "IpK " << IpK << std::endl;
 //    // std::cout << "INaL " << INaL << std::endl;
 //    // std::cout << "Istim " << Istim << std::endl;
 //    // std::cout << Iion << std::endl;

 //    // if(IKs < 0.0){std::cout << "IKs = " << IKs << std::endl;/* exit(0);*/}
}





void TNNP06MCIKr::get_optional_output(const double &Vm,
									const Vector<double> &CellVariables,
									const double &t,
									const unsigned &cell_type,
									const double &Istim,
									const Vector<double> &Other_Parameters,
									const Vector<double> &Other_Variables,
									Vector<double> &Out) const
{

}


} //end namespace
