#include "Explicit_TNNP_Vent.h"

namespace oomph{
	ExplicitTNNPVent::ExplicitTNNPVent() : ExplicitTimeStepCellModelBase(){
		this->Required_Storage = 25;

		intrinsic_dt = 0.02;

		//Overload the requests which are not the default value
		Requires_AB_Index = true;
		Requires_RV_Index = true;
		Requires_IS_Index = true;

		// constants
		TTCell_Ko = 5.4;
		TTCell_Cao = 2.0;
		TTCell_Nao = 140.0;
		TTCell_Vc = 0.016404;
		TTCell_Vsr = 0.001094;
		TTCell_Vss = 0.00005468;
		TTCell_Bufc = 0.2;
		TTCell_Kbufc = 0.001;
		TTCell_Bufsr = 10.;
		TTCell_Kbufsr = 0.3;
		TTCell_Bufss = 0.4;
		TTCell_Kbufss = 0.00025;
		TTCell_Vmaxup = 0.006375;
		TTCell_Kup = 0.00025;
		TTCell_Vrel = 0.102; //40.8;
		TTCell_k1_ = 0.15;
		TTCell_k2_ = 0.045;
		TTCell_k3 = 0.060;
		TTCell_k4 = 0.005; //0.000015;
		TTCell_EC = 1.5;
		TTCell_maxsr = 2.5;
		TTCell_minsr = 1.;
		TTCell_Vleak = 0.00036;
		TTCell_Vxfer = 0.0038;
		TTCell_R = 8314.472;
		TTCell_F = 96485.3415;
		TTCell_T = 310.0;
		TTCell_CAPACITANCE = 0.185;
		TTCell_pKNa = 0.03;
		TTCell_GbNa = 0.00029;
		TTCell_KmK = 1.0;
		TTCell_KmNa = 40.0;
		TTCell_knak = 2.724;
		TTCell_GCaL = 0.00003980;
		TTCell_GbCa = 0.000592;
		TTCell_knaca = 1000;
		TTCell_KmNai = 87.5;
		TTCell_KmCa = 1.38;
		TTCell_ksat = 0.1;
		TTCell_n = 0.35;
		TTCell_GpCa = 0.1238;
		TTCell_KpCa = 0.0005;
		TTCell_GpK = 0.0146;

		TTCell_RTONF = (TTCell_R *TTCell_T) / TTCell_F;
		TTCell_inverseVcF2 = 1.0 / (2.0 * TTCell_Vc *TTCell_F);
		TTCell_inverseVcF = 1.0 / (TTCell_Vc *TTCell_F);
		TTCell_inversevssF2 = 1.0 / (2.0 * TTCell_Vss *TTCell_F);
	}


	//Return the initial condition for the nth variable and cell_typeth cell type
	inline bool ExplicitTNNPVent::return_initial_value(const unsigned &n, double &v, const unsigned &cell_type){
		switch(n){
			case 0 : v = 0.00007;
					break;
			case 1 : v = 1.3;
					break;
			case 2 : v = 0.00007;
					break;
			case 3 : v = 7.67;
					break;
		    case 4 : v = 138.3;
					break;
		    case 5 : v = 0.0;
					break;
		    case 6 : v = 0.75;
					break;
		    case 7 : v = 0.75;
					break;
		    case 8 : v = 0.0;
					break;
		    case 9 : v = 1.0;
					break;
		    case 10 : v = 0.0;
					break;
		    case 11 : v = 1.0;
					break;
		    case 12 : v = 0.0;
					break;
		    case 13 : v = 0.0;
					break;
		    case 14 : v = 1.0;
					break;
		    case 15 : v = 1.0;
					break;
		    case 16 : v = 1.0;
					break;
		    case 17 : v = 1.0;
					break;
		    case 18 : v = 0.0;
					break;
		    case 19 : v = 0.0;
					break;
		    case 20 : v = 1.0;
					break;
		    case 21 : v = 0.0;
					break;
		    case 22 : v = 0.0;
					break;
		    case 23 : v = 0.0;
					break;
		    case 24 : v = 0.0;
					break;
			default : return false;
		}
		return true;
	}


	void ExplicitTNNPVent::explicit_timestep(CellState &Cellstate, Vector<double> &new_state){
		//make a vector of size two to avoid having to rewrite a load of code
		Vector<double> state;
		state.resize(2);
		state[0]  	  = Cellstate.vm(); // -86.2;
		state[1]      = new_state[0]; // 0.00007;
		double CaSR   = new_state[1]; // 1.3;
		double CaSS   = new_state[2]; // 0.00007;
		double Nai    = new_state[3]; // 7.67;
		double Ki     = new_state[4]; // 138.3;
		double sm     = new_state[5]; // 0.0;
		double sh     = new_state[6]; // 0.75;
		double sj     = new_state[7]; // 0.75;
		double sxr1   = new_state[8]; // 0.0;
		double sxr2   = new_state[9]; // 1.0;
		double sxs    = new_state[10]; // 0.0;
		double ss     = new_state[11]; // 1.0;
		double sr     = new_state[12]; // 0.0;
		double sd     = new_state[13]; // 0.0;
		double sf     = new_state[14]; // 1.0;
		double sf2    = new_state[15]; // 1.0;
		double sfcass = new_state[16]; // 1.0;
		double sRR    = new_state[17]; // 1.0;
		double mNaL   = new_state[18];
		double hNaL   = new_state[19];
		double wt_C3  = new_state[20];  // 1.0;
		double wt_O   = new_state[21];  // 0.0;
		double wt_C1  = new_state[22];  // 0.0;
		double wt_C2  = new_state[23];  // 0.0;
		double wt_I   = new_state[24];  // 0.0;

		double ABIndex = Cellstate.ab_index();
		unsigned m_celltype = Cellstate.cell_type();
		double IS_Index = Cellstate.is_index();
		double RV_Index = Cellstate.rv_index();
		double dt = Cellstate.dt();

		// std::cout << ABIndex << "\t" << m_celltype << "\t" << IS_Index << "\t" << RV_Index << "\t" << dt << std::endl;

		/*  default values for ABIndex = 0.0 (Value set 0 == Apical Cells)
		    and ABIndex ranges from 0.0 -> 1.0 */
		/*
		*   0.0 -- Apical Cells
		*   1.0 -- Basal Cells
		*/

		double GTo_ABh = 1.0;
		double InAcTo_Vhalf_ABh = 0.0;
		double InAcTo_Vk_ABh = 1.0;
		double GKs_ABh = 1.0;
		double TauKs_ABh = 1.0;
		double AcKs_Vhalf_ABh = 0.0;
		double AB_IKr_ActVhalf = 0.0;

		if (ABIndex < 0.0 or ABIndex > 1.0) {
			std::cerr << "Wrong Apcial-Basal Ratio, at ABIndex = " << ABIndex << std:: endl;
			std::cerr << "Program Existing... " << std::endl;
			std::exit(0);
		}
		/* V1, ABIndex   == 0.0 -> Apical Cells */
		/*GTo_ABh          = (1.0 + (1 - ABIndex) * (29.6 / 16.5 - 1));      // Apical cells: 29.6/16.5; Basal Cells: 1.0
		InAcTo_Vhalf_ABh = -(1 - ABIndex) * 4.0;                          // Apical Cells: -4; Basal Cells: 0.0;
		InAcTo_Vk_ABh    = 1.0 + (1 - ABIndex) * ((4.5 - 3.4) / 3.4);   // Apical cells: 4.5/3.4; Basal Cells: 1.0
		GKs_ABh          = 1.0 + (1 - ABIndex) * (5.6 / 2.1 - 1.0);      // Apical Cells: 5.6/2.1; Basal Cells: 1.0
		// TauKs_ABh        = 1.0 + (1 - ABIndex) * ((358.0 - 516.0) / 516.0);      // Apical Cells:358/516 ; Basal Cells: 1.0
		AB_IKr_ActVhalf = (1-ABIndex) * 4.9;									// Apical cells: 4.9; Basal cells 0.0
		*/


		// longer action potential

		/*GTo_ABh          = 1.0 -  ABIndex * (1.0 - 16.5 / 29.6); // Apical cells: 1.0; Basal Cells: 16.5/29.6
		InAcTo_Vhalf_ABh = ABIndex * 4.0;                          // Apical Cells: 0.0; Basal Cells: 4.0;
		InAcTo_Vk_ABh    = 1.0 -  ABIndex * (1.0 - 3.4 / 4.5); // Apical cells: 1.0; Basal Cells: 3.4/4.5
		GKs_ABh          = 1.0 -  ABIndex * (1.0 - 2.1 / 4.2);//1.0 -  ABIndex * (1.0 - 2.1 / 5.6);  // Apical Cells: 1.0; Basal Cells: 2.1/5.6
		TauKs_ABh = 1.0 + ABIndex * (516.0 / 318.0 - 1); // Apical Cells:1.0 ; Basal Cells: 516/358
		AB_IKr_ActVhalf = -ABIndex * 4.9; // Apical cells: 0.0; Basal cells -4.9*/



		//version three , find the middel value;
		/*{double x;
			x = (29.6-16.5)/(29.6+16.5);
			GTo_ABh          = 1.0 -  (ABIndex-0.5) * x; // Apical cells: 1.0; Basal Cells: 16.5/29.6
			InAcTo_Vhalf_ABh = (ABIndex-0.5) * 4.0;                          // Apical Cells: -2.0; Basal Cells: 2.0;
			x = (4.5-3.4)/(4.5+3.4);
			InAcTo_Vk_ABh    = 1.0 -  (ABIndex-0.5) * x; // Apical cells: 1.0; Basal Cells: 3.4/4.5
			x = (5.6-2.1)/(5.6+2.1);
			GKs_ABh          =  1.0 -  (ABIndex-0.5) * x+0.1;     // Apical Cells: 1.0; Basal Cells: 2.1/5.6
			x = (516-318)/(516+318);
			TauKs_ABh = 1.0 +  (ABIndex-0.5) * x; // Apical Cells:1.0 ; Basal Cells: 516/358
			// TauKs_ABh = 1.0 - ABIndex*(516.0/318.0-1);// Apical Cells:1.0 ; Basal Cells: 516/358
			AB_IKr_ActVhalf = -(ABIndex-0.5) * 4.9; // Apical cells: 2.45; Basal cells -2.45
		} */


		// AB_IKr_ActVhalf = -ABIndex * 4.9; // Apical cells: 0.0; Basal cells -4.9

		// version four, middle 2.
		/*GTo_ABh          = 1.0 -  ABIndex * (1.0 - 16.5 / 29.6); // Apical cells: 1.0; Basal Cells: 16.5/29.6
		InAcTo_Vhalf_ABh = ABIndex * 4.0;                          // Apical Cells: 0.0; Basal Cells: 4.0;
		InAcTo_Vk_ABh    = 1.0 -  ABIndex * (1.0 - 3.4 / 4.5); // Apical cells: 1.0; Basal Cells: 3.4/4.5
		GKs_ABh          = 1.5 -  ABIndex * (1.0 - 2.1 / 5.6);//1.0 -  ABIndex * (1.0 - 2.1 / 5.6);  // Apical Cells: 1.0; Basal Cells: 2.1/5.6
		TauKs_ABh = 1.0 + ABIndex * (516.0 / 318.0 - 1); // Apical Cells:1.0 ; Basal Cells: 516/358*/
		{
			double x;
			double A = 0.4;  // was  0.5 for the LV, here change to 0.35, to increase APD, see what happens.

			if (m_celltype == 105 or m_celltype == 104 or m_celltype == 103)
			{
				A = 0.4;
			}

			if (m_celltype == 103 or m_celltype == 100) {
				A = 0.2;
			}
			x = (29.6 / 16.5 - 1) / (A + 29.6 / 16.5 * (1 - A) );
			GTo_ABh          = 1.0 -  (ABIndex - A) * x; // Apical cells: 1.0; Basal Cells: 16.5/29.6
			InAcTo_Vhalf_ABh = (ABIndex - 1.0) * 4.0;                        // Apical Cells: -2.0; Basal Cells: 2.0;
			x = (4.5 / 3.4 - 1) / (A + 4.5 / 3.4 * (1 - A) );
			InAcTo_Vk_ABh    = 1.0 -  (ABIndex - A) * x; // Apical cells: 1.0; Basal Cells: 3.4/4.5
			x = (5.6 / 2.1 - 1) / (A + 5.6 / 2.1 * (1 - A) );
			GKs_ABh          =  1.0 -  (ABIndex - A) * x;   // Apical Cells: 1.0; Basal Cells: 2.1/5.6
			// x = (516-318)/(516+318);
			x = (358.0 / 516.0 - 1) / (A + 358.0 / 516.0 * (1 - A) );
			// TauKs_ABh = 1.0 -  (ABIndex-1.0) * x; // Apical Cells:1.0 ; Basal Cells: 516/358
			// TauKs_ABh = 1.0 - ABIndex*(516.0/318.0-1);// Apical Cells:1.0 ; Basal Cells: 516/358
			AB_IKr_ActVhalf = -(ABIndex - 1.0) * 4.9; // Apical cells: 2.45; Basal cells -2.45
		}





		/*{
			double x;
			double A = 0.5;

			if (m_celltype == 105 or m_celltype == 104 or m_celltype == 103)
			{
				A = 0.4;
			}
			x = (29.6 / 16.5 - 1) / (A + 29.6 / 16.5 * (1 - A) );
			GTo_ABh          = 1.0 -  (ABIndex - A) * x; // Apical cells: 1.0; Basal Cells: 16.5/29.6
			InAcTo_Vhalf_ABh = (ABIndex - 1.0) * 4.0;                        // Apical Cells: -2.0; Basal Cells: 2.0;
			x = (4.5 / 3.4 - 1) / (A + 4.5 / 3.4 * (1 - A) );
			InAcTo_Vk_ABh    = 1.0 -  (ABIndex - A) * x; // Apical cells: 1.0; Basal Cells: 3.4/4.5
			x = (5.6 / 2.1 - 1) / (A + 5.6 / 2.1 * (1 - A) );
			GKs_ABh          =  1.0 -  (ABIndex - A) * x;   // Apical Cells: 1.0; Basal Cells: 2.1/5.6
			// x = (516-318)/(516+318);
			x = (358.0 / 516.0 - 1) / (A + 358.0 / 516.0 * (1 - A) );
			// TauKs_ABh = 1.0 -  (ABIndex-1.0) * x; // Apical Cells:1.0 ; Basal Cells: 516/358
			// TauKs_ABh = 1.0 - ABIndex*(516.0/318.0-1);// Apical Cells:1.0 ; Basal Cells: 516/358
			AB_IKr_ActVhalf = -(ABIndex - 1.0) * 4.9; // Apical cells: 2.45; Basal cells -2.45
		}
		*/

		/* end of Apical-Basal heterogeneity settings */
		/*if (m_celltype == 105 or m_celltype == 104 or m_celltype == 103)
		{
			if (ABIndex < 0.5)
			{
				ABIndex = 0.5;
			}
			GTo_ABh          = (1.0 + (1 - ABIndex) * (29.6 / 16.5 - 1));      // Apical cells: 29.6/16.5; Basal Cells: 1.0
			InAcTo_Vhalf_ABh = -(1 - ABIndex) * 4.0;                          // Apical Cells: -4; Basal Cells: 0.0;

		}*/

		/*implementation of Ischemia effects: */

		/*reference:
		Ivan V. Kazbanov, et al. 2014. Plos One. Effect of Global Cardiac Ischemia on Human Ventricular
												 Fibrillation: Insights from a Multi-scale Mechanistic
												 Model of the Human Heart */
		double Ischemia_TTCell_Ko = TTCell_Ko + IS_Index * 5.4;
		double Ek                 = TTCell_RTONF * (log((Ischemia_TTCell_Ko / Ki)));
		float f_atp               = 0.0055 * IS_Index;
		double GKATP              = 155.0;
		double IKATP              = GKATP * f_atp * pow(Ischemia_TTCell_Ko / 5.4, 0.3) * (state[0] - Ek) / (40 + 3.5 * exp(0.025 * state[0]));
		float Acidosis_factor     = 1.0  - 0.2 * IS_Index; // Acidosis_factor is the scaling factor to multiply INa and ICaL, to represent the effect of Acidosis following ischemia...



		/*end of ischemia implementation */



		double Ena = TTCell_RTONF * (log((TTCell_Nao / Nai)));
		double Eks = TTCell_RTONF * (log((Ischemia_TTCell_Ko + TTCell_pKNa * TTCell_Nao) / (Ki + TTCell_pKNa * Nai)));
		double Eca = 0.5 * TTCell_RTONF * (log((TTCell_Cao / state[1])));
		double Ak1 = 0.1 / (1. + exp(0.06 * (state[0] - Ek - 200)));
		double Bk1 = (3.0 * exp(0.0002 * (state[0] - Ek + 100)) + exp(0.1 * (state[0] - Ek - 10))) / (1. + exp(-0.5 * (state[0] - Ek)));
		double rec_iK1 = Ak1 / (Ak1 + Bk1);
		double rec_iNaK = (1. / (1. + 0.1245 * exp(-0.1 * state[0] * TTCell_F / (TTCell_R * TTCell_T)) + 0.0353 * exp(-state[0] * TTCell_F / (TTCell_R * TTCell_T))));
		double rec_ipK = 1. / (1. + exp((25 - state[0]) / 5.98));

		double R_INF;
		double S_INF;
		double TAU_R;
		double TAU_S;
		double GNaL, Gks, Gto;

		double Mfactor = 2.0;

		if (m_celltype == 100) {
			GNaL = 0.0065;
			Gks = 0.392 * 1.4;  // was 0.392, changed in the same way as the ORd model did, to increase the APD difference in EPI and ENDO by haibo
			Gto = 0.294;
			R_INF = 1. / (1. + exp((20 - state[0]) / 6.));
			S_INF = 1. / (1. + exp((state[0] + 20 + InAcTo_Vhalf_ABh) / (5.0 * InAcTo_Vk_ABh)));
			TAU_R = 9.5 * exp(-(state[0] + 40.) * (state[0] + 40.) / 1800.) + 0.8;
			TAU_S = 85.*exp(-(state[0] + 45.) * (state[0] + 45.) / 320.) + 5. / (1. + exp((state[0] - 20.) / 5.)) + 3.;
		} else if (m_celltype == 102) {
			GNaL = 0.0065;
			Gks = 0.392;
			Gto = 0.073;
			R_INF = 1. / (1. + exp((20 - state[0]) / 6.));
			S_INF = 1. / (1. + exp((state[0] + 28 + InAcTo_Vhalf_ABh) / (5.0 * InAcTo_Vk_ABh)));
			TAU_R = 9.5 * exp(-(state[0] + 40.) * (state[0] + 40.) / 1800.) + 0.8;
			TAU_S = 1000.*exp(-(state[0] + 67) * (state[0] + 67) / 1000.) + 8.;
		} else if (m_celltype == 101) {
			GNaL = 0.0095;
			Gks = 0.098 * Mfactor;
			Gto = 0.294;
			R_INF = 1. / (1. + exp((20 - state[0]) / 6.));
			S_INF = 1. / (1. + exp((state[0] + 20 + InAcTo_Vhalf_ABh) / (5.0 * InAcTo_Vk_ABh)));
			TAU_R = 9.5 * exp(-(state[0] + 40.) * (state[0] + 40.) / 1800.) + 0.8;
			TAU_S = 85.0 * exp(-(state[0] + 45.) * (state[0] + 45.) / 320.) + 5. / (1. + exp((state[0] - 20.) / 5.)) + 3.;
		} else if (m_celltype == 105) {
			GNaL = 0.0065;
			Gks = 0.392 * (1.0 + 0.8 * RV_Index); // was 2.0
			Gto = 0.073 * (1.0 + 0.8 * RV_Index); // was 3.0
			R_INF = 1. / (1. + exp((20 - state[0]) / 6.));
			S_INF = 1. / (1. + exp((state[0] + 28 + InAcTo_Vhalf_ABh) / (5.0 * InAcTo_Vk_ABh)));
			TAU_R = 9.5 * exp(-(state[0] + 40.) * (state[0] + 40.) / 1800.) + 0.8;
			TAU_S = 1000.*exp(-(state[0] + 67) * (state[0] + 67) / 1000.) + 8.;
		} else if (m_celltype == 104) {
			GNaL = 0.0095;
			Gks = 0.098 * Mfactor *  (1.0 + 0.8 * RV_Index); // was  0.098 * 2.0
			Gto = 0.294 * (1.0 + 0.8 * RV_Index); // was 3.0;
			R_INF = 1. / (1. + exp((20 - state[0]) / 6.));
			S_INF = 1. / (1. + exp((state[0] + 20 + InAcTo_Vhalf_ABh) / (5.0 * InAcTo_Vk_ABh)));
			TAU_R = 9.5 * exp(-(state[0] + 40.) * (state[0] + 40.) / 1800.) + 0.8;
			TAU_S = 85.0 * exp(-(state[0] + 45.) * (state[0] + 45.) / 320.) + 5. / (1. + exp((state[0] - 20.) / 5.)) + 3.;
		} else if (m_celltype == 103) {
			GNaL = 0.0065;
			Gks = 0.392 * 1.4 * (1.0 + 0.8 * RV_Index); // was 0.392, changed in the same way as the ORd model did, by haibo // was 0.392 * 1.4 * 2.0
			Gto = 0.294 * (1.0 + 0.8 * RV_Index); // was 3.0
			R_INF = 1. / (1. + exp((20 - state[0]) / 6.));
			S_INF = 1. / (1. + exp((state[0] + 20 + InAcTo_Vhalf_ABh) / (5.0 * InAcTo_Vk_ABh)));
			TAU_R = 9.5 * exp(-(state[0] + 40.) * (state[0] + 40.) / 1800.) + 0.8;
			TAU_S = 85.*exp(-(state[0] + 45.) * (state[0] + 45.) / 320.) + 5. / (1. + exp((state[0] - 20.) / 5.)) + 3.;
		} else {
			std::cerr << "Cell type wrong!" << std::endl;
			// std::exit(0);
			throw OomphLibError("Cell type wrong!",
		                       OOMPH_CURRENT_FUNCTION,
		                       OOMPH_EXCEPTION_LOCATION);
		}
		// cout << Gto << endl;
		/* Ito */
		ss = S_INF - (S_INF - ss) * exp(-dt / TAU_S);
		sr = R_INF - (R_INF - sr) * exp(-dt / TAU_R);
		double Ito = GTo_ABh * Gto * sr * ss * (state[0] - Ek);

		/* INa */
		// double AM = 1. / (1. + exp((-60. - state[0]) / 5.));
		// double BM = 0.1 / (1. + exp((state[0] + 35.) / 5.)) + 0.10 / (1. + exp((state[0] - 50.) / 200.));
		double TAU_M = (1. / (1. + exp((-60. - state[0]) / 5.))) * ( 0.1 / (1. + exp((state[0] + 35.) / 5.)) + 0.10 / (1. + exp((state[0] - 50.) / 200.))) ;
		double TAU_H;
		double TAU_J;
		double M_INF = 1. / ((1. + exp((-56.86 - state[0]) / 9.03)) * (1. + exp((-56.86 - state[0]) / 9.03)));
		if (state[0] >= -40.) {
			double AH_1 = 0.;
			double BH_1 = (0.77 / (0.13 * (1. + exp(-(state[0] + 10.66) / 11.1))));
			TAU_H = 1.0 / (AH_1 + BH_1);
		} else {
			double AH_2 = (0.057 * exp(-(state[0] + 80.) / 6.8));
			double BH_2 = (2.7 * exp(0.079 * state[0]) + (3.1e5) * exp(0.3485 * state[0]));
			TAU_H = 1.0 / (AH_2 + BH_2);
		}
		double H_INF = 1. / ((1. + exp((state[0] + 71.55) / 7.43)) * (1. + exp((state[0] + 71.55) / 7.43)));
		if (state[0] >= -40.) {
			double AJ_1 = 0.;
			double BJ_1 = (0.6 * exp((0.057) * state[0]) / (1. + exp(-0.1 * (state[0] + 32.))));
			TAU_J = 1.0 / (AJ_1 + BJ_1);
		} else {
			double AJ_2 = (((-2.5428e4) * exp(0.2444 * state[0]) - (6.948e-6) *
			                exp(-0.04391 * state[0])) * (state[0] + 37.78) /
			               (1. + exp(0.311 * (state[0] + 79.23))));
			double BJ_2 = (0.02424 * exp(-0.01052 * state[0]) / (1. + exp(-0.1378 * (state[0] + 40.14))));
			TAU_J = 1.0 / (AJ_2 + BJ_2);
		}
		double J_INF = H_INF;

		sm = M_INF - (M_INF - sm) * exp(-dt / TAU_M);
		sh = H_INF - (H_INF - sh) * exp(-dt / TAU_H);
		sj = J_INF - (J_INF - sj) * exp(-dt / TAU_J);
		double GNa = 14.838;
		if (m_celltype == 104 or m_celltype == 101)
		{
			// GNa *= 1.5;  // tranmural hereterogeneity includes INa
		}
		double INa = Acidosis_factor * GNa * sm * sm * sm * sh * sj * (state[0] - Ena);

		/* IKr */  // CalculateTTIKr
		/*double Xr1_INF;
		double Xr2_INF;
		double TAU_Xr1;
		double TAU_Xr2;
		double Gkr = 0.153;
		Xr1_INF = 1. / (1. + exp((-26. - state[0]) / 7.));
		double axr1 = 450. / (1. + exp((-45. - state[0]) / 10.));
		double bxr1 = 6. / (1. + exp((state[0] - (-30.)) / 11.5));
		TAU_Xr1 = axr1 * bxr1;
		Xr2_INF = 1. / (1. + exp((state[0] - (-88.)) / 24.));
		double axr2 = 3. / (1. + exp((-60. - state[0]) / 20.));
		double bxr2 = 1.12 / (1. + exp((state[0] - 60.) / 20.));
		TAU_Xr2 = axr2 * bxr2;
		sxr1 = Xr1_INF - (Xr1_INF - sxr1) * exp(-dt / TAU_Xr1);
		sxr2 = Xr2_INF - (Xr2_INF - sxr2) * exp(-dt / TAU_Xr2);*/
		// double IKr =  Gkr * sqrt(Ischemia_TTCell_Ko / 5.4) * sxr1 * sxr2 * (state[0] - Ek);
		// CalculateMarkovIKr();

		/* Marcov IKr Model */
		double m_epiMidRatio = 1.5;
		double epi_factor   = 1.8 * m_epiMidRatio;
		double endo_factor  = 1.8;
		double mcell_factor = 1.8;
		double wt_a1 =  2.172;
		double wt_b1 =  1.077;
		double wt_a2 =  0.00655   * exp(0.5 * 0.05547153 * (state[0] - 36. - AB_IKr_ActVhalf));
		double wt_a  =  0.00555   * exp(0.05547153 * (state[0] - 12. - AB_IKr_ActVhalf));
		double wt_b  =  0.002357  * exp(-0.036588 * (state[0] - AB_IKr_ActVhalf));
		double wt_b2 = 0.65 * 0.0029357 * exp(0.69 * -0.02158 * (state[0] - AB_IKr_ActVhalf));
		double wt_ai = 0.11 * 0.439     * exp(1.7 *  -0.02352 * (state[0] + 25. - AB_IKr_ActVhalf)) * (4.5 / Ischemia_TTCell_Ko);
		double wt_bi = 0.4 *  0.656     * exp(      0.000942 * (state[0] - AB_IKr_ActVhalf)) * ((pow((4.5 / Ischemia_TTCell_Ko), 0.3)));
		double wt_mu = (wt_ai * wt_b2) / wt_bi;

		double Gkr = 0.153;
		if (m_celltype == 102 or m_celltype == 105)
			Gkr = 0.0135 * pow(Ischemia_TTCell_Ko, 0.59) * endo_factor;
		else if (m_celltype == 101 or m_celltype == 104)
			Gkr = 0.0135 * pow(Ischemia_TTCell_Ko, 0.59) * mcell_factor;
		else
			Gkr = 0.0135 * pow(Ischemia_TTCell_Ko, 0.59) * epi_factor;

		{	double wt_dC3 = (wt_b * wt_C2) - (wt_a * wt_C3);
			double wt_dC2 = -((wt_b + wt_a1) * wt_C2) + (wt_a * wt_C3) + (wt_b1 * wt_C1);
			double wt_dC1 = -((wt_b1 + wt_a2 + wt_a2) * wt_C1) + (wt_a1 * wt_C2) + (wt_b2 * wt_O) + (wt_mu * wt_I);
			double wt_dO  =  -((wt_b2 + wt_bi) * wt_O) + (wt_a2 * wt_C1) + (wt_ai * wt_I);
			double wt_dI  = -((wt_mu + wt_ai) * wt_I) + (wt_a2 * wt_C1) + (wt_bi * wt_O);

			wt_O  += dt * wt_dO;
			wt_C1 += dt * wt_dC1;
			wt_C2 += dt * wt_dC2;
			wt_C3 += dt * wt_dC3;
			wt_I  += dt * wt_dI;
		}
		double IKr = Gkr * wt_O * (state[0] - Ek);
		
		//An addition so we can output the currents
		Cellstate.set_ikr_current(IKr);
		//////////////////

		/* IKs */
		double Xs_INF = 1. / (1. + exp((-5. - state[0]) / 14.));
		// double Axs = (1400. / (sqrt(1. + exp((5. - state[0]) / 6))));
		// double Bxs = (1. / (1. + exp((state[0] - 35.) / 15.)));
		double TAU_Xs = (1400. / (sqrt(1. + exp((5. - state[0]) / 6)))) * (1. / (1. + exp((state[0] - 35.) / 15.))) + 80;
		TAU_Xs = TauKs_ABh * TAU_Xs;
		sxs = Xs_INF - (Xs_INF - sxs) * exp(-dt / TAU_Xs);

		double IKs = GKs_ABh * Gks * sxs * sxs * (state[0] - Eks);  // with apicalbasal heterogeneity wit direct scaling factor, by haibo
		
		//An addition so we can output the currents
		Cellstate.set_iks_current(IKs);
		//////////////////

		/* ICaL */
		double D_INF = 1. / (1. + exp((-8 - state[0]) / 7.5));
		double Ad = 1.4 / (1. + exp((-35 - state[0]) / 13)) + 0.25;
		double Bd = 1.4 / (1. + exp((state[0] + 5) / 5));
		double Cd = 1. / (1. + exp((50 - state[0]) / 20));
		double TAU_D = Ad * Bd + Cd;
		double F_INF = 1. / (1. + exp((state[0] + 20) / 7));
		double Af = 1102.5 * exp(-(state[0] + 27) * (state[0] + 27) / 225);
		double Bf = 200. / (1 + exp((13 - state[0]) / 10.));
		double Cf = (180. / (1 + exp((state[0] + 30) / 10))) + 20;
		double TAU_F = Af + Bf + Cf;

		double F2_INF = 0.67 / (1. + exp((state[0] + 35) / 7)) + 0.33;
		double Af2 = 600 * exp(-(state[0] + 25) * (state[0] + 25) / 170);
		double Bf2 = 31 / (1. + exp((25 - state[0]) / 10));
		double Cf2 = 16 / (1. + exp((state[0] + 30) / 10));
		double TAU_F2 = Af2 + Bf2 + Cf2;

		sd = D_INF - (D_INF - sd) * exp(-dt / TAU_D);
		sf = F_INF - (F_INF - sf) * exp(-dt / TAU_F);
		sf2 = F2_INF - (F2_INF - sf2) * exp(-dt / TAU_F2);

		double FCaSS_INF = 0.6 / (1 + (CaSS / 0.05) * (CaSS / 0.05)) + 0.4;
		double TAU_FCaSS = 80. / (1 + (CaSS / 0.05) * (CaSS / 0.05)) + 2.;
		sfcass = FCaSS_INF - (FCaSS_INF - sfcass) * exp(-dt / TAU_FCaSS);

		double ICaL = Acidosis_factor * TTCell_GCaL * sd * sf * sf2 * sfcass * 4 * (state[0] - 15) * (TTCell_F * TTCell_F / (TTCell_R * TTCell_T)) *
		              (0.25 * exp(2 * (state[0] - 15) * TTCell_F / (TTCell_R * TTCell_T)) * CaSS - TTCell_Cao) / (exp(2 * (state[0] - 15) * TTCell_F / (TTCell_R * TTCell_T)) - 1.);

		//INaL gates

		// double mNaL_INF, hNaL_INF, TAU_mNaL, TAU_hNaL;
		// need to update INaL later...
		// mNaL = mNaL_INF - (mNaL_INF - mNaL) * exp(-dt / TAU_mNaL);
		// hNaL = hNaL_INF - (hNaL_INF - hNaL) * exp(-dt / TAU_hNaL);
		// double INaL = GNaL * mNaL * mNaL * mNaL * hNaL * (state[0] - Ena);
		double INaL = 0.0;



		/*INaL was taken from ORd 2011 model. */
		double TMP_INF = 1.0 / (1.0 + exp((-(state[0] + 42.85)) / 5.264));

		mNaL = TMP_INF - (TMP_INF - mNaL) * exp(-dt / TAU_M);  // double tmL = TAU_M;
		TMP_INF = 1.0 / (1.0 + exp((state[0] + 87.61) / 7.488));

		hNaL = TMP_INF - (TMP_INF - hNaL) * exp(-dt / 200.0);  // double thL = 200.0;
		/*double hLssp = 1.0 / (1.0 + exp((state[0] + 93.81) / 7.488));
		double thLp = 3.0 * thL;
		hLp = hLssp - (hLssp - hLp) * exp(-dt / thLp);*/
		// double fINaLp = (1.0 / (1.0 + KmCaMK / CaMKa));

		INaL = (1 + 10 * IS_Index) * GNaL * mNaL * hNaL * (state[0] - Ena) ;






		double IbNa = TTCell_GbNa * (state[0] - Ena);
		double IbCa = TTCell_GbCa * (state[0] - Eca);
		double IpK = TTCell_GpK * rec_ipK * (state[0] - Ek);
		double IpCa = TTCell_GpCa * state[1] / (TTCell_KpCa + state[1]);
		double INaK = (1 - 0.35 * IS_Index) * TTCell_knak * (Ischemia_TTCell_Ko / (Ischemia_TTCell_Ko + TTCell_KmK)) * (Nai / (Nai + TTCell_KmNa)) * rec_iNaK;
		double INaCa = TTCell_knaca * (1. / (TTCell_KmNai * TTCell_KmNai * TTCell_KmNai + TTCell_Nao * TTCell_Nao * TTCell_Nao)) * (1. / (TTCell_KmCa + TTCell_Cao)) *
		               (1. / (1 + TTCell_ksat * exp((TTCell_n - 1) * state[0] * TTCell_F / (TTCell_R * TTCell_T)))) *
		               (exp(TTCell_n * state[0] * TTCell_F / (TTCell_R * TTCell_T)) * Nai * Nai * Nai * TTCell_Cao -
		                exp((TTCell_n - 1) * state[0] * TTCell_F / (TTCell_R * TTCell_T)) * TTCell_Nao * TTCell_Nao * TTCell_Nao * state[1] * 2.5);

		double GK1 = 5.405;
		if (m_celltype == 100)
			GK1 *= 1.2; // according to the ORd model supple, added by haibo
		double IK1 = GK1 * rec_iK1 * (state[0] - Ek);
		// double sItot = IKr + IKs + IK1 + Ito + INa + IbNa + ICaL + IbCa + INaK + INaCa + IpCa + IpK + INaL + Istim + IKATP;

		/* update concentrations */
		double kCaSR = TTCell_maxsr - ((TTCell_maxsr - TTCell_minsr) / (1 + (TTCell_EC / CaSR) * (TTCell_EC / CaSR)));
		double k1 = TTCell_k1_ / kCaSR;
		double k2 = TTCell_k2_ * kCaSR;
		double dRR = TTCell_k4 * (1 - sRR) - k2 * CaSS * sRR;
		sRR += dt * dRR;
		double sOO = k1 * CaSS * CaSS * sRR / (TTCell_k3 + k1 * CaSS * CaSS);

		double Irel = TTCell_Vrel * sOO * (CaSR - CaSS);
		double Ileak = TTCell_Vleak * (CaSR - state[1]);
		double Iup = TTCell_Vmaxup / (1. + ((TTCell_Kup * TTCell_Kup) / (state[1] * state[1])));
		double Ixfer = TTCell_Vxfer * (CaSS - state[1]);

		{	double CaCSQN = TTCell_Bufsr * CaSR / (CaSR + TTCell_Kbufsr);
			double dCaSR = dt * (Iup - Irel - Ileak);
			double bjsr = TTCell_Bufsr - CaCSQN - dCaSR - CaSR + TTCell_Kbufsr;
			double cjsr = TTCell_Kbufsr * (CaCSQN + dCaSR + CaSR);
			CaSR = (sqrt(bjsr * bjsr + 4 * cjsr) - bjsr) / 2;
		}


		{	double CaSSBuf = TTCell_Bufss * CaSS / (CaSS + TTCell_Kbufss);
			double dCaSS = dt * (-Ixfer * (TTCell_Vc / TTCell_Vss) + Irel * (TTCell_Vsr / TTCell_Vss) + (-ICaL * TTCell_inversevssF2 * TTCell_CAPACITANCE));
			double bcss = TTCell_Bufss - CaSSBuf - dCaSS - CaSS + TTCell_Kbufss;
			double ccss = TTCell_Kbufss * (CaSSBuf + dCaSS + CaSS);
			CaSS = (sqrt(bcss * bcss + 4 * ccss) - bcss) / 2;
		}


		double CaBuf = TTCell_Bufc * state[1] / (state[1] + TTCell_Kbufc);
		double dCai = dt * ((-(IbCa + IpCa - 2 * INaCa) * TTCell_inverseVcF2 * TTCell_CAPACITANCE) - (Iup - Ileak) * (TTCell_Vsr / TTCell_Vc) + Ixfer);
		double bc = TTCell_Bufc - CaBuf - dCai - state[1] + TTCell_Kbufc;
		double cc = TTCell_Kbufc * (CaBuf + dCai + state[1]);
		state[1] = (sqrt(bc * bc + 4 * cc) - bc) / 2;


		double dNai = -(INa + INaL + IbNa + 3 * INaK + 3 * INaCa) * TTCell_inverseVcF * TTCell_CAPACITANCE;
		Nai += dt * dNai;

		double dKi = -(/*Istim +*/ IK1 + Ito + IKr + IKs - 2 * INaK + IpK + IKATP) * TTCell_inverseVcF * TTCell_CAPACITANCE;
		Ki += dt * dKi;

		new_state[0]   = state[1]    ;
		new_state[1]   = CaSR   ;
		new_state[2]   = CaSS   ;
		new_state[3]   = Nai    ;
		new_state[4]   = Ki     ;
		new_state[5]   = sm     ;
		new_state[6]   = sh     ;
		new_state[7]   = sj     ;
		new_state[8]   = sxr1   ;
		new_state[9]  = sxr2   ;
		new_state[10]  = sxs    ;
		new_state[11]  = ss     ;
		new_state[12]  = sr     ;
		new_state[13]  = sd     ;
		new_state[14]  = sf     ;
		new_state[15]  = sf2    ;
		new_state[16]  = sfcass ;
		new_state[17]  = sRR    ;
		new_state[18]  = mNaL   ;
		new_state[19]  = hNaL   ;
		new_state[20]  = wt_C3  ;
		new_state[21]  = wt_O   ;
		new_state[22]  = wt_C1  ;
		new_state[23]  = wt_C2  ;
		new_state[24]  = wt_I   ;

		/*FILE * out = fopen("INaL.dat", "a+");
		fprintf(out, "%f %f\n", t, INaL);
		fclose(out);*/
		Cellstate.set_cell_membrane_current(IKr + IKs + IK1 + Ito + INa + IbNa + ICaL + IbCa + INaK + INaCa + IpCa + IpK + INaL /*+ Istim*/ + IKATP);
		// std::cout << "reached end" << std::endl;
	}

} //end namespace