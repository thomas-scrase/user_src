#include "TNNP_Vent.h"

namespace oomph{

	TNNPVent::TNNPVent()	:	CellModelBase()
	{
		this->Required_Storage = 25;

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

	//====================================================================
    //====================================================================
    // Define the channel currents
    //====================================================================
    //====================================================================
    void TNNPVent::IKr_current( CellState &state) 
    {
    	double m_epiMidRatio = 1.5;
		double epi_factor   = 1.8 * m_epiMidRatio;
		double endo_factor  = 1.8;
		double mcell_factor = 1.8;

		double Gkr = 0.153;
		if (state.cell_type() == 102 or state.cell_type() == 105)
			Gkr = 0.0135 * pow(get_Ischemia_TTCell_Ko(state), 0.59) * endo_factor;
		else if (state.cell_type() == 101 or state.cell_type() == 104)
			Gkr = 0.0135 * pow(get_Ischemia_TTCell_Ko(state), 0.59) * mcell_factor;
		else
			Gkr = 0.0135 * pow(get_Ischemia_TTCell_Ko(state), 0.59) * epi_factor;

    	state.set_ikr_current(Gkr * state.get_var(0,21) * (state.get_vm() - state.ek()));
    }

	void TNNPVent::IKs_current( CellState &state) 
	{
		double InAcTo_Vhalf_ABh = 0.0;
		double InAcTo_Vk_ABh = 1.0;

		double x;
		double A = 0.4;  // was  0.5 for the LV, here change to 0.35, to increase APD, see what happens.

		if (state.cell_type() == 105 or state.cell_type() == 104 or state.cell_type() == 103)
		{
			A = 0.4;
		}

		if (state.cell_type() == 103 or state.cell_type() == 100) {
			A = 0.2;
		}
		x = (29.6 / 16.5 - 1) / (A + 29.6 / 16.5 * (1 - A) );
		InAcTo_Vhalf_ABh = (state.ab_index() - 1.0) * 4.0;                        // Apical Cells: -2.0; Basal Cells: 2.0;
		x = (4.5 / 3.4 - 1) / (A + 4.5 / 3.4 * (1 - A) );
		InAcTo_Vk_ABh    = 1.0 -  (state.ab_index() - A) * x; // Apical cells: 1.0; Basal Cells: 3.4/4.5
		x = (5.6 / 2.1 - 1) / (A + 5.6 / 2.1 * (1 - A) );
		// x = (516-318)/(516+318);
		x = (358.0 / 516.0 - 1) / (A + 358.0 / 516.0 * (1 - A) );

		double GKs_ABh = 1.0;
		GKs_ABh          =  1.0 -  (state.ab_index() - A) * x;   // Apical Cells: 1.0; Basal Cells: 2.1/5.6
		state.set_iks_current(GKs_ABh * get_Gks(state) * state.get_var(0,10) * state.get_var(0,10) * (state.get_vm() - state.eks()));  // with apicalbasal heterogeneity wit direct scaling factor, by haibo
	}

	void TNNPVent::IK1_current( CellState &state) 
	{
		double GK1 = 5.405;
		double Ak1 = 0.1 / (1. + exp(0.06 * (state.get_vm() - state.ek() - 200)));
		double Bk1 = (3.0 * exp(0.0002 * (state.get_vm() - state.ek() + 100)) + exp(0.1 * (state.get_vm() - state.ek() - 10))) / (1. + exp(-0.5 * (state.get_vm() - state.ek())));
		double rec_iK1 = Ak1 / (Ak1 + Bk1);
		if (state.cell_type() == 100)
			GK1 *= 1.2; // according to the ORd model supple, added by haibo
		state.set_ik1_current(GK1 * rec_iK1 * (state.get_vm() - state.ek()));
	}

	void TNNPVent::Ito_current( CellState &state) 
	{
		double InAcTo_Vhalf_ABh = 0.0;
		double InAcTo_Vk_ABh = 1.0;

		double x;
		double A = 0.4;  // was  0.5 for the LV, here change to 0.35, to increase APD, see what happens.

		if (state.cell_type() == 105 or state.cell_type() == 104 or state.cell_type() == 103)
		{
			A = 0.4;
		}

		if (state.cell_type() == 103 or state.cell_type() == 100) {
			A = 0.2;
		}
		x = (29.6 / 16.5 - 1) / (A + 29.6 / 16.5 * (1 - A) );
		InAcTo_Vhalf_ABh = (state.ab_index() - 1.0) * 4.0;                        // Apical Cells: -2.0; Basal Cells: 2.0;
		x = (4.5 / 3.4 - 1) / (A + 4.5 / 3.4 * (1 - A) );
		InAcTo_Vk_ABh    = 1.0 -  (state.ab_index() - A) * x; // Apical cells: 1.0; Basal Cells: 3.4/4.5
		x = (5.6 / 2.1 - 1) / (A + 5.6 / 2.1 * (1 - A) );
		// x = (516-318)/(516+318);
		x = (358.0 / 516.0 - 1) / (A + 358.0 / 516.0 * (1 - A) );
		double GTo_ABh = 1.0;
		GTo_ABh          = 1.0 -  (state.ab_index() - A) * x; // Apical cells: 1.0; Basal Cells: 16.5/29.6
		state.set_ito_current(GTo_ABh * get_Gto(state) * state.get_var(0,12) * state.get_var(0,11) * (state.get_vm() - state.ek()));
	}

	void TNNPVent::Isus_current(CellState &state)
	{

	}

	void TNNPVent::INa_current(CellState &state) 
	{
		double GNa = 14.838;
		if (state.cell_type() == 104 or state.cell_type() == 101)
		{
			// GNa *= 1.5;  // tranmural hereterogeneity includes INa
		}
		state.set_ina_current(get_Acidosis_factor(state) * GNa * state.get_var(0,5) * state.get_var(0,5) * state.get_var(0,5) * state.get_var(0,6) * state.get_var(0,7) * (state.get_vm() - state.ena()));
	}

	void TNNPVent::IbNa_current( CellState &state) 
	{
		state.set_ibna_current(TTCell_GbNa * (state.get_vm() - state.ena()));
	}

	void TNNPVent::ICaL_current( CellState &state)
	{
		state.set_ical_current(get_Acidosis_factor(state) * TTCell_GCaL * state.get_var(0,13) * state.get_var(0,14) * state.get_var(0,15) * state.get_var(0,16) * 4 * (state.get_vm() - 15) * (TTCell_F * TTCell_F / (TTCell_R * TTCell_T)) *
		              (0.25 * exp(2 * (state.get_vm() - 15) * TTCell_F / (TTCell_R * TTCell_T)) * state.get_var(0,2) - TTCell_Cao) / (exp(2 * (state.get_vm() - 15) * TTCell_F / (TTCell_R * TTCell_T)) - 1.));
	}

	void TNNPVent::IbCa_current( CellState &state) 
	{
		state.set_ibca_current(TTCell_GbCa * (state.get_vm() - state.eca()));
	}

	void TNNPVent::INaK_current( CellState &state) 
	{
		double rec_iNaK = (1. / (1. + 0.1245 * exp(-0.1 * state.get_vm() * TTCell_F / (TTCell_R * TTCell_T)) + 0.0353 * exp(-state.get_vm() * TTCell_F / (TTCell_R * TTCell_T))));
		state.set_inak_current((1 - 0.35 * state.is_index()) * TTCell_knak * (get_Ischemia_TTCell_Ko(state) / (get_Ischemia_TTCell_Ko(state) + TTCell_KmK)) * (state.get_var(0,3) / (state.get_var(0,3) + TTCell_KmNa)) * rec_iNaK);
	}

	void TNNPVent::INaCa_current( CellState &state) 
	{
		state.set_inaca_current(TTCell_knaca * (1. / (TTCell_KmNai * TTCell_KmNai * TTCell_KmNai + TTCell_Nao * TTCell_Nao * TTCell_Nao)) * (1. / (TTCell_KmCa + TTCell_Cao)) *
	               (1. / (1 + TTCell_ksat * exp((TTCell_n - 1) * state.get_vm() * TTCell_F / (TTCell_R * TTCell_T)))) *
	               (exp(TTCell_n * state.get_vm() * TTCell_F / (TTCell_R * TTCell_T)) * state.get_var(0,3) * state.get_var(0,3) * state.get_var(0,3) * TTCell_Cao -
	                exp((TTCell_n - 1) * state.get_vm() * TTCell_F / (TTCell_R * TTCell_T)) * TTCell_Nao * TTCell_Nao * TTCell_Nao * state.get_var(0,0) * 2.5));
	}

	void TNNPVent::IpCa_current( CellState &state) 
	{
		state.set_ipca_current(TTCell_GpCa * state.get_var(0,0) / (TTCell_KpCa + state.get_var(0,0)));
	}

	void TNNPVent::IpK_current( CellState &state) 
	{
		double rec_ipK = 1. / (1. + exp((25 - state.get_vm()) / 5.98));
		state.set_ipk_current(TTCell_GpK * rec_ipK * (state.get_vm() - state.ek()));
	}

	void TNNPVent::INaL_current( CellState &state) 
	{
		state.set_inal_current(get_GNaL(state) * state.get_var(0,18) * state.get_var(0,19) * (state.get_vm() - state.ena()));
		// state.set_inal_current((1 + 10 * state.is_index()) * get_GNaL(state) * state.get_var(0,18) * state.get_var(0,19) * (state.get_vm() - state.ena()));
	}

	void TNNPVent::IKATP_current( CellState &state) 
	{
		double GKATP              = 155.0;
		float f_atp               = 0.0055 * state.is_index();
		state.set_ikatp_current(GKATP * f_atp * pow(get_Ischemia_TTCell_Ko(state) / 5.4, 0.3) * (state.get_vm() - state.ek()) / (40 + 3.5 * exp(0.025 * state.get_vm())));
	}

	void TNNPVent::If_current(CellState &state)
	{

	}

	//====================================================================
	// Define the reversal potentials
	//====================================================================

	void TNNPVent::ENa_reversal( CellState &state) 
	{
		state.set_ena(TTCell_RTONF * (log((TTCell_Nao / state.get_var(0,3)))));
	}

	void TNNPVent::EK_reversal( CellState &state) 
	{
		state.set_ek(TTCell_RTONF * (log((get_Ischemia_TTCell_Ko(state) / state.get_var(0,4)))));
	}

	void TNNPVent::EKs_reversal( CellState &state) 
	{
		state.set_eks(TTCell_RTONF * (log((get_Ischemia_TTCell_Ko(state) + TTCell_pKNa * TTCell_Nao) / (state.get_var(0,4) + TTCell_pKNa * state.get_var(0,3)))));
		// std::cout << "Detailing calc EKs" << std::endl;
		// std::cout << "TTCell_RTONF " << TTCell_RTONF << std::endl;
		// std::cout << "get_Ischemia_TTCell_Ko(state) " << get_Ischemia_TTCell_Ko(state) << std::endl;
		// std::cout << "TTCell_pKNa " << TTCell_pKNa << std::endl;
		// std::cout << "TTCell_Nao " << TTCell_Nao << std::endl;
		// std::cout << "state.get_var(0,4)" << state.get_var(0,4) << std::endl;
		// std::cout << "state.get_var(0,3)" << state.get_var(0,3) << std::endl;
		// std::cout << "End detailing calc EKs" << std::endl;
	}

	void TNNPVent::ECa_reversal( CellState &state) 
	{
		state.set_eca(0.5 * TTCell_RTONF * (log((TTCell_Cao / state.get_var(0,0)))));
	}

	//====================================================================
	//====================================================================
	// Fill in residual for cell get_variables -
	//	grouped into single channel get_variables where appropriate
	//====================================================================
	//====================================================================

	inline void TNNPVent::IKr_Markov_residual( CellState &state,
									Vector<double> &residuals){

		//From the non-Markov model for IKr
		residuals[8] = state.get_var(1,8);
		residuals[9] = state.get_var(1,9);

		double TauKs_ABh = 1.0;
		double AcKs_Vhalf_ABh = 0.0;
		double AB_IKr_ActVhalf = 0.0;

		
		// TauKs_ABh = 1.0 -  (state.ab_index()-1.0) * x; // Apical Cells:1.0 ; Basal Cells: 516/358
		// TauKs_ABh = 1.0 - state.ab_index()*(516.0/318.0-1);// Apical Cells:1.0 ; Basal Cells: 516/358
		AB_IKr_ActVhalf = -(state.ab_index() - 1.0) * 4.9; // Apical cells: 2.45; Basal cells -2.45

		double m_epiMidRatio = 1.5;
		double epi_factor   = 1.8 * m_epiMidRatio;
		double endo_factor  = 1.8;
		double mcell_factor = 1.8;
		double wt_a1 =  2.172;
		double wt_b1 =  1.077;
		double wt_a2 =  0.00655   * exp(0.5 * 0.05547153 * (state.get_vm() - 36. - AB_IKr_ActVhalf));
		double wt_a  =  0.00555   * exp(0.05547153 * (state.get_vm() - 12. - AB_IKr_ActVhalf));
		double wt_b  =  0.002357  * exp(-0.036588 * (state.get_vm() - AB_IKr_ActVhalf));
		double wt_b2 = 0.65 * 0.0029357 * exp(0.69 * -0.02158 * (state.get_vm() - AB_IKr_ActVhalf));
		double wt_ai = 0.11 * 0.439     * exp(1.7 *  -0.02352 * (state.get_vm() + 25. - AB_IKr_ActVhalf)) * (4.5 / get_Ischemia_TTCell_Ko(state));
		double wt_bi = 0.4 *  0.656     * exp(      0.000942 * (state.get_vm() - AB_IKr_ActVhalf)) * ((pow((4.5 / get_Ischemia_TTCell_Ko(state)), 0.3)));
		double wt_mu = (wt_ai * wt_b2) / wt_bi;

		double wt_dC3 = (wt_b * state.get_var(0,23)) - (wt_a * state.get_var(0,20));
		double wt_dC2 = -((wt_b + wt_a1) * state.get_var(0,23)) + (wt_a * state.get_var(0,20)) + (wt_b1 * state.get_var(0,22));
		double wt_dC1 = -((wt_b1 + wt_a2 + wt_a2) * state.get_var(0,22)) + (wt_a1 * state.get_var(0,23)) + (wt_b2 * state.get_var(0,21)) + (wt_mu * state.get_var(0,24));
		double wt_dO  =  -((wt_b2 + wt_bi) * state.get_var(0,21)) + (wt_a2 * state.get_var(0,22)) + (wt_ai * state.get_var(0,24));
		double wt_dI  = -((wt_mu + wt_ai) * state.get_var(0,24)) + (wt_a2 * state.get_var(0,22)) + (wt_bi * state.get_var(0,21));

		residuals[21]  -= state.get_var(1,21) - wt_dO;
		residuals[22] -= state.get_var(1,22) - wt_dC1;
		residuals[23] -= state.get_var(1,23) - wt_dC2;
		residuals[20] -= state.get_var(1,20) - wt_dC3;
		residuals[24]  -= state.get_var(1,24) - wt_dI;
	}

	inline void TNNPVent::IKs_residual( CellState &state,
							Vector<double> &residuals){
		double TauKs_ABh = 1.0;
		double Xs_INF = 1. / (1. + exp((-5. - state.get_vm()) / 14.));

		double TAU_Xs = (1400. / (sqrt(1. + exp((5. - state.get_vm()) / 6)))) * (1. / (1. + exp((state.get_vm() - 35.) / 15.))) + 80;
		TAU_Xs = TauKs_ABh * TAU_Xs;
		residuals[10] -= state.get_var(1,10) - (Xs_INF - state.get_var(0,10))/TAU_Xs;
	}

	inline void TNNPVent::IK1_residual( CellState &state,
							Vector<double> &residuals){
	}

	inline void TNNPVent::Ito_residual( CellState &state,
							Vector<double> &residuals){
		residuals[11] -= state.get_var(1,11) - (get_S_INF(state) - state.get_var(0,11))/get_TAU_S(state);
		residuals[12] -= state.get_var(1,12) - (get_R_INF(state) - state.get_var(0,12))/get_TAU_R(state);
	}

	inline void TNNPVent::INa_residual( CellState &state,
							Vector<double> &residuals){
		double TAU_M = (1. / (1. + exp((-60. - state.get_vm()) / 5.))) * ( 0.1 / (1. + exp((state.get_vm() + 35.) / 5.)) + 0.10 / (1. + exp((state.get_vm() - 50.) / 200.))) ;
		double TAU_H;
		double TAU_J;
		double M_INF = 1. / ((1. + exp((-56.86 - state.get_vm()) / 9.03)) * (1. + exp((-56.86 - state.get_vm()) / 9.03)));
		if (state.get_vm() >= -40.) {
			double AH_1 = 0.;
			double BH_1 = (0.77 / (0.13 * (1. + exp(-(state.get_vm() + 10.66) / 11.1))));
			TAU_H = 1.0 / (AH_1 + BH_1);
		} else {
			double AH_2 = (0.057 * exp(-(state.get_vm() + 80.) / 6.8));
			double BH_2 = (2.7 * exp(0.079 * state.get_vm()) + (3.1e5) * exp(0.3485 * state.get_vm()));
			TAU_H = 1.0 / (AH_2 + BH_2);
		}
		double H_INF = 1. / ((1. + exp((state.get_vm() + 71.55) / 7.43)) * (1. + exp((state.get_vm() + 71.55) / 7.43)));
		if (state.get_vm() >= -40.) {
			double AJ_1 = 0.;
			double BJ_1 = (0.6 * exp((0.057) * state.get_vm()) / (1. + exp(-0.1 * (state.get_vm() + 32.))));
			TAU_J = 1.0 / (AJ_1 + BJ_1);
		} else {
			double AJ_2 = (((-2.5428e4) * exp(0.2444 * state.get_vm()) - (6.948e-6) *
			                exp(-0.04391 * state.get_vm())) * (state.get_vm() + 37.78) /
			               (1. + exp(0.311 * (state.get_vm() + 79.23))));
			double BJ_2 = (0.02424 * exp(-0.01052 * state.get_vm()) / (1. + exp(-0.1378 * (state.get_vm() + 40.14))));
			TAU_J = 1.0 / (AJ_2 + BJ_2);
		}
		double J_INF = H_INF;

		residuals[5] -= state.get_var(1,5) - (M_INF - state.get_var(0,5))/TAU_M;
		residuals[6] -= state.get_var(1,6) - (H_INF - state.get_var(0,6))/TAU_H;
		residuals[7] -= state.get_var(1,7) - (J_INF - state.get_var(0,7))/TAU_J;
	}

	inline void TNNPVent::IbNa_residual( CellState &state,
							Vector<double> &residuals){
		
	}

	inline void TNNPVent::ICaL_residual( CellState &state,
							Vector<double> &residuals){
		double D_INF = 1. / (1. + exp((-8 - state.get_vm()) / 7.5));
		double Ad = 1.4 / (1. + exp((-35 - state.get_vm()) / 13)) + 0.25;
		double Bd = 1.4 / (1. + exp((state.get_vm() + 5) / 5));
		double Cd = 1. / (1. + exp((50 - state.get_vm()) / 20));
		double TAU_D = Ad * Bd + Cd;
		double F_INF = 1. / (1. + exp((state.get_vm() + 20) / 7));
		double Af = 1102.5 * exp(-(state.get_vm() + 27) * (state.get_vm() + 27) / 225);
		double Bf = 200. / (1 + exp((13 - state.get_vm()) / 10.));
		double Cf = (180. / (1 + exp((state.get_vm() + 30) / 10))) + 20;
		double TAU_F = Af + Bf + Cf;

		double F2_INF = 0.67 / (1. + exp((state.get_vm() + 35) / 7)) + 0.33;
		double Af2 = 600 * exp(-(state.get_vm() + 25) * (state.get_vm() + 25) / 170);
		double Bf2 = 31 / (1. + exp((25 - state.get_vm()) / 10));
		double Cf2 = 16 / (1. + exp((state.get_vm() + 30) / 10));
		double TAU_F2 = Af2 + Bf2 + Cf2;

		residuals[13] -= state.get_var(1,13) - (D_INF - state.get_var(0,13))/TAU_D;
		residuals[14] -= state.get_var(1,14) - (F_INF - state.get_var(0,14))/TAU_F;
		residuals[15] -= state.get_var(1,15) - (F2_INF - state.get_var(0,15))/TAU_F2;

		double FCaSS_INF = 0.6 / (1 + (state.get_var(0,2) / 0.05) * (state.get_var(0,2) / 0.05)) + 0.4;
		double TAU_FCaSS = 80. / (1 + (state.get_var(0,2) / 0.05) * (state.get_var(0,2) / 0.05)) + 2.;

		residuals[16] -= state.get_var(1,16) - (FCaSS_INF - state.get_var(0,16))/TAU_FCaSS;
	}

	inline void TNNPVent::IbCa_residual( CellState &state,
							Vector<double> &residuals){
		
	}

	inline void TNNPVent::INaK_residual( CellState &state,
							Vector<double> &residuals){
		
	}

	inline void TNNPVent::INaCa_residual( CellState &state,
							Vector<double> &residuals){
		
	}

	inline void TNNPVent::IpCa_residual( CellState &state,
							Vector<double> &residuals){
		
	}

	inline void TNNPVent::IpK_residual( CellState &state,
							Vector<double> &residuals){
		
	}

	inline void TNNPVent::INaL_residual( CellState &state,
							Vector<double> &residuals){
		double TMP_INF = 1.0 / (1.0 + exp((-(state.get_vm() + 42.85)) / 5.264));
		double TAU_M = (1. / (1. + exp((-60. - state.get_vm()) / 5.))) * ( 0.1 / (1. + exp((state.get_vm() + 35.) / 5.)) + 0.10 / (1. + exp((state.get_vm() - 50.) / 200.))) ;

		residuals[18] -= state.get_var(1,18) - (TMP_INF - state.get_var(0,18))/TAU_M;

		TMP_INF = 1.0 / (1.0 + exp((state.get_vm() + 87.61) / 7.488));

		residuals[19] -= state.get_var(1,19) - (TMP_INF - state.get_var(0,19))/200.0;

		//NO residuals for states 18 and 19 in original code
	}

	inline void TNNPVent::IKATP_residual( CellState &state,
							Vector<double> &residuals){
		
	}


	//====================================================================
	// Fill in residual for intracellular concentrations
	//====================================================================

	inline void TNNPVent::Ca_i_residual( CellState &state,
							Vector<double> &residuals){
		double kCaSR = TTCell_maxsr - ((TTCell_maxsr - TTCell_minsr) / (1 + (TTCell_EC / state.get_var(0,1)) * (TTCell_EC / state.get_var(0,1))));
		double k1 = TTCell_k1_ / kCaSR;
		double k2 = TTCell_k2_ * kCaSR;
		double dRR = TTCell_k4 * (1 - state.get_var(0,17)) - k2 * state.get_var(0,2) * state.get_var(0,17);
		residuals[17] -= state.get_var(1,17) - dRR;
		double sOO = k1 * state.get_var(0,2) * state.get_var(0,2) * state.get_var(0,17) / (TTCell_k3 + k1 * state.get_var(0,2) * state.get_var(0,2));

		double Irel = TTCell_Vrel * sOO * (state.get_var(0,1) - state.get_var(0,2));
		double Ileak = TTCell_Vleak * (state.get_var(0,1) - state.get_var(0,0));
		double Iup = TTCell_Vmaxup / (1. + ((TTCell_Kup * TTCell_Kup) / (state.get_var(0,0) * state.get_var(0,0))));
		double Ixfer = TTCell_Vxfer * (state.get_var(0,2) - state.get_var(0,0));

		//CaSR

		residuals[1] -= state.get_var(1,1)*(
											(state.get_var(0,1)+TTCell_Kbufsr)*TTCell_Bufsr
											- state.get_var(0,1)*TTCell_Bufsr
											+ pow((state.get_var(0,1)+TTCell_Kbufsr),2.0)
										)
						-(Iup-Ileak-Irel)*pow((state.get_var(0,1)+TTCell_Kbufsr),2.0);

		//CaSS
		residuals[2] -= state.get_var(1,2)*(
											(state.get_var(0,2)+TTCell_Kbufss)*TTCell_Bufss
											- state.get_var(0,2)*TTCell_Bufss
											+ pow((state.get_var(0,2)+TTCell_Kbufss),2.0)
										)
						-(-Ixfer * (TTCell_Vc / TTCell_Vss) + Irel * (TTCell_Vsr / TTCell_Vss) + (-state.ical() * TTCell_inversevssF2 * TTCell_CAPACITANCE))*pow((state.get_var(0,2)+TTCell_Kbufsr),2.0);

		//Cai
		residuals[0] -= state.get_var(1,0)*(
											(state.get_var(0,0)+TTCell_Kbufc)*TTCell_Bufc
											- state.get_var(0,0)*TTCell_Bufc
											+ pow((state.get_var(0,0)+TTCell_Kbufc),2.0)
										)
						-((-(state.ibca() + state.ipca() - 2.0 * state.inaca()) * TTCell_inverseVcF2 * TTCell_CAPACITANCE) - (Iup - Ileak) * (TTCell_Vsr / TTCell_Vc) + Ixfer)*pow((state.get_var(0,0)+TTCell_Kbufsr),2.0);
	}

	inline void TNNPVent::Na_i_residual( CellState &state,
							Vector<double> &residuals){
		// std::cout << state.get_var(1,3) << "\t";
		// std::cout << state.ina() << "\t";
		// std::cout << state.inal() << "\t";
		// std::cout << state.ibna() << "\t";
		// std::cout << 3 * state.inak() << "\t";
		// std::cout << 3 * state.inaca() << "\t";
		// std::cout << TTCell_inverseVcF << "\t";
		// std::cout << TTCell_CAPACITANCE << std::endl;

		residuals[3] -= state.get_var(1,3) + (state.ina() + state.inal() + state.ibna() + 3 * state.inak() + 3 * state.inaca()) * TTCell_inverseVcF * TTCell_CAPACITANCE;

		// std::cout << residuals[3] << std::endl;
		// throw OomphLibError("This is an intentional error",
		//                        OOMPH_CURRENT_FUNCTION,
		//                        OOMPH_EXCEPTION_LOCATION);
	}

	inline void TNNPVent::K_i_residual( CellState &state,
							Vector<double> &residuals){
		residuals[4] -= state.get_var(1,4) + (/*Istim +*/ state.ik1() + state.ito() + state.ikr() + state.iks() - 2 * state.inak() + state.ipk() + state.ikatp()) * TTCell_inverseVcF * TTCell_CAPACITANCE;
	}




	void TNNPVent::fill_in_generic_residual_contribution_cell_base( CellState &state,
																	Vector<double> &residuals,
																	DenseMatrix<double> &jacobian,
																	unsigned flag)
	{	
		//Call the reversal functions
		ENa_reversal(state);
		EK_reversal(state);
		EKs_reversal(state);
		ECa_reversal(state);

		//Call the channel current functions
		IKr_current(state);
		IKs_current(state);
		IK1_current(state);
		Ito_current(state);
		INa_current(state);
		IbNa_current(state);
		ICaL_current(state);
		IbCa_current(state);
		INaK_current(state);
		INaCa_current(state);
		IpCa_current(state);
		IpK_current(state);
		INaL_current(state);
		IKATP_current(state);

		//Call the channel residual functions
		IKr_Markov_residual(state, residuals);
		IKs_residual(state, residuals);
		IK1_residual(state, residuals);
		Ito_residual(state, residuals);
		INa_residual(state, residuals);
		IbNa_residual(state, residuals);
		ICaL_residual(state, residuals);
		IbCa_residual(state, residuals);
		INaK_residual(state, residuals);
		INaCa_residual(state, residuals);
		IpCa_residual(state, residuals);
		IpK_residual(state, residuals);
		INaL_residual(state, residuals);
		IKATP_residual(state, residuals);

		//Call the current residual functions
		Ca_i_residual(state, residuals);
		Na_i_residual(state, residuals);
		K_i_residual(state, residuals);
	}

	inline void TNNPVent::return_initial_membrane_potential(double &v, const unsigned &cell_type){
		v = -86.2;
	}

	//Return the initial condition for the nth get_variable and cell_typeth cell type
	inline bool TNNPVent::return_initial_value(const unsigned &n, double &v, const unsigned &cell_type){
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

}