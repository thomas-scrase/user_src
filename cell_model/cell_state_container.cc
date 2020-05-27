#include "cell_state_container.h"

namespace oomph{
	//Empty constructor and destructor
	CellState::CellState(){	}
	// CellState::~CellState(){	}
	
	//Access functions for the base class member data
	void CellState::set_vars(const DenseMatrix<double> &new_vars){
		Vars.resize(new_vars.nrow(), new_vars.ncol());
		for(unsigned i = 0; i < new_vars.nrow(); i++){
			for(unsigned j = 0; j < new_vars.ncol(); j++){
				Vars(i,j) = new_vars(i,j);
			}
		}
	}
	double CellState::var(const unsigned &d, const unsigned &n) const {return Vars(d,n);}

	//time_stepper_weights
	void CellState::set_time_stepper_weights(const DenseMatrix<double> &new_time_stepper_weights){
		Time_Stepper_Weights.resize(new_time_stepper_weights.nrow(), new_time_stepper_weights.ncol());
		for(unsigned i = 0; i < new_time_stepper_weights.nrow(); i++){
			for(unsigned j = 0; j < new_time_stepper_weights.ncol(); j++){
				Time_Stepper_Weights(i,j) = new_time_stepper_weights(i,j);
			}
		}
	}
	double CellState::time_stepper_weights(const unsigned &d, const unsigned &n) const {return Time_Stepper_Weights(d,n);}

	//Vm
	void CellState::set_vm(const double &new_vm){Vm = new_vm;}
	double CellState::vm() const {return Vm;}

	//stress
	void CellState::set_stress(const double &new_stress){Stress = new_stress;}
	double CellState::stress() const {return Stress;}

	//Na_o
	void CellState::set_nao(const double &new_nao){Na_o = new_nao;}
	double CellState::na_o() const {return Na_o;}

	//K_o
	void CellState::set_ko(const double &new_ko){K_o = new_ko;}
	double CellState::k_o() const {return K_o;}

	//Ca_o
	void CellState::set_cao(const double &new_cao){Ca_o = new_cao;}
	double CellState::ca_o() const {return Ca_o;}

	//cell_type
	void CellState::set_cell_type(const unsigned &new_cell_type){Cell_Type = new_cell_type;}
	unsigned CellState::cell_type() const {return Cell_Type;}

	//fibrosis_type
	void CellState::set_fibrosis_type(const unsigned &new_fibrosis_type){Fibrosis_Type = new_fibrosis_type;}
	unsigned CellState::fibrosis_type() const {return Fibrosis_Type;}

	//AB_index
	void CellState::set_ab_index(const double &new_ab_index){AB_index = new_ab_index;}
	double CellState::ab_index() const {return AB_index;}

	//RV_index
	void CellState::set_rv_index(const double &new_rv_index){RV_index = new_rv_index;}
	double CellState::rv_index() const {return RV_index;}

	//IS_index
	void CellState::set_is_index(const double &new_is_index){IS_index = new_is_index;}
	double CellState::is_index() const {return IS_index;}

	//Currents
	void CellState::set_ikr_current(const double &new_ikr_current){IKr_current = new_ikr_current;}
	double CellState::ikr() const {return IKr_current;}

	void CellState::set_iks_current(const double &new_iks_current){IKs_current = new_iks_current;}
	double CellState::iks() const {return IKs_current;}

	void CellState::set_ik1_current(const double &new_ik1_current){IK1_current = new_ik1_current;}
	double CellState::ik1() const {return IK1_current;}
	
	void CellState::set_ito_current(const double &new_ito_current){Ito_current = new_ito_current;}
	double CellState::ito() const {return Ito_current;}

	void CellState::set_ina_current(const double &new_ina_current){INa_current = new_ina_current;}
	double CellState::ina() const {return INa_current;}

	void CellState::set_ibna_current(const double &new_ibna_current){IbNa_current = new_ibna_current;}
	double CellState::ibna() const {return IbNa_current;}

	void CellState::set_ical_current(const double &new_ical_current){ICaL_current = new_ical_current;}
	double CellState::ical() const {return ICaL_current;}

	void CellState::set_ibca_current(const double &new_ibca_current){IbCa_current = new_ibca_current;}
	double CellState::ibca() const {return IbCa_current;}

	void CellState::set_inak_current(const double &new_inak_current){INaK_current = new_inak_current;}
	double CellState::inak() const {return INaK_current;}

	void CellState::set_inaca_current(const double &new_inaca_current){INaCa_current = new_inaca_current;}
	double CellState::inaca() const {return INaCa_current;}

	void CellState::set_ipca_current(const double &new_ipca_current){IpCa_current = new_ipca_current;}
	double CellState::ipca() const {return IpCa_current;}

	void CellState::set_ipk_current(const double &new_ipk_current){IpK_current = new_ipk_current;}
	double CellState::ipk() const {return IpK_current;}

	void CellState::set_inal_current(const double &new_inal_current){INaL_current = new_inal_current;}
	double CellState::inal() const {return INaL_current;}

	void CellState::set_ikatp_current(const double &new_ikatp_current){IKATP_current = new_ikatp_current;}
	double CellState::ikatp() const {return IKATP_current;}

	void CellState::set_icat_current(const double &new_icat_current){ICaT_current = new_icat_current;}
	double CellState::icat() const {return ICaT_current;}
	
	void CellState::set_isac_ca_current(const double &new_isac_ca_current){ISAC_Ca_current = new_isac_ca_current;}
	double CellState::isac_ca() const {return ISAC_Ca_current;}
	
	void CellState::set_ikv_current(const double &new_ikv_current){IKv_current = new_ikv_current;}
	double CellState::ikv() const {return IKv_current;}
	
	void CellState::set_ik1f_current(const double &new_ik1f_current){IK1f_current = new_ik1f_current;}
	double CellState::ik1f() const {return IK1f_current;}
	
	void CellState::set_inakf_current(const double &new_inakf_current){INaKf_current = new_inakf_current;}
	double CellState::inakf() const {return INaKf_current;}
	
	void CellState::set_ibnaf_current(const double &new_ibnaf_current){IbNaf_current = new_ibnaf_current;}
	double CellState::ibnaf() const {return IbNaf_current;}
	
	void CellState::set_igap_current(const double &new_igap_current){IGap_current = new_igap_current;}
	double CellState::igap() const {return IGap_current;}
	
	void CellState::set_ikur_current(const double &new_ikur_current){IKur_current = new_ikur_current;}
	double CellState::ikur() const {return IKur_current;}
	
	void CellState::set_ibk_current(const double &new_ibk_current){IbK_current = new_ibk_current;}
	double CellState::ibk() const {return IbK_current;}
	
	void CellState::set_isac_k_current(const double &new_isac_k_current){ISAC_K_current = new_isac_k_current;}
	double CellState::isac_k() const {return ISAC_K_current;}
	
	void CellState::set_icap_current(const double &new_icap_current){ICap_current = new_icap_current;}
	double CellState::icap() const {return ICap_current;}

	void CellState::set_isac_na_current(const double &new_isac_na_current){ISAC_Na_current = new_isac_na_current;}
	double CellState::isac_na() const {return ISAC_Na_current;}

	void CellState::set_iab_current(const double &new_iab_current){Iab_current = new_iab_current;}
	double CellState::iab() const {return Iab_current;}


	//Reversal potentials
	void CellState::set_ena(const double &new_ena){ENa = new_ena;}
	double CellState::ena() const {return ENa;}

	void CellState::set_ek(const double &new_ek){EK = new_ek;}
	double CellState::ek() const {return EK;}

	void CellState::set_eks(const double &new_eks){EKs = new_eks;}
	double CellState::eks() const {return EKs;}

	void CellState::set_eca(const double &new_eca){ECa = new_eca;}
	double CellState::eca() const {return ECa;}

	void CellState::set_ekf(const double &new_ekf){EKf = new_ekf;}
	double CellState::ekf() const {return EKf;}

	void CellState::set_enaf(const double &new_enaf){ENaf = new_enaf;}
	double CellState::enaf() const {return ENaf;}
	
	//Force build CellState
	class CellState;
}