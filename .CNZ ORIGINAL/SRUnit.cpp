#include <cmath>
#include "HumanAtriaConstants.h"
// #include "cvode_struct.h"
// #include <cvode/cvode_dense.h>
/* Header files with a description of contents used */
// #include <cvode/cvode.h>                  /* prototypes for CVODE fcts., consts.  */
// #include <nvector/nvector_serial.h>       /* serial N_Vector types, fcts., macros */
#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ    */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ    */


#ifndef SRUnit_HPP
#define SRUnit_HPP

class SRUnit
{
public:
	SRUnit() {Initialise(); };
	~SRUnit() {};

	void Initialise();
	double Cai, Caj, CaSRj, CaSRN;
	double Casl;
	double dt, ICa_tot_j, ICa_tot_nj;
	double RyR_O, RyR_R, RyR_I, RyR_RC;
	double J_SR_Rel, J_SRleak, J_SERCA, J_SR_NJ_Diff;
	double Buffer[10];
	double Y[20], dY[20];
	double CSQN_Ca;


	// static const int NEQ = 18;  // using C++11 standard.
	void GB_RyR(double dt);
	void SERCA_UP();
	void UpdateCalciumConcentrations(double dt, double ICa_tot_j, double ICa_tot_nj);
	void SRUnitSingleTimeStep(double dt, double ICa_tot_j, double ICa_tot_nj) {
		Cai       = Y[0] ;
		Caj       = Y[1] ;
		Casl      = Y[2] ;
		CaSRj     = Y[3] ;
		CaSRN     = Y[4] ;
		RyR_O     = Y[5] ;
		RyR_I     = Y[6] ;
		RyR_R     = Y[7] ;
		Buffer[0] = Y[8] ;
		Buffer[1] = Y[9] ;
		Buffer[2] = Y[10];
		Buffer[3] = Y[11];
		Buffer[4] = Y[12];
		Buffer[5] = Y[13];
		Buffer[6] = Y[14];
		CSQN_Ca   = Y[15];

		GB_RyR(dt);
		SERCA_UP();
		UpdateCalciumConcentrations(dt, ICa_tot_j, ICa_tot_nj);
	}

	void lsoda_f(double t, double *y, double *ydot, void *data) {

		int i;
		for (i = 0; i < 16; i++)
			Y[i] = y[i];

		SRUnitSingleTimeStep(dt, ICa_tot_j, ICa_tot_nj);

		for (i = 0; i < 16; i++)
			ydot[i] = dY[i];
	}
};




void SRUnit::Initialise() {
	Cai = 8.597401e-5;
	Caj = 1.737475e-4;
	Casl = 1.737475e-4;
	CaSRj = CaSRN = 0.5; // SR content from Grandi model.

	RyR_O = 3.81797e-07;
	RyR_I = 7.35896e-08;
	RyR_R = 0.838402;

	Buffer[0] = 2.911916e-4; // 3: CaM (mM) (in Cytosolic_Ca_Buffers)
	Buffer[1] = 1.298754e-3;  // 4: Myo_c (mM) (in Cytosolic_Ca_Buffers)
	Buffer[2] = 1.381982e-1;    // 5: Myo_m (mM) (in Cytosolic_Ca_Buffers)
	Buffer[3] = 2.143165e-3;  // 6: SRB (mM) (in Cytosolic_Ca_Buffers)
	Buffer[4] = 1.078283e-1;    // 7: Tn_CHc (mM) (in Cytosolic_Ca_Buffers)
	Buffer[5] = 1.524002e-2;   // 8: Tn_CHm (mM) (in Cytosolic_Ca_Buffers)
	Buffer[6] = 8.773191e-3;  // 9: Tn_CL (mM) (in Cytosolic_Ca_Buffers)
	CSQN_Ca = 1.242988; // buffer in junctional SR;
	J_SR_Rel = J_SRleak = J_SERCA = J_SR_NJ_Diff = 0.0;

	Y[0]  = Cai;
	Y[1]  = Caj;
	Y[2]  = Casl;
	Y[3]  = CaSRj;
	Y[4]  = CaSRN;
	Y[5]  = RyR_O;
	Y[6]  = RyR_I;
	Y[7]  = RyR_R;
	Y[8]  = Buffer[0];
	Y[9]  = Buffer[1];
	Y[10] = Buffer[2];
	Y[11] = Buffer[3];
	Y[12] = Buffer[4];
	Y[13] = Buffer[5];
	Y[14] = Buffer[6];
	Y[15] = CSQN_Ca;
}



/*solving A Markov RyR Model using FE solver */
void SRUnit::GB_RyR(double dt) {

	// SR release
	double kCaSR = MaxSR - (MaxSR - MinSR) / (1.0 + pow(ec50SR / CaSRj, 2.5));
	double koSRCa = koCa / kCaSR;
	double kiSRCa = kiCa * kCaSR;

	// double RI = 1.0 - RyR_R - RyR_O - RyR_I;
	// double Caj_Square = Caj * Caj;
	// double dRyR_R = kim * RI - kiSRCa * Caj * RyR_R - (koSRCa * Caj_Square/*pow(Caj, 2.0)*/ * RyR_R - kom * RyR_O);
	// double dRyR_O = koSRCa * Caj_Square/*pow(Caj, 2.0)*/ * RyR_R - kom * RyR_O - (kiSRCa * Caj * RyR_O - kim * RyR_I);
	// double dRyR_I = kiSRCa * Caj * RyR_O - kim * RyR_I - (kom * RyR_I - koSRCa * Caj_Square /*pow(Caj, 2.0)*/ * RI);
	// kCaSR = MaxSR - (MaxSR - MinSR) / (1.0 + pow(ec50SR / Y[33], 2.5));
	// koSRCa = koCa / kCaSR;
	// kiSRCa = kiCa * kCaSR;
	double RI = 1.0 - RyR_R - RyR_O - RyR_I;
	double dRyR_R = kim * RI - kiSRCa * Caj * RyR_R - (koSRCa * pow(Caj, 2.0) * RyR_R - kom * RyR_O);
	double dRyR_O = koSRCa * pow(Caj, 2.0) * RyR_R - kom * RyR_O - (kiSRCa * Caj * RyR_O - kim * RyR_I);
	double dRyR_I = kiSRCa * Caj * RyR_O - kim * RyR_I - (kom * RyR_I - koSRCa * pow(Caj, 2.0) * RI);


	J_SR_Rel = ks * RyR_O * (CaSRj - Caj);
	// SR leak
	// update the state variables of RyR model using FE.
	J_SRleak = 5.348e-6 * (1.0 + 0.25 * AF) * (CaSRj - Cai);

	/*RyR_O += dt * dRyR_O;
	RyR_R += dt * dRyR_R;
	RyR_I += dt * dRyR_I;*/

	dY[5]  = dRyR_O;
	dY[6]  = dRyR_I;
	dY[7]  = dRyR_R;
}


void SRUnit::SERCA_UP() {
	// SERCA uptake
	// std::cout << Cai << std::endl;
	J_SERCA = Vmax_SRCaP * (pow(Cai / Kmf, hillSRCaP) - pow(CaSRN / Kmr, hillSRCaP))
	          / (1.0 + pow(Cai / Kmf, hillSRCaP) + pow(CaSRN / Kmr, hillSRCaP));
}

void SRUnit::UpdateCalciumConcentrations(double dt, double ICa_tot_j, double ICa_tot_sl) {

	// double J_Diff = (Caj - Cai) / 0.2;//
	double J_Diff = (Caj - Cai) * 46; //


	double I_pca_junc = 0.1 * pow(Q10SLCaP, Qpow) * IbarSLCaP * pow(Caj, 1.6) / (pow(KmPCa, 1.6) + pow(Caj, 1.6));
	double I_pca_sl = 0.9 * pow(Q10SLCaP, Qpow) * IbarSLCaP * pow(Casl, 1.6) / (pow(KmPCa, 1.6) + pow(Casl, 1.6));


	// notes: in the original Grandi Model, J_ca_juncsl/Vjunc = 46.32435221117347; By Haibo
	double dCasl =  -(ICa_tot_sl + I_pca_sl) * Cmem / (Vsl * 2.0 * Frdy) +  J_ca_juncsl / Vsl * (Caj - Casl) + J_ca_slmyo / Vsl * (Cai - Casl);
	double dCaj = -(ICa_tot_j + I_pca_junc) * Cmem / (Vjunc * 2.0 * Frdy)   // ICa influx
	              +  J_ca_juncsl / Vjunc * (Casl - Caj)
	              + J_SR_Rel * Vsr_junc / Vjunc + J_SRleak * Vmyo / Vjunc;

	// cytolic  Ca buffers:
	double dBuffer[7];
	dBuffer[0] = kon_cam * Cai * (Bmax_CaM - Buffer[0]) - koff_cam * Buffer[0];
	dBuffer[1] = kon_myoca * Cai * (Bmax_myosin - Buffer[1] - Buffer[2]) - koff_myoca * Buffer[1];
	dBuffer[2] = kon_myomg * Mgi * (Bmax_myosin - Buffer[1] - Buffer[2]) - koff_myomg * Buffer[2];
	dBuffer[3] = kon_sr * Cai * (Bmax_SR - Buffer[3]) - koff_sr * Buffer[3];
	dBuffer[4] = kon_tnchca * Cai * (Bmax_TnChigh - Buffer[4] - Buffer[5]) - koff_tnchca * Buffer[4];
	dBuffer[5] = kon_tnchmg * Mgi * (Bmax_TnChigh - Buffer[4] - Buffer[5]) - koff_tnchmg * Buffer[5];
	dBuffer[6] = kon_tncl * Cai * (Bmax_TnClow - Buffer[6]) - koff_tncl * Buffer[6];


	double J_Ca_Bulk_Buffer = dBuffer[6] + dBuffer[5] + dBuffer[4] + dBuffer[3] + dBuffer[2] + dBuffer[1] + dBuffer[0];

	double dCai = -J_SERCA * Vsr_network / Vmyo - J_Ca_Bulk_Buffer + J_ca_slmyo / Vmyo * (Casl - Cai); // ;


	double dCSQN_Ca = kon_csqn * CaSRj * (Bmax_Csqn - CSQN_Ca) - koff_csqn * CSQN_Ca;

	double J_SR_NJ_Diff = (CaSRN - CaSRj) / 100.0; // diffusion of network SR to junctional SR, in ORd, 100 ms

	double dCaSRj = - J_SR_Rel - dCSQN_Ca + J_SR_NJ_Diff;
	double dCaSRN = J_SERCA - J_SRleak - J_SR_NJ_Diff * Vsr_junc / Vsr_network;


	// updated Buffer using FE.
	// Buffer[0] += dt * dBuffer[0];
	// Buffer[1] += dt * dBuffer[1];
	// Buffer[2] += dt * dBuffer[2];
	// Buffer[3] += dt * dBuffer[3];
	// Buffer[4] += dt * dBuffer[4];
	// Buffer[5] += dt * dBuffer[5];
	// Buffer[6] += dt * dBuffer[6];

	// CSQN_Ca += dCSQN_Ca * dt;
	// Cai     += dt * dCai;
	// CaSRj   += dt * dCaSRj;
	// CaSRN = CaSRj;
	// // CaSRN   += dt * dCaSRN;
	// // CaSRj = CaSRN;
	// Caj += dt * dCaj;
	// Casl += dt * dCasl;


	dY[0]  = dCai;
	dY[1]  = dCaj;
	dY[2]  = dCasl;
	dY[3]  = dCaSRj;
	dY[4]  =  dCaSRN;
	dY[8]  = dBuffer[0];
	dY[9]  = dBuffer[1];
	dY[10] = dBuffer[2];
	dY[11] = dBuffer[3];
	dY[12] = dBuffer[4];
	dY[13] = dBuffer[5];
	dY[14] = dBuffer[6];
	// dY[15] = dY[16] = 0.0;
	dY[15] = dCSQN_Ca;
}


////////////////////////////////////////////////////////////////////////////////////////////
//THIS IS COMMENTED BECAUSE N_VECTOR IS A PART OF SUNDAILS, THEY ARE NOT USED ANYWHERE ELSE THOUGH SO I DON'T KNOW WHY THEY'RE HERE
//IF THEY TURN OUT TO BE IMPORTANT THEN REWRITE TO USE STD::VECTOR
////////////////////////////////////////////////////////////////////////////////////////////


// static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data) {

// 	SRUnit * data = (SRUnit*) user_data;
// 	int i;
// 	for (i = 0; i < 16; i++)
// 		data->Y[i] = Ith(y, i + 1);
// 	data->SRUnitSingleTimeStep(data->dt, data->ICa_tot_j, data->ICa_tot_nj);
// 	for (i = 0; i < 16; i++)
// 		Ith(ydot, i + 1) = data->dY[i];
// 	return 0;

// }




// void frk2f( double dt, double t, double *y1, double *ydot, SRUnit & caunit ) {

// 	int i;
// 	for (i = 0; i < 16; i++)
// 		caunit.Y[i] = y1[i];

// 	caunit.SRUnitSingleTimeStep(dt, caunit.ICa_tot_j, caunit.ICa_tot_nj);

// 	for (i = 0; i < 16; i++)
// 		ydot[i] = caunit.dY[i];
// }


// void InitialiseStates(N_Vector y, SRUnit &data) {
// 	for (int i = 0; i < 16; i++)
// 		Ith(y, i + 1) = data.Y[i];
// }


// void InitialiseStates(double *y, SRUnit &data) {
// 	for (int i = 0; i < 16; i++)
// 		y[i+1] = data.Y[i];
// }



#endif