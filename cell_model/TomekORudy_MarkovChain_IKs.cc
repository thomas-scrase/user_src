#include "TomekORudy_MarkovChain_IKs.h"

namespace oomph{

TomekORudyVentMarkovChainIKs::TomekORudyVentMarkovChainIKs() : CellModelBase()
{
	Intrinsic_dt = 0.02;

	nao =140.0;
	ko = 5.0;
	cao = 1.8;


	R=8314.0;
	T=310.0;
	F=96485.0;

	//Added for efficiency
	FNORT = F/(R*T);

	L=0.01;
	rad=0.0011;
	vcell=1000*3.14*rad*rad*L;
	Ageo=2*3.14*rad*rad+2*3.14*rad*L;
	Acap=2*Ageo;
	vmyo=0.68*vcell;
	vnsr=0.0552*vcell;
	vjsr=0.0048*vcell;
	vss=0.02*vcell;
}

bool TomekORudyVentMarkovChainIKs::compatible_cell_types(const unsigned& cell_type)
{
	switch(cell_type){
		case 100 : return true;
		case 101 : return true;
		case 102 : return true;
	}
	return false;
}

inline void TomekORudyVentMarkovChainIKs::return_initial_membrane_potential(double &v, const unsigned &cell_type)
{
	// v = -88.7638;
	switch(cell_type){
		case 100:
			v = -90.7456332496277;
			break;
		case 101:
			v = -91.3391803967798;
			break;
		case 102:
			v = -89.7480828168753;
			break;
	}
}

inline bool TomekORudyVentMarkovChainIKs::return_initial_state_variable(const unsigned&n, double& v, const unsigned& cell_type)
{	
	double STATES[3][59] = 
	{
		{13.4006197127101,	13.4009366375605,	152.363889991411,	152.363844281265,	6.62181569329109e-05,	5.74992101186664e-05,	1.80679354108253,	1.80504703846161,	0.00052533142244124,	0.731365613931879,	0.864514804478946,	0.864457115015362,	0.864352715733040,	0.000111796946200753,	0.591653600903076,	0.347681236790268,	0.000832040787904695,	0.999724237643709,	0.999723490323258,	0.000423912054921203,	0.999724237767365,	0.999724054768924,	-2.48652723154831e-36,	0.999999995656968,	0.951060213798798,	0.999999995656976,	0.999937692440505,	0.999988620215189,	0.000304952292598877,	0.000527266836293694,	0.999999995656685,	0.999999995656399,
			1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			6.77882673136480e-25,	0.0127354136594547,	0.998473340797544,	0.000739304520671507,	0.000602907876760312,	0.000178778329437106,	5.67825481839631e-06,	-1.58194108230975e-23,	34.3172099045416,	34.3171879616709},

		{15.9486733631902,	15.9492201219199,	156.713067183465,	156.713012917889,	8.29757614119249e-05,	6.64218738239930e-05,	2.01222465016069,	2.01641461550468,	0.000461956534829239,	0.747897241513971,	0.873907657623256,	0.873784149746068,	0.873537507470071,	9.98770874860868e-05,	0.598611759466370,	0.333989914508018,	0.000799404188333773,	0.999751418648699,	0.570253781914665,	0.000407277681198651,	0.999751424044492,	0.635192659820007,	-8.33460419955820e-30,	0.999999996299523,	0.918358695116332,	0.999999996300839,	0.999754043234341,	0.999974289100453,	0.000533651965198686,	0.00125786136830066,	0.999999996207845,	0.999999996298728,
			1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			-1.30048587380675e-21,	0.0201882046405105,	0.998345129712001,	0.000708645966770289,	0.000591004733696651,	0.000344573276330018,	1.06423014130810e-05,	-7.61071425848620e-20,	48.9127681450109,	48.9127355777872},

		{12.3973564703237,	12.3976953330256,	147.711474407329,	147.711433360862,	7.45348093580908e-05,	6.49734059262396e-05,	1.52800059077496,	1.52569253436820,	0.000651715361520436,	0.701845427662345,	0.847326693399099,	0.847165668535608,	0.846901438940383,	0.000135120261194560,	0.556601650459093,	0.311549070692613,	0.000889925901945137,	0.999671576537361,	0.598890849399207,	0.000453416519958675,	0.999671583620741,	0.662069197871101,	1.58884121071311e-31,	0.999999994310640,	0.940179072110115,	0.999999994310551,	0.999901376748386,	0.999984619576572,	0.000489937844976504,	0.000832600925714172,	0.999999994343118,	0.999999994312043,
			1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			1.80824774078763e-22,	0.0109502564192484,	0.998251134783098,	0.000793602046257213,	0.000653214338491609,	0.000292244895058368,	9.80408293117642e-06,	4.35860784390795e-21,	29.2069793887423,	29.2069557375901}
	};

	if(compatible_cell_types(cell_type)){
		if(n<59){
			v = STATES[cell_type-100][n];
			return true;
		}
		else{std::cout << "Failed at variable index check " << n << " is not within the expected range" << std::endl;}
	}
	else{std::cout << "Failed at cell type check " << cell_type << " is not supported" << std::endl;}
	return true;
}

inline void TomekORudyVentMarkovChainIKs::extract_black_box_parameters_TomekORudyVent(CellState &Cellstate)
{

}

void TomekORudyVentMarkovChainIKs::explicit_timestep(CellState &Cellstate, Vector<double> &new_state)
{

	//Unpack the membrane potential, cell type, and time step
	double v         = Cellstate.get_vm();
	double dt = Cellstate.get_dt();
	unsigned celltype = Cellstate.get_cell_type();

	//give names to the state vector values
	// v=X(1);
	const double nai=new_state[0];
	const double nass=new_state[1];
	const double ki=new_state[2];
	const double kss=new_state[3];
	const double cai=new_state[4];
	const double cass=new_state[5];
	const double cansr=new_state[6];
	const double cajsr=new_state[7];
	const double m=new_state[8];
	const double hp=new_state[9];
	const double h=new_state[10];
	const double j=new_state[11];

	const double jp=new_state[12];
	const double mL=new_state[13];
	const double hL=new_state[14];
	const double hLp=new_state[15];
	const double a=new_state[16];
	const double iF=new_state[17];
	const double iS=new_state[18];
	const double ap=new_state[19];
	const double iFp=new_state[20];
	const double iSp=new_state[21];
	// ical
	const double d=new_state[22];
	const double ff=new_state[23];
	const double fs=new_state[24];
	const double fcaf=new_state[25];
	const double fcas=new_state[26];
	const double jca=new_state[27];
	const double nca=new_state[28];
	const double nca_i=new_state[29];
	const double ffp=new_state[30];
	const double fcafp=new_state[31];
	
	//Old IKs variables	
	// const double xs1=new_state[32];
	// const double xs2=new_state[33];

	//The new IKs variables
	const double iks_o2 = new_state[32];
	const double iks_o1 = new_state[33];
	const double iks_c1 = new_state[34];
	const double iks_c2 = new_state[35];
	const double iks_c3 = new_state[36];
	const double iks_c4 = new_state[37];
	const double iks_c5 = new_state[38];
	const double iks_c6 = new_state[39];
	const double iks_c7 = new_state[40];
	const double iks_c8 = new_state[41];
	const double iks_c9 = new_state[42];
	const double iks_c10 = new_state[43];
	const double iks_c11 = new_state[44];
	const double iks_c12 = new_state[45];
	const double iks_c13 = new_state[46];
	const double iks_c14 = new_state[47];
	const double iks_c15 = new_state[48];



	const double Jrel_np=new_state[49];
	const double CaMKt=new_state[50];
	// MM ICaL states
	const double ikr_c0 = new_state[51];
	const double ikr_c1 = new_state[52];
	const double ikr_c2 = new_state[53];
	const double ikr_o = new_state[54];
	const double ikr_i = new_state[55];
	const double Jrel_p=new_state[56];
	// const double cli = new_state[57];
	const double clss = new_state[58];
	const double cli = 24; // Intracellular Cl  [mM]
	const double clo = 150;  // Extracellular Cl  [mM]


	//Some parameters which are passed to the timestepper in the original function
	//	These probably represent mutation/disease/pathological effects
	const double INa_Multiplier = 1.0;
	const double INaL_Multiplier = 1.0;
	const double Ito_Multiplier = 1.0;
	const double ICaL_PCaMultiplier = 1.0;
	const double ICaL_fractionSS = 0.8;
	const double IKs_Multiplier = Cellstate.get_black_box_nodal_parameters(0);
	const double IKr_Multiplier = 1.0;
	const double IK1_Multiplier = 1.0;
	const double INaCa_Multiplier = 1.0;
	const double INaCa_fractionSS = 0.35;
	const double INaK_Multiplier = 1.0;
	const double IKb_Multiplier = 1.0;
	const double INab_Multiplier = 1.0;
	const double ICab_Multiplier = 1.0;
	const double IpCa_Multiplier = 1.0;
	const double IClCa_Multiplier = 1.0;
	const double IClb_Multiplier = 1.0;
	const double Jrel_Multiplier = 1.0;
	const double Jup_Multiplier = 1.0;


	// CaMK constants
	const double KmCaMK=0.15;

	const double aCaMK=0.05;
	const double bCaMK=0.00068;
	const double CaMKo=0.05;
	const double KmCaM=0.0015;
	// update CaMK
	const double CaMKb=CaMKo*(1.0-CaMKt)/(1.0+KmCaM/cass);
	const double CaMKa=CaMKb+CaMKt;
	const double dCaMKt=aCaMK*CaMKb*(CaMKb+CaMKt)-bCaMK*CaMKt;


	// reversal potentials
	const double ENa=(R*T/F)*log(nao/nai);
	const double EK=(R*T/F)*log(ko/ki);
	const double PKNa=0.01833;
	const double EKs=(R*T/F)*log((ko+PKNa*nao)/(ki+PKNa*nai));

	// convenient shorthand calculations
	const double vffrt=v*F*F/(R*T);
	const double vfrt=v*F/(R*T);
	const double frt = F/(R*T);

	
	const double fINap=(1.0/(1.0+KmCaMK/CaMKa));
	const double fINaLp=(1.0/(1.0+KmCaMK/CaMKa));
	const double fItop=(1.0/(1.0+KmCaMK/CaMKa));
	const double fICaLp=(1.0/(1.0+KmCaMK/CaMKa));

	// INa
	// [INa, dm, dh, dhp, dj, djp] = getINa_Grandi(v, m, h, hp, j, jp, fINap, ENa, INa_Multiplier);

	// The Grandi implementation updated with INa phosphorylation.
	// m gate
	const double mss = 1 / pow((1 + exp( -(56.86 + v) / 9.03 )), 2.0);
	const double taum = 0.1292 * exp(-pow((v+45.79)/15.54,2.0)) + 0.06487 * exp(-pow((v-4.823)/51.12,2.0));
	const double dm = (mss - m) / taum;

	// h gate
	const double ah = (v >= -40) * (0) + (v < -40) * (0.057 * exp( -(v + 80) / 6.8 ));
	const double bh = (v >= -40) * (0.77 / (0.13*(1 + exp( -(v + 10.66) / 11.1 )))) + (v < -40) * ((2.7 * exp( 0.079 * v) + 3.1*1e5 * exp(0.3485 * v)));
	const double tauh = 1 / (ah + bh);
	const double hss = 1 / (pow(1 + exp( (v + 71.55)/7.43 ),2.0));
	const double dh = (hss - h) / tauh;
	// j gate
	const double aj = (v >= -40) * (0) +(v < -40) * (((-2.5428 * 1e4*exp(0.2444*v) - 6.948*1e-6 * exp(-0.04391*v)) * (v + 37.78)) / (1 + exp( 0.311 * (v + 79.23) )));
	const double bj = (v >= -40) * ((0.6 * exp( 0.057 * v)) / (1 + exp( -0.1 * (v + 32) ))) + (v < -40) * ((0.02424 * exp( -0.01052 * v )) / (1 + exp( -0.1378 * (v + 40.14) )));
	const double tauj = 1 / (aj + bj);
	const double jss = 1 / (pow(1 + exp( (v + 71.55)/7.43 ),2.0));
	const double dj = (jss - j) / tauj;

	// h phosphorylated
	const double hssp = 1 / (pow(1 + exp( (v + 71.55 + 6)/7.43 ),2.0));
	const double dhp = (hssp - hp) / tauh;
	// j phosphorylated
	const double taujp = 1.46 * tauj;
	const double djp = (jss - jp) / taujp;

	const double GNa = 11.7802;
	const double INa=INa_Multiplier * GNa*(v-ENa)*pow(m,3.0)*((1.0-fINap)*h*j+fINap*hp*jp);




	// INaL
	// [INaL,dmL,dhL,dhLp] = getINaL_ORd2011(v, mL, hL, hLp, fINaLp, ENa, celltype, INaL_Multiplier);

	//calculate INaL
	const double mLss=1.0/(1.0+exp((-(v+42.85))/5.264));
	const double tm = 0.1292 * exp(-pow((v+45.79)/15.54,2.0)) + 0.06487 * exp(-pow((v-4.823)/51.12,2.0)); 
	const double tmL=tm;
	const double dmL=(mLss-mL)/tmL;
	const double hLss=1.0/(1.0+exp((v+87.61)/7.488));
	const double thL=200.0;
	const double dhL=(hLss-hL)/thL;
	const double hLssp=1.0/(1.0+exp((v+93.81)/7.488));
	const double thLp=3.0*thL;
	const double dhLp=(hLssp-hLp)/thLp;
	double GNaL=0.0279 * INaL_Multiplier;
	if(celltype==100){
	    GNaL=GNaL*0.6;
	}

	const double INaL=GNaL*(v-ENa)*mL*((1.0-fINaLp)*hL+fINaLp*hLp);



	// ITo
	// [Ito,da,diF,diS,dap,diFp, diSp] = getITo_ORd2011(v, a, iF, iS, ap, iFp, iSp, fItop, EK, celltype, Ito_Multiplier);

	//calculate Ito
	const double ass=1.0/(1.0+exp((-(v-14.34))/14.82));
	const double ta=1.0515/(1.0/(1.2089*(1.0+exp(-(v-18.4099)/29.3814)))+3.5/(1.0+exp((v+100.0)/29.3814)));
	const double da=(ass-a)/ta;
	const double iss=1.0/(1.0+exp((v+43.94)/5.711));
	double delta_epi;
	if(celltype==100){
	    delta_epi=1.0-(0.95/(1.0+exp((v+70.0)/5.0)));
	}
	else{
	    delta_epi=1.0;
	}
	double tiF=4.562+1/(0.3933*exp((-(v+100.0))/100.0)+0.08004*exp((v+50.0)/16.59));
	double tiS=23.62+1/(0.001416*exp((-(v+96.52))/59.05)+1.780e-8*exp((v+114.1)/8.079));
	tiF=tiF*delta_epi;
	tiS=tiS*delta_epi;
	const double AiF=1.0/(1.0+exp((v-213.6)/151.2));
	const double AiS=1.0-AiF;
	const double diF=(iss-iF)/tiF;
	const double diS=(iss-iS)/tiS;
	const double i=AiF*iF+AiS*iS;
	const double assp=1.0/(1.0+exp((-(v-24.34))/14.82));
	const double dap=(assp-ap)/ta;
	const double dti_develop=1.354+1.0e-4/(exp((v-167.4)/15.89)+exp(-(v-12.23)/0.2154));
	const double dti_recover=1.0-0.5/(1.0+exp((v+70.0)/20.0));
	const double tiFp=dti_develop*dti_recover*tiF;
	const double tiSp=dti_develop*dti_recover*tiS;
	const double diFp=(iss-iFp)/tiFp;
	const double diSp=(iss-iSp)/tiSp;
	const double ip=AiF*iFp+AiS*iSp;
	double Gto=0.16 * Ito_Multiplier;
	if(celltype==100){
	    Gto=Gto*2.0;
	}
	else if(celltype==101){
	    Gto=Gto*2.0;
	}

	const double Ito=Gto*(v-EK)*((1.0-fItop)*a*i+fItop*ap*ip);




	// ICaL
	// [ICaL_ss,ICaNa_ss,ICaK_ss,ICaL_i,ICaNa_i,ICaK_i,dd,dff,dfs,dfcaf,dfcas,djca,dnca,dnca_i,...
	//     dffp,dfcafp, PhiCaL_ss, PhiCaL_i, gammaCaoMyo, gammaCaiMyo] = getICaL_ORd2011_jt(v, d,ff,fs,fcaf,fcas,jca,nca,nca_i,ffp,fcafp,...
	//     fICaLp, cai, cass, cao, nai, nass, nao, ki, kss, ko, cli, clo, clss, celltype, ICaL_fractionSS, ICaL_Multiplier );

	//physical constants
	// const double vffrt=v*F*F/(R*T);
	// const double vfrt=v*F/(R*T);

	//calculate ICaL, ICaNa, ICaK

	double dss=1.0763*exp(-1.0070*exp(-0.0829*(v)));  // magyar
	if(v >31.4978){// activation cannot be greater than 1
	    dss = 1;
	}


	const double td= 0.6+1.0/(exp(-0.05*(v+6.0))+exp(0.09*(v+14.0)));

	const double dd=(dss-d)/td;
	const double fss=1.0/(1.0+exp((v+19.58)/3.696));
	const double tff=7.0+1.0/(0.0045*exp(-(v+20.0)/10.0)+0.0045*exp((v+20.0)/10.0));
	const double tfs=1000.0+1.0/(0.000035*exp(-(v+5.0)/4.0)+0.000035*exp((v+5.0)/6.0));
	const double Aff=0.6;
	const double Afs=1.0-Aff;
	const double dff=(fss-ff)/tff;
	const double dfs=(fss-fs)/tfs;
	const double f=Aff*ff+Afs*fs;
	const double fcass=fss;
	const double tfcaf=7.0+1.0/(0.04*exp(-(v-4.0)/7.0)+0.04*exp((v-4.0)/7.0));
	const double tfcas=100.0+1.0/(0.00012*exp(-v/3.0)+0.00012*exp(v/7.0));

	const double Afcaf=0.3+0.6/(1.0+exp((v-10.0)/10.0));

	const double Afcas=1.0-Afcaf;
	const double dfcaf=(fcass-fcaf)/tfcaf;
	const double dfcas=(fcass-fcas)/tfcas;
	const double fca=Afcaf*fcaf+Afcas*fcas;

	const double tjca = 72.5;
	const double jcass = 1.0/(1.0+exp((v+18.08)/(2.7916)));   
	const double djca=(jcass-jca)/tjca;
	const double tffp=2.5*tff;
	const double dffp=(fss-ffp)/tffp;
	const double fp=Aff*ffp+Afs*fs;
	const double tfcafp=2.5*tfcaf;
	const double dfcafp=(fcass-fcafp)/tfcafp;
	const double fcap=Afcaf*fcafp+Afcas*fcas;

	// SS nca
	const double Kmn=0.002;
	const double k2n=500.0;
	const double km2n=jca*1;
	const double anca=1.0/(k2n/km2n+pow(1.0+Kmn/cass,4.0));
	const double dnca=anca*k2n-nca*km2n;

	// myoplasmic nca
	const double anca_i = 1.0/(k2n/km2n+pow(1.0+Kmn/cai,4.0));
	const double dnca_i = anca_i*k2n-nca_i*km2n;

	// SS driving force
	// clo = 150; cli = 24;
	double Io = 0.5*(nao + ko + clo + 4*cao)/1000 ; // ionic strength outside. /1000 is for things being in micromolar
	double Ii = 0.5*(nass + kss + clss + 4*cass)/1000 ; // ionic strength outside. /1000 is for things being in micromolar
	// The ionic strength is too high for basic DebHuc. We'll use Davies
	double dielConstant = 74; // water at 37�.
	double constA = 1.82*1e6*pow(dielConstant*T,-1.5);

	double gamma_cai = exp(-constA * 4 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
	double gamma_cao = exp(-constA * 4 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));
	double gamma_nai = exp(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
	double gamma_nao = exp(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));
	double gamma_ki = exp(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
	double gamma_kao = exp(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));


	double PhiCaL_ss =  4.0*vffrt*(gamma_cai*cass*exp(2.0*vfrt)-gamma_cao*cao)/(exp(2.0*vfrt)-1.0);
	double PhiCaNa_ss =  1.0*vffrt*(gamma_nai*nass*exp(1.0*vfrt)-gamma_nao*nao)/(exp(1.0*vfrt)-1.0);
	double PhiCaK_ss =  1.0*vffrt*(gamma_ki*kss*exp(1.0*vfrt)-gamma_kao*ko)/(exp(1.0*vfrt)-1.0);

	// Myo driving force
	Io = 0.5*(nao + ko + clo + 4*cao)/1000 ; // ionic strength outside. /1000 is for things being in micromolar
	Ii = 0.5*(nai + ki + cli + 4*cai)/1000 ; // ionic strength outside. /1000 is for things being in micromolar
	// The ionic strength is too high for basic DebHuc. We'll use Davies
	dielConstant = 74; // water at 37�.
	constA = 1.82*1e6*pow(dielConstant*T,-1.5);

	gamma_cai = exp(-constA * 4 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
	gamma_cao = exp(-constA * 4 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));
	gamma_nai = exp(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
	gamma_nao = exp(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));
	gamma_ki = exp(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
	gamma_kao = exp(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));

	const double gammaCaoMyo = gamma_cao;
	const double gammaCaiMyo = gamma_cai;

	const double PhiCaL_i =  4.0*vffrt*(gamma_cai*cai*exp(2.0*vfrt)-gamma_cao*cao)/(exp(2.0*vfrt)-1.0);
	const double PhiCaNa_i =  1.0*vffrt*(gamma_nai*nai*exp(1.0*vfrt)-gamma_nao*nao)/(exp(1.0*vfrt)-1.0);
	const double PhiCaK_i =  1.0*vffrt*(gamma_ki*ki*exp(1.0*vfrt)-gamma_kao*ko)/(exp(1.0*vfrt)-1.0);
	// The rest
	double PCa=8.3757e-05 * ICaL_PCaMultiplier;

	if(celltype==100){
	    PCa=PCa*1.2;
	}
	else if(celltype==101){
	    PCa=PCa*2;
	}

	const double PCap=1.1*PCa;
	const double PCaNa=0.00125*PCa;
	const double PCaK=3.574e-4*PCa;
	const double PCaNap=0.00125*PCap;
	const double PCaKp=3.574e-4*PCap;

	double ICaL_ss=(1.0-fICaLp)*PCa*PhiCaL_ss*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCap*PhiCaL_ss*d*(fp*(1.0-nca)+jca*fcap*nca);
	double ICaNa_ss=(1.0-fICaLp)*PCaNa*PhiCaNa_ss*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaNap*PhiCaNa_ss*d*(fp*(1.0-nca)+jca*fcap*nca);
	double ICaK_ss=(1.0-fICaLp)*PCaK*PhiCaK_ss*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaKp*PhiCaK_ss*d*(fp*(1.0-nca)+jca*fcap*nca);

	double ICaL_i=(1.0-fICaLp)*PCa*PhiCaL_i*d*(f*(1.0-nca_i)+jca*fca*nca_i)+fICaLp*PCap*PhiCaL_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i);
	double ICaNa_i=(1.0-fICaLp)*PCaNa*PhiCaNa_i*d*(f*(1.0-nca_i)+jca*fca*nca_i)+fICaLp*PCaNap*PhiCaNa_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i);
	double ICaK_i=(1.0-fICaLp)*PCaK*PhiCaK_i*d*(f*(1.0-nca_i)+jca*fca*nca_i)+fICaLp*PCaKp*PhiCaK_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i);


	// And we weight ICaL (in ss) and ICaL_i
	ICaL_i = ICaL_i * (1-ICaL_fractionSS);
	ICaNa_i = ICaNa_i * (1-ICaL_fractionSS);
	ICaK_i = ICaK_i * (1-ICaL_fractionSS);
	ICaL_ss = ICaL_ss * ICaL_fractionSS;
	ICaNa_ss = ICaNa_ss * ICaL_fractionSS;
	ICaK_ss = ICaK_ss * ICaL_fractionSS;

	double ICaL = ICaL_ss + ICaL_i;
	double ICaNa = ICaNa_ss + ICaNa_i;
	double ICaK = ICaK_ss + ICaK_i;
	double ICaL_tot = ICaL + ICaNa + ICaK;
	
	// IKr
	// [IKr, dt_ikr_c0, dt_ikr_c1, dt_ikr_c2, dt_ikr_o, dt_ikr_i ] = getIKr_ORd2011_MM(v,ikr_c0,ikr_c1, ikr_c2, ikr_o, ikr_i,...
	//     ko, EK, celltype, IKr_Multiplier);


	// Extracting state vector
	// c3 = y(1);
	// c2 = y(2);
	// c1 = y(3);
	// o = y(4);
	// i = y(5);
	const double b = 0; // no channels blocked in via the mechanism of specific MM states
	// const double vfrt = v*F/(R*T);

	// transition rates
	// from c0 to c1 in l-v model,
	const double alpha = 0.1161 * exp(0.2990 * vfrt);
	// from c1 to c0 in l-v/
	const double beta =  0.2442 * exp(-1.604 * vfrt);

	// from c1 to c2 in l-v/
	const double alpha1 = 1.25 * 0.1235 ;
	// from c2 to c1 in l-v/
	const double beta1 =  0.1911;

	// from c2 to o/           c1 to o
	const double alpha2 =0.0578 * exp(0.9710 * vfrt); //
	// from o to c2/
	const double beta2 = 0.349e-3* exp(-1.062 * vfrt); //

	// from o to i
	const double alphai = 0.2533 * exp(0.5953 * vfrt); //
	// from i to o
	const double betai = 1.25* 0.0522 * exp(-0.8209 * vfrt); //

	// from c2 to i (from c1 in orig)
	const double alphac2ToI = 0.52e-4 * exp(1.525 * vfrt); //
	// from i to c2
	// betaItoC2 = 0.85e-8 * exp(-1.842 * vfrt); //
	const double betaItoC2 = (beta2 * betai * alphac2ToI)/(alpha2 * alphai); //
	// transitions themselves
	// for reason of backward compatibility of naming of an older version of a
	// MM IKr, c3 in code is c0 in article diagram, c2 is c1, c1 is c2.

	const double d_ikr_c0 = ikr_c1 * beta - ikr_c0 * alpha; // delta for c0
	const double d_ikr_c1 = ikr_c0 * alpha + ikr_c2*beta1 - ikr_c1*(beta+alpha1); // c1
	const double d_ikr_c2 = ikr_c1 * alpha1 + ikr_o*beta2 + ikr_i*betaItoC2 - ikr_c2 * (beta1 + alpha2 + alphac2ToI); // subtraction is into c2, to o, to i. // c2
	const double d_ikr_o = ikr_c2 * alpha2 + ikr_i*betai - ikr_o*(beta2+alphai);
	const double d_ikr_i = ikr_c2*alphac2ToI + ikr_o*alphai - ikr_i*(betaItoC2 + betai);

	double GKr = 0.0321 * sqrt(ko/5) * IKr_Multiplier; // 1st element compensates for change to ko (sqrt(5/5.4)* 0.0362)
	if(celltype==100){
	    GKr=GKr*1.3;
	}
	else if(celltype==101){
	    GKr=GKr*0.8;
	}

	const double IKr = GKr * ikr_o  * (v-EK);



	// IKs
	// [IKs,dxs1, dxs2] = getIKs_ORd2011(v,xs1, xs2, cai,  EKs,  celltype, IKs_Multiplier);


	//calculate IKs
	// const double xs1ss=1.0/(1.0+exp((-(v+11.60))/8.932));
	// const double txs1=817.3+1.0/(2.326e-4*exp((v+48.28)/17.80)+0.001292*exp((-(v+210.0))/230.0));
	// const double dxs1=(xs1ss-xs1)/txs1;
	// const double xs2ss=xs1ss;
	// const double txs2=1.0/(0.01*exp((v-50.0)/20.0)+0.0193*exp((-(v+66.54))/31.0));
	// const double dxs2=(xs2ss-xs2)/txs2;
	// const double KsCa=1.0+0.6/(1.0+pow(3.8e-5/cai,1.4));
	// double GKs= 0.0011 * IKs_Multiplier;
	// if(celltype==100){
	//     GKs=GKs*1.4;
	// }

	// double IKs=GKs*KsCa*xs1*xs2*(v-EKs);



	//Toms markov chain formulation of IKs

	const double alpha_Iks=3.72e-003*exp(v*FNORT*2.10e-001);
	const double beta_Iks=2.35e-004*exp(v*FNORT*-2.42e-001);
	const double gamma_Iks=7.25e-003*exp(v*FNORT*2.43e+000);
	const double delta_Iks=1.53e-003*exp(v*FNORT*-6.26e-001);
	const double theta_Iks=1.96e-003;
	const double eta_Iks=1.67e-002*exp(v*FNORT*-1.34e+000);
	const double psi_Iks=0.00852459080491*exp(v*FNORT*0.0124);
	const double omega_Iks=0.001631335712234*exp(v*FNORT*-0.457531296450783);

	const double diks_o2 = psi_Iks*iks_o1-omega_Iks*iks_o2;
	const double diks_o1 = theta_Iks*iks_c15+omega_Iks*iks_o2-(psi_Iks+eta_Iks)*iks_o1;
	const double diks_c15 = gamma_Iks*iks_c14+eta_Iks*iks_o1-(4*delta_Iks+theta_Iks)*iks_c15;
	const double diks_c14 = alpha_Iks*iks_c13+4*delta_Iks*iks_c15+2*gamma_Iks*iks_c12-(beta_Iks+3*delta_Iks+gamma_Iks)*iks_c14;
	const double diks_c13 = beta_Iks*iks_c14+gamma_Iks*iks_c11-(alpha_Iks+3*delta_Iks)*iks_c13;
	const double diks_c12 = alpha_Iks*iks_c11+3*delta_Iks*iks_c14+3*gamma_Iks*iks_c9-(2*beta_Iks+2*delta_Iks+2*gamma_Iks)*iks_c12;
	const double diks_c11 = 2*alpha_Iks*iks_c10+2*beta_Iks*iks_c12+2*gamma_Iks*iks_c8+3*delta_Iks*iks_c13-(beta_Iks+alpha_Iks+gamma_Iks+2*delta_Iks)*iks_c11;
	const double diks_c10 = beta_Iks*iks_c11+gamma_Iks*iks_c7-(2*alpha_Iks+2*delta_Iks)*iks_c10;
	const double diks_c9 = alpha_Iks*iks_c8+2*delta_Iks*iks_c12+4*gamma_Iks*iks_c5-(3*beta_Iks+delta_Iks+3*gamma_Iks)*iks_c9;
	const double diks_c8 = 2*alpha_Iks*iks_c7+3*beta_Iks*iks_c9+3*gamma_Iks*iks_c4+2*delta_Iks*iks_c11-(2*beta_Iks+alpha_Iks+2*gamma_Iks+delta_Iks)*iks_c8;
	const double diks_c7 = 3*alpha_Iks*iks_c6+2*beta_Iks*iks_c8+2*gamma_Iks*iks_c3+2*delta_Iks*iks_c10-(beta_Iks+2*alpha_Iks+delta_Iks+gamma_Iks)*iks_c7;
	const double diks_c6 = beta_Iks*iks_c7+gamma_Iks*iks_c2-(3*alpha_Iks+delta_Iks)*iks_c6;
	const double diks_c5 = alpha_Iks*iks_c4+delta_Iks*iks_c9-(4*beta_Iks+4*gamma_Iks)*iks_c5;
	const double diks_c4 = 2*alpha_Iks*iks_c3+4*beta_Iks*iks_c5+delta_Iks*iks_c8-(3*beta_Iks+alpha_Iks+3*gamma_Iks)*iks_c4;
	const double diks_c3 = 3*alpha_Iks*iks_c2+3*beta_Iks*iks_c4+delta_Iks*iks_c7-(2*beta_Iks+2*alpha_Iks+2*gamma_Iks)*iks_c3;
	const double diks_c2 = 4*alpha_Iks*iks_c1+2*beta_Iks*iks_c3+delta_Iks*iks_c6-(beta_Iks+3*alpha_Iks+gamma_Iks)*iks_c2;
	const double diks_c1 = beta_Iks*iks_c2-4*alpha_Iks*iks_c1;

	//Incremement the channel Markov states
	// iks_o2 += diks_o2*dt;
	// iks_o1 += diks_o1*dt;
	// iks_c1 += diks_c1*dt;
	// iks_c2 += diks_c2*dt;
	// iks_c3 += diks_c3*dt;
	// iks_c4 += diks_c4*dt;
	// iks_c5 += diks_c5*dt;
	// iks_c6 += diks_c6*dt;
	// iks_c7 += diks_c7*dt;
	// iks_c8 += diks_c8*dt;
	// iks_c9 += diks_c9*dt;
	// iks_c10 += diks_c10*dt;
	// iks_c11 += diks_c11*dt;
	// iks_c12 += diks_c12*dt;
	// iks_c13 += diks_c13*dt;
	// iks_c14 += diks_c14*dt;
	// iks_c15 += diks_c15*dt;

	// double GKs = 0.0034;
	double GKs=0.0011 * IKs_Multiplier;
	if (celltype==100){
		GKs*=1.4;
	}
	GKs*=(1.0+0.6/(1.0+pow(3.8e-5/cai,1.4)));

	const double IKs = GKs*(iks_o1 + iks_o2)*(v-EKs);

	//Record the IKs
	Cellstate.set_new_general_cell_model_data(IKs);








	// IK1
	// IK1 = getIK1_CRLP(v, ko, EK, celltype, IK1_Multiplier);

	// IK1
	const double aK1 = 4.094/(1+exp(0.1217*(v-EK-49.934)));
	const double bK1 = (15.72*exp(0.0674*(v-EK-3.257))+exp(0.0618*(v-EK-594.31)))/(1+exp(-0.1629*(v-EK+14.207)));
	const double K1ss = aK1/(aK1+bK1);

	double GK1=IK1_Multiplier  * 0.6992; //0.7266; //* sqrt(5/5.4))
	if(celltype==100){
	    GK1=GK1*1.2;
	}
	else if(celltype==101){
	    GK1=GK1*1.3;
	}

	const double IK1=GK1*sqrt(ko/5)*K1ss*(v-EK);

	// INaCa
	// [ INaCa_i, INaCa_ss] = getINaCa_ORd2011(v,F,R,T, nass, nai, nao, cass, cai, cao, celltype, INaCa_Multiplier, INaCa_fractionSS);

	double zca = 2.0;
	double kna1=15.0;
	double kna2=5.0;
	double kna3=88.12;
	double kasymm=12.5;
	double wna=6.0e4;
	double wca=6.0e4;
	double wnaca=5.0e3;
	double kcaon=1.5e6;
	double kcaoff=5.0e3;
	double qna=0.5224;
	double qca=0.1670;
	double hca=exp((qca*v*F)/(R*T));
	double hna=exp((qna*v*F)/(R*T));
	double h1=1+nai/kna3*(1+hna);
	double h2=(nai*hna)/(kna3*h1);
	double h3=1.0/h1;
	double h4=1.0+nai/kna1*(1+nai/kna2);
	double h5=nai*nai/(h4*kna1*kna2);
	double h6=1.0/h4;
	double h7=1.0+nao/kna3*(1.0+1.0/hna);
	double h8=nao/(kna3*hna*h7);
	double h9=1.0/h7;
	double h10=kasymm+1.0+nao/kna1*(1.0+nao/kna2);
	double h11=nao*nao/(h10*kna1*kna2);
	double h12=1.0/h10;
	double k1=h12*cao*kcaon;
	double k2=kcaoff;
	double k3p=h9*wca;
	double k3pp=h8*wnaca;
	double k3=k3p+k3pp;
	double k4p=h3*wca/hca;
	double k4pp=h2*wnaca;
	double k4=k4p+k4pp;
	double k5=kcaoff;
	double k6=h6*cai*kcaon;
	double k7=h5*h2*wna;
	double k8=h8*h11*wna;
	double x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
	double x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
	double x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
	double x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
	double E1=x1/(x1+x2+x3+x4);
	double E2=x2/(x1+x2+x3+x4);
	double E3=x3/(x1+x2+x3+x4);
	double E4=x4/(x1+x2+x3+x4);
	double KmCaAct=150.0e-6;
	double allo=1.0/(1.0+pow(KmCaAct/cai,2.0));
	double zna=1.0;
	double JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
	double JncxCa=E2*k2-E1*k1;
	double Gncx= 0.0034* INaCa_Multiplier;
	if(celltype==100){
	    Gncx=Gncx*1.1;
	}
	else if(celltype==101){
	    Gncx=Gncx*1.4;
	}

	const double INaCa_i=(1-INaCa_fractionSS)*Gncx*allo*(zna*JncxNa+zca*JncxCa);

	//calculate INaCa_ss
	h1=1+nass/kna3*(1+hna);
	h2=(nass*hna)/(kna3*h1);
	h3=1.0/h1;
	h4=1.0+nass/kna1*(1+nass/kna2);
	h5=nass*nass/(h4*kna1*kna2);
	h6=1.0/h4;
	h7=1.0+nao/kna3*(1.0+1.0/hna);
	h8=nao/(kna3*hna*h7);
	h9=1.0/h7;
	h10=kasymm+1.0+nao/kna1*(1+nao/kna2);
	h11=nao*nao/(h10*kna1*kna2);
	h12=1.0/h10;
	k1=h12*cao*kcaon;
	k2=kcaoff;
	k3p=h9*wca;
	k3pp=h8*wnaca;
	k3=k3p+k3pp;
	k4p=h3*wca/hca;
	k4pp=h2*wnaca;
	k4=k4p+k4pp;
	k5=kcaoff;
	k6=h6*cass*kcaon;
	k7=h5*h2*wna;
	k8=h8*h11*wna;
	x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
	x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
	x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
	x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
	E1=x1/(x1+x2+x3+x4);
	E2=x2/(x1+x2+x3+x4);
	E3=x3/(x1+x2+x3+x4);
	E4=x4/(x1+x2+x3+x4);
	KmCaAct=150.0e-6 ;
	allo=1.0/(1.0+pow(KmCaAct/cass,2.0));
	JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
	JncxCa=E2*k2-E1*k1;
	double INaCa_ss=INaCa_fractionSS*Gncx*allo*(zna*JncxNa+zca*JncxCa);

	// INaK
	// INaK = getINaK_ORd2011(v, F, R, T, nai, nao, ki, ko, celltype, INaK_Multiplier);

	//calculate INaK
	zna=1.0;
	const double k1p=949.5;
	const double k1m=182.4;
	const double k2p=687.2;
	const double k2m=39.4;
	k3p=1899.0;
	const double k3m=79300.0;
	k4p=639.0;
	const double k4m=40.0;
	const double Knai0=9.073;
	const double Knao0=27.78;
	const double delta=-0.1550;
	const double Knai=Knai0*exp((delta*v*F)/(3.0*R*T));
	const double Knao=Knao0*exp(((1.0-delta)*v*F)/(3.0*R*T));
	const double Kki=0.5;
	const double Kko=0.3582;
	const double MgADP=0.05;
	const double MgATP=9.8;
	const double Kmgatp=1.698e-7;
	const double H=1.0e-7;
	const double eP=4.2;
	const double Khp=1.698e-7;
	const double Knap=224.0;
	const double Kxkur=292.0;
	const double P=eP/(1.0+H/Khp+nai/Knap+ki/Kxkur);
	const double a1=(k1p*pow(nai/Knai,3.0))/(pow(1.0+nai/Knai,3.0)+pow(1.0+ki/Kki,2.0)-1.0);
	const double b1=k1m*MgADP;
	const double a2=k2p;
	const double b2=(k2m*pow(nao/Knao,3.0))/(pow(1.0+nao/Knao,3.0)+pow(1.0+ko/Kko,2.0)-1.0);
	const double a3=(k3p*pow(ko/Kko,2.0))/(pow(1.0+nao/Knao,3.0)+pow(1.0+ko/Kko,2.0)-1.0);
	const double b3=(k3m*P*H)/(1.0+MgATP/Kmgatp);
	const double a4=(k4p*MgATP/Kmgatp)/(1.0+MgATP/Kmgatp);
	const double b4=(k4m*pow(ki/Kki,2.0))/(pow(1.0+nai/Knai,3.0)+pow(1.0+ki/Kki,2.0)-1.0);
	x1=a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2;
	x2=b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4;
	x3=a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1;
	x4=b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1;
	E1=x1/(x1+x2+x3+x4);
	E2=x2/(x1+x2+x3+x4);
	E3=x3/(x1+x2+x3+x4);
	E4=x4/(x1+x2+x3+x4);
	const double zk=1.0;
	const double JnakNa=3.0*(E1*a3-E2*b3);
	const double JnakK=2.0*(E4*b1-E3*a1);
	double Pnak= 15.4509 * INaK_Multiplier;
	if(celltype==100){
	    Pnak=Pnak*0.9;
	}
	else if(celltype==101){
	    Pnak=Pnak*0.7;
	}

	const double INaK=Pnak*(zna*JnakNa+zk*JnakK);


	// Minor/background currents
	// calculate IKb
	const double xkb=1.0/(1.0+exp(-(v-10.8968)/(23.9871)));
	double GKb=0.0189*IKb_Multiplier;
	if(celltype==100){
	    GKb=GKb*0.6;
	}
	const double IKb=GKb*xkb*(v-EK);

	// calculate INab
	const double PNab=1.9239e-09*INab_Multiplier;
	const double INab=PNab*vffrt*(nai*exp(vfrt)-nao)/(exp(vfrt)-1.0);

	// calculate ICab
	const double PCab=5.9194e-08*ICab_Multiplier;
	
	const double ICab=PCab*4.0*vffrt*(gammaCaiMyo*cai*exp(2.0*vfrt)-gammaCaoMyo*cao)/(exp(2.0*vfrt)-1.0);

	// calculate IpCa
	const double GpCa=5e-04*IpCa_Multiplier;
	const double IpCa=GpCa*cai/(0.0005+cai);

	//  Chloride
	// I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current

	const double ecl = (R*T/F)*log(cli/clo);            // [mV]
	const double eclss = (R*T/F)*log(clss/clo);            // [mV]

	const double Fjunc = 1;
	const double Fsl = 1-Fjunc; // fraction in SS and in myoplasm - as per literature, I(Ca)Cl is in junctional subspace

	// Fsl = 1-Fjunc; // fraction in SS and in myoplasm
	const double GClCa = IClCa_Multiplier * 0.2843;   // [mS/uF]
	const double GClB = IClb_Multiplier * 1.98e-3;        // [mS/uF]
	const double KdClCa = 0.1;    // [mM]

	const double I_ClCa_junc = Fjunc*GClCa/(1+KdClCa/cass)*(v-eclss);
	const double I_ClCa_sl = Fsl*GClCa/(1+KdClCa/cai)*(v-ecl);

	const double I_ClCa = I_ClCa_junc+I_ClCa_sl;
	const double I_Clbk = GClB*(v-ecl);

	// Calcium handling
	// calculate ryanodione receptor calcium induced calcium release from the jsr
	const double fJrelp=(1.0/(1.0+KmCaMK/CaMKa));

	// Jrel
	// [Jrel, dJrel_np, dJrel_p] = getJrel_ORd2011(Jrel_np, Jrel_p, ICaL_ss,cass, cajsr, fJrelp, celltype, Jrel_Multiplier);

	const double jsrMidpoint = 1.7;

	const double bt=4.75;
	const double a_rel=0.5*bt;
	double Jrel_inf=a_rel*(-ICaL)/(1.0+pow(jsrMidpoint/cajsr,8.0));
	if(celltype==101){
	    Jrel_inf=Jrel_inf*1.7;
	}
	double tau_rel=bt/(1.0+0.0123/cajsr);

	if( tau_rel<0.001){
	    tau_rel=0.001;
	}

	const double dJrel_np=(Jrel_inf-Jrel_np)/tau_rel;
	const double btp=1.25*bt;
	const double a_relp=0.5*btp;
	double Jrel_infp=a_relp*(-ICaL)/(1.0+pow(jsrMidpoint/cajsr,8.0));
	if(celltype==101){
	    Jrel_infp=Jrel_infp*1.7;
	}
	double tau_relp=btp/(1.0+0.0123/cajsr);

	if( tau_relp<0.001){
	    tau_relp=0.001;
	}

	const double dJrel_p=(Jrel_infp-Jrel_p)/tau_relp;

	const double Jrel=Jrel_Multiplier * 1.5378 * ((1.0-fJrelp)*Jrel_np+fJrelp*Jrel_p);




	const double fJupp=(1.0/(1.0+KmCaMK/CaMKa));
	// [Jup, Jleak] = getJup_ORd2011(cai, cansr, fJupp, celltype, Jup_Multiplier);

	//calculate serca pump, ca uptake flux
	double Jupnp=Jup_Multiplier * 0.005425*cai/(cai+0.00092);
	double Jupp=Jup_Multiplier * 2.75*0.005425*cai/(cai+0.00092-0.00017);
	if(celltype==100){
	    Jupnp=Jupnp*1.3;
	    Jupp=Jupp*1.3;
	}

	const double Jleak=Jup_Multiplier* 0.0048825*cansr/15.0;
	const double Jup=(1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak;



	// calculate tranlocation flux
	const double Jtr=(cansr-cajsr)/60;

	

	// calculate the stimulus current, Istim
	// amp=stimAmp;
	// duration=stimDur;
	// if t<=duration
	//     Istim=amp;
	// else
	//     Istim=0.0;
	// end

	

	// update the membrane voltage

	// dv=-(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab+I_ClCa+I_Clbk/*+Istim*/);



	// calculate diffusion fluxes
	const double JdiffNa=(nass-nai)/2.0;
	const double JdiffK=(kss-ki)/2.0;
	const double JdiffCl=(clss-cli)/2.0;
	const double Jdiff=(cass-cai)/0.2;


	// calcium buffer constants 
	double cmdnmax= 0.05; 
	if(celltype==100){
	    cmdnmax=cmdnmax*1.3;
	}

	const double kmcmdn=0.00238; 
	const double trpnmax=0.07;
	const double kmtrpn=0.0005;
	const double BSRmax=0.047;
	const double KmBSR = 0.00087;
	const double BSLmax=1.124;
	const double KmBSL = 0.0087;
	const double csqnmax=10.0;
	const double kmcsqn=0.8;

	// update intracellular concentrations, using buffers for cai, cass, cajsr
	const double dnai=-(ICaNa_i+INa+INaL+3.0*INaCa_i+3.0*INaK+INab)*Acap/(F*vmyo)+JdiffNa*vss/vmyo;
	const double dnass=-(ICaNa_ss+3.0*INaCa_ss)*Acap/(F*vss)-JdiffNa;

	const double dki=-(ICaK_i+Ito+IKr+IKs+IK1+IKb/*+Istim*/-2.0*INaK)*Acap/(F*vmyo)+JdiffK*vss/vmyo;
	const double dkss=-(ICaK_ss)*Acap/(F*vss)-JdiffK;

	const double Bcai=1.0/(1.0+cmdnmax*kmcmdn/pow(kmcmdn+cai,2.0)+trpnmax*kmtrpn/pow(kmtrpn+cai,2.0));
	const double dcai=Bcai*(-(ICaL_i + IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo);


	const double Bcass=1.0/(1.0+BSRmax*KmBSR/pow(KmBSR+cass,2.0)+BSLmax*KmBSL/pow(KmBSL+cass,2.0));
	const double dcass=Bcass*(-(ICaL_ss-2.0*INaCa_ss)*Acap/(2.0*F*vss)+Jrel*vjsr/vss-Jdiff);

	const double dcli = - (I_Clbk + I_ClCa_sl)*Acap/(-1*F*vmyo)+JdiffCl*vss/vmyo;
	const double dclss = - I_ClCa_junc*Acap/(-1*F*vss)-JdiffCl;



	const double dcansr=Jup-Jtr*vjsr/vnsr;

	const double Bcajsr=1.0/(1.0+csqnmax*kmcsqn/pow(kmcsqn+cajsr,2.0));
	const double dcajsr=Bcajsr*(Jtr-Jrel);

	//Update the variables
	new_state[0] += dt * dnai;
	new_state[1] += dt * dnass;
	new_state[2] += dt * dki;
	new_state[3] += dt * dkss;
	new_state[4] += dt * dcai;
	new_state[5] += dt * dcass;
	new_state[6] += dt * dcansr;
	new_state[7] += dt * dcajsr;
	new_state[8] += dt * dm;
	new_state[9] += dt * dhp;
	new_state[10] += dt * dh;
	new_state[11] += dt * dj;

	new_state[12] += dt * djp;
	new_state[13] += dt * dmL;
	new_state[14] += dt * dhL;
	new_state[15] += dt * dhLp;
	new_state[16] += dt * da;
	new_state[17] += dt * diF;
	new_state[18] += dt * diS;
	new_state[19] += dt * dap;
	new_state[20] += dt * diFp;
	new_state[21] += dt * diSp;
	// ical
	new_state[22] += dt * dd;
	new_state[23] += dt * dff;
	new_state[24] += dt * dfs;
	new_state[25] += dt * dfcaf;
	new_state[26] += dt * dfcas;
	new_state[27] += dt * djca;
	new_state[28] += dt * dnca;
	new_state[29] += dt * dnca_i;
	new_state[30] += dt * dffp;
	new_state[31] += dt * dfcafp;

	//Old IKs variables
	// new_state[32] += dt * dxs1;
	// new_state[33] += dt * dxs2;

	//New IKs variables for markov chain
	new_state[32] += dt * diks_o2;
	new_state[33] += dt * diks_o1;
	new_state[34] += dt * diks_c1;
	new_state[35] += dt * diks_c2;
	new_state[36] += dt * diks_c3;
	new_state[37] += dt * diks_c4;
	new_state[38] += dt * diks_c5;
	new_state[39] += dt * diks_c6;
	new_state[40] += dt * diks_c7;
	new_state[41] += dt * diks_c8;
	new_state[42] += dt * diks_c9;
	new_state[43] += dt * diks_c10;
	new_state[44] += dt * diks_c11;
	new_state[45] += dt * diks_c12;
	new_state[46] += dt * diks_c13;
	new_state[47] += dt * diks_c14;
	new_state[48] += dt * diks_c15;


	new_state[49] += dt * dJrel_np;
	new_state[50] += dt * dCaMKt;
	// MM ICaL states
	new_state[51] += dt * d_ikr_c0;
	new_state[52] += dt * d_ikr_c1;
	new_state[53] += dt * d_ikr_c2;
	new_state[54] += dt * d_ikr_o;
	new_state[55] += dt * d_ikr_i;
	new_state[56] += dt * dJrel_p;
	new_state[57] += dt * dcli;
	new_state[58] += dt * dclss;

	//The active strain and ion membrane current
	Cellstate.set_active_strain(0.0);
	Cellstate.set_membrane_current((INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab+I_ClCa+I_Clbk/*+Istim*/));
}




}