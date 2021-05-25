#include "TomekORudy_Updated.h"

namespace oomph{

TomekORudyVentUpdated::TomekORudyVentUpdated()
{
	nao =140.0;
	ko = 5.0;
	cao = 1.8;


	R=8314.0;
	T=310.0;
	F=96485.0;

	L=0.01;
	rad=0.0011;
	vcell=1000*3.14*rad*rad*L;
	Ageo=2*3.14*rad*rad+2*3.14*rad*L;
	Acap=2*Ageo;
	vmyo=0.68*vcell;
	vnsr=0.0552*vcell;
	vjsr=0.0048*vcell;
	vss=0.02*vcell;

	Names_Of_Cell_Variables=
	{
		"nai",
		"nass",
		"ki",
		"kss",
		"cai",
		"cass",
		"cansr",
		"cajsr",
		"m",
		"hp",
		"h",
		"j",
		"jp",
		"mL",
		"hL",
		"hLp",
		"a",
		"iF",
		"iS",
		"ap",
		"iFp",
		"iSp",
		"d",
		"ff",
		"fs",
		"fcaf",
		"fcas",
		"jca",
		"nca",
		"nca_i",
		"ffp",
		"fcafp",
		"xs1",
		"xs2",
		"Jrel_np",
		"CaMKt",
		"ikr_c0",
		"ikr_c1",
		"ikr_c2",
		"ikr_o",
		"ikr_i",
		"Jrel_p",
		"cli",
		"clss"
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


double TomekORudyVentUpdated::return_initial_state_variable(const unsigned &v, const unsigned &cell_type)
{
	double STATES[3][44] = 
	{
		{13.4006197127101,	13.4009366375605,	152.363889991411,	152.363844281265,	6.62181569329109e-05,	5.74992101186664e-05,	1.80679354108253,	1.80504703846161,	0.000525323142244124,	0.731365613931879,	0.864514804478946,	0.864457115015362,	0.864352715733040,	0.000111796946200753,	0.591653600903076,	0.347681236790268,	0.000832040787904695,	0.999724237643709,	0.999723490323258,	0.000423912054921203,	0.999724237767365,	0.999724054768924,	-2.48652723154831e-36,	0.999999995656968,	0.951060213798798,	0.999999995656976,	0.999937692440505,	0.999988620215189,	0.000304952292598877,	0.000527266836293694,	0.999999995656685,	0.999999995656399,	0.223358448594285,	0.000141824743798310,	6.77882673136480e-25,	0.0127354136594547,	0.998473340797544,	0.000739304520671507,	0.000602907876760312,	0.000178778329437106,	5.67825481839631e-06,	-1.58194108230975e-23,	34.3172099045416,	34.3171879616709},
		
		{15.9486733631902,	15.9492201219199,	156.713067183465,	156.713012917889,	8.29757614119249e-05,	6.64218738239930e-05,	2.01222465016069,	2.01641461550468,	0.000461956534829239,	0.747897241513971,	0.873907657623256,	0.873784149746068,	0.873537507470071,	9.98770874860868e-05,	0.598611759466370,	0.333989914508018,	0.000799404188333773,	0.999751418648699,	0.570253781914665,	0.000407277681198651,	0.999751424044492,	0.635192659820007,	-8.33460419955820e-30,	0.999999996299523,	0.918358695116332,	0.999999996300839,	0.999754043234341,	0.999974289100453,	0.000533651965198686,	0.00125786136830066,	0.999999996207845,	0.999999996298728,	0.264229298875232,	0.000132734770983313,	-1.30048587380675e-21,	0.0201882046405105,	0.998345129712001,	0.000708645966770289,	0.000591004733696651,	0.000344573276330018,	1.06423014130810e-05,	-7.61071425848620e-20,	48.9127681450109,	48.9127355777872},
		
		{12.3973564703237,	12.3976953330256,	147.711474407329,	147.711433360862,	7.45348093580908e-05,	6.49734059262396e-05,	1.52800059077496,	1.52569253436820,	0.000651715361520436,	0.701845427662345,	0.847326693399099,	0.847165668535608,	0.846901438940383,	0.000135120261194560,	0.556601650459093,	0.311549070692613,	0.000889925901945137,	0.999671576537361,	0.598890849399207,	0.000453416519958675,	0.999671583620741,	0.662069197871101,	1.58884121071311e-31,	0.999999994310640,	0.940179072110115,	0.999999994310551,	0.999901376748386,	0.999984619576572,	0.000489937844976504,	0.000832600925714172,	0.999999994343118,	0.999999994312043,	0.243959016060676,	0.000158616694622612,	1.80824774078763e-22,	0.0109502564192484,	0.998251134783098,	0.000793602046257213,	0.000653214338491609,	0.000292244895058368,	9.80408293117642e-06,	4.35860784390795e-21,	29.2069793887423,	29.2069557375901}
	};

	if(v<44){
		return STATES[cell_type-100][v];
	}
	else{
		std::cout << "Failed at variable index check " << v << " is not within the expected range" << std::endl;
	}
}

double TomekORudyVentUpdated::return_initial_membrane_potential(const unsigned &cell_type)
{
	switch(cell_type){
		case 100:
			return -90.7456332496277;
		case 101:
			return -91.3391803967798;
		case 102:
			return -89.7480828168753;
	}
}

void TomekORudyVentUpdated::Calculate_Derivatives(const double &Vm,
											const Vector<double> &CellVariables,
											const double &t,
											const unsigned &cell_type,
											const double &Istim,
											const Vector<double> &Other_Parameters,
											const Vector<double> &Other_Variables,
											
											Vector<double> &Variable_Derivatives,
											double &Iion)
{
	//Unpack the membrane potential, cell type, and time step

	//give names to the state vector values
	// v=X(1);
	const double nai = CellVariables[nai_Tomek];
	const double nass = CellVariables[nass_Tomek];
	const double ki = CellVariables[ki_Tomek];
	const double kss = CellVariables[kss_Tomek];
	const double cai = CellVariables[cai_Tomek];
	const double cass = CellVariables[cass_Tomek];
	const double cansr = CellVariables[cansr_Tomek];
	const double cajsr = CellVariables[cajsr_Tomek];
	const double m = CellVariables[m_Tomek];
	const double hp = CellVariables[hp_Tomek];
	const double h = CellVariables[h_Tomek];
	const double j = CellVariables[j_Tomek];
	const double jp = CellVariables[jp_Tomek];
	const double mL = CellVariables[mL_Tomek];
	const double hL = CellVariables[hL_Tomek];
	const double hLp = CellVariables[hLp_Tomek];
	const double a = CellVariables[a_Tomek];
	const double iF = CellVariables[iF_Tomek];
	const double iS = CellVariables[iS_Tomek];
	const double ap = CellVariables[ap_Tomek];
	const double iFp = CellVariables[iFp_Tomek];
	const double iSp = CellVariables[iSp_Tomek];
	const double d = CellVariables[d_Tomek];
	const double ff = CellVariables[ff_Tomek];
	const double fs = CellVariables[fs_Tomek];
	const double fcaf = CellVariables[fcaf_Tomek];
	const double fcas = CellVariables[fcas_Tomek];
	const double jca = CellVariables[jca_Tomek];
	const double nca = CellVariables[nca_Tomek];
	const double nca_i = CellVariables[nca_i_Tomek];
	const double ffp = CellVariables[ffp_Tomek];
	const double fcafp = CellVariables[fcafp_Tomek];
	const double xs1 = CellVariables[xs1_Tomek];
	const double xs2 = CellVariables[xs2_Tomek];
	const double Jrel_np = CellVariables[Jrel_np_Tomek];
	const double CaMKt = CellVariables[CaMKt_Tomek];
	const double ikr_c0 = CellVariables[ikr_c0_Tomek];
	const double ikr_c1 = CellVariables[ikr_c1_Tomek];
	const double ikr_c2 = CellVariables[ikr_c2_Tomek];
	const double ikr_o = CellVariables[ikr_o_Tomek];
	const double ikr_i = CellVariables[ikr_i_Tomek];
	const double Jrel_p = CellVariables[Jrel_p_Tomek];
	const double cli = CellVariables[cli_Tomek];
	const double clss = CellVariables[clss_Tomek];
	// const double cli = 24; // Intracellular Cl  [mM]
	const double clo = 150;  // Extracellular Cl  [mM]


	//Some parameters which are passed to the timestepper in the original function
	//	These probably represent mutation/disease/pathological effects
	const double INa_Multiplier = 1.0;
	const double INaL_Multiplier = 1.0;
	const double Ito_Multiplier = 1.0;
	const double ICaL_PCaMultiplier = 1.0;
	const double ICaL_fractionSS = 0.8;
	const double IKs_Multiplier = 1.0;
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
	const double vffrt=Vm*F*F/(R*T);
	const double vfrt=Vm*F/(R*T);
	const double frt = F/(R*T);

	
	const double fINap=(1.0/(1.0+KmCaMK/CaMKa));
	const double fINaLp=(1.0/(1.0+KmCaMK/CaMKa));
	const double fItop=(1.0/(1.0+KmCaMK/CaMKa));
	const double fICaLp=(1.0/(1.0+KmCaMK/CaMKa));

	// INa
	// [INa, dm, dh, dhp, dj, djp] = getINa_Grandi(Vm, m, h, hp, j, jp, fINap, ENa, INa_Multiplier);

	// The Grandi implementation updated with INa phosphorylation.
	// m gate
	const double mss = 1 / pow((1 + exp( -(56.86 + Vm) / 9.03 )), 2.0);
	const double taum = 0.1292 * exp(-pow((Vm+45.79)/15.54,2.0)) + 0.06487 * exp(-pow((Vm-4.823)/51.12,2.0));
	const double dm = (mss - m) / taum;

	// h gate
	const double ah = (Vm >= -40) * (0) + (Vm < -40) * (0.057 * exp( -(Vm + 80) / 6.8 ));
	const double bh = (Vm >= -40) * (0.77 / (0.13*(1 + exp( -(Vm + 10.66) / 11.1 )))) + (Vm < -40) * ((2.7 * exp( 0.079 * Vm) + 3.1*1e5 * exp(0.3485 * Vm)));
	const double tauh = 1 / (ah + bh);
	const double hss = 1 / (pow(1 + exp( (Vm + 71.55)/7.43 ),2.0));
	const double dh = (hss - h) / tauh;
	// j gate
	const double aj = (Vm >= -40) * (0) +(Vm < -40) * (((-2.5428 * 1e4*exp(0.2444*Vm) - 6.948*1e-6 * exp(-0.04391*Vm)) * (Vm + 37.78)) / (1 + exp( 0.311 * (Vm + 79.23) )));
	const double bj = (Vm >= -40) * ((0.6 * exp( 0.057 * Vm)) / (1 + exp( -0.1 * (Vm + 32) ))) + (Vm < -40) * ((0.02424 * exp( -0.01052 * Vm )) / (1 + exp( -0.1378 * (Vm + 40.14) )));
	const double tauj = 1 / (aj + bj);
	const double jss = 1 / (pow(1 + exp( (Vm + 71.55)/7.43 ),2.0));
	const double dj = (jss - j) / tauj;

	// h phosphorylated
	const double hssp = 1 / (pow(1 + exp( (Vm + 71.55 + 6)/7.43 ),2.0));
	const double dhp = (hssp - hp) / tauh;
	// j phosphorylated
	const double taujp = 1.46 * tauj;
	const double djp = (jss - jp) / taujp;

	const double GNa = 11.7802;
	const double INa=INa_Multiplier * GNa*(Vm-ENa)*pow(m,3.0)*((1.0-fINap)*h*j+fINap*hp*jp);




	// INaL
	// [INaL,dmL,dhL,dhLp] = getINaL_ORd2011(Vm, mL, hL, hLp, fINaLp, ENa, cell_type, INaL_Multiplier);

	//calculate INaL
	const double mLss=1.0/(1.0+exp((-(Vm+42.85))/5.264));
	const double tm = 0.1292 * exp(-pow((Vm+45.79)/15.54,2.0)) + 0.06487 * exp(-pow((Vm-4.823)/51.12,2.0)); 
	const double tmL=tm;
	const double dmL=(mLss-mL)/tmL;
	const double hLss=1.0/(1.0+exp((Vm+87.61)/7.488));
	const double thL=200.0;
	const double dhL=(hLss-hL)/thL;
	const double hLssp=1.0/(1.0+exp((Vm+93.81)/7.488));
	const double thLp=3.0*thL;
	const double dhLp=(hLssp-hLp)/thLp;
	double GNaL=0.0279 * INaL_Multiplier;
	if(cell_type==100){
	    GNaL=GNaL*0.6;
	}

	const double INaL=GNaL*(Vm-ENa)*mL*((1.0-fINaLp)*hL+fINaLp*hLp);



	// ITo
	// [Ito,da,diF,diS,dap,diFp, diSp] = getITo_ORd2011(Vm, a, iF, iS, ap, iFp, iSp, fItop, EK, cell_type, Ito_Multiplier);

	//calculate Ito
	const double ass=1.0/(1.0+exp((-(Vm-14.34))/14.82));
	const double ta=1.0515/(1.0/(1.2089*(1.0+exp(-(Vm-18.4099)/29.3814)))+3.5/(1.0+exp((Vm+100.0)/29.3814)));
	const double da=(ass-a)/ta;
	const double iss=1.0/(1.0+exp((Vm+43.94)/5.711));
	double delta_epi;
	if(cell_type==100){
	    delta_epi=1.0-(0.95/(1.0+exp((Vm+70.0)/5.0)));
	}
	else{
	    delta_epi=1.0;
	}
	double tiF=4.562+1/(0.3933*exp((-(Vm+100.0))/100.0)+0.08004*exp((Vm+50.0)/16.59));
	double tiS=23.62+1/(0.001416*exp((-(Vm+96.52))/59.05)+1.780e-8*exp((Vm+114.1)/8.079));
	tiF=tiF*delta_epi;
	tiS=tiS*delta_epi;
	const double AiF=1.0/(1.0+exp((Vm-213.6)/151.2));
	const double AiS=1.0-AiF;
	const double diF=(iss-iF)/tiF;
	const double diS=(iss-iS)/tiS;
	const double i=AiF*iF+AiS*iS;
	const double assp=1.0/(1.0+exp((-(Vm-24.34))/14.82));
	const double dap=(assp-ap)/ta;
	const double dti_develop=1.354+1.0e-4/(exp((Vm-167.4)/15.89)+exp(-(Vm-12.23)/0.2154));
	const double dti_recover=1.0-0.5/(1.0+exp((Vm+70.0)/20.0));
	const double tiFp=dti_develop*dti_recover*tiF;
	const double tiSp=dti_develop*dti_recover*tiS;
	const double diFp=(iss-iFp)/tiFp;
	const double diSp=(iss-iSp)/tiSp;
	const double ip=AiF*iFp+AiS*iSp;
	double Gto=0.16 * Ito_Multiplier;
	if(cell_type==100){
	    Gto=Gto*2.0;
	}
	else if(cell_type==101){
	    Gto=Gto*2.0;
	}

	const double Ito=Gto*(Vm-EK)*((1.0-fItop)*a*i+fItop*ap*ip);




	// ICaL
	// [ICaL_ss,ICaNa_ss,ICaK_ss,ICaL_i,ICaNa_i,ICaK_i,dd,dff,dfs,dfcaf,dfcas,djca,dnca,dnca_i,...
	//     dffp,dfcafp, PhiCaL_ss, PhiCaL_i, gammaCaoMyo, gammaCaiMyo] = getICaL_ORd2011_jt(Vm, d,ff,fs,fcaf,fcas,jca,nca,nca_i,ffp,fcafp,...
	//     fICaLp, cai, cass, cao, nai, nass, nao, ki, kss, ko, cli, clo, clss, cell_type, ICaL_fractionSS, ICaL_Multiplier );

	//physical constants
	// const double vffrt=Vm*F*F/(R*T);
	// const double vfrt=Vm*F/(R*T);

	//calculate ICaL, ICaNa, ICaK

	double dss=1.0763*exp(-1.0070*exp(-0.0829*(Vm)));  // magyar
	if(Vm >31.4978){// activation cannot be greater than 1
	    dss = 1;
	}


	const double td= 0.6+1.0/(exp(-0.05*(Vm+6.0))+exp(0.09*(Vm+14.0)));

	const double dd=(dss-d)/td;
	const double fss=1.0/(1.0+exp((Vm+19.58)/3.696));
	const double tff=7.0+1.0/(0.0045*exp(-(Vm+20.0)/10.0)+0.0045*exp((Vm+20.0)/10.0));
	const double tfs=1000.0+1.0/(0.000035*exp(-(Vm+5.0)/4.0)+0.000035*exp((Vm+5.0)/6.0));
	const double Aff=0.6;
	const double Afs=1.0-Aff;
	const double dff=(fss-ff)/tff;
	const double dfs=(fss-fs)/tfs;
	const double f=Aff*ff+Afs*fs;
	const double fcass=fss;
	const double tfcaf=7.0+1.0/(0.04*exp(-(Vm-4.0)/7.0)+0.04*exp((Vm-4.0)/7.0));
	const double tfcas=100.0+1.0/(0.00012*exp(-Vm/3.0)+0.00012*exp(Vm/7.0));

	const double Afcaf=0.3+0.6/(1.0+exp((Vm-10.0)/10.0));

	const double Afcas=1.0-Afcaf;
	const double dfcaf=(fcass-fcaf)/tfcaf;
	const double dfcas=(fcass-fcas)/tfcas;
	const double fca=Afcaf*fcaf+Afcas*fcas;

	const double tjca = 72.5;
	const double jcass = 1.0/(1.0+exp((Vm+18.08)/(2.7916)));   
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

	if(cell_type==100){
	    PCa=PCa*1.2;
	}
	else if(cell_type==101){
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
	// [IKr, dt_ikr_c0, dt_ikr_c1, dt_ikr_c2, dt_ikr_o, dt_ikr_i ] = getIKr_ORd2011_MM(Vm,ikr_c0,ikr_c1, ikr_c2, ikr_o, ikr_i,...
	//     ko, EK, cell_type, IKr_Multiplier);


	// Extracting state vector
	// c3 = y(1);
	// c2 = y(2);
	// c1 = y(3);
	// o = y(4);
	// i = y(5);
	const double b = 0; // no channels blocked in via the mechanism of specific MM states
	// const double vfrt = Vm*F/(R*T);

	// transition rates
	// from c0 to c1 in l-Vm model,
	const double alpha = 0.1161 * exp(0.2990 * vfrt);
	// from c1 to c0 in l-Vm/
	const double beta =  0.2442 * exp(-1.604 * vfrt);

	// from c1 to c2 in l-Vm/
	const double alpha1 = 1.25 * 0.1235 ;
	// from c2 to c1 in l-Vm/
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
	if(cell_type==100){
	    GKr=GKr*1.3;
	}
	else if(cell_type==101){
	    GKr=GKr*0.8;
	}

	// GKr /= 90.0;

	const double IKr = GKr * ikr_o  * (Vm-EK);



	// IKs
	// [IKs,dxs1, dxs2] = getIKs_ORd2011(Vm,xs1, xs2, cai,  EKs,  cell_type, IKs_Multiplier);


	//calculate IKs
	const double xs1ss=1.0/(1.0+exp((-(Vm+11.60))/8.932));
	const double txs1=817.3+1.0/(2.326e-4*exp((Vm+48.28)/17.80)+0.001292*exp((-(Vm+210.0))/230.0));
	const double dxs1=(xs1ss-xs1)/txs1;
	const double xs2ss=xs1ss;
	const double txs2=1.0/(0.01*exp((Vm-50.0)/20.0)+0.0193*exp((-(Vm+66.54))/31.0));
	const double dxs2=(xs2ss-xs2)/txs2;
	const double KsCa=1.0+0.6/(1.0+pow(3.8e-5/cai,1.4));
	double GKs= 0.0011 * IKs_Multiplier;
	if(cell_type==100){
	    GKs=GKs*1.4;
	}
	// GKs *= 90.0;

	const double IKs=GKs*KsCa*xs1*xs2*(Vm-EKs);




	// IK1
	// IK1 = getIK1_CRLP(Vm, ko, EK, cell_type, IK1_Multiplier);

	// IK1
	const double aK1 = 4.094/(1+exp(0.1217*(Vm-EK-49.934)));
	const double bK1 = (15.72*exp(0.0674*(Vm-EK-3.257))+exp(0.0618*(Vm-EK-594.31)))/(1+exp(-0.1629*(Vm-EK+14.207)));
	const double K1ss = aK1/(aK1+bK1);

	double GK1=IK1_Multiplier  * 0.6992; //0.7266; //* sqrt(5/5.4))
	if(cell_type==100){
	    GK1=GK1*1.2;
	}
	else if(cell_type==101){
	    GK1=GK1*1.3;
	}

	const double IK1=GK1*sqrt(ko/5)*K1ss*(Vm-EK);

	// INaCa
	// [ INaCa_i, INaCa_ss] = getINaCa_ORd2011(Vm,F,R,T, nass, nai, nao, cass, cai, cao, cell_type, INaCa_Multiplier, INaCa_fractionSS);

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
	double hca=exp((qca*Vm*F)/(R*T));
	double hna=exp((qna*Vm*F)/(R*T));
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
	if(cell_type==100){
	    Gncx=Gncx*1.1;
	}
	else if(cell_type==101){
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
	// INaK = getINaK_ORd2011(Vm, F, R, T, nai, nao, ki, ko, cell_type, INaK_Multiplier);

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
	const double Knai=Knai0*exp((delta*Vm*F)/(3.0*R*T));
	const double Knao=Knao0*exp(((1.0-delta)*Vm*F)/(3.0*R*T));
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
	if(cell_type==100){
	    Pnak=Pnak*0.9;
	}
	else if(cell_type==101){
	    Pnak=Pnak*0.7;
	}

	const double INaK=Pnak*(zna*JnakNa+zk*JnakK);


	// Minor/background currents
	// calculate IKb
	const double xkb=1.0/(1.0+exp(-(Vm-10.8968)/(23.9871)));
	double GKb=0.0189*IKb_Multiplier;
	if(cell_type==100){
	    GKb=GKb*0.6;
	}
	const double IKb=GKb*xkb*(Vm-EK);

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

	const double I_ClCa_junc = Fjunc*GClCa/(1+KdClCa/cass)*(Vm-eclss);
	const double I_ClCa_sl = Fsl*GClCa/(1+KdClCa/cai)*(Vm-ecl);

	const double I_ClCa = I_ClCa_junc+I_ClCa_sl;
	const double I_Clbk = GClB*(Vm-ecl);

	// Calcium handling
	// calculate ryanodione receptor calcium induced calcium release from the jsr
	const double fJrelp=(1.0/(1.0+KmCaMK/CaMKa));

	// Jrel
	// [Jrel, dJrel_np, dJrel_p] = getJrel_ORd2011(Jrel_np, Jrel_p, ICaL_ss,cass, cajsr, fJrelp, cell_type, Jrel_Multiplier);

	const double jsrMidpoint = 1.7;

	const double bt=4.75;
	const double a_rel=0.5*bt;
	double Jrel_inf=a_rel*(-ICaL)/(1.0+pow(jsrMidpoint/cajsr,8.0));
	if(cell_type==101){
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
	if(cell_type==101){
	    Jrel_infp=Jrel_infp*1.7;
	}
	double tau_relp=btp/(1.0+0.0123/cajsr);

	if( tau_relp<0.001){
	    tau_relp=0.001;
	}

	const double dJrel_p=(Jrel_infp-Jrel_p)/tau_relp;

	const double Jrel=Jrel_Multiplier * 1.5378 * ((1.0-fJrelp)*Jrel_np+fJrelp*Jrel_p);




	const double fJupp=(1.0/(1.0+KmCaMK/CaMKa));
	// [Jup, Jleak] = getJup_ORd2011(cai, cansr, fJupp, cell_type, Jup_Multiplier);

	//calculate serca pump, ca uptake flux
	double Jupnp=Jup_Multiplier * 0.005425*cai/(cai+0.00092);
	double Jupp=Jup_Multiplier * 2.75*0.005425*cai/(cai+0.00092-0.00017);
	if(cell_type==100){
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
	if(cell_type==100){
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

	const double dki=-(ICaK_i+Ito+IKr+IKs+IK1+IKb+Istim-2.0*INaK)*Acap/(F*vmyo)+JdiffK*vss/vmyo;
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
	Variable_Derivatives[nai_Tomek] = dnai;
	Variable_Derivatives[nass_Tomek] = dnass;
	Variable_Derivatives[ki_Tomek] = dki;
	Variable_Derivatives[kss_Tomek] = dkss;
	Variable_Derivatives[cai_Tomek] = dcai;
	Variable_Derivatives[cass_Tomek] = dcass;
	Variable_Derivatives[cansr_Tomek] = dcansr;
	Variable_Derivatives[cajsr_Tomek] = dcajsr;
	Variable_Derivatives[m_Tomek] = dm;
	Variable_Derivatives[hp_Tomek] = dhp;
	Variable_Derivatives[h_Tomek] = dh;
	Variable_Derivatives[j_Tomek] = dj;
	Variable_Derivatives[jp_Tomek] = djp;
	Variable_Derivatives[mL_Tomek] = dmL;
	Variable_Derivatives[hL_Tomek] = dhL;
	Variable_Derivatives[hLp_Tomek] = dhLp;
	Variable_Derivatives[a_Tomek] = da;
	Variable_Derivatives[iF_Tomek] = diF;
	Variable_Derivatives[iS_Tomek] = diS;
	Variable_Derivatives[ap_Tomek] = dap;
	Variable_Derivatives[iFp_Tomek] = diFp;
	Variable_Derivatives[iSp_Tomek] = diSp;
	Variable_Derivatives[d_Tomek] = dd;
	Variable_Derivatives[ff_Tomek] = dff;
	Variable_Derivatives[fs_Tomek] = dfs;
	Variable_Derivatives[fcaf_Tomek] = dfcaf;
	Variable_Derivatives[fcas_Tomek] = dfcas;
	Variable_Derivatives[jca_Tomek] = djca;
	Variable_Derivatives[nca_Tomek] = dnca;
	Variable_Derivatives[nca_i_Tomek] = dnca_i;
	Variable_Derivatives[ffp_Tomek] = dffp;
	Variable_Derivatives[fcafp_Tomek] = dfcafp;
	Variable_Derivatives[xs1_Tomek] = dxs1;
	Variable_Derivatives[xs2_Tomek] = dxs2;
	Variable_Derivatives[Jrel_np_Tomek] = dJrel_np;
	Variable_Derivatives[CaMKt_Tomek] = dCaMKt;
	Variable_Derivatives[ikr_c0_Tomek] = d_ikr_c0;
	Variable_Derivatives[ikr_c1_Tomek] = d_ikr_c1;
	Variable_Derivatives[ikr_c2_Tomek] = d_ikr_c2;
	Variable_Derivatives[ikr_o_Tomek] = d_ikr_o;
	Variable_Derivatives[ikr_i_Tomek] = d_ikr_i;
	Variable_Derivatives[Jrel_p_Tomek] = dJrel_p;
	Variable_Derivatives[cli_Tomek] = dcli;
	Variable_Derivatives[clss_Tomek] = dclss;

	//The active strain and ion membrane current
	Iion = -(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab+I_ClCa+I_Clbk+Istim);
}

void TomekORudyVentUpdated::get_optional_output(const double &Vm,
							const Vector<double> &CellVariables,
							const double &t,
							const unsigned &cell_type,
							const double &Istim,
							const Vector<double> &Other_Parameters,
							const Vector<double> &Other_Variables,

							Vector<double> &Out) const
{

}




}