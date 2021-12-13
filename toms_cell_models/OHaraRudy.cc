// Copyright (c) 2011-2015 by Thomas O'Hara, Yoram Rudy, 
//                            Washington University in St. Louis.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright 
// notice, this list of conditions and the following disclaimer in the 
// documentation and/or other materials provided with the distribution.
// 
// 3. Neither the names of the copyright holders nor the names of its
// contributors may be used to endorse or promote products derived from 
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS 
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
// HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF 
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
// THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//

// C++ Implementation of the O'Hara-Rudy dynamic (ORd) model for the
// undiseased human ventricular action potential and calcium transient
//
// The ORd model is described in the article "Simulation of the Undiseased
// Human Cardiac Ventricular Action Potential: Model Formulation and
// Experimental Validation"
// by Thomas O'Hara, Laszlo Virag, Andras Varro, and Yoram Rudy
//
// The article and supplemental materails are freely available in the
// Open Access jounal PLoS Computational Biology
// Link to Article:
// http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1002061
// 
// Email: tom.ohara@gmail.com / rudy@wustl.edu
// Web: http://rudylab.wustl.edu
// 


#include "OHaraRudy.h"


namespace oomph{

OHaraRudy::OHaraRudy()
{
	//constants
	nao=140.0;//extracellular sodium in mM
	cao=1.8;//extracellular calcium in mM
	ko=5.4;//extracellular potassium in mM

	//buffer paramaters
	BSRmax=0.047;
	KmBSR=0.00087;
	BSLmax=1.124;
	KmBSL=0.0087;
	cmdnmax=0.05;
	kmcmdn=0.00238;
	trpnmax=0.07;
	kmtrpn=0.0005;
	csqnmax=10.0;
	kmcsqn=0.8;

	//CaMK paramaters
	aCaMK=0.05;
	bCaMK=0.00068;
	CaMKo=0.05;
	KmCaM=0.0015;
	KmCaMK=0.15;

	//physical constants
	R=8314.0;
	T=310.0;
	F=96485.0;

	//cell geometry
	L=0.01;
	rad=0.0011;
	vcell=1000*3.14*rad*rad*L;
	Ageo=2*3.14*rad*rad+2*3.14*rad*L;
	Acap=2*Ageo;
	vmyo=0.68*vcell;
	vmito=0.26*vcell;
	vsr=0.06*vcell;
	vnsr=0.0552*vcell;
	vjsr=0.0048*vcell;
	vss=0.02*vcell;

	//Assign the names of the variables used by this model
	Names_Of_Cell_Variables =
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
		"hf",
		"hs",
		"j",
		"hsp",
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
		"ffp",
		"fcafp",
		"xrf",
		"xrs",
		"xs1",
		"xs2",
		"xk1",
		"Jrelnp",
		"Jrelp",
		"CaMKt"
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

OHaraRudy::~OHaraRudy()
{

}

double OHaraRudy::return_initial_state_variable(const unsigned &v, const unsigned &cell_type)
{
	double Vars[40] = {
		7.0, 7.0, 145.0, 145.0, 1.0e-4, 1.0e-4, 1.2, 1.2, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0,
		1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0};

	if(v>=0 && v<=40){return Vars[v];}
}

double OHaraRudy::return_initial_membrane_potential(const unsigned &cell_type)
{
	return -87.5;
}

void OHaraRudy::Calculate_Derivatives(const Boost_State_Type &Variables,
									const double &t,
									const unsigned &cell_type,
									const double &Istim,
									const Vector<double> &Other_Parameters,
									const Vector<double> &Other_Variables,
									Vector<double> &Variable_Derivatives,
									double &Iion)
{
	//Extract the variables
	const double nai=Variables[0];
	const double nass=Variables[1];
	const double ki=Variables[2];
	const double kss=Variables[3];
	const double cai=Variables[4];
	const double cass=Variables[5];
	const double cansr=Variables[6];
	const double cajsr=Variables[7];
	const double m=Variables[8];
	const double hf=Variables[9];
	const double hs=Variables[10];
	const double j=Variables[11];
	const double hsp=Variables[12];
	const double jp=Variables[13];
	const double mL=Variables[14];
	const double hL=Variables[15];
	const double hLp=Variables[16];
	const double a=Variables[17];
	const double iF=Variables[18];
	const double iS=Variables[19];
	const double ap=Variables[20];
	const double iFp=Variables[21];
	const double iSp=Variables[22];
	const double d=Variables[23];
	const double ff=Variables[24];
	const double fs=Variables[25];
	const double fcaf=Variables[26];
	const double fcas=Variables[27];
	const double jca=Variables[28];
	const double nca=Variables[29];
	const double ffp=Variables[30];
	const double fcafp=Variables[31];
	const double xrf=Variables[32];
	const double xrs=Variables[33];
	const double xs1=Variables[34];
	const double xs2=Variables[35];
	const double xk1=Variables[36];
	const double Jrelnp=Variables[37];
	const double Jrelp=Variables[38];
	const double CaMKt=Variables[39];

	const double v=Variables[Num_Cell_Vars];

	//introduce varaibles for reversal potentials, currents, fluxes, and CaMK
	double ENa,EK,EKs;
	double INa,INaL,Ito,ICaL,ICaNa,ICaK,IKr,IKs,IK1,INaCa_i,INaCa_ss,INaCa,INaK,IKb,INab,IpCa,ICab;
	double Jrel,Jup,Jtr,Jdiff,JdiffNa,JdiffK,Jleak;
	double CaMKa,CaMKb;


	//Calculate reversal potentials
	ENa=(R*T/F)*log(nao/nai);
	EK=(R*T/F)*log(ko/ki);
	EKs=(R*T/F)*log((ko+0.01833*nao)/(ki+0.01833*nai));


	//Calculate the rates, gates, and currents
	CaMKb=CaMKo*(1.0-CaMKt)/(1.0+KmCaM/cass);
	CaMKa=CaMKb+CaMKt;
	double vffrt=v*F*F/(R*T);
	double vfrt=v*F/(R*T);

	double mss=1.0/(1.0+exp((-(v+39.57))/9.871));
	double tm=1.0/(6.765*exp((v+11.64)/34.77)+8.552*exp(-(v+77.42)/5.955));
	// m=mss-(mss-m)*exp(-dt/tm);
	double hss=1.0/(1+exp((v+82.90)/6.086));
	double thf=1.0/(1.432e-5*exp(-(v+1.196)/6.285)+6.149*exp((v+0.5096)/20.27));
	double ths=1.0/(0.009794*exp(-(v+17.95)/28.05)+0.3343*exp((v+5.730)/56.66));
	double Ahf=0.99;
	double Ahs=1.0-Ahf;
	// hf=hss-(hss-hf)*exp(-dt/thf);
	// hs=hss-(hss-hs)*exp(-dt/ths);
	double h=Ahf*hf+Ahs*hs;
	double jss=hss;
	double tj=2.038+1.0/(0.02136*exp(-(v+100.6)/8.281)+0.3052*exp((v+0.9941)/38.45));
	// j=jss-(jss-j)*exp(-dt/tj);
	double hssp=1.0/(1+exp((v+89.1)/6.086));
	double thsp=3.0*ths;
	// hsp=hssp-(hssp-hsp)*exp(-dt/thsp);
	double hp=Ahf*hf+Ahs*hsp;
	double tjp=1.46*tj;
	// jp=jss-(jss-jp)*exp(-dt/tjp);
	double GNa=75;
	double fINap=(1.0/(1.0+KmCaMK/CaMKa));
	INa=GNa*(v-ENa)*m*m*m*((1.0-fINap)*h*j+fINap*hp*jp);
	  
	double mLss=1.0/(1.0+exp((-(v+42.85))/5.264));
	double tmL=tm;
	// mL=mLss-(mLss-mL)*exp(-dt/tmL);
	double hLss=1.0/(1.0+exp((v+87.61)/7.488));
	double thL=200.0;
	// hL=hLss-(hLss-hL)*exp(-dt/thL);
	double hLssp=1.0/(1.0+exp((v+93.81)/7.488));
	double thLp=3.0*thL;
	// hLp=hLssp-(hLssp-hLp)*exp(-dt/thLp);
	double GNaL=0.0075;
	if (cell_type==100)
	{
	GNaL*=0.6;
	}
	double fINaLp=(1.0/(1.0+KmCaMK/CaMKa));
	INaL=GNaL*(v-ENa)*mL*((1.0-fINaLp)*hL+fINaLp*hLp);

	double ass=1.0/(1.0+exp((-(v-14.34))/14.82));
	double ta=1.0515/(1.0/(1.2089*(1.0+exp(-(v-18.4099)/29.3814)))+3.5/(1.0+exp((v+100.0)/29.3814)));
	// a=ass-(ass-a)*exp(-dt/ta);
	double iss=1.0/(1.0+exp((v+43.94)/5.711));
	double delta_epi;
	if (cell_type==100)
	{
	delta_epi=1.0-(0.95/(1.0+exp((v+70.0)/5.0)));
	}
	else
	{
	delta_epi=1.0;
	}
	double tiF=4.562+1/(0.3933*exp((-(v+100.0))/100.0)+0.08004*exp((v+50.0)/16.59));
	double tiS=23.62+1/(0.001416*exp((-(v+96.52))/59.05)+1.780e-8*exp((v+114.1)/8.079));
	tiF*=delta_epi;
	tiS*=delta_epi;
	double AiF=1.0/(1.0+exp((v-213.6)/151.2));
	double AiS=1.0-AiF;
	// iF=iss-(iss-iF)*exp(-dt/tiF);
	// iS=iss-(iss-iS)*exp(-dt/tiS);
	double i=AiF*iF+AiS*iS;
	double assp=1.0/(1.0+exp((-(v-24.34))/14.82));
	// ap=assp-(assp-ap)*exp(-dt/ta);
	double dti_develop=1.354+1.0e-4/(exp((v-167.4)/15.89)+exp(-(v-12.23)/0.2154));
	double dti_recover=1.0-0.5/(1.0+exp((v+70.0)/20.0));
	double tiFp=dti_develop*dti_recover*tiF;
	double tiSp=dti_develop*dti_recover*tiS;
	// ap=assp-(assp-ap)*exp(-dt/ta);
	// iSp=iss-(iss-iSp)*exp(-dt/tiSp);
	double ip=AiF*iFp+AiS*iSp;
	double Gto=0.02;
	if (cell_type==100)
	{
	Gto*=4.0;
	}
	if (cell_type==101)
	{
	Gto*=4.0;
	}
	double fItop=(1.0/(1.0+KmCaMK/CaMKa));
	Ito=Gto*(v-EK)*((1.0-fItop)*a*i+fItop*ap*ip);

	double dss=1.0/(1.0+exp((-(v+3.940))/4.230));
	double td=0.6+1.0/(exp(-0.05*(v+6.0))+exp(0.09*(v+14.0)));
	// d=dss-(dss-d)*exp(-dt/td);
	double fss=1.0/(1.0+exp((v+19.58)/3.696));
	double tff=7.0+1.0/(0.0045*exp(-(v+20.0)/10.0)+0.0045*exp((v+20.0)/10.0));
	double tfs=1000.0+1.0/(0.000035*exp(-(v+5.0)/4.0)+0.000035*exp((v+5.0)/6.0));
	double Aff=0.6;
	double Afs=1.0-Aff;
	// ff=fss-(fss-ff)*exp(-dt/tff);
	// fs=fss-(fss-fs)*exp(-dt/tfs);
	double f=Aff*ff+Afs*fs;
	double fcass=fss;
	double tfcaf=7.0+1.0/(0.04*exp(-(v-4.0)/7.0)+0.04*exp((v-4.0)/7.0));
	double tfcas=100.0+1.0/(0.00012*exp(-v/3.0)+0.00012*exp(v/7.0));
	double Afcaf=0.3+0.6/(1.0+exp((v-10.0)/10.0));
	double Afcas=1.0-Afcaf;
	// fcaf=fcass-(fcass-fcaf)*exp(-dt/tfcaf);
	// fcas=fcass-(fcass-fcas)*exp(-dt/tfcas);
	double fca=Afcaf*fcaf+Afcas*fcas;
	double tjca=75.0;
	// jca=fcass-(fcass-jca)*exp(-dt/tjca);
	double tffp=2.5*tff;
	// ffp=fss-(fss-ffp)*exp(-dt/tffp);
	double fp=Aff*ffp+Afs*fs;
	double tfcafp=2.5*tfcaf;
	// fcafp=fcass-(fcass-fcafp)*exp(-dt/tfcafp);
	double fcap=Afcaf*fcafp+Afcas*fcas;
	double Kmn=0.002;
	double k2n=1000.0;
	double km2n=jca*1.0;
	double anca=1.0/(k2n/km2n+pow(1.0+Kmn/cass,4.0));
	// nca=anca*k2n/km2n-(anca*k2n/km2n-nca)*exp(-km2n*dt);
	//These are undefined for vffrt = vfrt = v = 0 so we need to be careful, use L'Hopital
	double PhiCaL;
	if(abs(v)<1e-9){PhiCaL = 2*F*(cass-0.341*cao);}
	else{PhiCaL=4.0*vffrt*(cass*exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);}
	double PhiCaNa;
	if(abs(v)<1e-9){PhiCaNa = 0.75*F*(nass - nao);}
	else{PhiCaNa=0.75*vffrt*(nass*exp(vfrt)-nao)/(exp(vfrt)-1.0);}
	double PhiCaK;
	if(abs(v)<1e-9){PhiCaK = 0.75*F*(kss-ko);}
	else{PhiCaK=0.75*vffrt*(kss*exp(vfrt)-ko)/(exp(vfrt)-1.0);}
	double zca=2.0;
	double PCa=0.0001;
	if (cell_type==100)
	{
	PCa*=1.2;
	}
	if (cell_type==101)
	{
	PCa*=2.5;
	}
	double PCap=1.1*PCa;
	double PCaNa=0.00125*PCa;
	double PCaK=3.574e-4*PCa;
	double PCaNap=0.00125*PCap;
	double PCaKp=3.574e-4*PCap;
	double fICaLp=(1.0/(1.0+KmCaMK/CaMKa));
	ICaL=(1.0-fICaLp)*PCa*PhiCaL*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCap*PhiCaL*d*(fp*(1.0-nca)+jca*fcap*nca);
	ICaNa=(1.0-fICaLp)*PCaNa*PhiCaNa*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaNap*PhiCaNa*d*(fp*(1.0-nca)+jca*fcap*nca);
	ICaK=(1.0-fICaLp)*PCaK*PhiCaK*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaKp*PhiCaK*d*(fp*(1.0-nca)+jca*fcap*nca);

	double xrss=1.0/(1.0+exp((-(v+8.337))/6.789));
	double txrf=12.98+1.0/(0.3652*exp((v-31.66)/3.869)+4.123e-5*exp((-(v-47.78))/20.38));
	double txrs=1.865+1.0/(0.06629*exp((v-34.70)/7.355)+1.128e-5*exp((-(v-29.74))/25.94));
	double Axrf=1.0/(1.0+exp((v+54.81)/38.21));
	double Axrs=1.0-Axrf;
	// xrf=xrss-(xrss-xrf)*exp(-dt/txrf);
	// xrs=xrss-(xrss-xrs)*exp(-dt/txrs);
	double xr=Axrf*xrf+Axrs*xrs;
	double rkr=1.0/(1.0+exp((v+55.0)/75.0))*1.0/(1.0+exp((v-10.0)/30.0));
	double GKr=0.046;
	if (cell_type==100)
	{
	GKr*=1.3;
	}
	if (cell_type==101)
	{
	GKr*=0.8;
	}
	IKr=GKr*sqrt(ko/5.4)*xr*rkr*(v-EK);

	double xs1ss=1.0/(1.0+exp((-(v+11.60))/8.932));
	double txs1=817.3+1.0/(2.326e-4*exp((v+48.28)/17.80)+0.001292*exp((-(v+210.0))/230.0));
	// xs1=xs1ss-(xs1ss-xs1)*exp(-dt/txs1);
	double xs2ss=xs1ss;
	double txs2=1.0/(0.01*exp((v-50.0)/20.0)+0.0193*exp((-(v+66.54))/31.0));
	// xs2=xs2ss-(xs2ss-xs2)*exp(-dt/txs2);
	double KsCa=1.0+0.6/(1.0+pow(3.8e-5/cai,1.4));
	double GKs=0.0034;
	if (cell_type==100)
	{
	GKs*=1.4;
	}
	IKs=GKs*KsCa*xs1*xs2*(v-EKs);

	double xk1ss=1.0/(1.0+exp(-(v+2.5538*ko+144.59)/(1.5692*ko+3.8115)));
	double txk1=122.2/(exp((-(v+127.2))/20.36)+exp((v+236.8)/69.33));
	// xk1=xk1ss-(xk1ss-xk1)*exp(-dt/txk1);
	double rk1=1.0/(1.0+exp((v+105.8-2.6*ko)/9.493));
	double GK1=0.1908;
	if (cell_type==100)
	{
	GK1*=1.2;
	}
	if (cell_type==101)
	{
	GK1*=1.3;
	}
	IK1=GK1*sqrt(ko)*rk1*xk1*(v-EK);

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
	double Gncx=0.0008;
	if (cell_type==100)
	{
	Gncx*=1.1;
	}
	if (cell_type==101)
	{
	Gncx*=1.4;
	}
	INaCa_i=0.8*Gncx*allo*(zna*JncxNa+zca*JncxCa);

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
	KmCaAct=150.0e-6;
	allo=1.0/(1.0+pow(KmCaAct/cass,2.0));
	JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
	JncxCa=E2*k2-E1*k1;
	INaCa_ss=0.2*Gncx*allo*(zna*JncxNa+zca*JncxCa);

	INaCa=INaCa_i+INaCa_ss;

	double k1p=949.5;
	double k1m=182.4;
	double k2p=687.2;
	double k2m=39.4;
	k3p=1899.0;
	double k3m=79300.0;
	k4p=639.0;
	double k4m=40.0;
	double Knai0=9.073;
	double Knao0=27.78;
	double delta=-0.1550;
	double Knai=Knai0*exp((delta*v*F)/(3.0*R*T));
	double Knao=Knao0*exp(((1.0-delta)*v*F)/(3.0*R*T));
	double Kki=0.5;
	double Kko=0.3582;
	double MgADP=0.05;
	double MgATP=9.8;
	double Kmgatp=1.698e-7;
	double H=1.0e-7;
	double eP=4.2;
	double Khp=1.698e-7;
	double Knap=224.0;
	double Kxkur=292.0;
	double P=eP/(1.0+H/Khp+nai/Knap+ki/Kxkur);
	double a1=(k1p*pow(nai/Knai,3.0))/(pow(1.0+nai/Knai,3.0)+pow(1.0+ki/Kki,2.0)-1.0);
	double b1=k1m*MgADP;
	double a2=k2p;
	double b2=(k2m*pow(nao/Knao,3.0))/(pow(1.0+nao/Knao,3.0)+pow(1.0+ko/Kko,2.0)-1.0);
	double a3=(k3p*pow(ko/Kko,2.0))/(pow(1.0+nao/Knao,3.0)+pow(1.0+ko/Kko,2.0)-1.0);
	double b3=(k3m*P*H)/(1.0+MgATP/Kmgatp);
	double a4=(k4p*MgATP/Kmgatp)/(1.0+MgATP/Kmgatp);
	double b4=(k4m*pow(ki/Kki,2.0))/(pow(1.0+nai/Knai,3.0)+pow(1.0+ki/Kki,2.0)-1.0);
	x1=a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2;
	x2=b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4;
	x3=a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1;
	x4=b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1;
	E1=x1/(x1+x2+x3+x4);
	E2=x2/(x1+x2+x3+x4);
	E3=x3/(x1+x2+x3+x4);
	E4=x4/(x1+x2+x3+x4);
	double zk=1.0;
	double JnakNa=3.0*(E1*a3-E2*b3);
	double JnakK=2.0*(E4*b1-E3*a1);
	double Pnak=30;
	if (cell_type==100)
	{
	    Pnak*=0.9;
	}
	if (cell_type==101)
	{
	    Pnak*=0.7;
	}
	INaK=Pnak*(zna*JnakNa+zk*JnakK);

	double xkb=1.0/(1.0+exp(-(v-14.48)/18.34));
	double GKb=0.003;
	if (cell_type==100)
	{
	    GKb*=0.6;
	}
	IKb=GKb*xkb*(v-EK);

	double PNab=3.75e-10;
	//These are undefined for vffrt = vfrt = v = 0
	if(abs(v)<1e-9){INab = PNab*F*(nai-nao);}
	else{INab=PNab*vffrt*(nai*exp(vfrt)-nao)/(exp(vfrt)-1.0);}
	double PCab=2.5e-8;
	if(abs(v)<1e-9){ICab = 2.0*PCab*F*(cai - 0.341*cao);}
	else{ICab=PCab*4.0*vffrt*(cai*exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);}

	double GpCa=0.0005;
	IpCa=GpCa*cai/(0.0005+cai);


	//calculate fluxes, buffers, and concentrations

	CaMKb=CaMKo*(1.0-CaMKt)/(1.0+KmCaM/cass);
	CaMKa=CaMKb+CaMKt;
	// CaMKt+=dt*(aCaMK*CaMKb*(CaMKb+CaMKt)-bCaMK*CaMKt);

	JdiffNa=(nass-nai)/2.0;
	JdiffK=(kss-ki)/2.0;
	Jdiff=(cass-cai)/0.2;

	double bt=4.75;
	double a_rel=0.5*bt;
	double Jrel_inf=a_rel*(-ICaL)/(1.0+pow(1.5/cajsr,8.0));
	if (cell_type==101)
	{
	Jrel_inf*=1.7;
	}
	double tau_rel=bt/(1.0+0.0123/cajsr);
	if (tau_rel<0.005)
	{
	tau_rel=0.005;
	}
	// Jrelnp=Jrel_inf-(Jrel_inf-Jrelnp)*exp(-dt/tau_rel);
	double btp=1.25*bt;
	double a_relp=0.5*btp;
	double Jrel_infp=a_relp*(-ICaL)/(1.0+pow(1.5/cajsr,8.0));
	if (cell_type==101)
	{
	Jrel_infp*=1.7;
	}
	double tau_relp=btp/(1.0+0.0123/cajsr);
	if (tau_relp<0.005)
	{
	tau_relp=0.005;
	}
	// Jrelp=Jrel_infp-(Jrel_infp-Jrelp)*exp(-dt/tau_relp);
	double fJrelp=(1.0/(1.0+KmCaMK/CaMKa));
	Jrel=(1.0-fJrelp)*Jrelnp+fJrelp*Jrelp;

	double Jupnp=0.004375*cai/(cai+0.00092);
	double Jupp=2.75*0.004375*cai/(cai+0.00092-0.00017);
	if (cell_type==100)
	{
	Jupnp*=1.3;
	Jupp*=1.3;
	}
	double fJupp=(1.0/(1.0+KmCaMK/CaMKa));
	Jleak=0.0039375*cansr/15.0;
	Jup=(1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak;

	Jtr=(cansr-cajsr)/100.0;

	// nai+=dt*(-(INa+INaL+3.0*INaCa_i+3.0*INaK+INab)*Acap/(F*vmyo)+JdiffNa*vss/vmyo);
	// nass+=dt*(-(ICaNa+3.0*INaCa_ss)*Acap/(F*vss)-JdiffNa);

	// ki+=dt*(-(Ito+IKr+IKs+IK1+IKb+Istim-2.0*INaK)*Acap/(F*vmyo)+JdiffK*vss/vmyo);
	// kss+=dt*(-(ICaK)*Acap/(F*vss)-JdiffK);

	double Bcai;
	if (cell_type==100)
	{
	Bcai=1.0/(1.0+1.3*cmdnmax*kmcmdn/pow(kmcmdn+cai,2.0)+trpnmax*kmtrpn/pow(kmtrpn+cai,2.0));
	}
	else
	{
	Bcai=1.0/(1.0+cmdnmax*kmcmdn/pow(kmcmdn+cai,2.0)+trpnmax*kmtrpn/pow(kmtrpn+cai,2.0));
	}
	// cai+=dt*(Bcai*(-(IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo));

	double Bcass=1.0/(1.0+BSRmax*KmBSR/pow(KmBSR+cass,2.0)+BSLmax*KmBSL/pow(KmBSL+cass,2.0));
	// cass+=dt*(Bcass*(-(ICaL-2.0*INaCa_ss)*Acap/(2.0*F*vss)+Jrel*vjsr/vss-Jdiff));

	// cansr+=dt*(Jup-Jtr*vjsr/vnsr);

	double Bcajsr=1.0/(1.0+csqnmax*kmcsqn/pow(kmcsqn+cajsr,2.0));
	// cajsr+=dt*(Bcajsr*(Jtr-Jrel));


	//Set the new values of the variables
	Variable_Derivatives[0] = (-(INa+INaL+3.0*INaCa_i+3.0*INaK+INab)*Acap/(F*vmyo)+JdiffNa*vss/vmyo);
	Variable_Derivatives[1] = (-(ICaNa+3.0*INaCa_ss)*Acap/(F*vss)-JdiffNa);
	Variable_Derivatives[2] = (-(Ito+IKr+IKs+IK1+IKb+Istim-2.0*INaK)*Acap/(F*vmyo)+JdiffK*vss/vmyo);
	Variable_Derivatives[3] = (-(ICaK)*Acap/(F*vss)-JdiffK);
	Variable_Derivatives[4] = (Bcai*(-(IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo));
	Variable_Derivatives[5] = (Bcass*(-(ICaL-2.0*INaCa_ss)*Acap/(2.0*F*vss)+Jrel*vjsr/vss-Jdiff));
	Variable_Derivatives[6] = (Jup-Jtr*vjsr/vnsr);
	Variable_Derivatives[7] = (Bcajsr*(Jtr-Jrel));
	Variable_Derivatives[8] = (1.0/tm)*(mss - m);
	Variable_Derivatives[9] = (1.0/thf)*(hss-hf);
	Variable_Derivatives[10] = (1.0/ths)*(hss-hs);
	Variable_Derivatives[11] = (1.0/tj)*(jss-j);
	Variable_Derivatives[12] = (1.0/thsp)*(hssp-hsp);
	Variable_Derivatives[13] = (1.0/tjp)*(jss-jp);
	Variable_Derivatives[14] = (1.0/tmL)*(mLss-mL);
	Variable_Derivatives[15] = (1.0/thL)*(hLss-hL);
	Variable_Derivatives[16] = (1.0/thLp)*(hLssp - hLp);
	Variable_Derivatives[17] = (1.0/ta)*(ass-a);
	Variable_Derivatives[18] = (1.0/tiF)*(iss-iF);
	Variable_Derivatives[19] = (1.0/tiS)*(iss-iS);
	Variable_Derivatives[20] = (1.0/ta)*(assp-ap);
	Variable_Derivatives[21] = (1.0/ta)*(assp-ap);
	Variable_Derivatives[22] = (1.0/tiSp)*(iss-iSp);
	Variable_Derivatives[23] = (1.0/td)*(dss-d);
	Variable_Derivatives[24] = (1.0/tff)*(fss-ff);
	Variable_Derivatives[25] = (1.0/tfs)*(fss-fs);
	Variable_Derivatives[26] = (1.0/tfcaf)*(fcass-fcaf);
	Variable_Derivatives[27] = (1.0/tfcas)*(fcass-fcas);
	Variable_Derivatives[28] = (1.0/tjca)*(fcass-jca);
	Variable_Derivatives[29] = (1.0/km2n)*(anca*k2n/km2n-nca);
	Variable_Derivatives[30] = (1.0/tffp)*(fss-ffp);
	Variable_Derivatives[31] = (1.0/tfcafp)*(fcass-fcafp);
	Variable_Derivatives[32] = (1.0/txrf)*(xrss-xrf);
	Variable_Derivatives[33] = (1.0/txrs)*(xrss-xrs);
	Variable_Derivatives[34] = (1.0/txs1)*(xs1ss-xs1);
	Variable_Derivatives[35] = (1.0/txs2)*(xs2ss-xs2);
	Variable_Derivatives[36] = (1.0/txk1)*(xk1ss-xk1);
	Variable_Derivatives[37] = (1.0/tau_rel)*(Jrel_inf-Jrelnp);
	Variable_Derivatives[38] = (1.0/tau_relp)*(Jrel_infp-Jrelp);
	Variable_Derivatives[39] = (aCaMK*CaMKb*(CaMKb+CaMKt)-bCaMK*CaMKt);

	//The active strain and ion membrane current
	Iion = -(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa+INaK+INab+IKb+IpCa+ICab+Istim);

	// //If a derivative is non-finite then set it to some random value
	// for(unsigned i=0; i<Num_Cell_Vars; i++){
	// 	Variable_Derivatives[i] = (std::isfinite(Variable_Derivatives[i]))*Variable_Derivatives[i]
	// 							+ (!std::isfinite(Variable_Derivatives[i]))*((double)std::rand()/(RAND_MAX));
	// }
	// //Do the same for the ionic current
	// Iion = (std::isfinite(Iion))*Iion
	// 	 + (!std::isfinite(Iion))*((double)std::rand()/(RAND_MAX));
}


void OHaraRudy::get_optional_output(const Boost_State_Type &Variables,
								const double &t,
								const unsigned &cell_type,
								const double &Istim,
								const Vector<double> &Other_Parameters,
								const Vector<double> &Other_Variables,
								Vector<double> &Out)
{
	//Intentionally empty
}

}; //End namespace