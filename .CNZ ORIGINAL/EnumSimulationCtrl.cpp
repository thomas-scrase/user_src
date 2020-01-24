/* ENUMSIMULATIONCTRL_HPP
*	enum type incorporating Atrial regional cell types,
*	IKur mutation types
*	and AF types.

*   Haibo NI qiangzi.ni@gmail.com
*   Date: 21 January 2015 (Wednesday)
    15:04
*/



#include  <string>
#include  <iostream>
#include  <cstdlib>

#include "EnumSimulationCtrl.hpp"


AtriaCellType GetCellTypeFromLabel(int tissue) {
	if (tissue == 10) { // SAN
		return SAN_C;
	}
	else if (tissue == 11) { // CT
		return CT;

	}
	else if (tissue == 12) { // PM
		return PM;
	}
	else if (tissue == 202) { // RAA / APG
		return RAA;
	}
	else if (tissue == 201) { //LAA
		return LAA;
	}
	else if (tissue == 13 || tissue == 14 || tissue == 200) // RA (no change)11
	{
		return RA;
	}
	else if (tissue == 15) { // BB
		return BB;
	}
	else if (tissue == 16) { // LA
		return LA;
	}
	else if (tissue == 17) { // AS
		return AS;
	}
	else if (tissue == 18) { // AVR
		return AVR;
	}
	else if (tissue == 106) { // AVN - compact node
		return AVN;
	}
	else if (tissue == 104) { // AVN - inferior nodal extension
		return AVN;
	}

	else if (tissue == 105) { // AVN - Penetrating bundle
		return AVN;
	}
	else if (tissue == 102 || tissue == 103) { // AVN - transitional area
		return AVN;

	}
	else if (tissue == 101) { // PV
		return PV;
	}
	else {

		std::cout << ("wrong cell type!!!\n");
		std::exit(0);
	}
}



AtriaCellType GetCellTypeFromLabel(std::string tissue) {
	if (tissue == "SAN_C") { // SAN
		return SAN_C;
	}
	else if (tissue == "CT") { // CT
		return CT;
	}
	else if (tissue == "PM") { // PM
		return PM;
	}
	else if (tissue == "RAA") { // RAA / APG
		return RAA;
	}
	else if (tissue == "LAA") { //LAA
		return LAA;
	}
	else if (tissue == "RA" ) { // RA (no change)11
		return RA;
	}
	else if (tissue == "BB") { // BB
		return BB;
	}
	else if (tissue == "LA") { // LA
		return LA;
	}
	else if (tissue == "AS") { // AS
		return AS;
	}
	else if (tissue == "AVR") { // AVR
		return AVR;
	}
	else if (tissue == "AVN") { // AVN - compact node
		return AVN;
	}
	else if (tissue == "PV") { // PV
		return PV;
	}
	else {
		std::cout << ("wrong cell type!!!\n");
		std::exit(0);
	}
}

MuationType GetMutationTypeFromLabel(std::string Mut) {
	if (Mut == "WT") { // SAN
		return WT;
	}
	else if (Mut == "D322H") { // D322H
		return D322H;

	}
	else if (Mut == "E48G") { // E48G
		return E48G;
	}
	else if (Mut == "A305T") { // A305T / APG
		return A305T;
	}
	else if (Mut == "D469E") { //D469E
		return D469E;
	}
	else if (Mut == "Y155C" ) // Y155C (no change)11
	{
		return Y155C;
	}
	else if (Mut == "P488S") { // P488S
		return P488S;
	} else {

		std::cout << ("wrong Mutation type!!!\n");
		std::exit(0);
	}
}



std::string GetCellTypeToString(AtriaCellType tissue)  {
	if (tissue == SAN_C) { // SAN
		return "SAN_C";
	}
	else if (tissue == CT ) { // CT
		return "CT";

	}
	else if (tissue == PM ) { // PM
		return "PM";
	}
	else if (tissue == RAA ) { // RAA / APG
		return "RAA";
	}
	else if (tissue == LAA ) { //LAA
		return "LAA";
	}
	else if (tissue == RA  ) // RA (no change)11
	{
		return "RA";
	}
	else if (tissue == BB ) { // BB
		return "BB";
	}
	else if (tissue == LA ) { // LA
		return "LA";
	}
	else if (tissue == AS ) { // AS
		return "AS";
	}
	else if (tissue == AVR ) { // AVR
		return "AVR";
	}
	else if (tissue == AVN ) { // AVN - compact node
		return "AVN";
	}


	else if (tissue == PV ) { // PV
		return "PV";
	}
	else {

		std::cout << ( "wrong cell type!!!\n" );
		std::exit(0);
	}
}


std::string GetMutationTypeToString(MuationType Mut)  {
	if (Mut == WT) { // SAN
		return "WT";
	}
	else if (Mut == D322H) { // D322H
		return "D322H";

	}
	else if (Mut == E48G) { // E48G
		return "E48G";
	}
	else if (Mut == A305T) { // A305T / APG
		return "A305T";
	}
	else if (Mut == D469E) { //D469E
		return "D469E";
	}
	else if (Mut == Y155C ) // Y155C (no change)11
	{
		return "Y155C";
	}
	else if (Mut == P488S) { // P488S
		return "P488S";
	} else {

		std::cout << ("wrong Mutation type!!!\n");
		std::exit(0);
	}
}

