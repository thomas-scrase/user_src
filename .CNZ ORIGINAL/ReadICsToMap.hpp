


#include <iostream>
#include <map>
#include <tuple>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cmath>
#ifndef READICSTOMAP_HPP
#define READICSTOMAP_HPP

#include <map>
#include <tuple>
#include "EnumSimulationCtrl.hpp"
#include <fstream>

typedef std::tuple<AtriaCellType, MuationType> CellInfo;
typedef std::pair < CellInfo, std::vector<double> > Cell_ICs_Pair;  // integer to find the coordinate is typically easier and faster, by haibo Fri 30 Oct 2015 15:35:38 GMT
typedef std::map < CellInfo, std::vector<double> > Cell_ICs_Map;  // integer to find the coordinate is typically easier and faster, by haibo Fri 30 Oct 2015 15:35:38 GMT


Cell_ICs_Pair Create_CellICs_FromFile(AtriaCellType cell_type, MuationType mutation_type) {

	auto file_name = std::string("ICs/") + GetCellTypeToString(cell_type) + std::string("__") +  GetMutationTypeToString(mutation_type) + "_ICs.dat";
	// InitialiseElectricalStatesFromFile(outICs_name.c_str());


	std::fstream readin(file_name.c_str(), std::ios_base::in);
	if (!readin.is_open())
	{
		std::cerr << "Failed to open initial data file " << file_name << std::endl;
		std::exit(0);
	}
	// int i = 0;

	std::vector<double> vec;
	double temp;

	for (int i = 0; i < 50; ++i)
	{
		readin >> temp;
		vec.push_back(temp);
	}
	Cell_ICs_Pair T2;

	auto t1 = std::make_tuple(cell_type, mutation_type);
	T2 = std::make_pair(t1, vec);

	return T2;

}

void Create_Cell_ICs_Map_FromFile(Cell_ICs_Map &ICMap, MuationType mutation_type) {
	std::vector<AtriaCellType> celltypes = {RA, PM, CT, RAA, AS, AVR, BB, LA, LAA, PV};

	for(int i = 0 ;i < celltypes.size();i++) {


		ICMap.insert(Create_CellICs_FromFile(celltypes[i], mutation_type));
	}


}








#endif