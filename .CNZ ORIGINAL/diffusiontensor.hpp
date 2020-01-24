#ifndef CNZDIFFUSIONTENSOR_HPP
#define CNZDIFFUSIONTENSOR_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cctype>
#include <algorithm>
#include <vector>
#include <math.h>


using namespace std;

namespace ExtractorParameters{
	double xx, y, z;        // Integers used to cast coordinates
	double AF_factor = 0.4; //0.7;

	// SAN (centre and radii?)
	double xc = 35.0;
	double yc = 157.0;
	double zc = 100.0;
	double xa =(41 - 30) / 2.0;
	double yb = (161 - 151) / 2.0;
	double zz = (110 - 90) / 2.0;
}

// #define diracdelta(a,b) ((a) == (b) ? (1) : (0)) //delta function
// #define diracdeltavect(a,n) (((a)%(n+1)) == (0) ? (1.) : (0.)) //delta function wrt row length


//overload scalar*vector
template <class T, class Q>
std::vector <T> operator* (const Q c, std::vector <T> A){
  vector<T> B;
  for(int i=0;i<A.size();i++){
  	B.push_back(c*A.at(i));
  }
  return B;
}
//overload vector + vector
template <class T>
std::vector<T> operator+ (std::vector<T> A, std::vector<T> B){
	vector<T> C;
	for(int i=0; i < A.size(); i++){
		C.push_back(A.at(i) + B.at(i));
	}
	return C;
}


class CalcDiffusionTensor{
private:
	std::vector<double> coord;
	std::vector<double> D;			//D IS A VECTOR BECAUSE OOMPH LIB NODES STORE DATA AS A VECTOR
	std::vector<double> comps;
	double dfibre, dsheet;
	unsigned cell_type;

	bool m_AF;
	double Datria, Dsan, Dsanp, divi, Dpv;

public:
	// CalcDiffusionTensor(double x, double y, double z, unsigned _cell_type, double compx, double compy, double compz, bool AF) : coord(3), D(9), m_AF(AF)
	// {
	// 	//store coordinates and cell type
	// 	coord[0] = x;
	// 	coord[1] = y;
	// 	coord[2] = z;
	// 	cell_type = _cell_type;


	// 	//calculate AF dependent diffusion variables
	// 	if ( m_AF ) {
	// 		// AF
	// 		Datria = 2.0 * 0.18 * (7 * (2 - 1)) * ExtractorParameters::AF_factor;
	// 		Dsan   = ( 0.18 * (5 * (2 - 1)) / 3.0 ) * ExtractorParameters::AF_factor;
	// 		Dsanp  = Dsan;
	// 		Dpv    = Dsan;
	// 	}
	// 	else {
	// 		// Control
	// 		Datria = 0.62*2.5 * 0.18 * (7 * (2 - 1));
	// 		Dsan   = 0.18 * (5 * (2 - 1)) / 3.0;
	// 		Dsanp  = Dsan;
	// 		Dpv    = 1.0*Datria; //1.1 * 0.18 * (7 * (2 - 1));;
	// 		divi   = 3.0;
	// 	}
	// 	divi   = 3.0;


	// 	//Calculate AA^T in vector format
	// 	std::vector<double> AAT(9);
	// 	AAT[0] = compx*compx;
	// 	AAT[1] = compx*compy;
	// 	AAT[2] = compx*compz;

	// 	AAT[3] = compy*compx;
	// 	AAT[4] = compy*compy;
	// 	AAT[5] = compy*compz;

	// 	AAT[6] = compz*compx;
	// 	AAT[7] = compz*compy;
	// 	AAT[8] = compz*compz;
	// 	//define identity matrix in vector format
	// 	std::vector<double> I(9);
	// 	for(int i=0; i<9;i++){I[i] = 0;}
	// 	I[0] = I[4] = I[8] = 1;


	// 	//begin calculation of diffusion tensor as in Haibo-Ni KCNA5_IKur_Force_Study/3D/lib/diffusion_tensor.h

	// 	//set scaling factor
	// 	float scale = 1.0;
	// 	if ( cell_type == 102)
	// 	{
	// 		scale = 2.0;
	// 	}
	// 	if (cell_type == 72 or cell_type == 98 or cell_type == 103 or cell_type == 102)
	// 	{
	// 		scale = 2.5;
	// 	}

	// 	//compute D from AA^T and cell type variables
	// 	if ( cell_type == 34 ) {//SAN
	// 		if ( std::sqrt( (ExtractorParameters::xx - ExtractorParameters::xc) * (ExtractorParameters::xx - ExtractorParameters::xc) / (ExtractorParameters::xa * ExtractorParameters::xa) + (y - ExtractorParameters::yc) * (y - ExtractorParameters::yc) / (ExtractorParameters::yb * ExtractorParameters::yb) + (z - ExtractorParameters::zc) * (z - ExtractorParameters::zc) / (ExtractorParameters::zz * ExtractorParameters::zz) ) < 0.7 ) { // SANC
	// 			D = (Dsan / divi) * I + ( Dsan - (Dsan / divi) ) * AAT;
	// 		}
	// 		else // SANP
	// 			D = (Dsanp / divi) * I + ( Dsanp - (Dsanp / divi) ) * AAT;
	// 	}
	// 	else if ( cell_type == 110   )  // PV
	// 		D = (Dpv / divi) * I    + ( Dpv - (Dpv / divi) ) * AAT;
	// 	else
	// 		D = (Datria / divi) * I + ( Datria - (Datria / divi) ) * AAT;

	// 	D=scale*D;

	// }
	CalcDiffusionTensor(double x, double y, double z, unsigned _cell_type, double compx, double compy, double compz, bool AF) : coord(3), D(9), m_AF(AF), comps(3)
	{
		//NEED TO ADD FUNCTIONALITY TO ALTER dfibre AND dsheet DEPENDING ON THE CELL TYPE, LOCATION, AND ATRIAL FIBRILATION.
		dfibre = 0.066588;
		dsheet = 0.00183037;

		//store coordinates and cell type
		coord[0] = x;
		coord[1] = y;
		coord[2] = z;
		cell_type = _cell_type;

		comps[0] = compx;
		comps[1] = compy;
		comps[2] = compz;


		//Calculate AA^T in vector format
		std::vector<double> AAT(9);
		AAT[0] = compx*compx;
		AAT[1] = compx*compy;
		AAT[2] = compx*compz;

		AAT[3] = compy*compx;
		AAT[4] = compy*compy;
		AAT[5] = compy*compz;

		AAT[6] = compz*compx;
		AAT[7] = compz*compy;
		AAT[8] = compz*compz;
		//define identity matrix in vector format
		std::vector<double> I(9);
		for(int i=0; i<9;i++){I[i] = 0;}
		I[0] = I[4] = I[8] = 1;

		for(int i=0; i<9;i++){D[i] = dsheet*I[i] + (dfibre - dsheet)*AAT[i];}
	}

	~CalcDiffusionTensor(){}

	// void set_dfibre(double Dfibre){dfibre=Dfibre;}
	// void set_dnormal(double Dsheet){dsheet=Dsheet;}
	double Coord(int i){return coord[i];}
	double D_get(int i){return D[i];}
	double Comps(int i){return comps[i];}
	unsigned Cell_type(){return cell_type;}
};




class Extractor{
private:
	std::vector<CalcDiffusionTensor> nodes;
public:
	Extractor(std::string const& cell_data_name, bool AF){

		//Extract all nodes from .txt file

		//Method 1 FOR 4 OR SO VALUES STILL RESULTS IN e-314 AS EXTRACTED VALUE, NOT SURE WHY
		std::ifstream txtfile(cell_data_name);
		double locx, locy, locz, componentx, componenty, componentz;
		unsigned type;
		if(!txtfile.is_open()){
			std::cout << "File " << cell_data_name << " cannot be opened" << std::endl;
			exit(0);
		}
		unsigned lecounter = 0;
		while(txtfile.good()){
			txtfile >> locx;
			txtfile >> locy;
			txtfile >> locz;
			txtfile >> type;
			txtfile >> componentx;
			txtfile >> componenty;
			txtfile >> componentz;
			if(!txtfile.eof()){
				nodes.push_back(CalcDiffusionTensor(locx, locy, locz, type, componentx, componenty, componentz, AF));
			}
			lecounter++;
		}


		txtfile.close();
		//Report nodes contained in file (not all will be used)
		std::cout << "file nodes: " << nodes.size()<<std::endl;

	}

	~Extractor(){}
	
	std::vector<CalcDiffusionTensor> ExtractorNodes(){return nodes;}

};


#endif