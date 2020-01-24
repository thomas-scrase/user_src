/* ENUMSIMULATIONCTRL_HPP
*	enum type incorporating Atrial regional cell types,
*	IKur mutation types
*	and AF types.

*   Haibo NI qiangzi.ni@gmail.com
*   Date: 21 January 2015 (Wednesday)
    15:04
*/


#ifndef ENUMSIMULATIONCTRL_HPP
#define ENUMSIMULATIONCTRL_HPP

#include  <string>
#include  <iostream>
#include  <cstdlib>

typedef enum AtriaCellType {RA=0, PM, CT, RAA, AS, AVR, BB, LA, LAA, PV, SAN_C, SAN_P,AVN};
typedef enum MuationType {WT=0, D322H, E48G, A305T, Y155C, D469E, P488S};
typedef enum AFType {NONE=0, AF1, AF2, AF3, AF4};

AtriaCellType GetCellTypeFromLabel(int tissue) ;
AtriaCellType GetCellTypeFromLabel(std::string tissue) ;
MuationType GetMutationTypeFromLabel(std::string tissue) ;
std::string GetCellTypeToString(AtriaCellType tissue) ;
std::string GetMutationTypeToString(MuationType Mut);

#endif