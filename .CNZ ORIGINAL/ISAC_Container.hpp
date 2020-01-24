/* ISAC_CONTAINER_HPP
*  single class container to calculate stretch activated channel;
*  the permiability to Na+, K+ and Ca2+ are 1:1:1 in default;

*   Haibo NI qiangzi.ni@gmail.com
*   Date: 21 January 2015 (Wednesday)
    14:38
*/

#ifndef ISAC_CONTAINER_HPP
#define ISAC_CONTAINER_HPP
#include <iostream>
#include <exception>
#include <cmath>

class ISAC_Container
{
public:

    ISAC_Container(int pna= 1, int pk = 1, int pca = 1) {
        PK      = pk;
        PNa     = pna;
        PCa     = pca;
        Isac_Ca = 0.0;
        Isac_Na = 0.0;
        Isac_K  = 0.0;
        ISAC = 0.0;
    };

    ~ISAC_Container() {};

    void Compute_ISAC(double strain, double Vm, double Cm) {
        double GSac = 0.00485 * 1.25; // need to be updated  ms/uF
        double Pm;
        const double Ke = 0.026590;// new fitting, was 0.02; // Dr Ismail's thesis
        const double strain_half = 0.163;// using 0.233; // Modeling of Stretch-Activated Sarcolemmal Channels in Smooth Muscle Cells
        // http://link.springer.com/chapter/10.1007%2F978-3-642-03882-2_197#page-1
        // 0.163 Role of Stretch-activated Channels in the Heart: Action Potential and Ca2+ Transients
        // in http://www.ncbi.nlm.nih.gov/books/NBK7490/
        const double Esac = -1.0; //mv see  Dr Ismail's thesis
        Pm = 1.0 / (1 + exp(-(strain - strain_half) / Ke));
        ISAC = Cm * (GSac * Pm * (Vm - Esac));
        Isac_Na, Isac_K, Isac_Ca;
        int total_Sca = PNa + PCa + PK;
        Isac_Ca = PCa * ISAC / total_Sca;
        Isac_Na = PNa * ISAC / total_Sca;
        Isac_K = PK * ISAC / total_Sca;
    };


    double ISAC, Isac_Na, Isac_K, Isac_Ca;
    int PK, PNa, PCa;
};

#endif 