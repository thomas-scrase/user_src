/*
 * stimulus.h
 * Haibo Ni <qiangzini@gmail.com>
 *
 * generate stimulus according to S1 or S1S2 proctol
 *
 * measure ERP
 *
 */

#ifndef STIMULUS_H
#define STIMULUS_H

#include <math.h>
#include <stdio.h>
// #include <iostream>
// #include <fstream>
// #include "vector"

/*  all the times are in ms */
typedef struct 
{
	float current_time;
	float BCL;
	float stim_duration;
	float stim_start_time;
	float time_since_last_stim;
	float time_since_1st_stim;
} Stim_Control;

Stim_Control * creat_stim_control(float BCL, float stim_start_time, float stim_duration);
 /* all current in pA/pF */
float Stim_flag_S1(Stim_Control *stim_infor, double current_time, double dt);
double S1(double stim_start, double stim, double BCL, double current_time, double stim_duration);

double S1S2(double stim_start_time, double stim_current, double BCL_s1, int S1_num, double BCL_s2, double current_time, double stim_duration);
double S1S2_num(double stim_start_time, double stim_current, double BCL_s1, int S1_num, double BCL_s2, int S2_num, double current_time, double stim_duration) ;

double S1S2_num_stim(double stim_start_time, double BCL_s1, double stim_current_1, int S1_num, 
						double BCL_s2, double stim_current_2, int S2_num,
						double current_time, double stim_duration) ;


/*typedef struct
{
	double ERP;
	double BCL_s1, BCL_s2;
	double ap_amp_ratio;
	double BCL_adjust;

} ERP_measure_struct;
*/


#endif


 // /* END OF STIMULUS.H 