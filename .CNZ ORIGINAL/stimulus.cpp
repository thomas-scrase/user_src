/*
 * stimulus.h
 * Haibo Ni <qiangzini@gmail.com>
 *
 * generate stimulus according to S1 or S1S2 proctol
 *
 * measure ERP
 *
 */

#ifndef STIMULUS_C
#define STIMULUS_C

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "stimulus.h"
// #include <iostream>
// #include <fstream>
// #include "vector"
Stim_Control *creat_stim_control(float BCL, float stim_start_time, float stim_duration) {
    Stim_Control *stim_infor = (Stim_Control *) malloc (sizeof (Stim_Control));
    stim_infor->BCL = BCL;
    stim_infor->stim_start_time = stim_start_time;
    stim_infor->stim_duration = stim_duration;
    stim_infor->current_time = 0.0;
    stim_infor->time_since_last_stim = 0.0;
    stim_infor->time_since_1st_stim = -1;
    return (stim_infor);
}

float Stim_flag_S1(Stim_Control *stim_infor, double current_time, double dt) {
    float flag;
    // stim_infor->time_since_1st_stim =  (current_time - stim_infor->stim_start_time) >= 0 ? (current_time - stim_infor->stim_start_time) : -1;
    if (current_time - stim_infor->stim_start_time >= 0.0)
    {
        if (stim_infor ->time_since_last_stim >= 0.0 && stim_infor ->time_since_last_stim <= 2) flag = 1.0; else flag = 0.0;
        if (stim_infor ->time_since_last_stim > stim_infor ->BCL) stim_infor->time_since_last_stim = 0.0;
        stim_infor ->time_since_last_stim +=dt;
    } else
    flag = 0.0;
    return flag;
}


/*  all the times are in ms */

/* all current in pA/pF */
/*  implementation of S1  */
double S1(double stim_start_time, double stim_current, double BCL, double current_time, double stim_duration) {
    double Istim = 0.0;
    double remain;
    double time_elapsed;

    time_elapsed = (current_time - stim_start_time) >= 0 ? (current_time - stim_start_time) : 0;
    remain =  fmod(time_elapsed, BCL);
    if (remain >= 0 && remain < stim_duration) {
        Istim = stim_current;
    } else {
        Istim = 0.0;
    }


    return Istim;
}

double S1S2(double stim_start_time, double stim_current, double BCL_s1, int S1_num, double BCL_s2, double current_time, double stim_duration) {
    double Istim = 0.0;
    double remain;
    double time_elapsed;

    time_elapsed = (current_time - stim_start_time) > 0 ? (current_time - stim_start_time) : 0;
    if (time_elapsed <= BCL_s1 * (S1_num - 1) + stim_duration)
    {
        remain =  fmod(time_elapsed, BCL_s1);
        if (remain >= 0 && remain <= stim_duration) {
            Istim = stim_current;
        } else {
            Istim = 0.0;
        }
    } else {
        time_elapsed = time_elapsed - BCL_s1 * (S1_num - 1) ;
        remain =  fmod(time_elapsed, BCL_s2);
        if (remain >= 0 && remain <= stim_duration) {
            Istim = stim_current;
        } else {
            Istim = 0.0;
        }
    }
    return Istim;
}


double S1S2_num(double stim_start_time, double stim_current, double BCL_s1, int S1_num, double BCL_s2, int S2_num, double current_time, double stim_duration) {
    double Istim = 0.0;
    double remain;
    double time_elapsed;

    time_elapsed = (current_time - stim_start_time) > 0 ? (current_time - stim_start_time) : 0;
    if (time_elapsed <= BCL_s1 * (S1_num - 1) + stim_duration)
    {
        remain =  fmod(time_elapsed, BCL_s1);
        if (remain >= 0 && remain <= stim_duration) {
            Istim = stim_current;
        } else {
            Istim = 0.0;
        }
    }
    else if ((time_elapsed > BCL_s1 * (S1_num - 1) + stim_duration) && (time_elapsed <= BCL_s1 * (S1_num - 1) + BCL_s2 * S2_num + stim_duration) )
    {
        time_elapsed = time_elapsed - BCL_s1 * (S1_num - 1) ;
        remain =  fmod(time_elapsed, BCL_s2);
        if (remain >= 0 && remain <= stim_duration) {
            Istim = stim_current;
        } else {
            Istim = 0.0;
        };
    }
    else
    {
        Istim = 0.0;
    }
    return Istim;
}

double S1S2_num_stim(double stim_start_time, double BCL_s1, double stim_current_1, int S1_num,
                     double BCL_s2, double stim_current_2, int S2_num,
                     double current_time, double stim_duration) {

    double Istim = 0.0;
    double remain;
    double time_elapsed;

    time_elapsed = (current_time - stim_start_time) > 0 ? (current_time - stim_start_time) : 0;
    if (time_elapsed <= BCL_s1 * (S1_num - 1) + stim_duration)
    {
        remain =  fmod(time_elapsed, BCL_s1);
        if (remain >= 0 && remain <= stim_duration) {
            Istim = stim_current_1;
        } else {
            Istim = 0.0;
        }
    }
    else if ((time_elapsed > BCL_s1 * (S1_num - 1) + stim_duration) && (time_elapsed <= BCL_s1 * (S1_num - 1) + BCL_s2 * S2_num + stim_duration) )
    {
        time_elapsed = time_elapsed - BCL_s1 * (S1_num - 1) ;
        remain =  fmod(time_elapsed, BCL_s2);
        if (remain >= 0 && remain <= stim_duration) {
            Istim = stim_current_2;
        } else {
            Istim = 0.0;
        };
    }
    else
    {
        Istim = 0.0;
    }
    return Istim;
}


/*typedef struct
{
    double ERP;
    double BCL_s1, BCL_s2;
    double ap_amp_ratio;
    double BCL_adjust;

} ERP_measure_struct;
*/


#endif
// /* END OF STIMULUS.C