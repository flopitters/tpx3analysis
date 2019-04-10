/********************************************************************
* File: calibrate.h
* ------------------------
*
* Description:
* Calibrate testpulse and/or source data from tpx3 calibration.
*
* Version:
* Author: Florian Pitters
*
*******************************************************************/

#ifndef __CALIBRATE_H_INCLUDED__
#define __CALIBRATE_H_INCLUDED__

#include "common.h"

void cal_analysis(std::string dev_id, int thr_dac, int ik_dac);

void calibrate_testpulses(std::string dev_id, int thr_dac, int ik_dac);
void calibrate_sources(std::string dev_id, int thr_dac, int ik_dac, std::vector<std::string> srcs);


#endif
