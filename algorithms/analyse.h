/********************************************************************
* File: analyse.h
* ------------------------
*
* Description:
* Analyse data from tpx3 source acquisition.
*
* Version:
* Author: Florian Pitters
*
*******************************************************************/

#ifndef __ANALYSE_H_INCLUDED__
#define __ANALYSE_H_INCLUDED__

#include "common.h"

void calibrate_sources_prepare(std::string dev_id, int thr_dac, int ik_dac, std::string src);
void acq_analysis(std::string dev_id, int thr_dac, int ik_dac, std::string src);


#endif
