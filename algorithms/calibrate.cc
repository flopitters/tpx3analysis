/********************************************************************
* File: calibrate.cc
* ------------------------
*
* Description:
* Calibrate testpulse and/or source data from tpx3 calibration.
*
* Version:
* Author: Florian Pitters
*
*******************************************************************/

#include "calibrate.h"


void calibrate_sources(std::string dev_id, int thr_dac, int ik_dac, std::vector<std::string> srcs)
{
	gROOT->Reset();
	gROOT->SetBatch(kTRUE);

	styleCommon();
	styleMap();
	styleGraph();


	// Preperations
	// ---------------------------------------

	int ret, cnt;
	std::vector<float> ref_ene, ref_ene_err, x, y, x_err, y_err, peak_global, peak_global_err, a_list, b_list;
	std::vector<std::vector<float>> dat_tot, dat_volts, out_volts;
	std::vector<std::vector<std::vector<float>>> dat_src;
	std::string tmp_file, out_file, src;


	// load calibration data
	bool info = true;
	tmp_file = dat_path + "/lab/" + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_tot.txt";
	ret = read_file(tmp_file.c_str(), ' ', dat_tot, info);

	// load pixel by pixel data
	for (int k = 0; k < srcs.size(); k++){
		src = srcs.at(k);

		if (src == "fe"){
			ref_ene.push_back(6.40);
			ref_ene_err.push_back(0.01);
		}
		else if (src == "cu"){
			ref_ene.push_back(8.05);
			ref_ene_err.push_back(0.010);
		}
		else if (src == "in"){
			ref_ene.push_back(24.21);
			ref_ene_err.push_back(0.1);
		}

		tmp_file = dat_path + "/lab/" + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac)
			+ "/" + dev_id + "_" + src + ".txt";
		ret = read_file(tmp_file.c_str(), ' ', dat_volts, info);

		dat_src.push_back(dat_volts);
		dat_volts.clear();
	}



	// Analysis
	// ---------------------------------------

	// global analysis
	TF1 *linef = new TF1("linef", "[0]*x + [1]", 0, 800);
	linef->SetNpx(400);
	linef->SetParameter(0, 0.08);
	linef->SetParLimits(0, 0.01, 0.2);
	linef->SetParameter(1, 1);
	linef->SetParLimits(1, -5, 5);
	linef->SetLineColor(kRed);

   	// matrix analysis
	int fskip = 0;
	int nfailed = 0;
	float a, b, thr;
	std::vector<float> res_volts;
	
	TH1F *hist_a_lin = new TH1F("hist_a_lin", "Slope in keV/mV", 150, 0, 0.15);
	TH1F *hist_b_lin = new TH1F("hist_b_lin", "Offset in keV", 100, 0, 3);
	TH2F *map_a_lin = new TH2F("map_a_lin", "Slope in keV/mV", 256, 0, 256, 256, 0, 256);
	TH2F *map_b_lin = new TH2F("map_b_lin", "Offset in keV", 256, 0, 256, 256, 0, 256);
	TH1F *hist_thr = new TH1F("hist_thr", "Threshold in mV", 50, 0, 50);
	TH2F *map_thr = new TH2F("map_thr", "Threshold in mV", 256, 0, 256, 256, 0, 256);
	TH2F *map_scatter = new TH2F("map_scatter", "Threshold vs Offset", 100, 0, 3, 50, 0, 50);

	TH1F *hist_a_err_lin = new TH1F("hist_a_err_lin", "Slope in keV/mV", 100, 0, 0.1);
	TH1F *hist_b_err_lin = new TH1F("hist_b_err_lin", "Offset in keV", 100, 0, 3);
	TH2F *map_a_err_lin = new TH2F("map_a_err_lin", "Slope in keV/mV", 256, 0, 256, 256, 0, 256);
	TH2F *map_b_err_lin = new TH2F("map_b_err_lin", "Offset in keV", 256, 0, 256, 256, 0, 256);

	TH1F *hist_chi2fe = new TH1F("hist_chi2fe", "Threshold vs Offset", 100, 0, 100);
	TH1F *hist_chi2in = new TH1F("hist_chi2in", "Threshold vs Offset", 100, 0, 100);
	TH1F *hist_mpv_fe = new TH1F("hist_mpv_fe", "Threshold vs Offset", 100, 0, 100);
	TH1F *hist_mpv_in = new TH1F("hist_mpv_in", "Threshold vs Offset", 100, 0, 400);
	TH1F *hist_mpv_err_fe = new TH1F("hist_mpv_err_fe", "Threshold vs Offset", 50, 0, 50);
	TH1F *hist_mpv_err_in = new TH1F("hist_mpv_err_in", "Threshold vs Offset", 50, 0, 50);
	TH2F *map_chi2fe = new TH2F("map_chi2fe", "Threshold vs Offset", 256, 0, 256, 256, 0, 256);
	TH2F *map_chi2in = new TH2F("map_chi2in", "Threshold vs Offset", 256, 0, 256, 256, 0, 256);
	TH2F *map_mpv_fe = new TH2F("map_mpv_fe", "Threshold vs Offset", 256, 0, 256, 256, 0, 256);
	TH2F *map_mpv_in = new TH2F("map_mpv_in", "Threshold vs Offset", 256, 0, 256, 256, 0, 256);
	TH2F *map_mpv_err_fe = new TH2F("map_mpv_err_fe", "Threshold vs Offset", 256, 0, 256, 256, 0, 256);
	TH2F *map_mpv_err_in = new TH2F("map_mpv_err_in", "Threshold vs Offset", 256, 0, 256, 256, 0, 256);
	for (int j = 0; j < 256; j++) {
		for (int i = 0; i < 256; i++) {

			// reset
			a = b = thr = -1;
			x.clear();
			y.clear();
			x_err.clear();
			y_err.clear();
			res_volts.clear();

			res_volts.push_back(i);
			res_volts.push_back(j);

			for (int k = 0; k < srcs.size(); k++){
				x.push_back(dat_src.at(k).at(256*j + i).at(2));
				y.push_back(ref_ene.at(k));
				x_err.push_back(dat_src.at(k).at(256*j + i).at(3));
				y_err.push_back(ref_ene_err.at(k));

				if ( k == 0 ) {
					hist_mpv_fe->Fill(dat_src.at(k).at(256*j + i).at(2)); map_mpv_fe->Fill(i, j, dat_src.at(k).at(256*j + i).at(2));
					hist_mpv_err_fe->Fill(dat_src.at(k).at(256*j + i).at(3)); map_mpv_err_fe->Fill(i, j, dat_src.at(k).at(256*j + i).at(3));
					hist_chi2fe->Fill(dat_src.at(k).at(256*j + i).at(4)); map_chi2fe->Fill(i, j, dat_src.at(k).at(256*j + i).at(4));
				}
				if ( k == 1 ) {
					hist_mpv_in->Fill(dat_src.at(k).at(256*j + i).at(2)); map_mpv_in->Fill(i, j, dat_src.at(k).at(256*j + i).at(2));
					hist_mpv_err_in->Fill(dat_src.at(k).at(256*j + i).at(3)); map_mpv_err_in->Fill(i, j, dat_src.at(k).at(256*j + i).at(3));
					hist_chi2in->Fill(dat_src.at(k).at(256*j + i).at(4)); map_chi2in->Fill(i, j, dat_src.at(k).at(256*j + i).at(4));
				}
			}

			auto gr = new TGraphErrors(srcs.size(), &x[0], &y[0], &x_err[0], &y_err[0]);
			gr->Fit("linef", "q");

			fskip = 0;
			if (x.at(0) == -1 || x.at(1) == -1 || linef->GetParError(1) > 1.5) {
				fskip = 1;
			}

			if (fskip == 0) {
				a = linef->GetParameter(0);
				b = linef->GetParameter(1);

				float params[4] = {dat_tot.at(256*j + i).at(2), dat_tot.at(256*j + i).at(3), dat_tot.at(256*j + i).at(4), dat_tot.at(256*j + i).at(5)};
				thr = tot_to_volts(0, params);

				if (i == xpix && j == ypix){
					out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/example_volt_fit.pdf";
					printGraph(gr, out_file, "Single Pixel Example Linear Fit", "voltage (mV)", "energy (keV)", -1.0, -1.0, 1.3);
				}

				if (0) {
					if (linef->GetParError(1) > 0.8) {
						std::cout << i << "\t" << j << "\t";
						for (int k = 0; k < srcs.size(); k++){
							std::cout << x.at(k) << "\t" << x_err.at(k) << "\t";
						}
						std::cout << linef->GetParError(0) << "\t" << linef->GetParError(1) << "\t";
						std::cout << std::endl;
						out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/example_bad_volt_curve_" + std::to_string(i) + "_" + std::to_string(j) + "_" + dev_id + ".pdf";
						// printGraph(gr, out_file, "example volt", "voltage (mV)", "peak energy (keV)", -1.0, -1.0, 1.3);
					}
				}

				hist_a_lin->Fill(a);
				hist_b_lin->Fill(b);
				map_a_lin->Fill(i, j, a);
				map_b_lin->Fill(i, j, b);
				hist_thr->Fill(thr);
				map_thr->Fill(i, j, thr);
				map_scatter->Fill(b, thr);

				hist_a_err_lin->Fill(linef->GetParError(0));
				hist_b_err_lin->Fill(linef->GetParError(1));
				map_a_err_lin->Fill(i, j, linef->GetParError(0));
				map_b_err_lin->Fill(i, j, linef->GetParError(1));

				res_volts.push_back(a);
				res_volts.push_back(b);
				res_volts.push_back(linef->GetParError(0));
				res_volts.push_back(linef->GetParError(1));
				res_volts.push_back(linef->GetChisquare() / linef->GetNDF());
			}

			else {
				res_volts.push_back(-1);
				res_volts.push_back(-1);
				res_volts.push_back(-1);
				res_volts.push_back(-1);
				res_volts.push_back(-1);
				nfailed += 1;
			}

			out_volts.push_back(res_volts);
		}
	}


	// use average from row if single pixel failed
	for (int i = 0; i < 256; i++) {
		a_list.clear();
		b_list.clear();
		for (int j = 0; j < 256; j++) {
			if (out_volts.at(256*i + j).at(2) !=  -1){
				a_list.push_back(out_volts.at(256*i + j).at(2));
				b_list.push_back(out_volts.at(256*i + j).at(3));
			}
		}
		for (int j = 0; j < 256; j++) {
			if (out_volts.at(256*i + j).at(2) == -1){
				if(!std::isnan(mean(a_list))) out_volts.at(256*i + j).at(2) = mean(a_list);
				if(!std::isnan(mean(b_list))) out_volts.at(256*i + j).at(3) = mean(b_list);
			}
		}
	}




	// Printing
	// --------------------------------

	styleGraph();

	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_hist_a_lin.pdf";
	printHist(hist_a_lin, out_file, "hist_a_lin", "energy (keV)", "# entries", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_hist_b_lin.pdf";
	printHist(hist_b_lin, out_file, "hist_b_lin", "energy (keV)", "# entries", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_hist_thr.pdf";
	printHist(hist_thr, out_file, "hist_thr", "voltage (mV)", "# entries", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_hist_chi2fe.pdf";
	printHist(hist_chi2fe, out_file, "hist_chi2fe", "chi2/ndf (-)", "# entries", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_hist_chi2in.pdf";
	printHist(hist_chi2in, out_file, "hist_chi2in", "chi2/ndf (-)", "# entries", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_hist_mpv_fe.pdf";
	printHist(hist_mpv_fe, out_file, "hist_mpv_fe", "mpv_/ndf (-)", "# entries", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_hist_mpv_in.pdf";
	printHist(hist_mpv_in, out_file, "hist_mpv_in", "mpv_/ndf (-)", "# entries", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_hist_mpv_err_fe.pdf";
	printHist(hist_mpv_err_fe, out_file, "hist_mpv_err_fe", "mpv_err_/ndf (-)", "# entries", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_hist_mpv_err_in.pdf";
	printHist(hist_mpv_err_in, out_file, "hist_mpv_err_in", "mpv_err_/ndf (-)", "# entries", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_hist_a_err_lin.pdf";
	printHist(hist_a_err_lin, out_file, "hist_a_err_lin", "energy (keV)", "# entries", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_hist_b_err_lin.pdf";
	printHist(hist_b_err_lin, out_file, "hist_b_err_lin", "energy (keV)", "# entries", -1.0, -1.0, 1.3);


	styleMap();

	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_map_a_lin.pdf";
	printMap(map_a_lin, out_file, "map_a_lin", "column", "row", "slope (keV/mV)", 0, 0.15, 1.1, 1.3);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_map_b_lin.pdf";
	printMap(map_b_lin, out_file, "map_b_lin", "column", "row", "offset (keV)", 0, 3, 1.1, 1.1);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_map_thr.pdf";
	printMap(map_thr, out_file, "map_thr", "column", "row", "threshold (mV)", 0, 50, 1.1, 1.1);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_map_scatter.pdf";
	printMap(map_scatter, out_file, "map_thr", "offset (keV)", "threshold (mV)", "# entries", -1.0, -1.0, 1.1, 1.1);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_map_chi2fe.pdf";
	printMap(map_chi2fe, out_file, "map_chi2fe", "column", "row", "chi2/ndf (-)", 0, 10, 1.1, 1.1);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_map_chi2in.pdf";
	printMap(map_chi2in, out_file, "map_chi2in", "column", "row", "chi2/ndf (-)", 0, 10, 1.1, 1.1);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_map_mpv_fe.pdf";
	printMap(map_mpv_fe, out_file, "map_mpv_fe", "column", "row", "peak (mV)", 40, 80, 1.1, 1.1);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_map_mpv_in.pdf";
	printMap(map_mpv_in, out_file, "map_mpv_in", "column", "row", "peak (mV)", 200, 400, 1.1, 1.1);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_map_mpv_err_fe.pdf";
	printMap(map_mpv_err_fe, out_file, "map_mpv_err_fe", "column", "row", "peak (mV)", 0, 10, 1.1, 1.1);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_map_mpv_err_in.pdf";
	printMap(map_mpv_err_in, out_file, "map_mpv_err_in", "column", "row", "peak (mV)", 0, 30, 1.1, 1.1);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_map_a_err_lin.pdf";
	printMap(map_a_err_lin, out_file, "map_a_err_lin", "column", "row", "slope (keV/mV)", 0, 0.15, 1.1, 1.3);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_map_b_err_lin.pdf";
	printMap(map_b_err_lin, out_file, "map_b_err_lin", "column", "row", "offset (keV)", 0, 3, 1.1, 1.1);

	// save dat files
	std::string hd_volts = "# column | row | a_volts (eV/mV) | b_volts (eV) | err_a_volts (eV/mV) | err_b_volts (eV) | chi2/ndf";

	out_file = dat_path + "/lab/" + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_volts.txt";
	write_file(out_file, ' ', out_volts, hd_volts, info);

	std::string cmd = "cp -r " + dat_path + "/lab/" + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac)
		+ "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_volts.txt " + cfg_path + dev_id + "/";
	system(cmd.c_str());

	printf("---->  Failed so far (%d)\n", nfailed);
}





// ---------------------------------------
// Analyse Calibration
// ---------------------------------------

void cal_analysis(std::string dev_id, int thr_dac, int ik_dac)
{
	gROOT->Reset();
	gROOT->SetBatch(kTRUE);

	styleCommon();
	styleMap();
	styleGraph();



	// Preperations
	// ---------------------------------------

	// initialise
	int cnt, ret, nsteps;;
	std::vector<int> amps;
	std::vector<float> tmp, volts;
	std::vector<std::vector<float>> table, mat;
	std::vector<std::vector<float>> dat_tot, dat_toa, dat_volts;
	std::vector<std::vector<std::vector<float>>> dat_pc_raw, dat_tot_raw, dat_toa_raw, dat_tot_std_raw, dat_toa_std_raw;
	std::string tmp_file, out_file;

	bool info = false;
	tmp_file = dat_path + "/lab/" + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_tot.txt";
	ret = read_file(tmp_file.c_str(), ' ', dat_tot, info);
	tmp_file = dat_path + "/lab/" + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_toa.txt";
	ret = read_file(tmp_file.c_str(), ' ', dat_toa, info);

	tmp_file = dat_path + "/lab/" + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_volts.txt";
	ret = read_file(tmp_file.c_str(), ' ', dat_volts, info);

	// write calibration data
	int col, row;
	float thr, thr_err, tot, toa;
	float a_tot, b_tot, c_tot, t_tot, c_toa, t_toa, d_toa, chi2_tot, chi2_toa;



	// Input
	// ---------------------------------------

	// read energy table with voltage values
	tmp_file = dat_path + "/lab/" + dev_id + "/cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/Energy_table.dat";
	ret = read_file(tmp_file, '\t', table, info);
	nsteps = table.size();

	// read data files
	for (int i = 0; i < nsteps; i++){
		tmp = table.at(i);
		amps.push_back(tmp.at(0));
		volts.push_back(tmp.at(3));

		tmp_file = dat_path + "/lab/" + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac)
			+ "/PC_" + std::to_string(amps.at(i)) + ".map";
		ret = read_file(tmp_file, ' ', mat, info);
		dat_pc_raw.push_back(mat);

		tmp_file = dat_path + "/lab/" + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac)
			+ "/TOT_" + std::to_string(amps.at(i)) + ".map";
		ret = read_file(tmp_file, ' ', mat, info);
		dat_tot_raw.push_back(mat);

		tmp_file = dat_path + "/lab/" + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac)
			+ "/TOA_" + std::to_string(amps.at(i)) + ".map";
		ret = read_file(tmp_file, ' ', mat, info);
		dat_toa_raw.push_back(mat);

		tmp_file = dat_path + "/lab/" + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac)
			+ "/TOT_VAR_" + std::to_string(amps.at(i)) + ".map";
		ret = read_file(tmp_file, ' ', mat, info);
		dat_tot_std_raw.push_back(mat);

		tmp_file = dat_path + "/lab/" + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac)
			+ "/TOA_VAR_" + std::to_string(amps.at(i)) + ".map";
		ret = read_file(tmp_file, ' ', mat, info);
		dat_toa_std_raw.push_back(mat);
	}




	// Analysis
	// ---------------------------------------

	// Fitted
	auto hist_atot = new TH1F("hist_atot", "a_tot", 100, 0, 0.5);
	auto hist_btot = new TH1F("hist_btot", "b_tot", 100, -50, 50);
	auto hist_ctot = new TH1F("hist_ctot", "c_tot", 150, 0, 1500);
	auto hist_ttot = new TH1F("hist_ttot", "t_tot", 100, -50, 50);
	auto hist_ctoa = new TH1F("hist_ctoa", "c_toa", 100, 0, 1500);
	auto hist_ttoa = new TH1F("hist_ttoa", "t_toa", 100, -50, 50);
	auto hist_dtoa = new TH1F("hist_dtoa", "d_toa", 100, -50, 50);
	auto hist_thr = new TH1F("hist_thr", "thr", 20, 0, 60);

	auto map_atot = new TH2F("map_atot", "a_tot", 256, 0, 256, 256, 0, 256);
	auto map_btot = new TH2F("map_btot", "b_tot", 256, 0, 256, 256, 0, 256);
	auto map_ctot = new TH2F("map_ctot", "c_tot", 256, 0, 256, 256, 0, 256);
	auto map_ttot = new TH2F("map_ttot", "t_tot", 256, 0, 256, 256, 0, 256);
	auto map_ctoa = new TH2F("map_ctoa", "c_toa", 256, 0, 256, 256, 0, 256);
	auto map_ttoa = new TH2F("map_ttoa", "t_toa", 256, 0, 256, 256, 0, 256);
	auto map_dtoa = new TH2F("map_dtoa", "d_toa", 256, 0, 256, 256, 0, 256);
	auto map_thr = new TH2F("map_thr", "thr", 256, 0, 256, 256, 0, 256);
	auto map_surrogate = new TH2F("map_surrogate", "ToT Response Function", 80, 0, 800, 80, 0, 80);
	auto map_timewalk = new TH2F("map_timewalk", "ToA Response Function", 80, 0, 800, 50, 0, 50);
	auto hist_chi2tot = new TH1F("hist_chi2tot", "chi2_tot", 50, 0, 50);
	auto hist_chi2toa = new TH1F("hist_chi2toa", "chi2_toa", 50, 0, 50);
	auto map_chi2tot = new TH2F("map_chi2tot", "chi2_tot", 256, 0, 256, 256, 0, 256);
	auto map_chi2toa = new TH2F("map_chi2toa", "chi2_toa", 256, 0, 256, 256, 0, 256);
	for (int i = 0; i < dat_tot.size(); i++) {
		col = dat_tot.at(i).at(0);
		row = dat_tot.at(i).at(1);
		a_tot = dat_tot.at(i).at(2);
		b_tot = dat_tot.at(i).at(3);
		c_tot = dat_tot.at(i).at(4);
		t_tot = dat_tot.at(i).at(5);
		c_toa = dat_toa.at(i).at(2);
		t_toa = dat_toa.at(i).at(3);
		d_toa = dat_toa.at(i).at(4);

		thr = dat_tot.at(i).at(6);
		thr_err = dat_tot.at(i).at(7);

		chi2_tot = dat_tot.at(i).at(8);
		chi2_toa = dat_toa.at(i).at(5);

		hist_atot->Fill(a_tot);
		hist_btot->Fill(b_tot);
		hist_ctot->Fill(c_tot);
		hist_ttot->Fill(t_tot);
		hist_ctoa->Fill(c_toa);
		hist_ttoa->Fill(t_toa);
		hist_dtoa->Fill(d_toa);
		hist_thr->Fill(thr);
		hist_chi2tot->Fill(chi2_tot);
		hist_chi2toa->Fill(chi2_toa);
		map_atot->Fill(col, row, a_tot);
		map_btot->Fill(col, row, b_tot);
		map_ctot->Fill(col, row, c_tot);
		map_ttot->Fill(col, row, t_tot);
		map_ctoa->Fill(col, row, c_toa);
		map_ttoa->Fill(col, row, t_toa);
		map_dtoa->Fill(col, row, d_toa);
		map_thr->Fill(col, row, thr);
		map_chi2tot->Fill(col, row, chi2_tot);
		map_chi2toa->Fill(col, row, chi2_toa);

		float params_tot[4] = {a_tot, b_tot, c_tot, t_tot};
		float params_toa[3] = {c_toa, t_toa, d_toa};
		for (int a = 0; a < 800; a += 2) {
			tot = volts_to_tot(a, params_tot);
			map_surrogate->Fill(a, tot);
			toa = volts_to_toa(a, params_toa);
			map_timewalk->Fill(a, toa);
		}
	}

	// Raw
	cnt = 0;
	std::vector<float> electrons(nsteps, 0), zeros(nsteps, 0);
	std::vector<float> toa_err_per_pixel(nsteps, 0), tot_err_per_pixel(nsteps, 0), toa_per_pixel(nsteps, 0), tot_per_pixel(nsteps, 0);
	std::vector<float> cnt_toa(nsteps, 0), cnt_tot(nsteps, 0),  toa_rel_pixel(nsteps, 0), tot_rel_pixel(nsteps, 0);
	std::vector<float> toa_err_over_matrix(nsteps, 0), tot_err_over_matrix(nsteps, 0);
	auto map_surrogate_raw = new TH2F("map_surrogate_raw", "ToT Response Function", 80, 0, 800, 160, 0, 80);
	auto map_timewalk_raw = new TH2F("map_timewalk_raw", "ToA Response Function", 80, 0, 800, 100, 0, 50);
	auto map_offset_mean = new TH2F("", "", 256, 0, 256, 256, 0, 256);
	auto map_offset_std = new TH2F("", "", 256, 0, 256, 256, 0, 256);
	for (int j = 0; j < 256; j++) {
		for (int i = 0; i < 256; i++) {
			for (int k = 1; k < nsteps-1; k++) {
				// get all voltage values that recieved all pulses sent
				// skip first one that did so, to reduce noise influence
				if (dat_pc_raw[k][j][i] == npulses && dat_pc_raw[k][j][i] == dat_pc_raw[k-1][j][i]) {
					if (dat_tot_raw[k][j][i] != 0 && dat_tot_raw[k][j][i] != -1) {
						map_surrogate_raw->Fill(volts[k], dat_tot_raw[k][j][i]);
					}
					if (dat_toa_raw[k][j][i] != 0 && dat_toa_raw[k][j][i] != -1) {
						map_timewalk_raw->Fill(volts[k], dat_toa_raw[k][j][i]);
					}
				}
				if (k == nsteps-2){
					map_offset_mean->Fill(i, j, dat_toa_raw[k][j][i]);
					map_offset_std->Fill(i, j, dat_toa_std_raw[k][j][i]);
				}
			}
		}
	}
	for (int j = 0; j < 256; j++) {
		for (int i = 0; i < 256; i++) {
			for (int k = 1; k < nsteps-1; k++) {
				if (dat_pc_raw[k][j][i] == npulses && dat_pc_raw[k][j][i] == dat_pc_raw[k-1][j][i]) {

					if (dat_tot_raw[k][j][i] != 0 && dat_tot_raw[k-1][j][i] != 0 && dat_tot_raw[k][j][i] != -1 && dat_tot_std_raw[k][j][i] != -1 && 		dat_volts.at(256*j+i).at(2) != -1) {
						cnt_tot.at(k) += 1;
						tot_per_pixel.at(k) += dat_tot_raw[k][j][i];
						tot_err_per_pixel.at(k) += sqrt(dat_tot_std_raw[k][j][i]);
					}

					if (dat_toa_raw[k][j][i] != 0 && dat_toa_raw[k][j][i] != -1 && dat_toa_std_raw[k][j][i] != -1 && dat_volts.at(256*j+i).at(2) != -1) {
						cnt_toa.at(k) += 1;
						electrons.at(k) += abs((volts[k]*dat_volts.at(256*j+i).at(2) + dat_volts.at(256*j+i).at(3)) / 3.65 * 1000.);

						toa_per_pixel.at(k) += dat_toa_raw[k][j][i];
						toa_err_per_pixel.at(k) += sqrt(dat_toa_std_raw[k][j][i]);
					}
				}
			}
		}
	}

	for (int k = 0; k < nsteps; k++) {
		electrons.at(k) /= cnt_toa.at(k);
		tot_per_pixel.at(k) /= cnt_tot.at(k);
		tot_err_per_pixel.at(k) /= cnt_tot.at(k);
		toa_per_pixel.at(k) /= cnt_toa.at(k);
		toa_err_per_pixel.at(k) /= cnt_tot.at(k);

		tot_rel_pixel.at(k) = tot_err_per_pixel.at(k) / tot_per_pixel.at(k);
		toa_rel_pixel.at(k) = toa_err_per_pixel.at(k) / toa_per_pixel.at(k);

		std::cout << electrons.at(k) << " " << tot_per_pixel.at(k) << " " << tot_err_per_pixel.at(k)
									 << " " << toa_per_pixel.at(k) << " " << toa_err_per_pixel.at(k)  << std::endl;
	}
	for (int k = 0; k < nsteps; k++) {
		toa_per_pixel.at(k) -= 11; // centre at zero
		if (k < 12){
			toa_per_pixel.at(k) = -2;
		}
	}

	std::vector<float> time_sn_contribution(nsteps, 0), time_adc_contribution(nsteps, 0), time_tot_contribution(nsteps, 0);
	for (int k = 0; k < nsteps; k++) {
		time_sn_contribution.at(k) = 25. / (electrons.at(k) / 90.);
		time_adc_contribution.at(k) = 1.5625/sqrt(12);
		time_tot_contribution.at(k) = sqrt(pow(time_adc_contribution.at(k), 2) + pow(toa_err_per_pixel.at(k), 2));
	}

	out_file = res_path + "dev_time_contributions.pdf";
	TCanvas *canv_time_contributions = new TCanvas("hist_time_contributions", "hist");
	TMultiGraph *mg_time_contributions = new TMultiGraph();
	TLegend *lg_time_contributions = new TLegend(0.45, 0.745, 0.88, 0.88);
	addGraph(mg_time_contributions, lg_time_contributions, electrons.size(), &(electrons[0]), &(time_tot_contribution[0]),
		colours.at(2), markers.at(3), 0.9, styles.at(0), "p", " root sum of the squares");
	addGraph(mg_time_contributions, lg_time_contributions, electrons.size(), &(electrons[0]), &(toa_err_per_pixel[0]),
		colours.at(0), markers.at(0), 0.9, styles.at(0), "p", " jitter contribution");
	addGraph(mg_time_contributions, lg_time_contributions, electrons.size(), &(electrons[0]), &(time_adc_contribution[0]),
		colours.at(1), markers.at(2), 0.9, styles.at(0), "p", " clock contribution");
	// gPad->SetLogy();
	mg_time_contributions->Draw("AP");
	mg_time_contributions->GetXaxis()->SetTitle("pixel charge [e-]");
	mg_time_contributions->GetYaxis()->SetTitle("uncertainty contribution [ns]");
	mg_time_contributions->GetYaxis()->SetTitleOffset(1.1);
	mg_time_contributions->GetXaxis()->SetLimits(0, 12000);
	mg_time_contributions->SetMinimum(0);
	mg_time_contributions->SetMaximum(2.5);
	lg_time_contributions->Draw();
	// canv_time_contributions->SetLogy();
	canv_time_contributions->SaveAs(out_file.c_str());

	out_file = res_path + "dev_tot_seperates.pdf";
	TCanvas *canv_tot_seperates = new TCanvas("hist_tot_seperates", "hist");
	TMultiGraph *mg_tot_seperates = new TMultiGraph();
	TLegend *lg_tot_seperates = new TLegend(0.64, 0.79, 0.88, 0.88);
	addGraph(mg_tot_seperates, lg_tot_seperates, electrons.size(), &(electrons[0]), &(tot_per_pixel[0]),
		colours.at(1), markers.at(1), 0.9, styles.at(0), "p", " mean");
	addGraph(mg_tot_seperates, lg_tot_seperates, electrons.size(), &(electrons[0]), &(tot_err_per_pixel[0]),
		colours.at(0), markers.at(0), 0.9, styles.at(0), "p", " rms");
	mg_tot_seperates->Draw("AP");
	mg_tot_seperates->SetTitle("Title");
	mg_tot_seperates->GetXaxis()->SetTitle("charge [e-]");
	mg_tot_seperates->GetYaxis()->SetTitle("time over threshold [LSB]");
	mg_tot_seperates->GetYaxis()->SetTitleOffset(1.1);
	mg_tot_seperates->GetXaxis()->SetLimits(0, 10000);
	mg_tot_seperates->SetMinimum(0);
	mg_tot_seperates->SetMaximum(80);
	lg_tot_seperates->Draw();
	canv_tot_seperates->SaveAs(out_file.c_str());

	out_file = res_path + "dev_tot_relative.pdf";
	TCanvas *canv_tot_relative = new TCanvas("hist_tot_relative", "hist");
	TMultiGraph *mg_tot_relative = new TMultiGraph();
	TLegend *lg_tot_relative = new TLegend(0.64, 0.79, 0.88, 0.88);
	addGraph(mg_tot_relative, lg_tot_relative, electrons.size(), &(electrons[0]), &(tot_rel_pixel[0]),
		colours.at(0), markers.at(0), 0.9, styles.at(0), "p", " ");
	mg_tot_relative->Draw("AP");
	mg_tot_relative->SetTitle("Title");
	mg_tot_relative->GetXaxis()->SetTitle("charge [e-]");
	mg_tot_relative->GetYaxis()->SetTitle("#Delta ToT/ToT");
	mg_tot_relative->GetYaxis()->SetTitleOffset(1.1);
	mg_tot_relative->GetXaxis()->SetLimits(0, 10000);
	mg_tot_relative->SetMinimum(0);
	mg_tot_relative->SetMaximum(0.25);
	canv_tot_relative->SaveAs(out_file.c_str());

	out_file = res_path + "dev_tot_errorband.pdf";
	TCanvas *canv_tot_errorband = new TCanvas("hist_tot_errorband", "hist");
	TMultiGraph *mg_tot_errorband = new TMultiGraph();
	TLegend *lg_tot_errorband = new TLegend(0.64, 0.79, 0.88, 0.88);
	addErrGraph(mg_tot_errorband, lg_tot_errorband, electrons.size(), &(electrons[0]), &(tot_per_pixel[0]), &(zeros[0]), &(tot_err_per_pixel[0]),
		colours.at(0), markers.at(0), 0.9, styles.at(0), "p", " timewalk");
	mg_tot_errorband->Draw("AP4");
	mg_tot_errorband->SetTitle("Title");
	mg_tot_errorband->GetXaxis()->SetTitle("charge [e-]");
	mg_tot_errorband->GetYaxis()->SetTitle("time over threshold [LSB]");
	mg_tot_errorband->GetYaxis()->SetTitleOffset(1.1);
	mg_tot_errorband->GetXaxis()->SetLimits(0, 10000);
	mg_tot_errorband->SetMinimum(0);
	mg_tot_errorband->SetMaximum(80);
	canv_tot_errorband->SaveAs(out_file.c_str());

	out_file = res_path + "dev_toa_seperate.pdf";
	TCanvas *canv_toa_seperate = new TCanvas("hist_toa_seperate", "hist");
	TMultiGraph *mg_toa_seperate = new TMultiGraph();
	TLegend *lg_toa_seperate = new TLegend(0.60, 0.78, 0.88, 0.88);
	// addGraph(mg_toa_seperate, lg_toa_seperate, electrons.size(), &(electrons[0]), &(toa_err_per_pixel[0]),
	// 	colours.at(1), markers.at(1), 0.9, styles.at(0), "p", " constant");
	addGraph(mg_toa_seperate, lg_toa_seperate, electrons.size(), &(electrons[0]), &(toa_per_pixel[0]),
		colours.at(1), markers.at(1), 0.9, styles.at(0), "p", " timewalk mean");
	addGraph(mg_toa_seperate, lg_toa_seperate, electrons.size(), &(electrons[0]), &(toa_err_per_pixel[0]),
		colours.at(0), markers.at(0), 0.9, styles.at(0), "p", " timewalk rms");
	mg_toa_seperate->Draw("AP");
	mg_toa_seperate->SetTitle("Title");
	mg_toa_seperate->GetXaxis()->SetTitle("injected charge [e-]");
	mg_toa_seperate->GetYaxis()->SetTitle("time [ns]");
	mg_toa_seperate->GetYaxis()->SetTitleOffset(1.1);
	mg_toa_seperate->GetXaxis()->SetLimits(0, 12000);
	mg_toa_seperate->SetMinimum(-1);
	mg_toa_seperate->SetMaximum(30);
	lg_toa_seperate->Draw();
	// canv_toa_seperate->SetLogy();
	canv_toa_seperate->SaveAs(out_file.c_str());

	out_file = res_path + "dev_toa_relative.pdf";
	TCanvas *canv_toa_relative = new TCanvas("hist_toa_relative", "hist");
	TMultiGraph *mg_toa_relative = new TMultiGraph();
	TLegend *lg_toa_relative = new TLegend(0.64, 0.79, 0.88, 0.88);
	addGraph(mg_toa_relative, lg_toa_relative, electrons.size(), &(electrons[0]), &(toa_rel_pixel[0]),
		colours.at(0), markers.at(0), 0.9, styles.at(0), "p", " timewalk");
	mg_toa_relative->Draw("AP4");
	mg_toa_relative->SetTitle("Title");
	mg_toa_relative->GetXaxis()->SetTitle("charge [e-]");
	mg_toa_relative->GetYaxis()->SetTitle("time [ns]");
	mg_toa_relative->GetYaxis()->SetTitleOffset(1.1);
	mg_toa_relative->GetXaxis()->SetLimits(0, 10000);
	mg_toa_relative->SetMinimum(0);
	mg_toa_relative->SetMaximum(0.25);
	canv_toa_relative->SaveAs(out_file.c_str());

	for (int i=0; i<electrons.size(); i++){
		std::cout << electrons.at(i) << " " << toa_err_per_pixel.at(i) << std::endl;
	}

	out_file = res_path + "dev_toa_errorband.pdf";
	TCanvas *canv_toa_errorband = new TCanvas("hist_toa_errorband", "hist");
	TMultiGraph *mg_toa_errorband = new TMultiGraph();
	TLegend *lg_toa_errorband = new TLegend(0.64, 0.79, 0.88, 0.88);
	addErrGraph(mg_toa_errorband, lg_toa_errorband, electrons.size(), &(electrons[0]), &(toa_per_pixel[0]), &(zeros[0]), &(toa_err_per_pixel[0]),
		colours.at(0), markers.at(0), 0.9, styles.at(0), "p", " timewalk");
	mg_toa_errorband->Draw("AP4");
	mg_toa_errorband->SetTitle("Title");
	mg_toa_errorband->GetXaxis()->SetTitle("charge [e-]");
	mg_toa_errorband->GetYaxis()->SetTitle("time [ns]");
	mg_toa_errorband->GetYaxis()->SetTitleOffset(1.1);
	mg_toa_errorband->GetXaxis()->SetLimits(0, 10000);
	mg_toa_errorband->SetMinimum(-1);
	mg_toa_errorband->SetMaximum(20);
	canv_toa_errorband->SaveAs(out_file.c_str());



	// Printing
	// ---------------------------------------

	styleMap();

	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_surrogate_raw.pdf";
	printMap(map_surrogate_raw, out_file, "map_surrogate raw", "voltage (mV)", "time over threshold (ADC)", "# entries", -10, 4E4, 1.1, 1.5, "COLZ");
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_timewalk_raw.pdf";
	printMap(map_timewalk_raw, out_file, "map_timewalk_raw", "voltage (mV)", "time of arrival (ns)", "# entries", -10, 3E4,  1.1, 1.5, "COLZ");

	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_offset_mean.pdf";
	printMap(map_offset_mean, out_file, "map_offset_mean", "column", "row", "offset [ns]", 0, 20, 1.1, 1.2, "COLZ");
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_offset_std.pdf";
	printMap(map_offset_std, out_file, "map_offset_std", "column", "row", "offset rms [ns]", -0.01, 1, 1.1, 1.2, "COLZ");

	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_surrogate.pdf";
	printMap(map_surrogate, out_file, "map_surrogate", "voltage (mV)", "time over threshold (ADC)", "# entries", -10, 1.1E5, 1.1, 1.2, "COLZ");
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_timewalk.pdf";
	printMap(map_timewalk, out_file, "map_timewalk", "voltage (mV)", "time of arrival (ns)", "# entries", -10, 1.6E5, 1.1, 1.2, "COLZ");

	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_map_tot_a.pdf";
	printMap(map_atot, out_file, "map_tot_a", "column", "row", "a_tot paramater (ADC/mV)", 0.05, 0.15, 1.1, 1.2);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_map_tot_b.pdf";
	printMap(map_btot, out_file, "map_tot_b", "column", "row", "b_tot paramater (ADC)", -20, 20, 1.1, 1.2);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_map_tot_c.pdf";
	printMap(map_ctot, out_file, "map_tot_c", "column", "row", "c_tot paramater (ADC*mV)", 0, 200, 1.1, 1.2);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_map_tot_t.pdf";
	printMap(map_ttot, out_file, "map_tot_t", "column", "row", "t_tot paramater (mV)", 0, 40, 1.1, 1.2);

	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_map_toa_c.pdf";
	printMap(map_ctoa, out_file, "map_toa_c", "column", "row", "c_toa paramater (ns*mV)", 0, 800, 1.1, 1.2);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_map_toa_t.pdf";
	printMap(map_ttoa, out_file, "map_toa_t", "column", "row", "t_toa paramater (mV)", -20, 20, 1.1, 1.2);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_map_toa_d.pdf";
	printMap(map_dtoa, out_file, "map_toa_d", "column", "row", "d_toa paramater (ns)", 0, 20, 1.1, 1.2);

	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_map_thr.pdf";
	printMap(map_thr, out_file, "map_thr", "column", "row", "estimated thr (mV)", 0, 50, 1.1, 1.0);

	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_map_chi2_tot.pdf";
	printMap(map_chi2tot, out_file, "map_tot_chi2", "column", "row", "chi2/ndf", 0, 10, 1.1, 1.3);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_map_chi2_toa.pdf";
	printMap(map_chi2toa, out_file, "map_toa_chi2", "column", "row", "chi2/ndf", 0, 10, 1.1, 1.3);

	styleGraph();

	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_hist_tot_a.pdf";
	printHist(hist_atot, out_file, "hist_tot_a", "a_tot paramater (ADC/mV)", "# entries", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_hist_tot_b.pdf";
	printHist(hist_btot, out_file, "hist_tot_b", "b_tot paramater (ADC)", "# entries", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_hist_tot_c.pdf";
	printHist(hist_ctot, out_file, "hist_tot_c", "c_tot paramater (ADC*mV)", "# entries", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_hist_tot_t.pdf";
	printHist(hist_ttot, out_file, "hist_tot_t", "t_tot paramater (mV)", "# entries", -1.0, -1.0, 1.3);

	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_hist_toa_c.pdf";
	printHist(hist_ctoa, out_file, "hist_toa_c", "c_toa paramater (ns*mV)", "# entries", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_hist_toa_t.pdf";
	printHist(hist_ttoa, out_file, "hist_toa_t", "t_toa paramater (mV)", "# entries", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_hist_toa_d.pdf";
	printHist(hist_dtoa, out_file, "hist_toa_d", "d_toa paramater (ns)", "# entries", -1.0, -1.0, 1.3);

	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_hist_thr.pdf";
	printHist(hist_thr, out_file, "hist_thr", "estimated thr (mV)", "# entries", -1.0, -1.0, 1.3);

	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_hist_chi2_tot.pdf";
	printHist(hist_chi2tot, out_file, "hist_chi2_tot", "chi2/ndf", "# entries", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_hist_chi2_toa.pdf";
	printHist(hist_chi2toa, out_file, "hist_chi2_toa", "chi2/ndf", "# entries", -1.0, -1.0, 1.3);

}




// ---------------------------------------
// Fit Calibration Data
// ---------------------------------------

void calibrate_testpulses(std::string dev_id, int thr_dac, int ik_dac)
{
	gROOT->Reset();
	gROOT->SetBatch(kTRUE);

	styleCommon();
	styleMap();
	styleGraph();



	// Preperations
	// ---------------------------------------

	// initialise
	int ret, nsteps;
	float thr, thr_err, thr_min, thr_max;
	std::string tmp_file, out_file;

	// write calibration data
	Int_t col, row;
	Float_t tot, toa;
	Float_t a_tot, b_tot, c_tot, t_tot, c_toa, t_toa, d_toa, chi2_tot, chi2_toa;

	std::vector<int> amps;
	std::vector<float> tmp, volts;
	std::vector<float> codes_thr, codes_tot, codes_tot_lin, codes_toa, codes_thr_err, codes_tot_err, codes_toa_err;
	std::vector<float> vals_pc, vals_thr, vals_tot, vals_tot_lin, vals_toa, vals_thr_std, vals_tot_std, vals_toa_std;
	std::vector<std::vector<float>> table, mat;
	std::vector<std::vector<std::vector<float>>> dat_pc, dat_tot, dat_tot_std, dat_toa, dat_toa_std;



	// Input
	// ---------------------------------------

	// read energy table with voltage values
	bool info = 0;
	tmp_file = dat_path + "/lab/" + dev_id + "/cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/Energy_table.dat";
	ret = read_file(tmp_file, '\t', table, info);
	nsteps = table.size();

	// read data files
	for (Int_t i = 0; i < nsteps; i++){
		tmp = table.at(i);
		amps.push_back(tmp.at(0));
		volts.push_back(tmp.at(3));

		tmp_file = dat_path + "/lab/" + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac)
			+ "/PC_" + std::to_string(amps.at(i)) + ".map";
		ret = read_file(tmp_file, ' ', mat, info);
		dat_pc.push_back(mat);

		tmp_file = dat_path + "/lab/" + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac)
			+ "/TOT_" + std::to_string(amps.at(i)) + ".map";
		ret = read_file(tmp_file, ' ', mat, info);
		dat_tot.push_back(mat);

		tmp_file = dat_path + "/lab/" + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac)
			+ "/TOA_" + std::to_string(amps.at(i)) + ".map";
		ret = read_file(tmp_file, ' ', mat, info);
		dat_toa.push_back(mat);

		tmp_file = dat_path + "/lab/" + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac)
			+ "/TOT_VAR_" + std::to_string(amps.at(i)) + ".map";
		ret = read_file(tmp_file, ' ', mat, info);
		dat_tot_std.push_back(mat);

		tmp_file = dat_path + "/lab/" + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac)
			+ "/TOA_VAR_" + std::to_string(amps.at(i)) + ".map";
		ret = read_file(tmp_file, ' ', mat, info);
		dat_toa_std.push_back(mat);

		if (0) {
			std::cout << amps.at(i) << "\t" << volts.at(i) << "\t" << dat_pc[i][0][0] << std::endl;
		}
	}

	if (0){
		// tot_hd = tot + ftoa - 0.5
		for (Int_t k = 0; k < dat_tot.size(); k++) {
			for (int i = 0; i < dat_tot.at(0).size(); i++){
				for (int j = 0; j < dat_tot.at(0).at(0).size(); j++){
					dat_tot.at(k).at(i).at(j) = fmod(dat_tot.at(k).at(i).at(j) + 0.5, 1);
				}
			}
		}
	}



	// Fitting
	// ---------------------------------------

	auto *linef = new TF1("linef", "[0]*x + [1]", 0, 1000);
	linef->SetNpx(500);

	auto *totf = new TF1("totf", "[0]*x + [1] - [2]/(x-[3])", 0, 1000);
	totf->SetParameter(0, 0.1);
	totf->SetParameter(1, 10);
	totf->SetParameter(2, 100);
	totf->SetParameter(3, 30);
	totf->SetNpx(500);

	auto *toaf = new TF1("toaf", "[0]/(x-[1]) + [2]", 0, 1000);
	toaf->SetParameter(0, 800);
	toaf->SetParameter(1, 10);
	toaf->SetParameter(2, 10);
	toaf->SetNpx(500);

	auto errorf = new TF1("errorf", "[0] * 0.5 * (1 + erf((x-[1]) / (1.414*abs([2]))))", 0, 80);
	errorf->SetParameter(0, 100);
	errorf->SetParameter(1, 25);
	errorf->SetParameter(2, 1);
	errorf->SetParLimits(0, 95, 105);
	errorf->SetParLimits(1, 10, 40);
	errorf->SetParLimits(2, 0.5, 2);
	errorf->SetNpx(500);


	int nfail_tot = 0;
	int nfail_toa = 0;

	// timer
	double timer = 0;
	std::clock_t t0 = std::clock();

	std::vector<float> res_tot, res_toa;
	std::vector<std::vector<float>> out_tot, out_toa;
	for (int j = 0; j < 256; j++) {
		for (int i = 0; i < 256; i++) {

			// clear vectors
			res_tot.clear();
			res_toa.clear();
			codes_thr.clear();
			codes_tot.clear();
			codes_tot_lin.clear();
			codes_toa.clear();
			codes_tot_err.clear();
			codes_toa_err.clear();
			vals_pc.clear();
			vals_thr.clear();
			vals_tot.clear();
			vals_tot_lin.clear();
			vals_toa.clear();
			vals_tot_std.clear();
			vals_toa_std.clear();

			// get next pixel data
			thr_min = 0;
			thr_max = 0;
			vals_pc.push_back(0);
			for (Int_t k = 1; k < nsteps -1; k++) {
				vals_pc.push_back(dat_pc[k][j][i]);

				// determine scurve region
				if (dat_pc[k][j][i] > 0 && dat_pc[k-1][j][i] == 0 && dat_pc[k+1][j][i] > 0 && thr_min == 0) {
					thr_min = volts[k];
				}
				if (dat_pc[k][j][i] == npulses && dat_pc[k-1][j][i] < npulses && dat_pc[k+1][j][i] == npulses && thr_max == 0) {
					thr_max = volts[k-1];
				}

				// add threshold point from scurve region
				// set strict t_tot limits around the threshold
				if (thr_min != 0 && thr_max != 0) {
					codes_tot.push_back(thr_min + (thr_max-thr_min)/2.);
					codes_tot_err.push_back((thr_max-thr_min)/2.);
					vals_tot.push_back(0);
					vals_tot_std.push_back(0);
					totf->SetParameter(3, (thr_min+thr_max-15)/2.);
					totf->SetParLimits(3, thr_min-15, thr_max);
				}

				// get all voltage values that recieved all pulses sent
				// skip first one that did so, to reduce noise influence
				if (dat_pc[k][j][i] == npulses && dat_pc[k][j][i] == dat_pc[k-1][j][i]) {
					if (dat_tot[k][j][i] != 0 && dat_tot[k][j][i] != -1) {
						codes_tot.push_back(volts[k]);
						codes_tot_err.push_back(0);
						vals_tot.push_back(dat_tot[k][j][i]);
						vals_tot_std.push_back(sqrt(dat_tot_std[k][j][i] / float(npulses - 1) + 0.05));
					}
					if (dat_toa[k][j][i] != 0 && dat_toa[k][j][i] != -1) {
						codes_toa.push_back(volts[k]);
						codes_toa_err.push_back(0);
						vals_toa.push_back(dat_toa[k][j][i]);
						vals_toa_std.push_back(sqrt(dat_toa_std[k][j][i] / float(npulses - 1) + 0.05));
					}
				}
			}

			res_tot.push_back(i);
			res_tot.push_back(j);
			res_toa.push_back(i);
			res_toa.push_back(j);

			// Working pixels
			if (codes_tot.size() > 25) {

				// pc data
				auto *graph = new TGraph(volts.size(), &(volts[0]), &(vals_pc[0]));

				// thr data
				for (Int_t k = 1; k < volts.size()-5; k++) {
					if (volts[k] < thr_min) {
						codes_thr.push_back(volts[k]);
						vals_thr.push_back(0);
					}
					else if (volts[k] >= thr_min && volts[k] <= thr_max){
						codes_thr.push_back(volts[k]);
						vals_thr.push_back(vals_pc[k]);
					}
					else if (volts[k-5] < thr_max && volts[k] > thr_max){
						codes_thr.push_back(volts[k]);
						vals_thr.push_back(100);
					}
					codes_thr_err.push_back(2);
					vals_thr_std.push_back(0);
				}
				auto *graph_thr = new TGraphErrors(codes_thr.size(), &(codes_thr[0]), &(vals_thr[0]), &(codes_thr_err[0]), &(vals_thr_std[0]));
				errorf->SetParameter(1, (thr_max+thr_min)/2);
				errorf->SetParLimits(1, (thr_max+thr_min)/2 - 3, (thr_max+thr_min)/2 + 3);
				graph_thr->Fit("errorf", "QM");

				// set tot limits for a and b according to linear region
				for (Int_t k = 0; k < codes_tot.size(); k++) {
					if (codes_tot[k] >= 300){
						codes_tot_lin.push_back(codes_tot[k]);
						vals_tot_lin.push_back(vals_tot[k]);
					}
				}
				auto *graph_tot_lin = new TGraph(codes_tot_lin.size(), &(codes_tot_lin[0]), &(vals_tot_lin[0]));
				graph_tot_lin->Fit("linef", "Q");
				totf->SetParameter(0, linef->GetParameter(0));
				totf->SetParameter(1, linef->GetParameter(1));
				totf->SetParLimits(0, linef->GetParameter(0)*0.9, linef->GetParameter(0)*1.1);
				totf->SetParLimits(1, linef->GetParameter(1)*0.9, linef->GetParameter(1)*1.1);

				// set toa limits
				toaf->SetParameters(0, 600);
				toaf->SetParameters(1, (-20 + codes_toa.at(0))/2.);
				toaf->SetParameters(2, (vals_toa.at(vals_toa.size()-1) + vals_toa.at(vals_toa.size()-2))/4.);
				toaf->SetParLimits(0, 0, 1300);
				toaf->SetParLimits(1, -20, codes_toa.at(0));
				toaf->SetParLimits(2, 0, (vals_toa.at(vals_toa.size()-1) + vals_toa.at(vals_toa.size()-2))/2.);

				// fit the full range tot data
				auto *graph_tot = new TGraphErrors(codes_tot.size(), &(codes_tot[0]), &(vals_tot[0]), &(codes_tot_err[0]), &(vals_tot_std[0]));
				graph_tot->Fit("totf", "QM");

				// fit full range toa data
				auto *graph_toa = new TGraphErrors(codes_toa.size(), &(codes_toa[0]), &(vals_toa[0]), &(codes_toa_err[0]), &(vals_toa_std[0]));
				graph_toa->Fit("toaf", "QM");

				if (totf->GetChisquare()/totf->GetNDF() < 50){
					res_tot.push_back(totf->GetParameter(0));
					res_tot.push_back(totf->GetParameter(1));
					res_tot.push_back(totf->GetParameter(2));
					res_tot.push_back(totf->GetParameter(3));
					res_tot.push_back(errorf->GetParameter(1)); // codes_tot.at(0) or errorf->GetParameter(1)
					res_tot.push_back(errorf->GetParError(1)); // codes_tot_err.at(0) or errorf->GetParError(1)
				}
				else {
					res_tot.push_back(-1);
					res_tot.push_back(-1);
					res_tot.push_back(-1);
					res_tot.push_back(-1);
					res_tot.push_back(-1);
					res_tot.push_back(-1);
					nfail_tot++;
				}
				res_tot.push_back(totf->GetChisquare()/totf->GetNDF());

				if (toaf->GetChisquare()/toaf->GetNDF() < 50){
					res_toa.push_back(toaf->GetParameter(0));
					res_toa.push_back(toaf->GetParameter(1));
					res_toa.push_back(toaf->GetParameter(2));
				}
				else {
					res_toa.push_back(-1);
					res_toa.push_back(-1);
					res_toa.push_back(-1);
					nfail_toa++;
				}
				res_toa.push_back(toaf->GetChisquare()/toaf->GetNDF());

				if ((i == 2 && j == 2) || (i == 1 && j == 1) || (i == 13 && j == 13)) {
					out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/example_a_tot_curve_" + std::to_string(i) + "_" + std::to_string(j) + "_" + dev_id + ".pdf";
					printGraph(graph_tot, out_file, "example tot", "voltage (mV)", "tot", -1.0, -1.0, 1.3);
					out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/example_a_toa_curve_" + std::to_string(i) + "_" + std::to_string(j) + "_" + dev_id + ".pdf";
					printGraph(graph_toa, out_file, "example toa", "voltage (mV)", "time (ns)", -1.0, -1.0, 1.3);
					out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/example_a_thr_curve_" + std::to_string(i) + "_" + std::to_string(j) + "_" + dev_id + ".pdf";
					printGraph(graph_thr, out_file, "example thr", "voltage (mV)", "counts (-)", 0, 150, 1.3);
				}

				if (1) {
					if (totf->GetChisquare()/totf->GetNDF() > 50) {
						out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/example_bad_tot_curve_" + std::to_string(i) + "_" + std::to_string(j) + "_" + dev_id + ".pdf";
						printGraph(graph_tot, out_file, "example tot", "voltage (mV)", "tot", -1.0, -1.0, 1.3);
					}
					if (toaf->GetChisquare()/toaf->GetNDF() > 50) {
						out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/example_bad_toa_curve_" + std::to_string(i) + "_" + std::to_string(j) + "_" + dev_id + ".pdf";
						printGraph(graph_toa, out_file, "example toa", "voltage (mV)", "time (ns)", -1.0, -1.0, 1.3);
					}
				}
				else {
					if (totf->GetChisquare()/totf->GetNDF() > 50)
						std::cout << "Failed TOT fit for " << i << " " << j << std::endl;
					if (toaf->GetChisquare()/toaf->GetNDF() > 50)
						std::cout << "Failed TOA fit for " << i << " " << j << std::endl;
				}
			}

			// Masked pixels
			else {
				res_tot.push_back(-1);
				res_tot.push_back(-1);
				res_tot.push_back(-1);
				res_tot.push_back(-1);
				res_tot.push_back(-1);
				res_tot.push_back(-1);
				res_tot.push_back(-1);

				res_toa.push_back(-1);
				res_toa.push_back(-1);
				res_toa.push_back(-1);
				res_toa.push_back(-1);
			}

			// push back results
			out_tot.push_back(res_tot);
			out_toa.push_back(res_toa);

			// print info
			timer = (std::clock() - t0) / (double)CLOCKS_PER_SEC;
			std::cout << std::fixed << std::setprecision(2);
			if (i == 0) {
				std::cout << " ---> Fitting row " << j
						  << "\tFailed so far " << nfail_tot << " " << nfail_toa
						  << "\tTime since start: " << timer << " sec " <<std::endl;
			}
		}
	}



	// Analysis
	// ---------------------------------------

	// Fitted
	auto hist_atot = new TH1F("hist_atot", "a_tot", 100, 0, 0.5);
	auto hist_btot = new TH1F("hist_btot", "b_tot", 100, -50, 50);
	auto hist_ctot = new TH1F("hist_ctot", "c_tot", 150, 0, 1500);
	auto hist_ttot = new TH1F("hist_ttot", "t_tot", 100, -50, 50);
	auto hist_ctoa = new TH1F("hist_ctoa", "c_toa", 100, 0, 1500);
	auto hist_ttoa = new TH1F("hist_ttoa", "t_toa", 100, -50, 50);
	auto hist_dtoa = new TH1F("hist_dtoa", "d_toa", 100, -50, 50);
	auto hist_thr = new TH1F("hist_thr", "thr", 20, 0, 60);
	auto map_atot = new TH2F("map_atot", "a_tot", 256, 0, 256, 256, 0, 256);
	auto map_btot = new TH2F("map_btot", "b_tot", 256, 0, 256, 256, 0, 256);
	auto map_ctot = new TH2F("map_ctot", "c_tot", 256, 0, 256, 256, 0, 256);
	auto map_ttot = new TH2F("map_ttot", "t_tot", 256, 0, 256, 256, 0, 256);
	auto map_ctoa = new TH2F("map_ctoa", "c_toa", 256, 0, 256, 256, 0, 256);
	auto map_ttoa = new TH2F("map_ttoa", "t_toa", 256, 0, 256, 256, 0, 256);
	auto map_dtoa = new TH2F("map_dtoa", "d_toa", 256, 0, 256, 256, 0, 256);
	auto map_thr = new TH2F("map_thr", "thr", 256, 0, 256, 256, 0, 256);
	auto map_surrogate = new TH2F("map_surrogate", "ToT Response Function", 80, 0, 800, 80, 0, 80);
	auto map_timewalk = new TH2F("map_timewalk", "ToA Response Function", 80, 0, 800, 50, 0, 50);
	auto hist_chi2tot = new TH1F("hist_chi2tot", "chi2_tot", 50, 0, 50);
	auto hist_chi2toa = new TH1F("hist_chi2toa", "chi2_toa", 50, 0, 50);
	auto map_chi2tot = new TH2F("map_chi2tot", "chi2_tot", 256, 0, 256, 256, 0, 256);
	auto map_chi2toa = new TH2F("map_chi2toa", "chi2_toa", 256, 0, 256, 256, 0, 256);
	for (int i = 0; i < out_tot.size(); i++) {
		col = out_tot.at(i).at(0);
		row = out_tot.at(i).at(1);
		a_tot = out_tot.at(i).at(2);
		b_tot = out_tot.at(i).at(3);
		c_tot = out_tot.at(i).at(4);
		t_tot = out_tot.at(i).at(5);
		c_toa = out_toa.at(i).at(2);
		t_toa = out_toa.at(i).at(3);
		d_toa = out_toa.at(i).at(4);

		thr = out_tot.at(i).at(6);
		thr_err = out_tot.at(i).at(7);

		chi2_tot = out_tot.at(i).at(8);
		chi2_toa = out_toa.at(i).at(5);

		hist_atot->Fill(a_tot);
		hist_btot->Fill(b_tot);
		hist_ctot->Fill(c_tot);
		hist_ttot->Fill(t_tot);
		hist_ctoa->Fill(c_toa);
		hist_ttoa->Fill(t_toa);
		hist_dtoa->Fill(d_toa);
		hist_thr->Fill(thr);
		hist_chi2tot->Fill(chi2_tot);
		hist_chi2toa->Fill(chi2_toa);
		map_atot->Fill(col, row, a_tot);
		map_btot->Fill(col, row, b_tot);
		map_ctot->Fill(col, row, c_tot);
		map_ttot->Fill(col, row, t_tot);
		map_ctoa->Fill(col, row, c_toa);
		map_ttoa->Fill(col, row, t_toa);
		map_dtoa->Fill(col, row, d_toa);
		map_thr->Fill(col, row, thr);
		map_chi2tot->Fill(col, row, chi2_tot);
		map_chi2toa->Fill(col, row, chi2_toa);

		float params_tot[4] = {a_tot, b_tot, c_tot, t_tot};
		float params_toa[3] = {c_toa, t_toa, d_toa};
		for (int a = 0; a < 800; a += 2) {
			tot = volts_to_tot(a, params_tot);
			map_surrogate->Fill(a, tot);
			toa = volts_to_toa(a, params_toa);
			map_timewalk->Fill(a, toa);
		}
	}

	// Raw
	TH2F *map_surrogate_raw = new TH2F("map_surrogate_raw", "ToT Response Function", 80, 0, 800, 160, 0, 80);
	TH2F *map_timewalk_raw = new TH2F("map_timewalk_raw", "ToA Response Function", 80, 0, 800, 100, 0, 50);
	for (Int_t j = 0; j < 256; j++) {
		for (Int_t i = 0; i < 256; i++) {
			for (Int_t k = 1; k < nsteps-1; k++) {

				// get all voltage values that recieved all pulses sent
				// skip first one that did so, to reduce noise influence
				if (dat_pc[k][j][i] == npulses && dat_pc[k][j][i] == dat_pc[k-1][j][i]) {
					if (dat_tot[k][j][i] != 0 && dat_tot[k][j][i] != -1) {
						map_surrogate_raw->Fill(volts[k], dat_tot[k][j][i]);
					}
					if (dat_toa[k][j][i] != 0 && dat_toa[k][j][i] != -1) {
						map_timewalk_raw->Fill(volts[k], dat_toa[k][j][i]);
					}
				}
			}
		}
	}



	// Printing
	// ---------------------------------------

	styleMap();

	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_surrogate_raw.pdf";
	printMap(map_surrogate_raw, out_file, "map_surrogate raw", "voltage (mV)", "time over threshold (ADC)", "# entries", -1.0, -1.0, 1.1, 1.3);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_timewalk_raw.pdf";
	printMap(map_timewalk_raw, out_file, "map_timewalk_raw", "voltage (mV)", "time of arrival (ns)", "# entries", -1.0, -1.0, 1.1, 1.3);

	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_surrogate.pdf";
	printMap(map_surrogate, out_file, "map_surrogate", "voltage (mV)", "time over threshold (ADC)", "# entries", -1.0, -1.0, 1.1, 1.3);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_timewalk.pdf";
	printMap(map_timewalk, out_file, "map_timewalk", "voltage (mV)", "time of arrival (ns)", "# entries", -1.0, -1.0, 1.1, 1.3);

	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_map_tot_a.pdf";
	printMap(map_atot, out_file, "map_tot_a", "column", "row", "a_tot paramater (ADC/mV)", 0.05, 0.15, 1.1, 1.2);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_map_tot_b.pdf";
	printMap(map_btot, out_file, "map_tot_b", "column", "row", "b_tot paramater (ADC)", -20, 20, 1.1, 1.2);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_map_tot_c.pdf";
	printMap(map_ctot, out_file, "map_tot_c", "column", "row", "c_tot paramater (ADC*mV)", 0, 500, 1.1, 1.2);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_map_tot_t.pdf";
	printMap(map_ttot, out_file, "map_tot_t", "column", "row", "t_tot paramater (mV)", 0, 40, 1.1, 1.2);

	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_map_toa_c.pdf";
	printMap(map_ctoa, out_file, "map_toa_c", "column", "row", "c_toa paramater (ns*mV)", 0, 1500, 1.1, 1.2);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_map_toa_t.pdf";
	printMap(map_ttoa, out_file, "map_toa_t", "column", "row", "t_toa paramater (mV)", -20, 20, 1.1, 1.2);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_map_toa_d.pdf";
	printMap(map_dtoa, out_file, "map_toa_d", "column", "row", "d_toa paramater (ns)", -20, 20, 1.1, 1.2);

	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_map_thr.pdf";
	printMap(map_thr, out_file, "map_thr", "column", "row", "estimated thr (mV)", 0, 50, 1.1, 1.2);

	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_map_chi2_tot.pdf";
	printMap(map_chi2tot, out_file, "map_tot_chi2", "column", "row", "chi2/ndf", 0, 50, 1.1, 1.2);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_map_chi2_toa.pdf";
	printMap(map_chi2toa, out_file, "map_toa_chi2", "column", "row", "chi2/ndf", 0, 50, 1.1, 1.2);

	styleGraph();

	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_hist_tot_a.pdf";
	printHist(hist_atot, out_file, "hist_tot_a", "a_tot paramater (ADC/mV)", "# entries", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_hist_tot_b.pdf";
	printHist(hist_btot, out_file, "hist_tot_b", "b_tot paramater (ADC)", "# entries", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_hist_tot_c.pdf";
	printHist(hist_ctot, out_file, "hist_tot_c", "c_tot paramater (ADC*mV)", "# entries", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_hist_tot_t.pdf";
	printHist(hist_ttot, out_file, "hist_tot_t", "t_tot paramater (mV)", "# entries", -1.0, -1.0, 1.3);

	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_hist_toa_c.pdf";
	printHist(hist_ctoa, out_file, "hist_toa_c", "c_toa paramater (ns*mV)", "# entries", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_hist_toa_t.pdf";
	printHist(hist_ttoa, out_file, "hist_toa_t", "t_toa paramater (mV)", "# entries", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_hist_toa_d.pdf";
	printHist(hist_dtoa, out_file, "hist_toa_d", "d_toa paramater (ns)", "# entries", -1.0, -1.0, 1.3);

	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_hist_thr.pdf";
	printHist(hist_thr, out_file, "hist_thr", "estimated thr (mV)", "# entries", -1.0, -1.0, 1.3);

	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id +"_cal_hist_chi2_tot.pdf";
	printHist(hist_chi2tot, out_file, "hist_chi2_tot", "chi2/ndf", "# entries", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id +"_cal_hist_chi2_toa.pdf";
	printHist(hist_chi2toa, out_file, "hist_chi2_toa", "chi2/ndf", "# entries", -1.0, -1.0, 1.3);




	// Output
	// ---------------------------------------

	std::string hd_tot = "# column | row | a_tot (ADC/mV) | b_tot (ADC) | c_tot (ADC*mV) | t_tot (mV) | thr (mV) | thr_err (mV) | chi2/ndf";
	std::string hd_toa = "# column | row | c_toa (ns*mV) | t_toa (mV) | d_toa (ns) | chi2/ndf";

	out_file = dat_path + "/lab/" + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_tot.txt";
	write_file(out_file, ' ', out_tot, hd_tot, info);

	out_file = dat_path + "/lab/" + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_toa.txt";
	write_file(out_file, ' ', out_toa, hd_toa, info);

	if (1) {
		std::vector<std::vector<float>> a_tot;
		std::vector<std::vector<float>> b_tot;
		std::vector<std::vector<float>> c_tot;
		std::vector<std::vector<float>> t_tot;
		std::vector<std::vector<float>> c_toa;
		std::vector<std::vector<float>> t_toa;
		std::vector<std::vector<float>> d_toa;

		for (Int_t j = 0; j < 256; j++) {
			std::vector<float> ra_tot; std::vector<float> rb_tot; std::vector<float> rc_tot; std::vector<float> rt_tot;
		    std::vector<float> rc_toa; std::vector<float> rt_toa; std::vector<float> rd_toa;

			for (Int_t i = 0; i < 256; i++) {
				ra_tot.push_back(out_tot[256*j + i][2]);
				rb_tot.push_back(out_tot[256*j + i][3]);
				rc_tot.push_back(out_tot[256*j + i][4]);
				rt_tot.push_back(out_tot[256*j + i][5]);
				rc_toa.push_back(out_toa[256*j + i][2]);
				rt_toa.push_back(out_toa[256*j + i][3]);
				rd_toa.push_back(out_toa[256*j + i][4]);
			}
			a_tot.push_back(ra_tot);
			b_tot.push_back(rb_tot);
			c_tot.push_back(rc_tot);
			t_tot.push_back(rt_tot);
			c_toa.push_back(rc_toa);
			t_toa.push_back(rt_toa);
			d_toa.push_back(rd_toa);
		}

		std::string hd = "";
		out_file = dat_path + "/lab/" + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/a_tot.dat";
		write_file(out_file, ' ', a_tot, hd, info);
		out_file = dat_path + "/lab/" + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/b_tot.dat";
		write_file(out_file, ' ', b_tot, hd, info);
		out_file = dat_path + "/lab/" + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/c_tot.dat";
		write_file(out_file, ' ', c_tot, hd, info);
		out_file = dat_path + "/lab/" + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/t_tot.dat";
		write_file(out_file, ' ', t_tot, hd, info);

		out_file = dat_path + "/lab/" + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/c_toa.dat";
		write_file(out_file, ' ', c_toa, hd, info);
		out_file = dat_path + "/lab/" + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/t_toa.dat";
		write_file(out_file, ' ', t_toa, hd, info);
		out_file = dat_path + "/lab/" + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/d_toa.dat";
		write_file(out_file, ' ', d_toa, hd, info);
	}

	std::string cmd = "cp -r " + dat_path + "/lab/" + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac)
		+ "_ik_" + std::to_string(ik_dac) + "/*[.dat,.txt] " + cfg_path + dev_id + "/";
	std::cout << cmd << std::endl;
	system(cmd.c_str());

	std::cout << "Finished.\nFailed Pixels ToT: " << nfail_tot << " ToA: " << nfail_toa << std::endl;
}
