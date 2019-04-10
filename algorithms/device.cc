/********************************************************************
* File: device.cc
* ------------------------
*
* Description:
* Analyse tpx device data.
*
*
* Version:
* Author: Florian Pitters
*
*******************************************************************/

#include "device.h"

#define IN 1
#define OUT 0


void dev_analysis(std::string dev_id, int pol)
{
	gROOT->Reset();
	gROOT->SetBatch(kTRUE);

	styleCommon();
	styleMap();
	styleGraph();



	// Preperations
	// ---------------------------------------

	// initalise
	int cnt, est;
	float eq_dat[3][ncol][nrow];
	float rms[ncol][nrow], bl[ncol][nrow], mask[ncol][nrow], trim[ncol][nrow], rms_before[ncol][nrow], bl_before[ncol][nrow], rms_gray[ncol][nrow];
	std::string temp_file, out_file;

	std::vector<float> x, x_err, y, y_err;
	std::vector<std::vector<float>> dat_noise, dat_gain;
	std::vector<std::vector<float>> gain_codes0, gain_codes1, gain_codes2, gain_codes3;
	std::vector<std::vector<float>> gain_vals0, gain_vals1, gain_vals2, gain_vals3;
	std::vector<std::vector<float>> noise_vals, noise_codes;

	// get short name
	std::string dev_id_short = "";
	for (Int_t i = 0; i < strlen(dev_id.c_str()); i++) {
		if (dev_id.at(i) != '0' || i == (strlen(dev_id.c_str())-1)) {
			dev_id_short.push_back(dev_id.at(i));
		}
	}

	// load eqalisation data
	bool info = false;
	temp_file = dat_path + "/lab/" + dev_id + "/" + cond + "/dev_equal/" + "dacs0_bl_thr.dat";
	read_data_mat(temp_file.c_str(), eq_dat[0], info);	// noise edge at trim dac 0
	temp_file = dat_path + "/lab/" + dev_id + "/" + cond + "/dev_equal/" + "dacsF_bl_thr.dat";
	read_data_mat(temp_file.c_str(), eq_dat[1], info); // noise edge at trim dac F

	temp_file = dat_path + "/lab/" + dev_id + "/" + cond + "/dev_equal/" + "eq_bl_measured_finestep_" + dev_id_short + ".dat";
	read_data_mat(temp_file.c_str(), eq_dat[2], info); // noise edge after equalisation
	temp_file = dat_path + "/lab/" + dev_id + "/" + cond + "/dev_equal/" + "eq_codes.dat";
	read_data_mat(temp_file.c_str(), trim, info); // trim dacs after equalisation
	temp_file = dat_path + "/lab/" + dev_id + "/" + cond + "/dev_equal/" + "eq_mask.dat";
	read_data_mat(temp_file.c_str(), mask, info); // mask after equalisation

	// load noise data
	temp_file = dat_path + "/lab/" + dev_id + "/" + cond + "/dev_noise/" + "rms.map";
	read_data_mat(temp_file.c_str(), rms, info); // width of gaussian after threshold scan over noise floor
	temp_file = dat_path + "/lab/" + dev_id + "/" + cond + "/dev_noise/" + "bl.map";
	read_data_mat(temp_file.c_str(), bl, info); // mean of gaussian after threshold scan over noise floor

	temp_file = dat_path + "/lab/" + dev_id + "/" + cond + "/dev_noise_gray_code_on/" + "rms.map";
	read_data_mat(temp_file.c_str(), rms_gray, info); // width of gaussian after threshold scan over noise floor
	temp_file = dat_path + "/lab/" + dev_id + "/" + cond + "/wafer_probing/test07_noise_scan/" + "rms.map";
	read_data_mat(temp_file.c_str(), rms_before, info); // width of gaussian after threshold scan over noise floor
	temp_file = dat_path + "/lab/" + dev_id + "/" + cond + "/wafer_probing/test07_noise_scan/" + "bl.map";
	read_data_mat(temp_file.c_str(), bl_before, info); // mean of gaussian after threshold scan over noise floor

	// load gain data
	temp_file = dat_path + "/lab/" + dev_id + "/" + cond + "/dev_scurves/details/amp00";
	read_gain_data(temp_file.c_str(), ' ', gain_vals0, gain_codes0, info);
	temp_file = dat_path + "/lab/" + dev_id + "/" + cond + "/dev_scurves/details/amp01";
	read_gain_data(temp_file.c_str(), ' ', gain_vals1, gain_codes1, info);
	temp_file = dat_path + "/lab/" + dev_id + "/" + cond + "/dev_scurves/details/amp02";
	read_gain_data(temp_file.c_str(), ' ', gain_vals2, gain_codes2, info);
	temp_file = dat_path + "/lab/" + dev_id + "/" + cond + "/dev_scurves/details/amp03";
	read_gain_data(temp_file.c_str(), ' ', gain_vals3, gain_codes3, info);

	temp_file = dat_path + "/lab/" + dev_id + "/" + cond + "/dev_noise/details";
	read_noise_data(temp_file.c_str(), ' ', noise_vals, noise_codes, info);

	float LSB_to_e_conv = 1./0.098;

	auto gaussf = new TF1("gaussf", "gaus", 800, 1200);
	gaussf->SetNpx(500);
	gaussf->SetParameter(0, 500);
	gaussf->SetParameter(1, 1000);
	gaussf->SetParameter(2, 10);

	auto errorf = new TF1("errorf", "[0] * 0.5 * (1 - erf((x-[1]) / (1.414*abs([2]))))", 0, 2000);
	errorf->SetNpx(500);
	errorf->SetParameter(0, 200);
	errorf->SetParameter(1, 1300);
	errorf->SetParameter(2, 10);

	auto linef = new TF1("linef", "[0] + x * [1]", 0, 1000);
	linef->SetNpx(500);
    linef->SetParameter(0, 1100);
    linef->SetParameter(1, 0);



	// Analysis
	// ---------------------------------------

	// noise & equalisation
	TH1F *hist_rms = new TH1F("hist_rms", "Noise in thr_dac", 50, 0, 100);
	TH1F *hist_bl = new TH1F("hist_bl", "Baseline in thr_dac", 200, 1000, 1200);
	TH1F *hist_mask = new TH1F("hist_mask", "Masked Pixels", 2, 0, 1);
	TH1F *hist_trim = new TH1F("hist_trim", "Trim Dacs", 16, 0, 16);

	TH2F *map_rms = new TH2F("map_rms", "Noise in thr_dac", 256, 0, 256, 256, 0, 256);
	TH2F *map_bl = new TH2F("map_bl", "Baseline in thr_dac", 256, 0, 256, 256, 0, 256);
	TH2F *map_mask = new TH2F("map_mask", "Masked Pixels", 256, 0, 256, 256, 0, 256);
	TH2F *map_trim = new TH2F("map_trim", "Trim Dacs", 256, 0, 256, 256, 0, 256);

	TH1F *hist_bleq = new TH1F("hist_bleq", "Baseline after Equalisation", 200, 900, 1300);
	TH1F *hist_bl0 = new TH1F("hist_bl0", "Baseline before Equalisation Trim DACs 0x0", 200, 900, 1300);
	TH1F *hist_blF = new TH1F("hist_blF", "Baseline before Equalisation Trim DACs 0xF", 200, 900, 1300);
	TH2F *map_bleq = new TH2F("map_bleq", "Baseline after Equalisation", 256, 0, 256, 256, 0, 256);
	TH2F *map_bl0 = new TH2F("map_bl0", "Baseline before Equalisation Trim DACs 0x0", 256, 0, 256, 256, 0, 256);
	TH2F *map_blF = new TH2F("map_blF", "Baseline before Equalisation Trim DACs 0xF", 256, 0, 256, 256, 0, 256);

	TH1F *hist_rms_before = new TH1F("hist_rms_before", "Noise in thr_dac", 50, 0, 100);
	TH2F *map_rms_before = new TH2F("map_rms_before", "Noise in thr_dac", 256, 0, 256, 256, 0, 256);

	for (Int_t i = 0; i < 256; i++) {
		for (Int_t j = 0; j < 256; j++) {
			hist_bl->Fill(bl[j][i]);
			hist_rms->Fill(rms[j][i]);
			hist_mask->Fill(mask[j][i]);
			hist_trim->Fill(trim[j][i]);

			map_bl->Fill(i, j, bl[j][i]);
			map_rms->Fill(i, j, rms[j][i]);
			map_mask->Fill(i, j, mask[j][i]);
			map_trim->Fill(i, j, trim[j][i]);

			hist_bleq->Fill(eq_dat[2][j][i]);
			hist_bl0->Fill(eq_dat[0][j][i]);
			hist_blF->Fill(eq_dat[1][j][i]);

			map_bleq->Fill(i, j, eq_dat[2][j][i]);
			map_bl0->Fill(i, j, eq_dat[0][j][i]);
			map_blF->Fill(i, j, eq_dat[1][j][i]);

			hist_rms_before->Fill(rms_before[j][i]);
			map_rms_before->Fill(i, j, rms_before[j][i]);

			dat_noise.push_back({static_cast<float>(i), static_cast<float>(j), bl[j][i], rms[j][i], bl_before[j][i], rms_before[j][i], rms_gray[j][i]});
		}
	}


	// gain
	TH1F *hist_gain = new TH1F("hist_gain", "Gain", 150, 1, 2.5);
	TH2F *map_gain = new TH2F("map_gain", "Gain", 256, 0, 256 , 1, 0, 1);
	TH1F *hist_offset = new TH1F("hist_offset", "Gain", 50, 0, 50);

	std::vector<float> p, m, w, me;
	std::vector<std::vector<float>> pixels, means, widths, means_err;
	std::vector<std::vector<std::vector<float>>> codes = {gain_codes0, gain_codes1, gain_codes2, gain_codes3};
	std::vector<std::vector<std::vector<float>>> vals = {gain_vals0, gain_vals1, gain_vals2, gain_vals3};

	// polarity 0 is electrons, 1 is holes
	for (int g = 0; g < 4; g++){
		if (codes.at(g).size() != 256){
			std::cout << g << " " << codes.at(g).size() << std::endl;
		}
		for (int pix = 0; pix < 256; pix++){
			std::vector<float> codes_err(codes.at(g).at(pix).size(), 0);
			std::vector<float> vals_err(codes.at(g).at(pix).size(), 10);

			if (codes.at(g).at(pix).at(0) == -1){
				continue;
			}

			auto egr = new TGraphErrors(codes.at(g).at(pix).size(), &codes.at(g).at(pix)[0], &vals.at(g).at(pix)[0],
					&codes_err[0], &vals_err[0]);

			if (g == 0) {
				est = codes.at(g).at(pix).at(find_extremum(vals.at(g).at(pix), "max"));

				gaussf->SetParameter(1, est);
				egr->Fit("gaussf", "MQ", "", est-100, est+100);
				if (gaussf->GetChisquare()/gaussf->GetNDF() < 20){
					p.push_back(pix);
					m.push_back(gaussf->GetParameter(1));
					w.push_back(gaussf->GetParameter(2));
					me.push_back(gaussf->GetParError(1));
				}
			}
			else {
				est = codes.at(g).at(pix).at(find_first_above_x(vals.at(g).at(pix), 100, pow(-1, pol)));

				errorf->SetParameter(1, est);
				egr->Fit("errorf", "MQ", "", est-40, est+40);
				if (errorf->GetChisquare()/errorf->GetNDF() < 20 && errorf->GetParError(1) < 2.5){
					p.push_back(pix);
					m.push_back(errorf->GetParameter(1));
					w.push_back(errorf->GetParameter(2));
					me.push_back(errorf->GetParError(1));

				}
				if (0) {
					std::cout << est << " " << errorf->GetParameter(1) << " " << errorf->GetParError(1) << " "
							  << errorf->GetChisquare()/errorf->GetNDF() << std::endl;
				}
				if (1) {
					if (abs(errorf->GetParameter(1) - est) > 10 || 1) {
						out_file = res_path + dev_id + "/dev/example_bad_fit_curve_" + std::to_string(g) + "_" + std::to_string(pix) + "_" + dev_id + ".pdf";
						printGraph(egr, out_file, "example volt", "threshold (LSB)", "# counts", -1.0, -1.0, 1.3);
					}
				}
			}
			delete egr;
		}
		pixels.push_back(p);
		means.push_back(m);
		means_err.push_back(me);
		widths.push_back(w);
		p.clear();
		m.clear();
		w.clear();
		me.clear();
	}

	int mlen = std::min(means.at(0).size(), std::min(means.at(1).size(), std::min(means.at(2).size(), means.at(3).size())));

	if (dev_id == "W0005_E02") { x = {49.5, 148.9, 296.8}; }
	else if (dev_id == "W0005_F1") { x = {49.1, 148.3, 297.5}; }
	else if (dev_id == "W0019_C07") { x = {49.2, 150.7, 296.8}; }
	else if (dev_id == "W0019_F07") { x = {50.9, 148.4, 298.0}; }
	else if (dev_id == "W0019_G07") { x = {48.5, 148.0, 297.1}; }
	else { x = {50, 150, 300}; }

	for (int pix = 0; pix < mlen; pix++){
		x_err = {3, 3, 3};
		y = {means.at(1).at(pix), means.at(2).at(pix), means.at(3).at(pix)};
		y_err = {means_err.at(1).at(pix), means_err.at(2).at(pix), means_err.at(3).at(pix)};

		auto eg = std::make_shared<TGraphErrors> (x.size(), &x[0], &y[0], &x_err[0], &y_err[0]);
		eg->Fit("linef", "meq", 0, 400);
		dat_gain.push_back({ static_cast<float>(pix), static_cast<float>(pix),
			static_cast<float>(means.at(0).at(pix)), static_cast<float>(widths.at(0).at(pix)),
			static_cast<float>(linef->GetParameter(0)), static_cast<float>(linef->GetParameter(1)) });

		hist_offset->Fill( abs(dat_gain.back().at(4) - dat_gain.back().at(2)) / abs(dat_gain.back().at(5)) );
		hist_gain->Fill( abs(dat_gain.back().at(5)) );

		map_gain->Fill(pixels.at(0).at(pix), 0.5, abs(dat_gain.back().at(5)) );
	}

	std::cout << "Baseline Dac: " << hist_bl->GetMean() << std::endl;
	std::cout << "Inj Point 1: Mean " << mean(means.at(0)) << " Number of successful fits: " << means.at(0).size() << std::endl;
	std::cout << "Inj Point 2: Mean " << mean(means.at(1)) << " Number of successful fits: " << means.at(1).size() << std::endl;
	std::cout << "Inj Point 3: Mean " << mean(means.at(2)) << " Number of successful fits: " << means.at(2).size() << std::endl;
	std::cout << "Inj Point 4: Mean " << mean(means.at(3)) << " Number of successful fits: " << means.at(3).size() << std::endl;

	std::string hd_noise = "# column | row | baseline (ADC) | noise (ADC) | baseline WP (ADC) | noise WP (ADC) ";
	std::string hd_gain = "# column | row | baseline (ADC) | noise (ADC) | offset (ADC) | gain (ADC/mV) ";

	out_file = dat_path + "/lab/" + dev_id + "/" + cond + "/dev_noise/" + dev_id + "_dev_noise.txt";
	write_file(out_file, ' ', dat_noise, hd_noise, info);

	out_file = dat_path + "/lab/" + dev_id + "/" + cond + "/dev_scurves/" + dev_id + "_dev_gain.txt";
	write_file(out_file, ' ', dat_gain, hd_gain, info);



	// investigate baseline shift with testpulses 
	TH1F *hist_bl_diag = new TH1F("hist_bl_diag", "Baseline in thr_dac", 200, 1000, 1200);
	TH1F *hist_bl_diag_with_tp = new TH1F("hist_bl_diag_with_tp", "Baseline in thr_dac", 200, 1000, 1200);

	for (int pix = 0; pix < 256; pix++) {
		est = noise_codes.at(pix).at(find_extremum(noise_vals.at(pix), "max"));
		auto gr = new TGraph(noise_codes.at(pix).size(), &noise_codes.at(pix)[0], &noise_vals.at(pix)[0]);
		gaussf->SetParameter(1, est);
		gr->Fit("gaussf", "MQ", "", est-100, est+100);
		hist_bl_diag->Fill(gaussf->GetParameter(1));

		if (pix == 100) {
			out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_bl_single.pdf";
			printGraph(gr, out_file, "Singles BL", "threshold (LSB)", "# counts", -1.0, -1.0, 1.3);
		}

		est = codes.at(0).at(pix).at(find_extremum(vals.at(0).at(pix), "max"));
		gr = new TGraphErrors(codes.at(0).at(pix).size(), &codes.at(0).at(pix)[0], &vals.at(0).at(pix)[0]);
		gaussf->SetParameter(1, est);
		gr->Fit("gaussf", "MQ", "", est-100, est+100);
		hist_bl_diag_with_tp->Fill(gaussf->GetParameter(1));

		if (pix == 100) {
			out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_bl_single_with_tp.pdf";
			printGraph(gr, out_file, "Singles BL", "threshold (LSB)", "# counts", -1.0, -1.0, 1.3);
		}
	}




	// Printing
	// ---------------------------------------

	printf("---->  Finalising analysis\n");

	auto t = new TLatex();
	t->SetNDC();
	t->SetTextFont(42);
	t->SetTextAlign(12);
	t->SetTextSize(0.040);
	std::stringstream tmp_stream;

	// Noise & Baseline
	styleMap();
	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_rms_matrix.pdf";
	printMap(map_rms, out_file, "map_rms", "column", "row", "noise [LSB]", 3, 8, 1.1, 1.0);
	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_bl_matrix.pdf";
	printMap(map_bl, out_file, "map_bl", "column", "row", "baseline [LSB]", 1080, 1180, 1.1, 1.3);
	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_trim_matrix.pdf";
	printMap(map_trim, out_file, "map_trim", "column", "row", "trim dac [LSB]", 0, 15, 1.1, 1.1);
	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_mask_matrix.pdf";
	printMap(map_mask, out_file, "map_mask", "column", "row", "mask", 0, 5, 1.1, 1.1);

	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_rms_matrix_before.pdf";
	printMap(map_rms_before, out_file, "map_rms_before", "column", "row", "noise [LSB]", 3, 8, 1.1, 1.0);

	styleGraph();
	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_rms_hist.pdf";
	printHist(hist_rms, out_file, "hist_rms", "noise [LSB]", "# pixels", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_bl_hist.pdf";
	printHist(hist_bl, out_file, "hist_bl", "baseline [LSB]", "# pixels", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_trim_hist.pdf";
	printHist(hist_trim, out_file, "hist_trim", "trim dac [LSB]", "# pixels", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_mask_hist.pdf";
	printHist(hist_mask, out_file, "hist_mask", "mask]", "# pixels", -1.0, -1.0, 1.3);

	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_rms_hist_before.pdf";
	printHist(hist_rms_before, out_file, "hist_rms_before", "noise [LSB]", "# pixels", -1.0, -1.0, 1.3);

	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_bl_hist_diag_with_tp.pdf";
	printHist(hist_bl_diag_with_tp, out_file, "hist_bl2", "baseline [LSB]", "# pixels", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_bl_hist_diag.pdf";
	printHist(hist_bl_diag, out_file, "hist_bl2", "baseline [LSB]", "# pixels", -1.0, -1.0, 1.3);



	// Equalisation
	styleGraph();
	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_bl_hist_vals0.pdf";
	printHist(hist_bl0, out_file, "hist_blvals0", "baseline [LSB]", "# pixels", -1.0, -1.0, 1.1);
	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_bl_hist_valsF.pdf";
	printHist(hist_blF, out_file, "hist_blvalsF", "baseline [LSB]", "# pixels", -1.0, -1.0, 1.1);
	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_bl_hist_valseq.pdf";
	printHist(hist_bleq, out_file, "hist_blvalseq", "baseline [LSB]", "# pixels", -1.0, -1.0, 1.1);

	styleMap();
	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_bl_map_vals0.pdf";
	printMap(map_bl0, out_file, "map_blvals0", "column", "row", "noise [LSB]", 900, 1400, 1.1, 1.1);
	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_bl_map_valsF.pdf";
	printMap(map_blF, out_file, "map_blvalsF", "column", "row", "baseline [LSB]", 980, 1400, 1.1, 1.1);
	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_bl_map_valseq.pdf";
	printMap(map_bleq, out_file, "map_blvalseq", "column", "row", "trim dac [LSB]", 900, 1400, 1.1, 1.1);

	styleGraph();
	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_eq_plot.pdf";
	TCanvas *canv_eq = new TCanvas("hist", "hist", 800, 400);
	THStack *hs_eq = new THStack();
	TLegend *lg_eq = new TLegend(0.68, 0.735, 0.88, 0.88);
	addHist(hs_eq, lg_eq, hist_bl0, colours.at(2), markers.at(2), 1, "l", " trim dacs 0", 1);
	addHist(hs_eq, lg_eq, hist_blF, colours.at(1), markers.at(1), 1, "l", " trim dacs F", 1);
	addHist(hs_eq, lg_eq, hist_bleq, colours.at(0), markers.at(0), 1, "l", " equalised", 1);
	hs_eq->Draw("nostack");
	hs_eq->SetTitle("Title");
	hs_eq->GetXaxis()->SetTitle("threshold [LSB]");
	hs_eq->GetYaxis()->SetTitle("# pixels");
	hs_eq->GetYaxis()->SetTitleOffset(1.1);
	hs_eq->GetXaxis()->SetLimits(950, 1250);
	lg_eq->Draw();
	t->SetTextColor(colours.at(1));
	tmp_stream.str(std::string());
	tmp_stream << "rms = " << std::fixed << std::setprecision(2) << hist_blF->GetRMS() << " LSB";
	t->DrawLatex(0.65, 0.43, tmp_stream.str().c_str());
	t->SetTextColor(colours.at(2));
	tmp_stream.str(std::string());
	tmp_stream << "rms = " << std::fixed << std::setprecision(2) << hist_bl0->GetRMS() << " LSB";
	t->DrawLatex(0.28, 0.43, tmp_stream.str().c_str());
	t->SetTextColor(colours.at(0));
	tmp_stream.str(std::string());
	tmp_stream << "rms = " << std::fixed << std::setprecision(2) << hist_bleq->GetRMS() << " LSB";
	t->DrawLatex(0.36, 0.79, tmp_stream.str().c_str());
	canv_eq->SaveAs(out_file.c_str());


	// Gain
	styleMap();
	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_gain_map.pdf";
	printMap(map_gain, out_file, "map_gain", "diagonal pixel", " ", "gain [LSB/mV]", 1.5, 2, 1.1, 1.1);

	styleGraph();
	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_gain_hist.pdf";
	printHist(hist_gain, out_file, "hist_gain", "gain [LSB/mV]", "# pixels", -1.0, -1.0, 1.1);
	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_gain_offset.pdf";
	printHist(hist_offset, out_file, "hist_offset", "offset [mV]", "# pixels", -1.0, -1.0, 1.1);

	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_gain_single.pdf";
	TCanvas *canv_gain_single = new TCanvas("gain_single", "gain_single");
	TMultiGraph *mg_gain_single = new TMultiGraph();
	TLegend *lg_gain_single = new TLegend(0.65, 0.79, 0.88, 0.88);
	for (Int_t i = 120; i < 121; i++) {
		addGraph(mg_gain_single, lg_gain_single, gain_codes0.at(i).size(), &gain_codes0.at(i)[0], &gain_vals0.at(i)[0],
					colours.at(3), markers.at(1), 0.5, styles.at(1), "p", "  gain0");
		addGraph(mg_gain_single, lg_gain_single, gain_codes1.at(i).size(), &gain_codes1.at(i)[0], &gain_vals1.at(i)[0],
					colours.at(2), markers.at(1), 0.5, styles.at(1), "p", "  gain1");
		addGraph(mg_gain_single, lg_gain_single, gain_codes2.at(i).size(), &gain_codes2.at(i)[0], &gain_vals2.at(i)[0],
					colours.at(1), markers.at(1), 0.5, styles.at(1), "p", "  gain2");
		addGraph(mg_gain_single, lg_gain_single, gain_codes3.at(i).size(), &gain_codes3.at(i)[0], &gain_vals3.at(i)[0],
					colours.at(0), markers.at(1), 0.5, styles.at(1), "p", "  gain3");
	}
	mg_gain_single->Draw("AP");
	mg_gain_single->GetXaxis()->SetTitle("theshold [LSB]");
	mg_gain_single->GetYaxis()->SetTitle("# above threshold");
	mg_gain_single->GetYaxis()->SetTitleOffset(1.3);
	if (pol == 0) mg_gain_single->GetXaxis()->SetLimits(1000, 2000);
	else mg_gain_single->GetXaxis()->SetLimits(400, 1400);
	mg_gain_single->SetMinimum(0);
	mg_gain_single->SetMaximum(1200);
	canv_gain_single->SaveAs(out_file.c_str());

	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_gain_all.pdf";
	TCanvas *canv_gain_all = new TCanvas("gain_all", "gain_all");
	TMultiGraph *mg_gain_all = new TMultiGraph();
	TLegend *lg_gain_all = new TLegend(0.55, 0.70, 0.88, 0.88);
	for (Int_t i = 0; i < 1; i++) {
		addGraph(mg_gain_all, lg_gain_all, gain_codes0.at(i).size(), &gain_codes0.at(i)[0], &gain_vals0.at(i)[0],
					colours.at(3), markers.at(1), 0.2, styles.at(1), "p", "  noise");
		addGraph(mg_gain_all, lg_gain_all, gain_codes1.at(i).size(), &gain_codes1.at(i)[0], &gain_vals1.at(i)[0],
					colours.at(2), markers.at(1), 0.2, styles.at(1), "p", "  low amplitude");
		addGraph(mg_gain_all, lg_gain_all, gain_codes2.at(i).size(), &gain_codes2.at(i)[0], &gain_vals2.at(i)[0],
					colours.at(1), markers.at(1), 0.2, styles.at(1), "p", "  mid amplitude");
		addGraph(mg_gain_all, lg_gain_all, gain_codes3.at(i).size(), &gain_codes3.at(i)[0], &gain_vals3.at(i)[0],
					colours.at(0), markers.at(1), 0.2, styles.at(1), "p", "  high amplitude");
	}
	for (Int_t i = 1; i < mlen; i++) {
		addGraph(mg_gain_all, lg_gain_all, gain_codes0.at(i).size(), &gain_codes0.at(i)[0], &gain_vals0.at(i)[0],
					colours.at(3), markers.at(1), 0.2, styles.at(1), "p", "");
		addGraph(mg_gain_all, lg_gain_all, gain_codes1.at(i).size(), &gain_codes1.at(i)[0], &gain_vals1.at(i)[0],
					colours.at(2), markers.at(1), 0.2, styles.at(1), "p", "");
		addGraph(mg_gain_all, lg_gain_all, gain_codes2.at(i).size(), &gain_codes2.at(i)[0], &gain_vals2.at(i)[0],
					colours.at(1), markers.at(1), 0.2, styles.at(1), "p", "");
		addGraph(mg_gain_all, lg_gain_all, gain_codes3.at(i).size(), &gain_codes3.at(i)[0], &gain_vals3.at(i)[0],
					colours.at(0), markers.at(1), 0.2, styles.at(1), "p", "");
	}
	mg_gain_all->Draw("AP");
	mg_gain_all->SetTitle("Gain Calibration");
	mg_gain_all->GetXaxis()->SetTitle("threshold [LSB]");
	mg_gain_all->GetYaxis()->SetTitle("# above threshold");
	mg_gain_all->GetYaxis()->SetTitleOffset(1.3);
	if (pol == 0) mg_gain_all->GetXaxis()->SetLimits(1050, 1800);
	else mg_gain_all->GetXaxis()->SetLimits(400, 1400);
	mg_gain_all->SetMinimum(-2);
	mg_gain_all->SetMaximum(600);
	lg_gain_all->Draw();
	canv_gain_all->SaveAs(out_file.c_str());


	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_gain_fit_single2.pdf";
	TCanvas *canv_gain_fit_single2 = new TCanvas("gain_fit_single2", "gain_fit_single2");
	TMultiGraph *mg_gain_fit_single2 = new TMultiGraph();
	TLegend *lg_gain_fit_single2 = new TLegend(0.65, 0.79, 0.88, 0.88);
	if (dev_id == "W0005_E02") { x = {49.5, 148.9, 296.8}; }
	else if (dev_id == "W0005_F1") { x = {49.1, 148.3, 297.5}; }
	else if (dev_id == "W0019_C07") { x = {49.2, 150.7, 296.8}; }
	else if (dev_id == "W0019_F07") { x = {50.9, 148.4, 298.0}; }
	else if (dev_id == "W0019_G07") { x = {48.5, 148.0, 297.1}; }
	else { x = {50, 150, 300}; }
	x_err = {1,1,1};
	y = {means.at(1).at(120), means.at(2).at(120), means.at(3).at(120)};
	y_err = {means_err.at(1).at(120), means_err.at(2).at(120), means_err.at(3).at(120)};
	auto fit1 = addFittedErrGraph(mg_gain_fit_single2, lg_gain_fit_single2, x.size(), &x[0], &y[0], &x_err[0], &y_err[0],
		colours.at(0), markers.at(0), 1, styles.at(1), "p", " injection data", "linef", 0, 400);
	x = {0};
	x_err = {1};
	y = {means.at(0).at(120)};
	y_err = {means_err.at(0).at(120)};
	addErrGraph(mg_gain_fit_single2, lg_gain_fit_single2, x.size(), &x[0], &y[0], &x_err[0], &y_err[0],
		colours.at(1), markers.at(1), 1, styles.at(1), "p", " noise data");
	mg_gain_fit_single2->Draw("AP");
	mg_gain_fit_single2->GetXaxis()->SetTitle("injected voltage [mV]");
	mg_gain_fit_single2->GetYaxis()->SetTitle("threshold [LSB]");
	mg_gain_fit_single2->GetYaxis()->SetTitleOffset(1.3);
	mg_gain_fit_single2->GetXaxis()->SetLimits(-10, 350);
	if (pol == 0) {
		mg_gain_fit_single2->SetMinimum(1000);
		mg_gain_fit_single2->SetMaximum(2000);
	}
	else {
		mg_gain_fit_single2->SetMinimum(400);
		mg_gain_fit_single2->SetMaximum(1400);
	}
	addFitResults(linef, lg_gain_fit_single2, fit1, colours.at(0), 1, " injection fit");
	t->SetTextColor(colours.at(0));
	tmp_stream.str(std::string());
	tmp_stream << "slope = " << std::fixed << std::setprecision(2) << fit1->Parameter(1);
	tmp_stream << " +/- " << std::fixed << std::setprecision(2) << fit1->ParError(1) << " LSB/mV";
	t->DrawLatex(0.30, 0.69, tmp_stream.str().c_str());
	tmp_stream.str(std::string());
	tmp_stream << "offset = " << std::fixed << std::setprecision(1) << fit1->Parameter(0);
	tmp_stream << " +/- " << std::fixed << std::setprecision(1) << fit1->ParError(0) << " LSB";
	t->DrawLatex(0.30, 0.66, tmp_stream.str().c_str());
	tmp_stream.str(std::string());
	tmp_stream << "baseline = " << std::fixed << std::setprecision(1) << means.at(0).at(120);
	tmp_stream << " +/- " << std::fixed << std::setprecision(1) << means_err.at(0).at(120) << " LSB";
	t->DrawLatex(0.30, 0.63, tmp_stream.str().c_str());
	tmp_stream.str(std::string());
	lg_gain_fit_single2->Draw();
	canv_gain_fit_single2->SaveAs(out_file.c_str());

	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_gain_fit_all2.pdf";
	TCanvas *canv_gain_fit_all2 = new TCanvas("gain_fit_all2", "gain_fit_all2");
	TMultiGraph *mg_gain_fit_all2 = new TMultiGraph();
	TLegend *lg_gain_fit_all2 = new TLegend(0.65, 0.79, 0.88, 0.88);
	x = {0};
	x_err = {1};
	y = {mean(means.at(0))};
	y_err = {mean(means_err.at(0))};
	addErrGraph(mg_gain_fit_all2, lg_gain_fit_all2, x.size(), &x[0], &y[0], &x_err[0], &y_err[0],
		colours.at(1), markers.at(1), 1, styles.at(1), "p", " noise data");
	if (dev_id == "W0005_E02") { x = {49.5, 148.9, 296.8}; }
	else if (dev_id == "W0005_F1") { x = {49.1, 148.3, 297.5}; }
	else if (dev_id == "W0019_C07") { x = {49.2, 150.7, 296.8}; }
	else if (dev_id == "W0019_F07") { x = {50.9, 148.4, 298.0}; }
	else if (dev_id == "W0019_G07") { x = {48.5, 148.0, 297.1}; }
	else { x = {50, 150, 300}; }
	x_err = {1,1,1};
	y = {mean(means.at(1)), mean(means.at(2)), mean(means.at(3))};
	y_err = {mean(means_err.at(1)), mean(means_err.at(2)), mean(means_err.at(3))};
	fit1 = addFittedErrGraph(mg_gain_fit_all2, lg_gain_fit_all2, x.size(), &x[0], &y[0], &x_err[0], &y_err[0],
		colours.at(0), markers.at(0), 1, styles.at(1), "p", " injection data", "linef", 0, 400);
	mg_gain_fit_all2->Draw("AP");
	mg_gain_fit_all2->GetXaxis()->SetTitle("injected voltage [mV]");
	mg_gain_fit_all2->GetYaxis()->SetTitle("threshold [LSB]");
	mg_gain_fit_all2->GetYaxis()->SetTitleOffset(1.3);
	mg_gain_fit_all2->GetXaxis()->SetLimits(-10, 350);
	if (pol == 0) {
		mg_gain_fit_all2->SetMinimum(1000);
		mg_gain_fit_all2->SetMaximum(2000);
	}
	else {
		mg_gain_fit_all2->SetMinimum(400);
		mg_gain_fit_all2->SetMaximum(1400);
	}
	addFitResults(linef, lg_gain_fit_all2, fit1, colours.at(0), 1, " injection fit");
	t->SetTextColor(colours.at(0));
	tmp_stream.str(std::string());
	tmp_stream << "slope = " << std::fixed << std::setprecision(2) << fit1->Parameter(1);
	tmp_stream << " +/- " << std::fixed << std::setprecision(2) << fit1->ParError(1) << " LSB/mV";
	t->DrawLatex(0.27, 0.72, tmp_stream.str().c_str());
	tmp_stream.str(std::string());
	tmp_stream << "offset = " << std::fixed << std::setprecision(1) << fit1->Parameter(0);
	tmp_stream << " +/- " << std::fixed << std::setprecision(1) << fit1->ParError(0) << " LSB";
	t->DrawLatex(0.27, 0.685, tmp_stream.str().c_str());
	tmp_stream.str(std::string());
	t->SetTextColor(colours.at(1));
	tmp_stream << "baseline = " << std::fixed << std::setprecision(1) << mean(means.at(0));
	tmp_stream << " +/- " << std::fixed << std::setprecision(1) << mean(means_err.at(0)) << " LSB";
	t->DrawLatex(0.27, 0.60, tmp_stream.str().c_str());
	tmp_stream.str(std::string());
	lg_gain_fit_all2->Draw();
	canv_gain_fit_all2->SaveAs(out_file.c_str());



	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_gain_fit_single.pdf";
	TCanvas *canv_gain_fit_single = new TCanvas("gain_fit_single", "gain_fit_single");
	TMultiGraph *mg_gain_fit_single = new TMultiGraph();
	TLegend *lg_gain_fit_single = new TLegend(0.53, 0.79, 0.88, 0.88);
	x = {0, 50, 150, 300};
	x_err = {2, 2, 2, 2};
	y = {means.at(0).at(120), means.at(1).at(120), means.at(2).at(120), means.at(3).at(120)};
	y_err = {means_err.at(0).at(120), means_err.at(1).at(120), means_err.at(2).at(120), means_err.at(3).at(120)};
	fit1 = addFittedErrGraph(mg_gain_fit_single, lg_gain_fit_single, x.size(), &x[0], &y[0], &x_err[0], &y_err[0],
		colours.at(0), markers.at(1), 1, styles.at(1), "p", " injection & noise data", "linef", 0, 400);
	mg_gain_fit_single->Draw("AP");
	mg_gain_fit_single->GetXaxis()->SetTitle("injected voltage [mV]");
	mg_gain_fit_single->GetYaxis()->SetTitle("threshold [LSB]");
	mg_gain_fit_single->GetYaxis()->SetTitleOffset(1.3);
	mg_gain_fit_single->GetXaxis()->SetLimits(-10, 350);
	if (pol == 0) {
		mg_gain_fit_single->SetMinimum(1000);
		mg_gain_fit_single->SetMaximum(2000);
	}
	else {
		mg_gain_fit_single->SetMinimum(400);
		mg_gain_fit_single->SetMaximum(1400);
	}
	addFitResults(linef, lg_gain_fit_single, fit1, colours.at(0), 1, " injection & noise fit");
	t->SetTextColor(colours.at(0));
	tmp_stream.str(std::string());
	tmp_stream << "slope = " << std::fixed << std::setprecision(2) << fit1->Parameter(1);
	tmp_stream << " +/- " << std::fixed << std::setprecision(2) << fit1->ParError(1) << " LSB/mV";
	t->DrawLatex(0.30, 0.69, tmp_stream.str().c_str());
	tmp_stream.str(std::string());
	tmp_stream << "offset = " << std::fixed << std::setprecision(1) << fit1->Parameter(0);
	tmp_stream << " +/- " << std::fixed << std::setprecision(1) << fit1->ParError(0) << " LSB";
	t->DrawLatex(0.30, 0.66, tmp_stream.str().c_str());
	tmp_stream.str(std::string());
	lg_gain_fit_single->Draw();
	canv_gain_fit_single->SaveAs(out_file.c_str());


	out_file = res_path + dev_id + "/dev/" + dev_id + "_dev_gain_fit_all.pdf";
	TCanvas *canv_gain_fit_all = new TCanvas("gain_fit_all", "gain_fit_all");
	TMultiGraph *mg_gain_fit_all = new TMultiGraph();
	TLegend *lg_gain_fit_all = new TLegend(0.53, 0.79, 0.88, 0.88);
	x = {0, 50, 150, 300};
	x_err = {2, 2, 2, 2};
	y = {mean(means.at(0)), mean(means.at(1)), mean(means.at(2)), mean(means.at(3))};
	y_err = {mean(means_err.at(0)), mean(means_err.at(1)), mean(means_err.at(2)), mean(means_err.at(3))};
	fit1 = addFittedErrGraph(mg_gain_fit_all, lg_gain_fit_all, x.size(), &x[0], &y[0], &x_err[0], &y_err[0],
		colours.at(0), markers.at(1), 1, styles.at(1), "p", " injection & noise data", "linef", 0, 400);
	mg_gain_fit_all->Draw("AP");
	mg_gain_fit_all->GetXaxis()->SetTitle("injected voltage [mV]");
	mg_gain_fit_all->GetYaxis()->SetTitle("threshold [LSB]");
	mg_gain_fit_all->GetYaxis()->SetTitleOffset(1.3);
	mg_gain_fit_all->GetXaxis()->SetLimits(-10, 350);
	if (pol == 0) {
		mg_gain_fit_all->SetMinimum(1000);
		mg_gain_fit_all->SetMaximum(2000);
	}
	else {
		mg_gain_fit_all->SetMinimum(400);
		mg_gain_fit_all->SetMaximum(1400);
	}
	addFitResults(linef, lg_gain_fit_all, fit1, colours.at(0), 1, " injection & noise fit");
	t->SetTextColor(colours.at(0));
	tmp_stream.str(std::string());
	tmp_stream << "slope = " << std::fixed << std::setprecision(2) << fit1->Parameter(1);
	tmp_stream << " +/- " << std::fixed << std::setprecision(2) << fit1->ParError(1) << " LSB/mV";
	t->DrawLatex(0.30, 0.69, tmp_stream.str().c_str());
	tmp_stream.str(std::string());
	tmp_stream << "offset = " << std::fixed << std::setprecision(1) << fit1->Parameter(0);
	tmp_stream << " +/- " << std::fixed << std::setprecision(1) << fit1->ParError(0) << " LSB";
	t->DrawLatex(0.30, 0.66, tmp_stream.str().c_str());
	tmp_stream.str(std::string());
	lg_gain_fit_all->Draw();
	canv_gain_fit_all->SaveAs(out_file.c_str());
}
