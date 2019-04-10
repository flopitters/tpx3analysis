/********************************************************************
* File: acq_analysis.cc
* ------------------------
*
* Description:
* Analyse data from tpx3 source acquisition.
*
* Version:
* Author: Florian Pitters
*
*******************************************************************/

#include "analyse.h"


void calibrate_sources_prepare(std::string dev_id, int thr_dac, int ik_dac, std::string src)
{
	gROOT->Reset();
	gROOT->SetBatch(kTRUE);

	styleCommon();
	styleMap();
	styleGraph();



	// Preperations
	// ---------------------------------------

	// initialisation

	int ret;
	std::vector<std::vector<float>> dat_tot, dat_toa;
	std::string temp_file, out_file;

	Int_t x, y, bin; // place holder for tbranch or histogramme data
	Int_t vmax_volts, vmax_tot, nbins_volts, nbins_tot; // histogramme parameters
	Float_t fit_scale, volts, tot; // place holder for tbranch data
	Float_t peak, val, lcut_tot, hcut_tot, lcut_volts, hcut_volts, lcut_fit, hcut_fit; // values for spectral analysis
	Float_t mpv_volts, mpv_err, chi2ndf; // most probable values
	Double_t rho; // value for kernel density estimation

	auto res = new fitResult();

	// configure paths
	std::string dat_file = dat_path + "/xraytube/" + dev_id + "/" + src + "/" + src + "_thr_"
			+ std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "_data.root";

	nbins_volts = nbins_tot = 50;
	fit_scale = 0.3;

	// configure source
	if (src == "fe"){
		lcut_volts = 40; hcut_volts = 120;
		vmax_volts = 150; nbins_volts = 30;
		fit_scale = 0.65;
		rho = pow(1000 * 3/4., (-1./5)); // silverman rule for bandwidth
	}
	else if (src == "in"){
		lcut_volts = 200; hcut_volts = 400;
		vmax_volts = 500; nbins_volts = 30;
		fit_scale = 0.25;
		rho = pow(1000 * 3/4., (-1./5)); // silverman rule for bandwidth
	}
	else {
		peak = 0; lcut_tot = 0; hcut_tot = 0;
		printf("---->  No source given!!\n\n\n\n\n");
	}

	printf("---->  Analyzing data from %s\n", dat_file.c_str());

	// load calibration data
	bool info = true;
	temp_file = dat_path + "/lab/"   + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_tot.txt";
	ret = read_file(temp_file.c_str(), ' ', dat_tot, info);
	temp_file = dat_path + "/lab/"   + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_toa.txt";
	ret = read_file(temp_file.c_str(), ' ', dat_toa, info);


	// open file, get tree, set adresses and add new branch
	TFile *file = new TFile(dat_file.c_str(), "UPDATE");
	file->ls();

	TTree *singles = (TTree*)file->Get("singleHits");

	TBranch *voltbr = singles->Branch("volts", &volts, "volts/F");
	TBranch *totbr = singles->GetBranch("tot_hd");
	TBranch *colbr = singles->GetBranch("col");
	TBranch *rowbr = singles->GetBranch("row");

	singles->SetBranchAddress("col", &x);
	singles->SetBranchAddress("row", &y);
	singles->SetBranchAddress("tot_hd", &tot);

	// fill voltage branch
	ULong64_t nevnt = singles->GetEntries();
	nevnt = int(nevnt*0.99);
	for (ULong64_t i = 0; i < nevnt; i++) {
		colbr->GetEntry(i); // x
		rowbr->GetEntry(i); // y
		totbr->GetEntry(i);

		float params[4] = {dat_tot[256*y + x][2], dat_tot[256*y + x][3], dat_tot[256*y + x][4], dat_tot[256*y + x][5]};
		volts = tot_to_volts(tot, params);

		voltbr->Fill();
	}



	// Analysis
	// ---------------------------------------

	printf("---->  Starting global analysis\n");

	// create histogramm for tot and energy
	auto hist_volts = new TH1F("hist_volts", "Volts Global", int(vmax_volts/2.), 0, vmax_volts);
	auto hist_svolts = new TH1F("hist_svolts", "Single Pixel", nbins_volts, 0, vmax_volts);
	for (ULong64_t i = 0; i < nevnt; i++) {
		rowbr->GetEntry(i);
		colbr->GetEntry(i);
		totbr->GetEntry(i);
		voltbr->GetEntry(i);

		hist_volts->Fill(volts);
		if (x == xpix && y == ypix) {
			hist_svolts->Fill(volts);
		}
	}

	res = fitGauss(hist_volts, fit_scale, lcut_volts, hcut_volts);
	std::cout<< "MPV Estimate Volts [mV]: " << res->mean << std::endl;


	// check beam profile
	TH2F *profile = new TH2F("pr", "Event Count Distribution", 256, 0, 256, 256, 0, 256);
	for (ULong64_t i = 0; i < nevnt; i++) {
		colbr->GetEntry(i);
		rowbr->GetEntry(i);
		profile->Fill(x, y);
	}
	printf("---->  Highest event number per pixel is %.0f\n", profile->GetMaximum());
	printf("---->  Lowest event number per pixel is %.0f\n", profile->GetMinimum());


	// check matrix correlations
	Int_t ntemp1, ntemp2, tempcnt;
	TH1D *htemp1, *htemp2;
	TKDE *kde;
	TH3F *map_fullvolts = new TH3F("map_fullvolts", "voltage Full", 256, 0, 256, 256, 0, 256, nbins_volts, 0, vmax_volts);
	for (ULong64_t i = 0; i < nevnt; i++) {
		colbr->GetEntry(i);
		rowbr->GetEntry(i);
		voltbr->GetEntry(i);
		map_fullvolts->Fill(x, y, volts);
	}

	// across matrix
	auto hist_mpv_volts = new TH1F("hist_mpv_volts", "Volts Global", int(vmax_volts/2.), 0, vmax_volts);
	auto hist_mpv_err_volts = new TH1F("hist_mpv_err_volts", "Volts Global", 100, 0, 25);

	int nfailed = 0;
	std::vector<float> res_volts;
	std::vector<std::vector<float>> out_volts;
	for (int j = 0; j < 256; j++) {
		for (int i = 0; i < 256; i++) {

			// clear vectors
			res_volts.clear();

			// reset
			mpv_volts = 0;
			mpv_err = 0;
			chi2ndf = 0;

			res_volts.push_back(i);
			res_volts.push_back(j);
			if (profile->GetBinContent(i, j) > 800){

				// fill voltage histogramms
				htemp2 = map_fullvolts->ProjectionZ("_pz", i+1, i+1, j+1, j+1);

				// kde option
				//if (fkde == true && fspec == false) {
				if (i == xpix && j%16 == 0 && 0) {
					htemp2->GetXaxis()->SetRangeUser(lcut_volts, hcut_volts);
					ntemp1 = htemp2->GetSize(); // number of bins
					ntemp2 = htemp2->GetEntries(); // number of events

					if (ntemp2 > 10 && ntemp2 <= profile->GetMaximum()){
						Double_t *dat = new Double_t[ntemp2];
						tempcnt = 0;
						for (Int_t k = 1; k < ntemp1-1; k++){
							bin = htemp2->GetBinCenter(k);
							val = htemp2->GetBinContent(k);
							for (Int_t l = 0; l < val; l++){
								dat[tempcnt+l] = bin;
							}
							tempcnt += val;
						}

						kde = new TKDE(ntemp2, &dat[0], 0, vmax_volts, "", rho);
						mpv_volts = kde->GetFunction(5000)->GetMaximumX();

						if (i == xpix && j%16 == 0){
							TF1 *kdef = kde->GetFunction(5000);
							TF1 *kdelow = kde->GetLowerFunction();
							TF1 *kdeup = kde->GetUpperFunction();

							out_file = res_path + dev_id + "/" +  src + "/" + src + "_single_voltage_kde_" + i + "_" + j + ".pdf";
							TCanvas *canvsingleex = new TCanvas("singleex", "singleex");
							htemp2->Scale(1./htemp2->Integral(), "width");
							htemp2->Draw();
							htemp2->GetXaxis()->SetTitle("voltage (mV)");
							htemp2->GetYaxis()->SetTitle("# events");
							htemp2->GetYaxis()->SetTitleOffset(1.0);
							kdef->SetLineColor(kRed);
							kdelow->SetLineColor(kRed-5);
							kdelow->SetLineStyle(8);
							kdeup->SetLineColor(kRed-5);
							kdeup->SetLineStyle(8);
							kdef->Draw("SAME");
							kdelow->Draw("SAME");
							kdeup->Draw("SAME");

							TLegend *legend = new TLegend(0.7,0.8,0.93,0.97);
							legend->AddEntry(htemp2, "Data");
							legend->AddEntry(kde->GetDrawnFunction(), "KDE");
							legend->AddEntry(kde->GetDrawnLowerFunction(), "KDE - 2 sigma");
							legend->AddEntry(kde->GetDrawnUpperFunction(), "KDE + 2 sigma");

							canvsingleex->SaveAs(out_file.c_str());
						}
					}

					else if (ntemp2 > profile->GetMaximum())
						printf("---->  Overflow pixel (%3d,%3d), number of bins %6d, number of entries %8d, peak at ()%f)\n", i, j, ntemp1, ntemp2, mpv_volts);

					else if (ntemp2 < 0)
						printf("---->  Zero event pixel (%3d,%3d), number of bins %6d, number of entries %8d\n", i, j, ntemp1, ntemp2);

					if (i == 0 && j == 0)
						printf("Using kde method\n");

				}

				// spectrum option
				else if (fspec == true && fkde == false) {
					res = fitGauss(htemp2, fit_scale, lcut_volts, hcut_volts);
					mpv_volts = res->mean;
					mpv_err = res->mean_err;
					chi2ndf = res->chi2/res->ndf;
				}

				// maximum option
				else {
					mpv_volts = htemp2->GetBinCenter(htemp2->GetMaximumBin());
				}

				if (mpv_volts == -1){
					nfailed += 1;
				}

				// fill histogramms
				res_volts.push_back(mpv_volts);
				res_volts.push_back(mpv_err);
				res_volts.push_back(chi2ndf);

				hist_mpv_volts->Fill(mpv_volts);
				hist_mpv_err_volts->Fill(mpv_err);

				if (mpv_err > 0.5*mpv_volts && 0){
					out_file = res_path + dev_id + "/" +  src + "/bad_example_single_volts_gauss_" + i + "_" + j + ".pdf";
					printHist(htemp2, out_file, "Single Pixel Example with Gaussian", "voltage (mV)", "# events", -1.0, -1.0, 1.1);
				}
			}

			else {
				res_volts.push_back(-1);
				res_volts.push_back(-1);
				res_volts.push_back(-1);
				nfailed += 1;
			}

			// push back results
			out_volts.push_back(res_volts);

			if ((i == 2 && j == 2) || (i == 1 && j == 1) || (i == 13 && j == 13)) {
				out_file = res_path + dev_id + "/" +  src + "/example_" + src + "_single_volts_"
									+ std::to_string(i) + "_" + std::to_string(j) + "_" +  dev_id + ".pdf";
				printHist(htemp2, out_file, "Single Pixel Example with Gaussian", "voltage (mV)", "# events", -1.0, -1.0, 1.1);
			}

			if (i == j)
				printf("---->  Progressed to (%3d,%3d), peak at (%.2f), failed so far (%d)\n", i, j, mpv_volts, nfailed);
		}
	}

	styleMap();
	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_profile.pdf";
	printMap(profile, out_file, "map_trim", "column", "row", "# events", 0, 1000, 1.1, 1.2);

	styleGraph();
	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_hist_cal_mpv_volts.pdf";
	printHist(hist_mpv_volts, out_file, "Single Voltage Distribution", "voltage (mV)", "# events", -1.0, -1.0, 1.4);
	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_hist_cal_mpv_err_volts.pdf";
	printHist(hist_mpv_err_volts, out_file, "Single Voltage Distribution", "voltage (mV)", "# events", -1.0, -1.0, 1.4);

	// save dat files
	if (fsave == true){
		std::string hd = "# column | row | mpv (mV) | mpv_err (mV) | chi2/ndf (-) ";
		out_file = dat_path + "lab/"  + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac)
			+ "/" + dev_id + "_" + src + ".txt";
		write_file(out_file, ' ', out_volts, hd, info);
	}

	delete map_fullvolts;
}





void acq_analysis(std::string dev_id, int thr_dac, int ik_dac, std::string src)
{

	gROOT->Reset();
	gROOT->SetBatch(kTRUE);

	styleCommon();
	styleMap();
	styleGraph();



	// Preperations
	// ---------------------------------------

	// initialisation
	int ret, cnt;
	std::vector<std::vector<float>> dat_tot, dat_volts;
	std::string temp_file, out_file;

	Int_t x, y, bin; // place holder for tbranch or histogramme data
	Int_t vmax_ene, vmax_volts, vmax_tot, nbins_ene, nbins_volts, nbins_tot; // histogramme parameters
	Float_t ene, tot, volts, fit_scale, fit_scale_global; // place holder for tbranch data
	Float_t peak, val, lcut_tot, hcut_tot, lcut_volts, hcut_volts, lcut_ene, hcut_ene, lcut_fit, hcut_fit; // values for spectral analysis
	Float_t mpv_tot, mpv_volts, mpv_ene, err_tot, err_volts, err_ene, width_tot, width_volts, width_ene; // most probable values
	Double_t rho; // value for kernel density estimation

	auto res = new fitResult();

	// configure paths
	std::string dat_file = dat_path + "/xraytube/" + dev_id + "/" + src + "/" + src + "_thr_"
			+ std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "_data.root";


	// configure source
	if (src == "fe"){
		peak = 6.40; // 6.40, 6.395
		lcut_tot = 3; hcut_tot = 20; lcut_volts = 50; hcut_volts = 120; lcut_ene = 5; hcut_ene = 8;
		vmax_ene = 12; vmax_volts = 150; vmax_tot = 20;
		nbins_ene = 200; nbins_volts = 150; nbins_tot = 40;
		fit_scale = 0.5; fit_scale_global = 0.7;
	}

	else if (src == "mn"){
		peak = 5.90; // 5.90, 5.958
		lcut_tot = 3; hcut_tot = 20; lcut_volts = 50; hcut_volts = 120; lcut_ene = 4; hcut_ene = 8;
		vmax_ene = 12; vmax_volts = 150; vmax_tot = 20;
		nbins_ene = 200; nbins_volts = 150; nbins_tot = 40;
		fit_scale = 0.5; fit_scale_global = 0.7;
	}

	else if (src == "cu"){
		peak = 8.05; // 8.05, 8.040
		lcut_tot = 10; hcut_tot = 30; lcut_volts = 50; hcut_volts = 150; lcut_ene = 5; hcut_ene = 12;
		vmax_ene = 15; vmax_volts = 200; vmax_tot = 30;
		nbins_ene = 300; nbins_volts = 200; nbins_tot = 120;
		fit_scale = 0.5; fit_scale_global = 0.7;
	}

	else if (src == "in"){
		peak = 24.21; // 24.21, 24.100
		lcut_tot = 35; hcut_tot = 50; lcut_volts = 250; hcut_volts = 400; lcut_ene = 20; hcut_ene = 30;
		vmax_ene = 30; vmax_volts = 500; vmax_tot = 60;
		nbins_ene = 600; nbins_volts = 500; nbins_tot = 240;
		fit_scale = 0.3; fit_scale_global = 0.5;
	}

	else if (src == "cd"){
		peak = 23.17; // 23.17, 23.100
		lcut_tot = 20; hcut_tot = 60; lcut_volts = 250; hcut_volts = 500; lcut_ene = 15; hcut_ene = 30;
		vmax_ene = 30; vmax_volts = 500; vmax_tot = 60;
		nbins_ene = 600; nbins_volts = 500; nbins_tot = 240;
		fit_scale = 0.35; fit_scale_global = 0.7;
	}

	else if (src == "unknown"){
		peak = 8.64; // 8.64
		lcut_tot = 10; hcut_tot = 30; lcut_volts = 50; hcut_volts = 150; lcut_ene = 5; hcut_ene = 12;
		vmax_ene = 15; vmax_volts = 200; vmax_tot = 30;
		nbins_ene = 300; nbins_volts = 200; nbins_tot = 120;
		fit_scale = 0.5; fit_scale_global = 0.8;
	}

	else if (src == "pb"){
		peak = 10.55; // 10.55, 10.5
		lcut_tot = 16; hcut_tot = 20; lcut_volts = 125; hcut_volts = 200; lcut_ene = 9.8; hcut_ene = 12;
		vmax_ene = 20; vmax_volts = 300; vmax_tot = 40;
		nbins_ene = 400; nbins_volts = 300; nbins_tot = 640;
		fit_scale = 0.35; fit_scale_global = 0.85;
	}

	else if (src == "zr"){
		peak = 15.77; // 15.77, 15.746
		lcut_tot = 24; hcut_tot = 40; lcut_volts = 160; hcut_volts = 300; lcut_ene = 10; hcut_ene = 20;
		vmax_ene = 22; vmax_volts = 400; vmax_tot = 50;
		nbins_ene = 440; nbins_volts = 400; nbins_tot = 200;
		fit_scale = 0.45; fit_scale_global = 0.7;
	}

	else {
		peak = 0; lcut_tot = 0; hcut_tot = 0;
		printf("---->  No source given!!\n\n\n\n\n");
	}

	printf("---->  Analyzing data from %s\n", dat_file.c_str());

	// load calibration data
	bool info = true;
	temp_file = dat_path + "/lab/"   + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_tot.txt";
	ret = read_file(temp_file.c_str(), ' ', dat_tot, info);
	temp_file = dat_path + "/lab/"   + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac) + "/" + dev_id + "_cal_volts.txt";
	ret = read_file(temp_file.c_str(), ' ', dat_volts, info);

	// open file, get tree, set adresses and add new branch
	TFile *file = new TFile(dat_file.c_str(), "UPDATE");
	file->ls();

	TTree *singles = (TTree*)file->Get("singleHits");

	TBranch *enebr = singles->Branch("ene", &ene, "ene/F");
	TBranch *voltbr = singles->Branch("volts", &volts, "volts/F");
	TBranch *totbr = singles->GetBranch("tot_hd");
	TBranch *colbr = singles->GetBranch("col");
	TBranch *rowbr = singles->GetBranch("row");

	singles->SetBranchAddress("col", &x);
	singles->SetBranchAddress("row", &y);
	singles->SetBranchAddress("tot_hd", &tot);

	// fill energy branch
	ULong64_t nevnt = singles->GetEntries();
	nevnt = int(0.8*nevnt);
	for (ULong64_t i = 0; i < nevnt; i++) {
		colbr->GetEntry(i);
		rowbr->GetEntry(i);
		totbr->GetEntry(i);

		float params[4] = {dat_tot[256*y + x][2], dat_tot[256*y + x][3], dat_tot[256*y + x][4], dat_tot[256*y + x][5]};
		volts = tot_to_volts(tot, params);

		float params2[2] = {dat_volts[256*y + x][2], dat_volts[256*y + x][3]};
		ene = volts_to_charge(volts, params2);

		voltbr->Fill();
		enebr->Fill();
	}

	// Get global fit values
	cnt = 0;
	std::vector<float> dat_tot_global = {0, 0, 0, 0};
	std::vector<float> dat_volts_global = {0.0659, 1.395};
	for (int i = 0; i < dat_tot.size(); i++){
		if (dat_tot.at(i).at(2) == -1) {
			continue;
		}
		dat_tot_global.at(0) += dat_tot.at(i).at(2);
		dat_tot_global.at(1) += dat_tot.at(i).at(3);
		dat_tot_global.at(2) += dat_tot.at(i).at(4);
		dat_tot_global.at(3) += dat_tot.at(i).at(5);
		cnt++;
	}
	for (int i = 0; i < dat_tot_global.size(); i++){
		dat_tot_global.at(i) /= cnt;
	}



	printf("---->  Starting global analysis\n");

	// create histogramm for tot and energy
	auto hist_tot = new TH1F("hist_tot", "ToT Global", vmax_tot, 0, vmax_tot);
	auto hist_volts = new TH1F("hist_volts", "Volts Global", 100, 0, vmax_volts);
	auto hist_ene = new TH1F("hist_ene", "Energy Global", vmax_ene*10, 0, vmax_ene);
	auto hist_stot = new TH1F("hist_stot", "Single Pixel", 30, 0, vmax_tot);
	auto hist_svolts = new TH1F("hist_svolts", "Single Pixel", 30, 0, vmax_volts);
	auto hist_sene = new TH1F("hist_sene", "Single Pixel", 30, 0, vmax_ene);
	auto hist_ftoa = new TH1F("hist_ftoa", "ToT Global", 16, 0, 1);
	for (ULong64_t i = 0; i < nevnt; i++) {
		rowbr->GetEntry(i);
		colbr->GetEntry(i);
		totbr->GetEntry(i);
		voltbr->GetEntry(i);
		enebr->GetEntry(i);

		hist_tot->Fill(tot+0.5);
		hist_volts->Fill(volts);
		hist_ene->Fill(ene);
		if (x == xpix && y == ypix) {
			hist_stot->Fill(tot);
			hist_svolts->Fill(volts);
			hist_sene->Fill(ene);
			hist_ftoa->Fill(fmod(tot+0.5, 1));
		}
	}


	out_file = res_path + dev_id + "/" +  src + "/" + src + "_global_ftoa_" + dev_id + ".pdf";
	printHist(hist_ftoa, out_file, "FTOA Global", "ftoa (-)", "# events", -1.0, -1.0, 1.1);


	float mpv_tot_global = 10.5;
	std::cout<< "MPV Estimate ToT [ADC]: " << mpv_tot_global << std::endl;
	res = fitGauss(hist_volts, fit_scale_global, lcut_volts, hcut_volts);
	float mpv_volts_global = res->mean;
	std::cout<< "MPV Estimate Volts [mV]: " << mpv_volts_global << std::endl;
	res = fitGauss(hist_ene, fit_scale_global, lcut_ene, hcut_ene);
	float mpv_ene_global = res->mean;
	std::cout<< "MPV Estimate Energy [keV]: " << mpv_ene_global << std::endl;


	// check beam profile
	auto profile = new TH2F("pr", "Event Count Distribution", 256, 0, 256, 256, 0, 256);
	for (ULong64_t i = 0; i < nevnt; i++) {
		colbr->GetEntry(i);
		rowbr->GetEntry(i);
		profile->Fill(x, y);
	}
	printf("---->  Highest event number per pixel is %.0f\n", profile->GetMaximum());
	printf("---->  Lowest event number per pixel is %.0f\n", profile->GetMinimum());


	// check matrix correlations
	Int_t ntemp1, ntemp2, tempcnt;
	TH1D *htemp1, *htemp2, *htemp3;
	auto map_fulltot = new TH3F("map_fulltot", "TOT Full", 256, 0, 256, 256, 0, 256, 30, 0, vmax_tot);
	auto map_fullvolts = new TH3F("map_fullvolts", "voltage Full", 256, 0, 256, 256, 0, 256, 30, 0, vmax_volts);
	auto map_fullene = new TH3F("map_fullene", "Energy Full", 256, 0, 256, 256, 0, 256, 30, 0, vmax_ene);
	for (ULong64_t i = 0; i < nevnt; i++) {
		colbr->GetEntry(i);
		rowbr->GetEntry(i);
		totbr->GetEntry(i);
		voltbr->GetEntry(i);
		enebr->GetEntry(i);
		map_fulltot->Fill(x, y, tot);
		map_fullvolts->Fill(x, y, volts);
		map_fullene->Fill(x, y, ene);
	}

	printf("---->  Starting col/row analysis\n");

	// along columns
	auto hist_coltot = new TH1F("hist_coltot", "TOT MPV along Columns", nbins_tot, 0, vmax_tot);
	auto hist_colvolts = new TH1F("hist_colvolts", "voltsry MPV along Columns", nbins_volts, 0, vmax_volts);
	auto hist_colene = new TH1F("hist_colene", "Enery MPV along Columns", nbins_ene, 0, vmax_ene);
	for (int i = 0; i < 256; i++) {
		htemp1 = map_fulltot->ProjectionZ("_pz", i, i, 0, 256);
		res = fitGauss(htemp1, fit_scale_global, lcut_tot, hcut_tot);
		mpv_tot = res->mean;

		htemp2 = map_fullvolts->ProjectionZ("_pz", i, i, 0, 256);
		res = fitGauss(htemp2, fit_scale_global, lcut_volts, hcut_volts);
		mpv_volts = res->mean;

		htemp3 = map_fullene->ProjectionZ("_pz", i, i, 0, 256);
		res = fitGauss(htemp3, fit_scale_global, lcut_ene, hcut_ene);
		mpv_ene = res->mean;

		hist_coltot->Fill(mpv_tot);
		hist_colvolts->Fill(mpv_volts);
		hist_colene->Fill(mpv_ene);
	}


	// along rows
	auto hist_rowtot = new TH1F("hist_rowtot", "TOT MPV along Rows", nbins_tot, 0, vmax_tot);
	auto hist_rowvolts = new TH1F("hist_rowvolts", "Volts MPV along Rows", nbins_volts, 0, vmax_volts);
	auto hist_rowene = new TH1F("hist_rowene", "Enery MPV along Rows", nbins_ene, 0, vmax_ene);
	for (int j = 0; j < 256; j++) {
		htemp1 = map_fulltot->ProjectionZ("_pz", 0, 256, j, j);
		res = fitGauss(htemp1, fit_scale_global, lcut_tot, hcut_tot);
		mpv_tot = res->mean;

		htemp2 = map_fullvolts->ProjectionZ("_pz", 0, 256, j, j);
		res = fitGauss(htemp2, fit_scale_global, lcut_volts, hcut_volts);
		mpv_volts = res->mean;

		htemp3 = map_fullene->ProjectionZ("_pz", 0, 256, j, j);
		res = fitGauss(htemp3, fit_scale_global, lcut_ene, hcut_ene);
		mpv_ene = res->mean;

		hist_rowtot->Fill(mpv_tot);
		hist_rowvolts->Fill(mpv_volts);
		hist_rowene->Fill(mpv_ene);
	}

	printf("---->  Starting matrix analysis\n");

	// across matrix
	auto hist_mpv_tot = new TH1F("hist_mpv_tot", "TOT MPV", nbins_tot, 0, vmax_tot);
	auto hist_mpv_volts = new TH1F("hist_mpv_volts", "voltage MPV", nbins_volts, 0, vmax_volts);
	auto hist_mpv_ene = new TH1F("hist_mpv_ene", "Energy MPV", nbins_ene, 0, vmax_ene);
	auto hist_width_tot = new TH1F("hist_width_tot", "ToT Peak Width", nbins_tot, 0, vmax_tot);
	auto hist_width_volts = new TH1F("hist_width_volts", "voltage Peak Width", nbins_volts, 0, vmax_volts);
	auto hist_width_ene = new TH1F("hist_width_ene", "Energy Peak Width", nbins_ene, 0, vmax_ene);

	auto hist_mpv_voltscol = new TH1F("hist_mpv_voltscol", "Voltage MPV along one Column", nbins_volts, 0, vmax_volts);
	auto hist_mpv_voltsrow = new TH1F("hist_mpv_voltsrow", "Voltage MPV along one Row", nbins_volts, 0, vmax_volts);

	auto hist_mpv_enecol = new TH1F("hist_mpv_enecol", "Energy MPV along one Column", nbins_ene, 0, vmax_ene);
	auto hist_mpv_enerow = new TH1F("hist_mpv_enerow", "Energy MPV along one Row", nbins_ene, 0, vmax_ene);

	auto map_mpv_tot = new TH2F("map_mpv_tot", "TOT MPV", 256, 0, 256, 256, 0, 256);
	auto map_mpv_volts = new TH2F("map_mpv_volts", "voltage MPV in mV", 256, 0, 256, 256, 0, 256);
	auto map_mpv_ene = new TH2F("map_mpv_ene", "Energy MPV in mV", 256, 0, 256, 256, 0, 256);

	auto map_width_tot = new TH2F("map_width_tot", "Energy Peak Width in ToT", 256, 0, 256, 256, 0, 256);
	auto map_width_volts = new TH2F("map_width_volts", "voltage Peak Width in mV", 256, 0, 256, 256, 0, 256);
	auto map_width_ene = new TH2F("map_width_ene", "Energy Peak Width in mV", 256, 0, 256, 256, 0, 256);

	auto map_mpv_ene_block = new TH2F("hist_mpv_ene_block", "Energy MPV along 16x2 sub-matrix", 2, 0, 2, 16, 0, 16);
	auto map_rmsene_block = new TH2F("hist_rmsene_block", "Energy MPV RMS along 16x2 sub-matrix", 2, 0, 2, 16, 0, 16);
	auto map_widthene_block = new TH2F("hist_widthene_block", "Energy Peak Width along 16x2 sub-matrix", 2, 0, 2, 16, 0, 16);

	auto scattercor_volts  = new TH2F("scattercor_volts", "MPV TOT vs. Voltage", nbins_tot, 0, vmax_tot, nbins_volts, 0, vmax_volts);
	auto scattercor_ene  = new TH2F("scattercor_ene", "MPV TOT vs. Energy", nbins_tot, 0, vmax_tot, nbins_ene, 0, vmax_ene);

	std::vector<float> res_volts;
	std::vector<std::vector<float>> out_volts;
	for (int j = 0; j < 256; j++) {
		for (int i = 0; i < 256; i++) {

			// clear vectors
			res_volts.clear();

			res_volts.push_back(i);
			res_volts.push_back(j);

			// reset
			mpv_tot = 0;
			mpv_ene = 0;
			mpv_volts = 0;
			err_tot = 0;
			err_ene = 0;
			err_volts = 0;
			width_tot = 0;
			width_ene = 0;
			width_volts = 0;

			// fill tot histogramms
			// 0 bin is underflow, n is nth bin, n+1 is overflow
			htemp1 = map_fulltot->ProjectionZ("_pz", i+1, i+1, j+1, j+1);

			// fill voltage histogramms
			htemp2 = map_fullvolts->ProjectionZ("_pz", i+1, i+1, j+1, j+1);

			// fill energy histogramms
			htemp3 = map_fullene->ProjectionZ("_pz", i+1, i+1, j+1, j+1);

			// spectrum option
			if (fspec == true && fkde == false) {
				res = fitGauss(htemp1, fit_scale, lcut_tot, hcut_tot);
				mpv_tot = res->mean;
				err_tot = res->mean_err;
				width_tot = res->std;

				res = fitGauss(htemp2, fit_scale, lcut_volts, hcut_volts);
				mpv_volts = res->mean;
				err_volts = res->mean_err;
				width_volts = res->std;

				res = fitGauss(htemp3, fit_scale, lcut_ene, hcut_ene);
				mpv_ene = res->mean;
				err_ene = res->mean_err;
				width_ene = res->std;
			}

			// maximum option
			else {
				mpv_tot = htemp1->GetBinCenter(htemp1->GetMaximumBin());
				mpv_volts = htemp2->GetBinCenter(htemp2->GetMaximumBin());
				mpv_ene = htemp3->GetBinCenter(htemp3->GetMaximumBin());
			}

			if (profile->GetBinContent(i, j) > 1000){

				// fill histogramms
				hist_mpv_tot->Fill(mpv_tot);
				hist_mpv_volts->Fill(mpv_volts);
				hist_mpv_ene->Fill(mpv_ene);

				hist_width_tot->Fill(width_tot);
				hist_width_volts->Fill(width_volts);
				hist_width_ene->Fill(width_ene);

				map_mpv_tot->Fill(i, j, mpv_tot);
				map_mpv_volts->Fill(i, j, mpv_volts);
				map_mpv_ene->Fill(i, j, mpv_ene);

				map_width_tot->Fill(i, j, width_tot);
				map_width_volts->Fill(i, j, width_volts);
				map_width_ene->Fill(i, j, width_ene);

				scattercor_volts->Fill(mpv_tot, mpv_volts);
				scattercor_ene->Fill(mpv_tot, mpv_ene);

				if (i == xpix) {
					hist_mpv_voltscol->Fill(mpv_volts);
					hist_mpv_enecol->Fill(mpv_ene);
				}
				if (j == ypix) {
					hist_mpv_voltsrow->Fill(mpv_volts);
					hist_mpv_enerow->Fill(mpv_ene);
				}

				// fill histogramms
				res_volts.push_back(mpv_tot);
				res_volts.push_back(err_tot);
				res_volts.push_back(width_tot);
				res_volts.push_back(mpv_volts);
				res_volts.push_back(err_volts);
				res_volts.push_back(width_volts);
				res_volts.push_back(mpv_ene);
				res_volts.push_back(err_ene);
				res_volts.push_back(width_ene);
			}

			else {
				res_volts.push_back(-1);
				res_volts.push_back(-1);
				res_volts.push_back(-1);
				res_volts.push_back(-1);
				res_volts.push_back(-1);
				res_volts.push_back(-1);
				res_volts.push_back(-1);
				res_volts.push_back(-1);
				res_volts.push_back(-1);
			}

			// push back results
			out_volts.push_back(res_volts);

			if (i == xpix && j == ypix){
				out_file = res_path + dev_id + "/" +  src + "/" + src + "_single_volts_gauss_" + dev_id + ".pdf";
				printHist(htemp2, out_file, "Single Pixel Example with Gaussian", "voltage (mV)", "# events", -1.0, -1.0, 1.1);

				out_file = res_path + dev_id + "/" +  src + "/" + src + "_single_energy_gauss_" + dev_id + ".pdf";
				printHist(htemp3, out_file, "Single Pixel Example with Gaussian", "energy (keV)", "# events", -1.0, -1.0, 1.1);
			}

			if (i == j)
				printf("---->  Progressed to (%3d,%3d), peak at (%.2f)\n", i, j, mpv_ene);
		}
	}



	// Printing
	// ---------------------------------------

	styleGraph();

	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_single_tot.pdf";
	printHist(hist_stot, out_file, "Singlel ToT Distribution", "time over threshold (adc)", "# events", 0, vmax_tot, 1.3);
	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_single_volts.pdf";
	printHist(hist_svolts, out_file, "Single Voltage Distribution", "voltage (mV)", "# events", 0, vmax_volts, 1.3);
	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_single_energy.pdf";
	printHist(hist_sene, out_file, "Single Energy Distribution", "energy (keV)", "# events", 0, vmax_ene, 1.3);

	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_global_tot.pdf";
	printHist(hist_tot, out_file, "Global ToT Distribution", "time over threshold (adc)", "# events", 0, vmax_tot, 1.3);
	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_global_voltage.pdf";
	printHist(hist_volts, out_file, "Global Voltage Distribution", "voltage (mV)", "# events", 0, vmax_volts, 1.3);
	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_global_energy.pdf";
	printHist(hist_ene, out_file, "Global Energy Distribution", "energy (keV)", "# events", 0, vmax_ene, 1.3);

	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_hist_mpv_tot.pdf";
	printHist(hist_mpv_tot, out_file, "Singlel ToT Distribution", "time over threshold (adc)", "# events", -1.0, -1.0, 1.4);
	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_hist_mpv_volts.pdf";
	printHist(hist_mpv_volts, out_file, "Single Voltage Distribution", "voltage (mV)", "# events", -1.0, -1.0, 1.4);
	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_hist_mpv_energy.pdf";
	printHist(hist_mpv_ene, out_file, "Single Energy Distribution", "energy (keV)", "# events", -1.0, -1.0, 1.4);

	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_hist_width_energy.pdf";
	printHist(hist_width_ene, out_file, "map_trim", "energy (keV)", "# entries", -1.0, -1.0, 1.4);
	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_hist_width_volts.pdf";
	printHist(hist_width_volts, out_file, "map_trim", "voltage (mV)", "# entries", -1.0, -1.0, 1.4);
	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_hist_width_tot.pdf";
	printHist(hist_width_tot, out_file, "map_trim", "time over threshold (adc)", "# entries", -1.0, -1.0, 1.4);

	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_hist_mpv_single_row_volts.pdf";
	printHist(hist_mpv_voltsrow, out_file, "Single Voltage Distribution", "voltage (mV)", "# events", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_hist_mpv_single_col_volts.pdf";
	printHist(hist_mpv_voltscol, out_file, "Single Voltage Distribution", "voltage (mV)", "# events", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_hist_mpv_single_col_energy.pdf";
	printHist(hist_mpv_enecol, out_file, "Single Energy Distribution", "energy (keV)", "# events", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_hist_mpv_single_row_energy.pdf";
	printHist(hist_mpv_enerow, out_file, "Single Energy Distribution", "energy (keV)", "# events", -1.0, -1.0, 1.3);

	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_hist_tot_mpv_col.pdf";
	printHist(hist_coltot, out_file, "Single ToT Distribution", "time over threshold (adc)", "# events", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_hist_volts_mpv_col.pdf";
	printHist(hist_colvolts, out_file, "Single Voltage Distribution", "voltage (mV)", "# events", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_hist_energy_mpv_col.pdf";
	printHist(hist_colene, out_file, "Single Energy Distribution", "energy (keV)", "# events", -1.0, -1.0, 1.3);

	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_hist_tot_mpv_row.pdf";
	printHist(hist_rowtot, out_file, "Singlel ToT Distribution", "time over threshold (adc)", "# events", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_hist_volts_mpv_row.pdf";
	printHist(hist_rowvolts, out_file, "Single Voltage Distribution", "voltage (mV)", "# events", -1.0, -1.0, 1.3);
	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_hist_energy_mpv_row.pdf";
	printHist(hist_rowene, out_file, "Single Energy Distribution", "energy (keV)", "# events", -1.0, -1.0, 1.3);


	styleMap();

	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_profile.pdf";
	printMap(profile, out_file, "map_trim", "column", "row", "# events", 0, 1000, 1.1, 1.2);
	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_tot_vs_volts.pdf";
	printMap(scattercor_volts, out_file, "map_trim", "time over threshold (adc)", "voltage (mV)", "# entries", -1.0, -1.0, 1.1, 1.1);
	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_tot_vs_energy.pdf";
	printMap(scattercor_ene, out_file, "map_trim", "time over threshold (adc)", "energy (keV)", "# entries", -1.0, -1.0, 1.1, 1.1);

	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_map_mpv_tot.pdf";
	printMap(map_mpv_tot, out_file, "map_trim", "column", "row", "time over threshold (adc)", 0.9*mpv_tot_global, 1.1*mpv_tot_global, 1.1, 1.1);
	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_map_mpv_volts.pdf";
	printMap(map_mpv_volts, out_file, "map_trim", "column", "row", "voltage (mV)", 0.9*mpv_volts_global, 1.1*mpv_volts_global, 1.1, 1.1);
	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_map_mpv_energy.pdf";
	printMap(map_mpv_ene, out_file, "map_trim", "column", "row", "energy (keV)", 0.9*mpv_ene_global, 1.1*mpv_ene_global, 1.1, 1.1);

	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_map_width_energy.pdf";
	printMap(map_width_ene, out_file, "map_trim", "column", "row", "energy (keV)", -1.0, -1.0, 1.1, 1.1);
	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_map_width_volts.pdf";
	printMap(map_width_volts, out_file, "map_trim", "column", "row", "voltage (mV)", -1.0, -1.0, 1.1, 1.1);
	out_file = res_path + dev_id + "/" +  src + "/" + dev_id + "_" + src + "_map_width_tot.pdf";
	printMap(map_width_tot, out_file, "map_trim", "column", "row", "time over threshold (adc)", -1.0, -1.0, 1.1, 1.1);

	// save dat files
	if (fsave == true){
		std::string hd = "# column | row | mpv_tot (tot) | err_tot (tot) | width_tot (tot) | mpv_volts (mV) | err_volts (mV) | width_volts (mV) | mpv_ene (keV) | 	err_ene (keV) | width_ene (keV) ";
		out_file = dat_path + "lab/"  + dev_id + "/" + "cal_thr_" + std::to_string(thr_dac) + "_ik_" + std::to_string(ik_dac)
			+ "/" + dev_id + "_" + src + "_corrected.txt";
		write_file(out_file, ' ', out_volts, hd, info);
	}

	delete map_fulltot;
	delete map_fullvolts;
	delete map_fullene;
	delete res;
}
