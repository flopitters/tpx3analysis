/********************************************************************
* File: common.cc
* ------------------------
*
* Description:
* Common functions.
*
* Version:
* Author: Florian Pitters
*
*******************************************************************/

#include "common.h"


void print_line(std::string line) {
	std::cout << " ---> " << line << std::endl;
}

// -----------------------------------------------------------------------------
// Mathematical functions
// -----------------------------------------------------------------------------

/* Gauss function */
float gauss(float *x, float *par)
{
	float A, mu, sigma;
	static float g;

	A = par[0];
	mu = par[1];
	sigma = par[2];

    g = A * exp(-pow(x[0]-mu, 2) / (2. * pow(sigma, 2)));

	return g;
}


/* Inverse surrogate function to convert tot to charge */
float tot_to_volts(float tot, float *par)
{
	float a, b, c, t;
	static float q;

	a = par[0];
	b = par[1];
	c = par[2];
	t = par[3];

	q = (t*a + tot - b + pow(pow(b+t*a-tot, 2) + 4*a*c, 0.5)) / (2*a);

	return q;
}


/* Inverse surrogate function to convert tot to charge */
float toa_to_volts(float toa, float *par)
{
	float a, b, c, t;
	static float q;

	a = par[0];
	b = par[1];
	c = par[2];
	t = par[3];

	q = (t * a + toa - b + pow(pow(b+t*a-toa, 2) + 4*a*c, 0.5)) / (2 * a);

	return q;
}

/* Surrogate function to convert charge to tot */
float volts_to_tot(float v, float *par)
{
    float a, b, c, t;
	static float tot;

	a = par[0];
	b = par[1];
	c = par[2];
	t = par[3];

	tot = a * v + b - c/(v-t);

    return tot;
}

/* Timewalk function to concert charge to toa */
float volts_to_toa(float v, float *par)
{
    float c, t, d;
	static float toa;

	c = par[0];
	t = par[1];
	d = par[2];

	toa = c/(v-t) + d;

    return toa;
}


/* Linear function to convert volts to charge */
float volts_to_charge(float volt, float *par)
{
	float a, b;
	static float q;

	a = par[0];
	b = par[1];

	q = a * volt + b;

	return q;
}

float mean(std::vector<float> vals)
{
	float m = 0;
	float n = 0;
	for (int i = 0; i < vals.size(); i++){
		if (!isnan(vals.at(i))){
			m += vals.at(i);
			n += 1;
		}
	}
	m /= n;

	return m;
}

float var(std::vector<float> vals)
{
	float v = 0;
	float m = 0;
	float n = 0;
	for (int i = 0; i < vals.size(); i++){
		if (!isnan(vals.at(i))){
			m += vals.at(i);
			n += 1;
		}
	}
	m /= float(n);

	for (int i = 0; i < vals.size(); i++){
		if (!isnan(vals.at(i))){
			v += pow((m - vals.at(i)), 2);
		}
	}

	v /= float(n);

	return v;
}

float median(std::vector<float> vals)
{
	float m = 0;

	std::sort(vals.begin(), vals.end());
	if (vals.size() % 2 == 0) {
  		m = (vals.at(vals.size()/2 - 1) + vals.at(vals.size()/2) / 2.);
	}
	else {
  		m = vals.at(vals.size()/2);
	}

	return m;
}

float mean_along_column(std::vector<std::vector<float>> vals, int col)
{
	float m = 0;
	float n = 0;

	for (int i = 0; i < vals.size(); i++){
		if (!isnan(vals.at(i).at(col)) && vals.at(i).at(col) != -1) {
			m += vals.at(i).at(col);
			n += 1;
		}
	}
	m /= n;

	return m;
}

float var_along_column(std::vector<std::vector<float>> vals, int col)
{
	float v = 0;
	float m = 0;
	float n = 0;
	for (int i = 0; i < vals.size(); i++){
		if (!isnan(vals.at(i).at(col)) && vals.at(i).at(col) != -1) {
			m += vals.at(i).at(col);
			n += 1;
		}
	}
	m /= float(n);

	for (int i = 0; i < vals.size(); i++){
		if (!isnan(vals.at(i).at(col))){
			v += pow((m - vals.at(i).at(col)), 2);
		}
	}

	v /= float(n);

	return v;
}


int find_first_above_x(std::vector<float> vals, float thr, int reverse)
{
	int it = 0;
	if (reverse == 0 || reverse == -1){
		for (int i = 0; i < vals.size(); i++){
			it++;
			if (vals.at(i) > thr){ break; }
		}
	}
	else {
		for (int i = vals.size()-2; i > 0; i--){
			//std::cout << it << " " << i << " " << vals.at(i) << std::endl;
			it++;
			if (vals.at(i) > thr){ break; }
		}
		it = vals.size()-1 - it;
	}

	return it;
}

int find_extremum(std::vector<float> vals, std::string option)
{
	int it = 0;
	float vmax = 0;
	for (int i = 0; i < vals.size(); i++){
		if (option == "max"){
			if (vals.at(i) > vmax) { vmax = vals.at(i); it = i; }
		}
		else if (option == "min"){
			if (vals.at(i) < vmax) { vmax = vals.at(i); it = i; }
		}
		else {
			std::cout << "Not a valid option." << std::endl;
			return -1;
		}
	}

	return it;
}



// -----------------------------------------------------------------------------
// Reading functions
// -----------------------------------------------------------------------------

/* Read data */
int read_data_mat(const char *file_name, float data[ncol][nrow], bool info = false)
{
	int i, j, ret;

	FILE *file;
	file = fopen(file_name, "r");
	if (file == NULL) {
		std::cout << "Cannot open input file: " << file_name << std::endl;
		return -1;
	}
	rewind(file);


	// Reading if each line is of format 'value1 value2 ...'
	for (i=0; i<ncol; ++i){
		for (j=0; j<nrow; ++j)
			ret = fscanf(file, "%f", &data[i][j]);
	};

 	fclose(file);

	if (info == true){
		printf("%8.5lf\n", data[83][142]);
		printf("%8.5lf\n", data[13][182]);
		printf("%8.5lf\n", data[77][44]);
	}

	return 0;
}


int read_data_list(const char *file_name, float data[][ncol*nrow], int length, bool info = false)
{
	int i, j, ret;

	FILE *file;
	file = fopen(file_name, "r");
	if (file == NULL) {
		std::cout << "Cannot open input file: " << file_name << std::endl;
		return -1;
	}
	rewind(file);

	// Dumping header
	char str[80];
	fscanf(file, "%[^\n]\n", str);

	// Reading if each line is of format 'value1 value2 ...'
	for (i=0; i<ncol*nrow; ++i){
		for(j=0; j<length; ++j)
			ret = fscanf(file, "%f", &data[j][i]);
	};

 	fclose(file);

	if (info == true){
		printf("%f\n", data[0][100]);
		printf("%f\n", data[1][100]);
		printf("%f\n", data[2][100]);
		printf("%f\n", data[3][100]);
	}

	return 0;
}

/* Write data */
int write_data_mat(const char *file_name, float data[ncol][nrow], bool info = false)
{
	int i, j;

	FILE *file;
	file = fopen(file_name, "w");
	rewind(file);

	// Reading if each line is of format 'value1 value2 ...'
	for (j=0; j<nrow; ++j) {
		for (i=0; i<ncol; ++i) {
			fprintf(file, "%3.4lf ", data[j][i]);
		}
		fprintf(file, "\n");
	}

 	fclose(file);

	if (info == true){
		printf("%8.5lf\n", data[83][142]);
		printf("%8.5lf\n", data[13][182]);
		printf("%8.5lf\n", data[77][44]);
	}

	return 0;
}

int write_data_list(const char *file_name, float data[][ncol*nrow], int length, bool info = false)
{
	int i, j;

	FILE *file;
	file = fopen(file_name, "w");
	rewind(file);

	// Write header
	//char header[80] = "# test";
	//fprintf("%s", header);

	// Write each line of format 'value1 value2 ...'
	for (i=0; i<ncol*nrow; ++i){
		for(j=0; j<length; ++j)
			fprintf(file, "%5.3lf ", data[j][i]);
		fprintf(file, "\n");
	};

 	fclose(file);

	if (info == true){
		printf("%f\n", data[0][100]);
		printf("%f\n", data[1][100]);
		printf("%f\n", data[2][100]);
		printf("%f\n", data[3][100]);
	}

	return 0;
}


int read_data_mat_to_tree(TTree *tr, const char *fn, int ncol, bool info = false)
{
	int i, j, x, y;
	float val;

	TBranch *colbr = tr->Branch("col", &x, "col/I");
	TBranch *rowbr = tr->Branch("row", &y, "row/I");
	TBranch *valbr = tr->Branch("val", &val, "val/F");

	std::ifstream f;
	f.open(fn);

	if (!f.is_open()) {
		std::cout << "Cannot open input file: " << fn << std::endl;
		return -1;
	}

	i = 0;
	j = 0;
	val = 0;

	std::string line;
	while (getline(f, line)) {

		// Ignore lines which contain a # symbol
		if (line.find("#") != std::string::npos)
			continue;

		// Ignore empty lines
		if (strlen(line.c_str()) == 0)
			continue;

		std::stringstream iss(line);

		std::string word;
		while (iss >> word) {
			x = i;
			y = j;
			val = atof(word.c_str());

			tr->Fill();

			j++;
			if (j == ncol) {
				i += 1;
				j = 0;
			}
		}
	}

	f.close();

	return 0;
}


int read_file(std::string fn, char delim, std::vector<std::vector<float>> &dat, bool info)
{
	std::ifstream f;
	f.open(fn);
	dat.clear();

	// check if file is open
	if (!f.is_open()) {
		std::cout << "Cannot open input file:\n\t" << fn << std::endl;
		return 1;
	}
	else {
		if (info == 1) {
			std::cout << "Reading data from file:\n\t" << fn << std::endl;
		}
	}

	// read file line by line
	int i = 0;
	std::string line;
	while (!f.eof()) {
	 	std::getline(f, line);

		// check if line is empty or a comment
		// if not write to output vector
		if (line.size() > 0 && isdigit(line.at(0))) {
			std::stringstream ss(line);
			std::string word;
			std::vector<float> row;
			while (std::getline(ss, word, delim)) {
				i += 1;
				row.push_back(stof(word));
			}
			dat.push_back(row);
		}
    }

	// debug info
	if (info == 1) {
		std::cout << "Number of elements found: " << i << " Number of Rows found: " << dat.size() << std::endl;
	}

	f.close();

	return 0;
}


int read_gain_data(std::string fn, char delim, std::vector<std::vector<float>> &vals, std::vector<std::vector<float>> &codes, bool info)
{
	int ret, cnt;
	std::ostringstream file_name;
	std::vector<float> c, v;
	std::vector<std::vector<float>> tmp;

	for (int i = 0; i < 256; i++){
		file_name << fn << "/";
		file_name << std::setw(3) << std::setfill('0') << std::to_string(i);
		file_name << "_";
		file_name << std::setw(3) << std::setfill('0') << std::to_string(i);
		file_name << ".dat";
		ret = read_file(file_name.str(), delim, tmp, info);
		if (ret == 0){
			for (auto element : tmp) {
				c.push_back(element.at(0));
				v.push_back(element.at(1));
			}
		} else {
			c.push_back(-1);
			v.push_back(-1);
		}
		codes.push_back(c);
		vals.push_back(v);
		c.clear();
		v.clear();
		tmp.clear();
		file_name.str(std::string());
	}

	return 0;
}


int read_noise_data(std::string fn, char delim, std::vector<std::vector<float>> &vals, std::vector<std::vector<float>> &codes, bool info)
{
	int ret, cnt;
	std::ostringstream file_name;
	std::vector<float> c, v;
	std::vector<std::vector<float>> tmp;

	for (int i = 0; i < 256; i++){
		file_name << fn << "/";
		file_name << std::setw(3) << std::setfill('0') << std::to_string(i);
		file_name << "/";
		file_name << std::setw(3) << std::setfill('0') << std::to_string(i);
		file_name << "_";
		file_name << std::setw(3) << std::setfill('0') << std::to_string(i);
		file_name << ".dat";
		ret = read_file(file_name.str(), delim, tmp, info);
		if (ret == 0){
			for (auto element : tmp) {
				c.push_back(element.at(0));
				v.push_back(element.at(1));
			}
		} else {
			c.push_back(-1);
			v.push_back(-1);
		}
		codes.push_back(c);
		vals.push_back(v);
		c.clear();
		v.clear();
		tmp.clear();
		file_name.str(std::string());
	}

	return 0;

}


int write_file(std::string fn, char delim, std::vector<std::vector<float>> &dat, std::string hd, bool info)
{
	std::ofstream f;
	f.open(fn, std::ofstream::trunc);

	// check if file is open
	if (!f.is_open()) {
		std::cout << "Cannot open output file:\n\t" << fn << std::endl;
		exit(1);
	}
	else {
		std::cout << "Writing data to file:\n\t" << fn << std::endl;
	}

	// write header
	if (hd != "") {
		f << hd << std::endl;
	}
	f << std::setprecision(6);

	// write data line by line
	int k = 0;
	int nrow = dat.size();
	int ncol = dat.at(0).size();
	for (int j = 0; j < nrow; j++) {
        for (int i = 0; i < ncol; i++) {
            f << dat.at(j).at(i) << delim;
			k += 1;
        }
		if (j != (nrow-1)){
			f << std::endl;
		}
    }

	// debug info
	if (info == 1) {
		std::cout << "Number of elements wrote: " << k << " Number of Rows wrote: " << dat.size() << std::endl;
	}

	f.close();

	return 0;
}


fitResult *fitGauss(TH1 *hist1d, float fit_scale = 0.9, float low = 9, float high = 11)
{
	int mbin, lbin, nbins;
	float mval, mpos;

	TF1 *func = new TF1("func", "gaus", low, high);
	func->SetNpx(400);

	nbins = hist1d->GetNbinsX();
	lbin = hist1d->GetBinCenter(nbins);
	hist1d->GetXaxis()->SetRangeUser(low, high);

	mbin = hist1d->GetMaximumBin();
	mpos = hist1d->GetBinCenter(hist1d->GetMaximumBin());
	mval = hist1d->GetMaximum();

	// find last bin above fit_scale * max
	int last_bin_above = mbin;
    for(int i = 1; i < (nbins - mbin); i++) {
		if(hist1d->GetBinContent(mbin+i) < fit_scale*mval) {
            break;
        }
		last_bin_above = mbin+i;
		if (hist1d->GetBinCenter(last_bin_above) > high){
			last_bin_above -= 1;
			break;
		}
    }
	//last_bin_above += 1; // Add anoter bin

	// find first value above fit_scale * max
	int first_bin_above = mbin;
    for(int i = 1; i < mbin; i++) {
		if(hist1d->GetBinContent(mbin-i) < fit_scale*mval) {
            break;
        }
		first_bin_above = mbin-i;
		if (hist1d->GetBinCenter(first_bin_above) < low){
			first_bin_above += 1;
			break;
		}
    }

	// std::cout << "\tHi: " << high << "\tFit Hi: " << hist1d->GetBinCenter(last_bin_above) << std::endl;
	// std::cout << "\tLo: " << low << "\tFit Lo: " << hist1d->GetBinCenter(first_bin_above) << std::endl;
	// std::cout << "\tMax Val: " << mval << "\tMax Pos: " << mpos << "\tMax Bin: " << mbin << std::endl;

	fitResult *res = new fitResult();
	if (mval == 0){
		res->mean = -1;
		res->mean_err = -1;
		res->std = -1;
		res->std_err = -1;
		res->chi2 = -1;
		res->ndf = -1;
	}
	else{
		func->SetRange(hist1d->GetBinCenter(first_bin_above), hist1d->GetBinCenter(last_bin_above));
		func->SetParameter(1, mpos);
		hist1d->Fit("func", "qrp");

		res->mean = func->GetParameter(1);
		res->mean_err = func->GetParError(1);
		res->std = func->GetParameter(2);
		res->std_err = func->GetParError(2);
		res->chi2 = func->GetChisquare();
		res->ndf = func->GetNDF();
	}

	hist1d->GetXaxis()->SetRangeUser(0, lbin);

	return res;
}




// ------------------------------------------------------------------------------
// Printing functions
// ------------------------------------------------------------------------------

void printMap(TH2* hist2d, std::string fn, std::string title, std::string titleX, std::string titleY, std::string titleZ, float limX, float limY, float offsetY,
              float offsetZ, std::string option) {
    TCanvas* canv_tmp = new TCanvas(fn.c_str(), fn.c_str());
    hist2d->Draw(option.c_str());
    hist2d->SetTitle(title.c_str());
    hist2d->GetXaxis()->SetTitle(titleX.c_str());
    hist2d->GetYaxis()->SetTitle(titleY.c_str());
    hist2d->GetZaxis()->SetTitle(titleZ.c_str());
    hist2d->GetXaxis()->SetTitleOffset(1.1);
    hist2d->GetYaxis()->SetTitleOffset(offsetY);
    hist2d->GetZaxis()->SetTitleOffset(offsetZ);
    hist2d->SetMarkerSize(0.1);
    if(limX != -1) hist2d->SetMinimum(limX);
    if(limY != -1) hist2d->SetMaximum(limY);
	// if(limX != -1) hist2d->GetXaxis()->SetRangeUser(150, 256);
	// if(limX != -1) hist2d->GetYaxis()->SetRangeUser(70, 256);
    canv_tmp->SaveAs(fn.c_str());
	delete canv_tmp;
}

void printHist(TH1* hist1d, std::string fn, std::string title, std::string titleX, std::string titleY, float limX, float limY, float offset, std::string option) {
    TCanvas* canv_tmp = new TCanvas(fn.c_str(), fn.c_str());
    hist1d->Draw(option.c_str());
    hist1d->SetTitle(title.c_str());
    hist1d->GetXaxis()->SetTitle(titleX.c_str());
    hist1d->GetYaxis()->SetTitle(titleY.c_str());
    hist1d->GetXaxis()->SetTitleOffset(1.1);
    hist1d->GetYaxis()->SetTitleOffset(offset);
    if(limX != -1 || limY != -1) hist1d->SetAxisRange(limX, limY);
    canv_tmp->SaveAs(fn.c_str());
	delete canv_tmp;
}

void printGraph(TGraphErrors* gr, std::string fn, std::string title, std::string titleX, std::string titleY, float limX, float limY, float offset, std::string option) {
    TCanvas* canv_tmp = new TCanvas(fn.c_str(), fn.c_str());
    gr->Draw(option.c_str());
    gr->SetTitle(title.c_str());
    gr->GetXaxis()->SetTitle(titleX.c_str());
    gr->GetYaxis()->SetTitle(titleY.c_str());
    gr->GetYaxis()->SetTitleOffset(offset);
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.8);
	gr->SetFillColor(kOneBlue);
   	gr->SetFillStyle(3001);
    if(limX != -1.0) gr->SetMinimum(limX);
    if(limY != -1.0) gr->SetMaximum(limY);
    canv_tmp->SaveAs(fn.c_str());
	delete canv_tmp;
}

void printGraph(TGraph* gr, std::string fn, std::string title, std::string titleX, std::string titleY, float limX, float limY, float offset, std::string option) {
    TCanvas* canv_tmp = new TCanvas(fn.c_str(), fn.c_str());
    gr->Draw(option.c_str());
    gr->SetTitle(title.c_str());
    gr->GetXaxis()->SetTitle(titleX.c_str());
    gr->GetYaxis()->SetTitle(titleY.c_str());
    gr->GetYaxis()->SetTitleOffset(offset);
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.8);
    if(limX != -1.0) gr->SetMinimum(limX);
    if(limY != -1.0) gr->SetMaximum(limY);
    canv_tmp->SaveAs(fn.c_str());
	delete canv_tmp;
}

void addGraph(TMultiGraph* mg, TLegend* lg, int nentries, float* x, float* y,
        Color_t c, Style_t m, Size_t ms, Style_t ls, std::string style, std::string label) {
    TGraph* gr = new TGraph(nentries, x, y);
    gr->SetLineColor(c);
    gr->SetLineStyle(ls);
    gr->SetMarkerColor(c);
    gr->SetMarkerStyle(m);
    gr->SetMarkerSize(ms);
    mg->Add(gr);
    if (label != "") lg->AddEntry(gr, label.c_str(), style.c_str());
}

void addErrGraph(TMultiGraph* mg, TLegend* lg, int nentries, float* x, float* y, float* x_err, float* y_err,
        Color_t c, Style_t m, Size_t ms, Style_t ls, std::string style, std::string label) {
    TGraphErrors* gr = new TGraphErrors(nentries, x, y, x_err, y_err);
    gr->SetLineColor(c);
    gr->SetLineStyle(ls);
    gr->SetMarkerColor(c);
    gr->SetMarkerStyle(m);
    gr->SetMarkerSize(ms);
    gr->SetFillStyle(3002);
    gr->SetFillColor(c);
    mg->Add(gr);
    lg->AddEntry(gr, label.c_str(), style.c_str());
}

void addFittedGraph(TMultiGraph* mg, TLegend* lg, int nentries, float* x, float* y,
        Color_t c, Style_t m, Size_t ms, Style_t ls, std::string style, std::string label, std::string func, Float_t fit_low, Float_t fit_up) {
    TGraph* gr = new TGraph(nentries, x, y);
    gr->SetLineColor(c);
    gr->SetLineStyle(ls);
    gr->SetMarkerColor(c);
    gr->SetMarkerStyle(m);
    gr->SetMarkerSize(ms);
    gr->Fit(func.c_str(), "me", "", fit_low, fit_up);
    mg->Add(gr);
    lg->AddEntry(gr, label.c_str(), style.c_str());
}


TFitResultPtr addFittedErrGraph(TMultiGraph* mg, TLegend* lg, int nentries, float* x, float* y, float* x_err, float* y_err,
		Color_t c, Style_t m, Size_t ms, Style_t ls, std::string style, std::string label, std::string func, Float_t fit_low, Float_t fit_up) {

	TGraphErrors* gr = new TGraphErrors(nentries, x, y, x_err, y_err);
    TFitResultPtr fitr = gr->Fit(func.c_str(), "MNS", "", fit_low, fit_up);
    //gr = new TGraphErrors(nentries, x, y, x_err, y_err);
    gr->SetLineColor(c);
    gr->SetLineStyle(ls);
    gr->SetMarkerColor(c);
    gr->SetMarkerStyle(m);
    gr->SetMarkerSize(ms);

    mg->Add(gr);
    lg->AddEntry(gr, label.c_str(), style.c_str());

    return fitr;
}

void addFitResults(TF1* func, TLegend* lg, TFitResultPtr fitr, Color_t c, Style_t ls, std::string label) {
    func->SetParameter(0, fitr->Parameter(0));
    func->SetParameter(1, fitr->Parameter(1));
    func->SetLineColor(c);
    func->SetLineStyle(ls);
    // lg->AddEntry(func, label.c_str(), " ");
    func->Draw("SAME");
}

void addScatterGraph(TLegend* lg, TH2* scatter,
        Color_t c, Style_t m, Size_t ms, std::string style, std::string label) {
    scatter->SetLineColor(c);
    scatter->SetMarkerStyle(m);
    scatter->SetMarkerColor(c);
    scatter->SetMarkerSize(ms);
    scatter->Draw("SAME");
    lg->AddEntry(scatter, label.c_str(), style.c_str());
}

void addHist(THStack* hs, TLegend* lg, TH1F* hist1d,
        Color_t c, Style_t m, Style_t ls, std::string style, std::string label, Float_t scale) {
    hist1d->SetLineColor(c);
    hist1d->SetLineStyle(ls);
    hist1d->SetStats(0);
    hist1d->Scale(scale, "nosw2");
    // hist1d->SetFillColorAlpha(kBlue - 3, 0.3);
    // hist1d->SetFillStyle(1001);
    // hist1d->Draw(style);
    hs->Add(hist1d);
    lg->AddEntry(hist1d, label.c_str(), style.c_str());
}


void addHist(THStack* hs, TLegend* lg, TH1D* hist1d,
        Color_t c, Style_t m, Style_t ls, std::string style, std::string label, Float_t scale) {
    hist1d->SetLineColor(c);
    hist1d->SetLineStyle(ls);
    hist1d->SetStats(0);
    hist1d->Scale(scale, "nosw2");
    // hist1d->SetFillColorAlpha(kBlue - 3, 0.3);
    // hist1d->SetFillStyle(1001);
    // hist1d->Draw(style);
    hs->Add(hist1d);
    lg->AddEntry(hist1d, label.c_str(), style.c_str());
}

void addHist(THStack* hs, TLegend* lg, TH2F* hist2d,
        Color_t c, Style_t m, Style_t ls, std::string style, std::string label, Float_t scale) {
    hist2d->SetLineColor(c);
    hist2d->SetLineStyle(ls);
    hist2d->SetStats(0);
    //hist2d->Scale(scale, "nosw2");
    // hist1d->SetFillColorAlpha(kBlue - 3, 0.3);
    // hist1d->SetFillStyle(1001);
    // hist1d->Draw(style);
    hs->Add(hist2d);
    lg->AddEntry(hist2d, label.c_str(), style.c_str());
}


TGraphErrors* getProjectionFrom2DAlongY(TH2F* scatter_plot, int step, double vmin, double vmax, int nmin, std::string option) {
    std::vector<double> valX;
    std::vector<double> errX;
    std::vector<double> valY;
    std::vector<double> errY;

	double mval, mbin, lo, hi;

	double x = 0;
	int nbins = scatter_plot->GetNbinsX();
	double vbin = scatter_plot->GetXaxis()->GetBinWidth(1);

	TF1* gaussf = new TF1("gaussf", "gaus", -100, +100);
    gaussf->SetNpx(400);

    TH1D* profileX = new TH1D();
    for(int i = 1; i < nbins+1; i = i+step) {
        profileX = scatter_plot->ProjectionY("_y", i, i+step-1, "");
		mval = profileX->GetMaximum();
		mbin = profileX->GetBinCenter(profileX->GetMaximumBin());
		lo = profileX->GetBinCenter(profileX->FindFirstBinAbove(0.05 * mval));
		hi = profileX->GetBinCenter(profileX->FindLastBinAbove(0.05 * mval));
		profileX->SetAxisRange(lo-5, hi+5);
		x = i * vbin + (vbin*step) / 2;
        if(profileX->GetEntries() > nmin && x > vmin && x < vmax) {
			valX.push_back(x);
			errX.push_back(0);

			// get max
            if(option == "0") {
				valY.push_back(profileX->GetMean());
				errY.push_back(profileX->GetMeanError());
            }

			if(option == "1") {
				valY.push_back(profileX->GetRMS());
				errY.push_back(profileX->GetRMSError());
			}

            // get gaussian mean
            if(option == "2") {
				gaussf->SetParameter(1, mbin);
                profileX->Fit("gaussf", "qre");
				valY.push_back(gaussf->GetParameter(1));
				errY.push_back(gaussf->GetParError(1));
            }

			// get gaussian mean
            if(option == "3") {
				gaussf->SetParameter(1, mbin);
                profileX->Fit("gaussf", "qre");
				valY.push_back(gaussf->GetParameter(2));
				errY.push_back(gaussf->GetParError(2));
            }
        }
    }
    TGraphErrors* graph = new TGraphErrors(valX.size(), &(valX[0]), &(valY[0]), &(errX[0]), &(errY[0]));

    return graph;
}

TH2F* getMapFrom3DAlongZ(TH3F* scatter_plot, int nbinsX, double vminX, double vmaxX, int nbinsY,
                                double vminY, double vmaxY, int nmin, std::string option) {

	std::vector<double> valX;
	std::vector<double> errX;
	std::vector<double> valY;
	std::vector<double> errY;

	double x, y;
	double mval, mbin, lo, hi;

    auto hist2d = new TH2F("sca_tmp", "sca_tmp", nbinsX, vminX, vmaxX, nbinsY, vminY, vmaxY);
    double vbinX = (vmaxX - vminX) / nbinsX;
    double vbinY = (vmaxY - vminY) / nbinsY;

    TF1* gaussf = new TF1("gaussf", "gaus", -100, +100);
    gaussf->SetNpx(400);

    TH1D* profileZ = new TH1D();
    for(int i = 0; i < nbinsX; i++) {
        for(int j = 0; j < nbinsY; j++) {
            profileZ = scatter_plot->ProjectionZ("_z", i+1, i+1, j+1, j+1, "");
			mval = profileZ->GetMaximum();
			mbin = profileZ->GetBinCenter(profileZ->GetMaximumBin());
			lo = profileZ->GetBinCenter(profileZ->FindFirstBinAbove(0.05 * mval));
			hi = profileZ->GetBinCenter(profileZ->FindLastBinAbove(0.05 * mval));
			profileZ->SetAxisRange(-10, +10);
	        if(profileZ->GetEntries() > nmin) {
				x = i * vbinX + vminX - 0.5*vbinX;
				y = j * vbinY + vminY - 0.5*vbinY;
				gaussf->SetParameter(1, mbin);
				profileZ->Fit("gaussf", "qre");

	            if(option == "0") hist2d->Fill(x, y, profileZ->GetMean());
				if(option == "1") hist2d->Fill(x, y, profileZ->GetRMS());
				if(option == "2") hist2d->Fill(x, y, gaussf->GetParameter(1));
				if(option == "3") hist2d->Fill(x, y, gaussf->GetParameter(2));
			}
        }
    }

    return hist2d;
}


std::vector<std::vector<double>> getRawFrom3DAlongZ(TH3F* scatter_plot, int nbinsX, double vminX, double vmaxX, int nbinsY,
                                double vminY, double vmaxY, int nmin, std::string option) {

	double x, y;
	double mval, mbin, lo, hi;
	double vbinX = (vmaxX - vminX) / nbinsX;
	double vbinY = (vmaxY - vminY) / nbinsY;

	std::vector<std::vector<double>> dat(256, std::vector<double>(0));

    TF1* gaussf = new TF1("gaussf", "gaus", -100, +100);
    gaussf->SetNpx(400);

    TH1D* profileZ = new TH1D();
    for(int i = 0; i < nbinsX; i++) {
        for(int j = 0; j < nbinsY; j++) {
            profileZ = scatter_plot->ProjectionZ("_z", i+1, i+1, j+1, j+1, "");
			mval = profileZ->GetMaximum();
			mbin = profileZ->GetBinCenter(profileZ->GetMaximumBin());
			lo = profileZ->GetBinCenter(profileZ->FindFirstBinAbove(0.05 * mval));
			hi = profileZ->GetBinCenter(profileZ->FindLastBinAbove(0.05 * mval));
			profileZ->SetAxisRange(-10, +10);
	        if(profileZ->GetEntries() > nmin) {

				x = i * vbinX + vminX - 0.5*vbinX;
				y = j * vbinY + vminY - 0.5*vbinY;
				gaussf->SetParameter(1, mbin);
				profileZ->Fit("gaussf", "qre");

	            if(option == "0") dat.at(i).at(j) = profileZ->GetMean();
				if(option == "1") dat.at(i).at(j) = profileZ->GetRMS();
				if(option == "2") dat.at(i).at(j) = gaussf->GetParameter(1);
				if(option == "3") dat.at(i).at(j) = gaussf->GetParameter(2);
			}
        }
    }

    return dat;
}


void prettyStats(float posx, float posy, float width, float height, bool frame, Color_t colour, int option){
	gStyle->SetStatFont(62);
	gStyle->SetStatX(posx);
	gStyle->SetStatY(posy);
	gStyle->SetStatW(width);
	gStyle->SetStatH(height);
	gStyle->SetStatTextColor(colour);
	gStyle->SetStatFormat("10.4f");
	gStyle->SetFitFormat("10.4f");

	gStyle->SetStatBorderSize(0);
	if (frame){
		gStyle->SetStatBorderSize(2);
	}

	// stats
	if (option == 0){
		gStyle->SetOptStat(1110);
		gStyle->SetOptFit(0);
	}

	// fit
	else if (option == 1){
		gStyle->SetOptStat(1110);
		// gStyle->SetOptFit(0011);
		gStyle->SetOptFit(111);
	}
}


void styleMap() {
	/* canvas */
    gStyle->SetPaperSize(20, 25);
    gStyle->SetCanvasDefH(400);
    gStyle->SetCanvasDefW(520);

    /* set the margins */
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadRightMargin(0.22);
    gStyle->SetPadLeftMargin(0.14);

    /* stats */
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(111);

    /* palette */
    gStyle->SetPalette(113);  // 51 = Jet, 56 = Temperature, 57 = Bird, 112 = Viridis, 113 = Cividis

}


void styleGraph() {
	/* canvas */
    gStyle->SetPaperSize(20, 25);
    gStyle->SetCanvasDefH(400);
    gStyle->SetCanvasDefW(520);

    /* margins */
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadRightMargin(0.10);
    gStyle->SetPadLeftMargin(0.18);

    /* font size */
	gStyle->SetTextSize(0.05);
	gStyle->SetTitleSize(0.055, "xyz");
	gStyle->SetLabelSize(0.050, "xyz");
	gStyle->SetStatFontSize(0.045);
	gStyle->SetLegendTextSize(0.045);

    /* stats */
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(111);

    /* palette */
    gStyle->SetPalette(113);  // 51 = Jet, 56 = Temperature, 57 = Bird, 112 = Viridis, 113 = Cividis
}



void styleCommon() {

    /* default white background for all plots */
    gROOT->SetStyle("Plain");

    /* canvas */
    gStyle->SetPaperSize(20, 25);
    gStyle->SetCanvasDefH(400);
    gStyle->SetCanvasDefW(520);

    // gStyle->SetCanvasColor(kWhite);
    // gStyle->SetFrameFillColor(kWhite);
    // gStyle->SetStatColor(kWhite);
    // gStyle->SetPadColor(kWhite);
    // gStyle->SetFillColor(kWhite);
    // gStyle->SetTitleFillColor(kWhite);
    gStyle->SetDrawBorder(0); // no yellow border around histogram
    gStyle->SetCanvasBorderMode(0); // remove border of canvas

    /* offsets */
    gStyle->SetLabelOffset(0.015, "xyz");
    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(1.1, "yz");

    /* default text size */
    gStyle->SetTextSize(0.05);
    gStyle->SetTitleSize(0.055, "xyz");
    gStyle->SetLabelSize(0.050, "xyz");
    gStyle->SetStatFontSize(0.042);
    gStyle->SetLegendTextSize(0.042);

    /* fonts */
    int font = 42;
    gStyle->SetTitleFont(font);
    gStyle->SetStatFont(font);
    gStyle->SetTextFont(font);
    gStyle->SetLabelFont(font, "xyz");
    gStyle->SetTitleFont(font, "xyz");
    gStyle->SetLegendFont(font);

	/* colours */
	gStyle->SetTextColor(kBlack);

    /* title and stuff */
    gStyle->SetOptTitle(0);  // remove histogramm title
    gStyle->SetOptDate(0); // print date on canvas

    /* stat box */
    gStyle->SetStatBorderSize(2);
    // gStyle->SetStatFormat("6.2e");
    // gStyle->SetFitFormat("6.2e");
    //gStyle->SetStatH(0.20);
    //gStyle->SetStatW(0.18);
	gStyle->SetStatX(0.87);
	gStyle->SetStatY(0.87);

    /* legend box */
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(kWhite);

    /* frame */
    gStyle->SetFrameLineWidth(2);

    /* histogrammes */
    gStyle->SetHistLineWidth(2);
    gStyle->SetHistLineColor(kBlack);
    gStyle->SetHistLineStyle(0);

    /* graphs */
    gStyle->SetLineStyle(1);
    gStyle->SetLineWidth(2);

    /* fits */
    gStyle->SetFuncWidth(2);
    gStyle->SetFuncColor(kRed+1);

    /* errorbars */
	gStyle->SetFillStyle(3002);
	gStyle->SetFillColor(4);
    gStyle->SetEndErrorSize(0);
    gStyle->SetTickLength(0.02, "xyz");

    /* markers */
    gStyle->SetMarkerStyle(kFullCircle);
    gStyle->SetMarkerSize(1);

    /* grid */
    gStyle->SetPadGridX(0);
    gStyle->SetPadGridY(0);
    gStyle->SetGridStyle(3);

    /* ticks */
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetNdivisions(506, "xyz");

    gROOT->ForceStyle();
}
