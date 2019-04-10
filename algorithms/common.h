/********************************************************************
* File: common.h
* ------------------------
*
* Description:
* Common functions.
*
* Version:
* Author: Florian Pitters
*
*******************************************************************/

#ifndef __COMMON_H_INCLUDED__
#define __COMMON_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <ctime>

#include "TTree.h"
#include "TBranch.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TList.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "THStack.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TGaxis.h"
#include "TKDE.h"
#include "TText.h"
#include "TLatex.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TPaveStats.h"
#include "TStopwatch.h"

#define IN 1
#define OUT 0



// --------------------------------------------------------
// Configuration
// --------------------------------------------------------

// Files and directories
const std::string home_path = getenv("HOME");

const std::string cfg_path = home_path + "/Documents/Works/tpx3/config/";
const std::string res_path = home_path + "/Documents/Works/tpx3/results/";
const std::string dat_path = home_path + "/Documents/Works/tpx3/data/";

// Colors and styles
const Color_t kOneBlue = (new TColor(0/256., 43/256., 128/256.))->GetNumber();
const Color_t kOneMagenta = (new TColor(199/256., 68/256., 24/256.))->GetNumber();
const Color_t kOneOrange = (new TColor(222/256., 143/256., 5/256.))->GetNumber();
const Color_t kOneGreen = (new TColor(81/256., 153/256., 28/256.))->GetNumber();
const Color_t kOneCyan = (new TColor(24/256., 156/256., 199/256.))->GetNumber();
const Color_t kOneRed = (new TColor(214/256., 42/256., 46/256.))->GetNumber();

const std::vector<Marker_t> markers = {kFullSquare, kFullCircle, kFullTriangleUp, kFullTriangleDown, kFullDiamond, kFullCrossX, kFullDoubleDiamond};
const std::vector<Marker_t> markers_open = {kOpenSquare, kOpenCircle, kOpenTriangleUp, kOpenTriangleDown};
const std::vector<Style_t> styles = {1, 2, 4, 6, 7, 10};
const std::vector<Color_t> colours = {kOneBlue, kOneMagenta, kOneOrange, kOneGreen, kOneCyan, kOneRed};



// Chip parameters
const int ncol = 256;
const int nrow = 256;

// Calibration parameters
const int npulses = 100;

// Parameters needed for special scripts
const int nsrc = 2;
const std::string cond = "";

// Example pixels
const int xpix = 65;
const int ypix = 142;

// Flags
const bool fsave = true;
const bool fkde = false;
const bool fspec = true;




// --------------------------------------------------------
// Methods
// --------------------------------------------------------

struct fitResult {
	Float_t mean;
	Float_t mean_err;
	Float_t std;
	Float_t std_err;
	Float_t chi2;
	Float_t ndf;
};

void print_line(std::string line);

float gauss(float *x, float *par);
float mean(std::vector<float> vals);
float var(std::vector<float> vals);
float median(std::vector<float> vals);
float iqr(std::vector<float> vals);
float mean_along_column(std::vector<std::vector<float>> vals, int col);
float var_along_column(std::vector<std::vector<float>> vals, int col);
int find_first_above_x(std::vector<float> vals, float thr, int reverse = 0);
int find_extremum(std::vector<float> vals, std::string option = "max");
float volts_to_tot(float volts, float *par);
float volts_to_toa(float volts, float *par);
float tot_to_volts(float tot, float *par);
float toa_to_volts(float toa, float *par);
float volts_to_charge(float volt, float *par);
int read_data_mat(const char *file_name, float data[ncol][nrow], bool info);
int read_data_list(const char *file_name, float data[][ncol*nrow], int length, bool info);
int write_data_mat(const char *file_name, float data[ncol][nrow], bool info);
int write_data_list(const char *file_name, float data[][ncol*nrow], int length, bool info);
int read_data_mat_to_tree(TTree *tr, const char *file_name, int ncol, bool info);
int read_file(std::string fn, char delim, std::vector<std::vector<float>> &dat, bool info);
int read_gain_data(std::string fn, char delim, std::vector<std::vector<float>> &vals, std::vector<std::vector<float>> &codes, bool info);
int read_noise_data(std::string fn, char delim, std::vector<std::vector<float>> &vals, std::vector<std::vector<float>> &codes, bool info);
int write_file(std::string fn, char delim, std::vector<std::vector<float>> &dat, std::string hd, bool info);
int create_directories(std::string path, mode_t mode);
fitResult *fitGauss(TH1 *hist1d, float fit_scale, float low, float high);


void prettyStats(float posx, float posy, float width, float height, bool frame, Color_t colour, int option);

void styleGraph();
void styleMap();
void styleCommon();


void printGraph(TGraph* gr, std::string fn, std::string title, std::string titleX, std::string titleY, float limX = -1.0, float limY = -1.0, float offset = 1.1,
                std::string option = "AP");
void printGraph(TGraphErrors* gr, std::string fn, std::string title, std::string titleX, std::string titleY, float limX = -1.0, float limY = -1.0,
			float offset = 1.1,std::string option = "APE");
void printMap(TH2* hist2d, std::string fn, std::string title, std::string titleX, std::string titleY, std::string titleZ, float limX = -1.0, float limY = -1.0, float offsetY = 1.0,
              float offsetZ = 1.0, std::string option = "COLZ1");
void printHist(TH1* hist1d, std::string fn, std::string title, std::string titleX, std::string titleY, float limX = -1.0, float limY = -1.0, float offset = 1.1,
               std::string option = "");

void addGraph(TMultiGraph* mg, TLegend* lg, int nentries, float* x, float* y,
        Color_t c, Style_t m, Size_t ms, Style_t ls, std::string style, std::string label);
void addErrGraph(TMultiGraph* mg, TLegend* lg, int nentries, float* x, float* y, float* x_err, float* y_err,
        Color_t c, Style_t m, Size_t ms, Style_t ls, std::string style, std::string label);
void addFittedGraph(TMultiGraph* mg, TLegend* lg, int nentries, float* x, float* y,
        Color_t c, Style_t m, Size_t ms, Style_t ls, std::string style, std::string label, std::string func, Float_t fit_low, Float_t fit_up);
TFitResultPtr addFittedErrGraph(TMultiGraph* mg, TLegend* lg, int nentries, float* x, float* y, float* x_err, float* y_err,
		Color_t c, Style_t m, Size_t ms, Style_t ls, std::string style, std::string label, std::string func, Float_t fit_low, Float_t fit_up);
void addFitResults(TF1* func, TLegend* lg, TFitResultPtr fitr, Color_t c, Style_t ls, std::string label);
void addScatterGraph(TLegend* lg, TH2* scatter,
        Color_t c, Style_t m, Size_t ms, std::string style, std::string label);
void addHist(THStack* hs, TLegend* lg, TH1F* hist1d,
        Color_t c, Style_t m, Style_t ls, std::string style, std::string label, Float_t scale);
void addHist(THStack* hs, TLegend* lg, TH1D* hist1d,
		Color_t c, Style_t m, Style_t ls, std::string style, std::string label, Float_t scale);
void addHist(THStack* hs, TLegend* lg, TH2F* hist2d,
		Color_t c, Style_t m, Style_t ls, std::string style, std::string label, Float_t scale);

TGraphErrors* getProjectionFrom2DAlongY(TH2F* scatter_plot, int step, double vmin, double vmax, int nmin, std::string option);
TGraphErrors* getProjectionFrom3DAlongZ(TH3F* scatter_plot, int nbinsX, double vminX, double vmaxX, int binZ, int nmin, std::string option);
TH2F* getMapFrom3DAlongZ(TH3F* scatter_plot, int nbinsX, double vminX, double vmaxX, int nbinsY,
    	double vminY, double vmaxY, int nmin, std::string option);

std::vector<std::vector<double>> getRawFrom3DAlongZ(TH3F* scatter_plot, int nbinsX, double vminX, double vmaxX, int nbinsY,
        double vminY, double vmaxY, int nmin, std::string option);

#endif
