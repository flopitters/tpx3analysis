#include <string>

#include "hello.h"
#include "common.h"
#include "device.h"
#include "calibrate.h"
#include "analyse.h"



int main(int argc,  char **argv)
{
	bool ret = 0;
	bool fReadConfig = false;
	bool fPrintHelp = false;
	bool fAlgorithm = false;
	int thr_dac = 1190;
	int ik_dac = 10;
	int pol = 0;
	std::string src = "fe";
	std::string dev_id = "";

	std::string algorithm = "";
	std::string configFile = "/Users/Home/Documents/Code/Works/tpx3Analysis/config/example.cfg";


	// Read arguments and determine what to do
	if (argc == 1){
		fPrintHelp = true;
	}
	else {
		for (int i = 1; i < argc; i += 2) {
			std::string option(argv[i]);
			if (option.compare("-h") == 0) fPrintHelp = true;
			if (option.compare("-?") == 0) fPrintHelp = true;
			if (option.compare("-a") == 0) fAlgorithm = true;
			if (option.compare("-c") == 0) fReadConfig = true;

			std::string argument(argv[i+1]);
			if (option.compare("-a") == 0) algorithm = argv[i+1];
			if (option.compare("-c") == 0) configFile = argv[i+1];
			if (option.compare("--dev") == 0) dev_id = argv[i+1];
			if (option.compare("--thr_dac") == 0) thr_dac = std::stoi(argv[i+1]);
			if (option.compare("--ik_dac") == 0) ik_dac = std::stoi(argv[i+1]);
			if (option.compare("--pol") == 0) pol = std::stoi(argv[i+1]);
			if (option.compare("--src") == 0) src = argv[i+1];

			if (option.compare("-h") == 1 && option.compare("-?") == 1 && option.compare("-a") == 1 &&
				option.compare("-c") == 1 && option.compare("--dev") == 1 && option.compare("--thr_dac") == 1 &&
				option.compare("--ik_dac") == 1 && option.compare("--pol") == 1 && option.compare("--src") == 1){
					std::cout << " Unknown paramater. " << std::endl;
					return -1;
				}
		}
	}

	// Print help
	if (fPrintHelp) {
		std::cout << "Usage: " << argv[0] << " -a <algorithm> " << std::endl;
		std::cout << "\n" << std::endl;

		std::cout << "Avalable options:" << std::endl;
		std::cout << " -h           Display this help" << std::endl;
		std::cout << " -?           Display this help" << std::endl;
		std::cout << " -a           Determine algorithm to use" << std::endl;
		std::cout << " -c           Determine config file to read" << std::endl;
		std::cout << " --dev		Device ID, e.g. W0005_E02" << std::endl;
		std::cout << " --pol        Polarity of majority charge carriers, electrons is 0" << std::endl;
		std::cout << " --src        Source used for Xray characterisation" << std::endl;
		std::cout << " --thr_dac	Used threshold" << std::endl;
		std::cout << " --ik_dac     Used discharge current" << std::endl;
		std::cout << "\n" << std::endl;

		std::cout << "Example commands:" << std::endl;
		std::cout << "./bin/run -a hello_world" << std::endl;
		std::cout << "./bin/run -a dev_analysis --dev W0019_C07 --pol 0" << std::endl;
		std::cout << "./bin/run -a cal_analysis --dev W0019_C07 --thr_dac 1190 --ik_dac 10" << std::endl;
		std::cout << "./bin/run -a calibrate_testpulses --dev W0019_C07 --thr_dac 1190 --ik_dac 10" << std::endl;
		std::cout << "./bin/run -a calibrate_sources --dev W0019_C07 --thr_dac 1190 --ik_dac 10" << std::endl;
		std::cout << "\n" << std::endl;
		std::cout << "Full Calibration Example:" << std::endl;
		std::cout << "./bin/run -a calibrate_testpulses --dev W0019_C07 --thr_dac 1190 --ik_dac 10" << std::endl;
		std::cout << "./bin/run -a prepare_sources --dev W0019_C07 --thr_dac 1190 --ik_dac 10" << std::endl;
		std::cout << "./bin/run -a calibrate_sources --dev W0019_C07 --thr_dac 1190 --ik_dac 10" << std::endl;
		std::cout << "./bin/run -a acq_analysis --dev W0019_C07 --thr_dac 1190 --ik_dac 10 --src fe" << std::endl;

		std::cout << "\n\n" << std::endl;
	}

	//ret = readConfig(configFile);
	std::vector<std::string> srcs = {"fe", "in"};

	// Process arguments
	if (fAlgorithm) {
		if (algorithm == "hello_world")
		 	hello_world();
		else if (algorithm == "calibrate_testpulses") {
			calibrate_testpulses(dev_id, thr_dac, ik_dac);
		}
		else if (algorithm == "prepare_sources") {
			calibrate_sources_prepare(dev_id, thr_dac, ik_dac, "fe");
			calibrate_sources_prepare(dev_id, thr_dac, ik_dac, "in");
		}
		else if (algorithm == "calibrate_sources") {
			calibrate_sources(dev_id, thr_dac, ik_dac, srcs);
		}
		else if (algorithm == "acq_analysis"){
			acq_analysis(dev_id, thr_dac, ik_dac, src);
		}
		else if (algorithm == "cal_analysis") {
			cal_analysis(dev_id, thr_dac, ik_dac);
		}
		else if (algorithm == "dev_analysis") {
			dev_analysis(dev_id, pol);
		}

		else
			std::cout << "not a valid algorithm" << std::endl;

		printf("---->  Success!\n\n");
	}

	return 0;
}
