# Install
At the folder's root type
```
make  
make test
```

# Requirements
root


# Usage

Configure by editing 'algorithms/common.h'. Fix data paths etc. Recompile.

To calibrate testpulse data use
```
./bin/run -a calibrate_testpulses --dev [id] --thr_dac [value] --ik_dac [value]
```

To analyse calibration use
```
./bin/run -a cal_analysis --dev [id] --thr_dac [value] --ik_dac [value]
```

To analyse equalisation and noise use
```
./bin/run -a dev_analysis --dev [id] --thr_dac [value] --pol [value]
```

# Examples

Running full calibration chain
```
./bin/run -a hello_world" << std::endl;
./bin/run -a dev_analysis --dev W0019_C07 --pol 0
./bin/run -a cal_analysis --dev W0019_C07 --thr_dac 1190 --ik_dac 10
./bin/run -a calibrate_testpulses --dev W0019_C07 --thr_dac 1190 --ik_dac 10
./bin/run -a calibrate_sources --dev W0019_C07 --thr_dac 1190 --ik_dac 10
```

Running full calibration chain
```
./bin/run -a calibrate_testpulses --dev W0019_C07 --thr_dac 1190 --ik_dac 10
./bin/run -a prepare_sources --dev W0019_C07 --thr_dac 1190 --ik_dac 10
./bin/run -a calibrate_sources --dev W0019_C07 --thr_dac 1190 --ik_dac 10
./bin/run -a acq_analysis --dev W0019_C07 --thr_dac 1190 --ik_dac 10 --src fe
```