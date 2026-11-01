source /home/t2k/nowak/MicroBooNE/production/nuwro_prod/setup_env.sh >/dev/null 2>&1
H=/home/t2k/nowak/MicroBooNE/production/nuwro_prod/genie_gen/hypfix
cd $H/test
g++ -O2 $(root-config --cflags) -I$H -I$GENIE_INC/GENIE -I$LOG4CPP_INC port_check.cxx -o port_check \
   -L$H/lib -lGHypFix $(genie-config --libs) -L$LOG4CPP_LIB -llog4cpp -L$LIBXML2_FQ_DIR/lib -lxml2 \
   -L$PYTHIA_LIB -lPythia6 -L$LHAPDF_LIB -lLHAPDF $(root-config --libs) -lGeom -lEG -lEGPythia6 -lMathMore -lGenVector 2>&1 | grep -i error | head -5
LD_LIBRARY_PATH=$H/lib:$LD_LIBRARY_PATH ./port_check
