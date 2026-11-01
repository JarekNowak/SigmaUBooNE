source /home/t2k/nowak/MicroBooNE/production/nuwro_prod/setup_env.sh >/dev/null 2>&1
H=/home/t2k/nowak/MicroBooNE/production/nuwro_prod/genie_gen/hypfix
source $H/env_fix.sh
cd $H/test
time gmkspl -p -14 -t 1000010010,1000000010,1000180400 -n 40 -e 20 --event-generator-list LambdaCCQE \
     --tune G18_10a_02_11a -o spl_test.xml > spl_test.log 2>&1
echo rc=$?; grep -c '<spline' spl_test.xml; grep 'spline name' spl_test.xml
grep -i 'HypQEL\|Pais\|error\|fatal' spl_test.log | head -5
