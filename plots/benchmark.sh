#!/bin/bash
set -e
DIR=`dirname "${BASH_SOURCE[0]}"`
BENCH=$DIR/benchmark

CORES=4
TRIALS=50
PARTICLES=(16 32 64 128 256 512 1024 2048 4096)

# Create the results directory structure
if [ ! -d $BENCH ]; then
    mkdir $BENCH
    mkdir $BENCH/smc
    mkdir $BENCH/sis
else
    read -p "The output directory '${BENCH}/' already exists. Do you want to replace the existing benchmark with a new one? [y/n] " yn
    case $yn in
	[yY]* ) rm -rf $BENCH;
	      mkdir $BENCH;
	      mkdir $BENCH/smc;
	      mkdir $BENCH/sis;;
	* ) exit 1;;
    esac
fi

echo -e "\nStarting benchmark, current configuration is:"
echo -e "\t max cores: $CORES"
echo -e "\t trials per number of particles: $TRIALS"
echo -e "\t number of particles: ${PARTICLES[@]}"
echo -e "\n(You can change these settings by editting the file ${BASH_SOURCE[0]})\n"

# Run the benchmarks in parallel
for i in ${PARTICLES[@]}; do
    echo -n "Running $TRIALS trials for $i particles. SMC: "
    seq 1 $TRIALS | xargs -P $CORES -I % python $DIR/run_approximation.py --method SMC --particles $i --filename $BENCH/smc/approx_${i}_%
    echo -n "done, SIS: "
    seq 1 $TRIALS | xargs -P $CORES -I % python $DIR/run_approximation.py --method SIS --particles $i --filename $BENCH/sis/approx_${i}_%
    echo "done."
done
