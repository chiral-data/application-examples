#!/bin/bash
#
# a free benchmark set of Gromacs
# by Max Planck Institute
# link: https://www.mpinat.mpg.de/grubmueller/bench
# paper: https://pubs.acs.org/doi/10.1021/acs.jcim.2c00044

BENCHMARK_INPUTS="
https://www.mpinat.mpg.de/3925268/cmet_eq.zip
https://www.mpinat.mpg.de/3925290/hif2a_eq.zip
https://www.mpinat.mpg.de/3925312/ligand_cmet_eq.zip
https://www.mpinat.mpg.de/3925334/shp2_eq.zip
"
GMX=/usr/local/gromacs/avx2_256/bin/gmx

for tuple_string in ${BENCHMARK_INPUTS}; do
    prefix=$(echo $tuple_string | cut -d '/' -f 4)
    benchmark=$(echo $tuple_string | cut -d '/' -f 5 | cut -d '.' -f 1)

    echo "Downloading input files for benchmark $benchmark ..."
    curl -O https://www.mpinat.mpg.de/$prefix/$benchmark.zip 
    unzip $benchmark.zip

    echo "Running the benchmark $benchmark ..."
    # $GMX mdrun -s $benchmark.tpr -nsteps 10000 -ntomp 2
    $GMX mdrun -s $benchmark.tpr -nsteps 1000 
    echo "Benchmark $benchmark done ..."
    echo ""
done

# Benchmarks of "Binding affinity study benchmarks"
#
# MSI laptop: Intel Core i7 4-core, Nvidia GTX 3060
#
# Instance            System              Steps      ns/day
# MSI laptop          CMET EQ             1000        6.530 
# MSI laptop          HIF2A EQ            1000        8.990
# MSI laptop          CMET Ligand EQ      1000       68.806 
# MSI laptop          SHP2 EQ             1000        4.517

# Benchmark results on Sakura Internet, DOK & Server
#                                 steps       ns/day    core t(s)   wall t(s)
# V100-DOK    1-MPI,2-OpenMP      1,000       14.524    23.814      11.910 
# V100-DOK    1-MPI,3-OpenMP      1,000       14.339    36.181      12.063
# V100-DOK    1-MPI,4-OpenMP      1,000       13.558    51.020      12.758 
# V100-DOK    1-MPI,2-OpenMP     10,000       20.655
# V100-DOK    1-MPI,3-OpenMP     10,000       20.447 
# V100-DOK    1-MPI,4-OpenMP     10,000       18.266
# V100-DOK    1-MPI,3-OpenMP    100,000       23.567 
# H100-DOK    1-MPI,10-OpenMP   100,000      105.692
# H100-DOK    1-MPI,10-OpenMP 1,000,000      108.374
# V100-Server 1-MPI,2-OpenMP      1,000       24.351    14.205       7.103 
# V100-Server 1-MPI,3-OpenMP      1,000       41.236    12.584       4.195 
# V100-Server 1-MPI,4-OpenMP      1,000       51.143    13.528       3.382
# V100-Server 1-MPI,2-OpenMP     10,000       27.727   124.653      62.329
# V100-Server 1-MPI,3-OpenMP     10,000       38.862   133.405      44.470 
# V100-Server 1-MPI,4-OpenMP     10,000       51.118 
