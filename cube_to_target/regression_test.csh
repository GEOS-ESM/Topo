#!/bin/tcsh



set case = "test3"

mkdir -p ../cases/${case}/output
cp *.F90 ../cases/${case}
cp Makefile ../cases/${case}
cp ../machine_settings.make ../cases/${case}
#cp clean_topo_files.csh ../cases/${case}
#cp -r analysis ../${case}


cd ../cases/${case}


# Assumes you are in the right directory, i.e, the one with F90 files and namelists
#----------------------------------------------------------------------------------
mkdir -p output
mkdir -p output/raw
mkdir -p output/clean


module load compiler/gnu/default
gmake clean
gmake

set ogrid = "ne30pg3"
set Co = "12"
set Fi = "1"
set Nsw = "8"



# This is now used for all. Doesn't matter, will eliminate
set Nrs=00



./cube_to_target --grid_descriptor_file='../../regression-test-data/ne30pg3.nc' --intermediate_cs_name='../../regression-test-data/gmted2010_bedmachine-ncube0540.nc' --output_grid=$ogrid --coarse_radius=$Co --fine_radius=$Fi -p -r -2 1

#./cube_to_target --grid_descriptor_file='../../regression-test-data/ne30pg3.nc' --intermediate_cs_name='../../regression-test-data/gmted2010_bedmachine-ncube0540.nc' --output_grid=$ogrid --coarse_radius=$Co --fine_radius=$Fi --smooth_topo_file='../regress-14/output/topo_smooth_nc0540_Co012_Fi001.dat' -p -r -2 1


exit
