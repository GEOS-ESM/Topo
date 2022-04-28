#!/bin/tcsh
./cube_to_target --grid_descriptor_file='../regression-test-data/ne30pg3.nc' --intermediate_cs_name='../regression-test-data/gmted2010_bedmachine-ncube0540.nc' --output_grid='ne30pg3' --coarse_radius=012 -p -r -u 'Peter Hjort Lauritzen, pel@ucar' -q 'output/' --grid_descriptor_file_gll='../regression-test-data/ne30np4.nc' --j
#
# CESM2 version
#
# ./cube_to_target --grid_descriptor_file='../regression-test-data/ne30pg3.nc' --intermediate_cs_name='../regression-test-data/gmted2010_bedmachine-ncube0540.nc' --output_grid='ne30pg3' --coarse_radius=012 -p -r -u 'Peter Hjort Lauritzen, pel@ucar' -q 'output/' --j
set file=`ls -t1 output/*.nc |  head -n 1`
echo $file
/glade/p/cesm/cseg/tools/cprnc/cprnc -m $file ../regression-test-data/ne30pg3_gmted2010_bedmachine_nc0540_Co012_20220311.nc
#source plot.sh
