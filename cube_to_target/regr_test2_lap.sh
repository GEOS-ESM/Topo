#!/bin/tcsh
make
./cube_to_target --grid_descriptor_file='../regression-test-data/ne0_15x2.nc' --intermediate_cs_name='../regression-test-data/gmted2010_bedmachine-ncube0540-220518.nc' --output_grid='ne0_15x2_SA' -y 2 -d -u "Peter Hjort Lauritzen, pel@ucar.edu" -m --smoothing_scale=200.0 -v
set file=`ls -t1 output/*.nc |  head -n 1`
echo $file
/glade/p/cesm/cseg/tools/cprnc/cprnc -m $file ../regression-test-data/ne0_15x2_SA_gmted2010_bedmachine_nc0540_Laplace0200_20220606.nc
#/fs/cgd/csm/tools/cprnc/cprnc  -m $file ../regression-test-data/ne0_15x2_SA_gmted2010_bedmachine_nc0540_Laplace0200_20220511.nc
#source plot.sh
