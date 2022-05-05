#!/bin/tcsh
make
#./cube_to_target --grid_descriptor_file='../regression-test-data/ne30pg3.nc' --intermediate_cs_name='../bin_to_cube/gmted2010_bedmachine-ncube540.nc' --output_grid='ne30pg3' --coarse_radius=012 -p -r -u 'Peter Hjort Lauritzen, pel@ucar' -q 'output/' -x

#./cube_to_target --grid_descriptor_file='/glade/p/cesmdata/cseg/inputdata/share/scripgrids/fv0.9x1.25_141008.nc' --intermediate_cs_name='../bin_to_cube/gmted2010_bedmachine-ncube540.nc' --output_grid='fv0.9x1.25' -u 'Peter Hjort Lauritzen, pel@ucar' -q 'output/' -b 20E7 -l 60 -m
#
# Laplacian smoother
#
#./cube_to_target --grid_descriptor_file='/glade/p/cesmdata/cseg/inputdata/share/scripgrids/fv0.9x1.25_141008.nc' --intermediate_cs_name='../bin_to_cube/gmted2010_bedmachine-ncube540.nc' --output_grid='fv0.9x1.25' -u 'Peter Hjort Lauritzen, pel@ucar' -q 'output/' -b 20E7 -l 60 -m
#
# Distance weighted smoother
#
./cube_to_target --grid_descriptor_file='/glade/p/cesmdata/cseg/inputdata/share/scripgrids/fv0.9x1.25_141008.nc' --intermediate_cs_name='../bin_to_cube/gmted2010_bedmachine-ncube540.nc' --output_grid='fv0.9x1.25' -u 'Peter Hjort Lauritzen, pel@ucar' -q 'output/'  --coarse_radius=011 -p -m

#
# Command
#
# qcmd -l walltime=01:30:00 -- source tmp.sh >& out &
#
#./cube_to_target --grid_descriptor_file='/glade/p/cesmdata/cseg/inputdata/share/scripgrids/fv0.9x1.25_141008.nc' --intermediate_cs_name='../bin_to_cube/gmted2010_modis-ncube3000.nc' --output_grid='fv0.9x1.25' -u 'Peter Hjort Lauritzen, pel@ucar' -q 'output/' -b 20E7 -l 1800 -m 
