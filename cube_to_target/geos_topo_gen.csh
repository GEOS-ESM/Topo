#!/bin/tcsh

if ( "$#argv" != 4) then
  echo "Wrong number of arguments specified:"
  set n = 1
  echo "ogrid = argv[1]"
  echo "Co    = argv[2]"
  echo "Fi    = argv[3]"
  echo "tag   = argv[4]"
  echo "     "
  echo "possible ogrid values: fv_0.9x1.25, ne30pg3 ..."
  exit
endif

set case = "Reg_"$argv[1]"_co"$argv[2]"_fi"$argv[3]"_"$argv[4]
echo $case

mkdir -p ../cases/${case}/output
cp *.F90 ../cases/${case}
cp Makefile ../cases/${case}
cp ../machine_settings.make ../cases/${case}

cd ../cases/${case}


# Assumes you are in the right directory, i.e, the one with F90 files and namelists
#----------------------------------------------------------------------------------
mkdir -p output
mkdir -p output/raw
mkdir -p output/clean


module load compiler/gnu/default
gmake clean
gmake

set n = 1
set ogrid = "$argv[$n]"
set n = 2
set Co = "$argv[$n]"
set n = 3
set Fi = "$argv[$n]"

if ( $ogrid == 'geos_fv_c12' ) then
   set scrip='/project/amp/juliob/Topo-generate-devel/Topo/inputdata/grid-descriptor-file/PE12x72-CF.nc4'
endif
if ( $ogrid == 'geos_fv_c24' ) then
   set scrip='/project/amp/juliob/Topo-generate-devel/Topo/inputdata/grid-descriptor-file/PE24x144-CF.nc4'
endif
if ( $ogrid == 'geos_fv_c48' ) then
   set scrip='/project/amp/juliob/Topo-generate-devel/Topo/inputdata/grid-descriptor-file/PE48x288-CF.nc4'
endif
if ( $ogrid == 'geos_fv_c90' ) then
   set scrip='/project/amp/juliob/Topo-generate-devel/Topo/inputdata/grid-descriptor-file/PE90x540-CF.nc4'
endif



./cube_to_target --grid_descriptor_file=$scrip --intermediate_cs_name='../../regression-test-data/gmted2010_bedmachine-ncube0540.nc' --output_grid=$ogrid --coarse_radius=$Co --fine_radius=$Fi -r -u 'userid@ucar.edu' -q 'output/' -z -a 2

# Variable res
#./cube_to_target --grid_descriptor_file='../../regression-test-data/'$scrip --intermediate_cs_name='../../regression-test-data/gmted2010_bedmachine-ncube0540.nc' --output_grid='ne0_15x2_SA' --coarse_radius=$Co --fine_radius=001 -r -y $Yfac -u 'userid@ucar.edu' -q 'output/'



exit
