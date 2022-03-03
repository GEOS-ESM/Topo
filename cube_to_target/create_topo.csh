#!/bin/tcsh


if ( "$#argv" != 4) then
  echo "Wrong number of arguments specified:"
set n = 1
  echo "ogrid = argv[1]"
  echo "Co    = argv[2]"
  echo "Fi    = argv[3]"
  ##echo "Nsw   = argv[4]"
  echo "tag  = argv[4]"
  echo "     "
  echo "possible ogrid values: fv_0.9x1.25, ne30pg3 ..."
  exit
endif

set n = 4
set case = $argv[1]"_co"$argv[2]"_fi"$argv[3]"_"$argv[4]

echo $case


mkdir -p ../cases/${case}/output
cp *.F90 ../cases/${case}
cp Makefile ../cases/${case}
cp ../machine_settings.make ../cases/${case}
#cp clean_topo_files.csh ../cases/${case}
#cp -r analysis ../${case}


cd ../cases/${case}


#cp output/topo_smooth* ../${case}/output/
#cd /project/amp/juliob/Topo-generate-devel/Topo/smooth_topo/
#set tops=`ls -1 *`
#cd /project/amp/juliob/Topo-generate-devel/Topo/${case}
#foreach foo ($tops)
#   echo $foo
#   ln -sf /project/amp/juliob/Topo-generate-devel/Topo/smooth_topo/${foo} output/${foo}
#end


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

# This is now used for all. Doesn't matter, will eliminate
#set Nrs=00

if ( $ogrid == 'geos_fv_c48' ) then
   set scrip='PE48x288-CF.nc4'
endif
if ( $ogrid == 'geos_fv_c90' ) then
   set scrip='PE90x540-CF.nc4'
endif
if ( $ogrid == 'geos_fv_c180' ) then
   set scrip='PE180x1080-CF.nc4'
endif
if ( $ogrid == 'geos_fv_c360' ) then
   set scrip='PE360x2160-CF.nc4'
endif
if ( $ogrid == 'geos_fv_c720' ) then
   set scrip='PE720x4320-CF.nc4'
endif
if ( $ogrid == 'geos_fv_c1440' ) then
   set scrip='PE1440x8640-CF.nc4'
endif
if ( $ogrid == 'fv_0.9x1.25' ) then
   set scrip='fv_0.9x1.25.nc'
endif
if ( $ogrid == 'ne120pg3' ) then
   set scrip='ne120pg3.nc'
endif
if ( $ogrid == 'ne30pg3' ) then
   set scrip='ne30pg3.nc'
endif
if ( $ogrid == 'Arctic' ) then
   set scrip='ne0ARCTICne30x4_scrip_c191212.nc'
endif

set scrip = '/project/amp/juliob/Topo-generate-devel/Topo/inputdata/grid-descriptor-file/'${scrip}
set cstopo = '/project/amp/juliob/Topo-generate-devel/Topo/inputdata/cubed-sphere-topo/gmted2010_modis-ncube3000-stitch.nc'
#set smtopo = '/project/amp/juliob/Topo-generate-devel/Topo/smooth_topo/topo_smooth_nc3000_Co060_Fi001.dat'
set smtopo = '/project/amp/juliob/Topo-generate-devel/Topo/Topo.git/cases/fv_co60_fi01_nothin/output/topo_smooth_nc3000_Co060_Fi001.dat'
echo $smtopo

./cube_to_target --grid_descriptor_file=$scrip --intermediate_cs_name=$cstopo --smooth_topo_file=$smtopo --output_grid=$ogrid --coarse_radius=$Co --fine_radius=$Fi -r

# --rrfac_max=4 --use_prefilter --find_ridges --precomputed_smooth_topo  --regional_refinement


exit
