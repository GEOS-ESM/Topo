ecbuild_add_executable (TARGET cube_to_target.x SOURCES reconstruct.F90  ridge_utils.F90  shr_kind_mod.F90  cube_to_target.F90 remap.F90 rot.F90 smooth_topo_cube.F90 f90getopt.F90 ridge_ana.F90 shared_vars.F90  subgrid_topo_ana.F90)
target_link_libraries (cube_to_target.x PRIVATE MPI::MPI_Fortran NetCDF::NetCDF_Fortran)
