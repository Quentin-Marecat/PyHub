set(MY_MODS   mods.f90
              CACHE INTERNAL "")

add_executable (main.x mods.f90 mainf.f90 diagmat.f90 basis.f90 hprod.f90 lanczos.f90 solve.f90 static_rq.f90 exact_diag.f90 boltzmann.f90 spgf_ed.f90 spgf_bandlanczos.f90 )
target_link_libraries(main.x ${BLAS_LIBRARIES} hdf5::hdf5_fortran)
