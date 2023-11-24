set(MY_MODS   hubmod.f90
              basismod.f90 
              opemod.f90
              funcmod.f90
              CACHE INTERNAL "")

add_executable (basis.x basismod.f90 funcmod.f90 compute_basis.f90 basis.f90 )
target_link_libraries(basis.x  hdf5::hdf5_fortran)

add_executable (operators.x basismod.f90 funcmod.f90 opemod.f90 compute_operators.f90 operators.f90 basis.f90)
target_link_libraries(operators.x  hdf5::hdf5_fortran)

add_executable (hubbard.x hubmod.f90 funcmod.f90 basismod.f90 hubbard.f90 diagmat.f90 basis.f90 hprod.f90 lanczos.f90 solve_hubbard.f90 static_rq.f90 exact_diag.f90 boltzmann.f90 spgf_ed.f90 spgf_bandlanczos.f90 )
target_link_libraries(hubbard.x ${BLAS_LIBRARIES} hdf5::hdf5_fortran)
