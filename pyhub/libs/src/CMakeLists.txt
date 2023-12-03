set(MY_MODS   hubmod.f90
              basismof.f90
              opemod.f90
              funcmod.f90
              CACHE INTERNAL "")

add_executable (basis.x basismod.f90 compute_basis.f90)
target_link_libraries(basis.x ${HDF5_Fortran_LIBRARIES})
# target_link_libraries(basis.x hdf5::hdf5_fortran)

add_executable (operators.x funcmod.f90 opemod.f90 basismod.f90 compute_operators.f90 operators.f90 operators_basis.f90)
target_link_libraries(operators.x  ${HDF5_Fortran_LIBRARIES})

add_executable (hubbard.x basismod.f90 hubmod.f90 funcmod.f90 compute_solver.f90 diagmat.f90 hprod_hubbard.f90 lanczos.f90 solve.f90 static_rq.f90 exact_diag.f90 boltzmann.f90 spgf_ed.f90 spgf_bandlanczos.f90 )
target_link_libraries(hubbard.x ${BLAS_LIBRARIES} ${HDF5_Fortran_LIBRARIES})

add_executable (heisenberg.x basismod.f90 hubmod.f90 funcmod.f90 compute_solver.f90 diagmat.f90 hprod_heis.f90 lanczos.f90 solve.f90 static_rq_heis.f90 exact_diag.f90 boltzmann.f90 spgf_ed.f90 spgf_bandlanczos.f90 )
target_link_libraries(heisenberg.x ${BLAS_LIBRARIES} ${HDF5_Fortran_LIBRARIES})

