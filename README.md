Educational Single Column Model
------


License: GPLv3

A simple single column model (SCM) for educational purposes. It uses the e-ε (k-ε) 1.5-order turbulence closure for a horizontally homogeneous column of atmosphere concentrating on the atmospheric boundary layer. So far, the moisture and moisture fluxes are ignored.

The implementation of the turbulent closure follows Zhang, Wang & Xue (2020) Monthly Weather Review 148, https://10.1175/MWR-D-19-0084.1

Configuration of the model is currently made by changing the variables inside the source code.

Compilation:

    gfortran escm.f90 -o escm -I$NETCDF_FORTRAN_INC  -L$NETCDF_FORTRAN_LIB -lnetcdff

Contact: vladimir.fuka@mff.cuni.cz

