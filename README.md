# 1D-Tight-Binding-Model (Under Construction)
Simple 1D TB model to obtain the dispersion curve and z component of spin-spin correlation.

 ---------------------- 1D Tight Binding Model --------------------

 Author: João Augusto Sobral da Silva / Github: @joaosds

 The code is commented, and each subroutine has a brief description (B.D.)
 of its functionality. For any discussions you can contact me at Github!

 History:

  v1.0 (04/30/20) - First implementation
  v1.1 (08/20/20 - Comments, cleaning code, and python script for plots 
  
  Future intended changes:
  - implementation of modern fortran (>=2008);
  - hamiltonian matrix determination from a neighbor's lattice;
  - python script for plots;
  - python interface to link both the fortran file and the plots automatically;

# File information

 File 'disp.dat' prints the eigenvalues ordered by increasing order
 and the original order from the LAPACK function DSYEV;


 File 'parameters.dat' informs the initial parameters N,t, BC: number of
 sites, value of hopping parameter t and lattice boundary condition. 
 Type 0 for periodic, 1 for anti-periodic and 2 for open conditions.

 File 'dispanaly.dat' displays the expected eigenvalues by the analytical
 expression. It can be compared to the 2nd column of file 'disp.dat'.
 
 File 'spincor.dat' informs the z component of the spin-spin correlation
 as obtained from the eigenfunctions of the DSEYV function;
 
 File 'spincor.dat' shows the expected spin-spin correlation by the discrete
 analytical expression.
 
 Note: For theoretical background see pdf "Tight Binding 1D implementation.pdf"

