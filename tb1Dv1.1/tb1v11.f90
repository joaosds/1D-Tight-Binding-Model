! ---------------------- Variational Monte Carlo 1D --------------------
!
! Author: João Augusto Sobral da Silva / Github: @joaosds
!
! The code is commented, and each subroutine has a brief description (B.D.)
! of its functionality. For any discussions you can contact me at Github!
!
! History:
!
!  v1.0 (04/30/20) - First implementation
!  v1.1 (08/20/20 - Comments, cleaning code, and python script for plots 
!
! File 'disp.dat' prints the eigenvalues ordered by increasing order
! and the original order from the LAPACK function DSYEV;
!
!
! File 'parameters.dat' informs the initial parameters N,t, BC: number of
! sites, value of hopping parameter t and lattice boundary condition. 
! Type 0 for periodic, 1 for anti-periodic and 2 for open conditions.
!
! File 'dispanaly.dat' displays the expected eigenvalues by the analytical
! expression. It can be compared to the 2nd column of file 'disp.dat'.
! 
! File 'spincor.dat' informs the z component of the spin-spin correlation
! as obtained from the eigenfunctions of the DSEYV function;
! 
! File 'spincor.dat' shows the expected spin-spin correlation by the discrete
! analytical expression.
! 
! Note: For theoretical background see pdf "Tight Binding 1D implementation.pdf"
! ---------------------- Informations ---------------------------------

module variables
 use ISO_FORTRAN_ENV
 implicit none
 real(real64), allocatable :: M(:,:) !real symmetric matrix
 real(real64), allocatable :: E(:)   ! matrix of analytical eigenvalues
 integer :: N, BC !Number of sites and order of matrix A / Boundary Condition option
 integer :: i,j !loop index
 real(real64) :: t,GS!Hopping (n.n.) and 
 real(real64) :: temp1,temp2! pointer variables

! Variables from external function
 integer :: LDA, INFO
 integer ::   LWORK
 character, external :: JOBZ, UPLO
 real(real64), allocatable :: W(:)
 real(real64), allocatable :: P(:)
 real(real64), allocatable :: k(:)
 real(real64), allocatable :: WORK(:)
 real(real64),parameter :: pi = 4.0d0*DATAN(1.0d0) ! pi in double precision * error at last digit (1, not 2)
end module variables

program tb1d

!       ##  Tight Binding 1D for first neighbors ##
!
! DSYEV routine computes all eigenvalues and, optionally, 
! eigenvectors of an n-by-n real symmetric matrix A.

 ! B.D.: main structure of the algorithm. First we read the initial
 ! parameters, construct the upper (or lower) triangle symmetric matrix,
 ! calculate the analytical eigenvalues (for further comparison) and call
 ! the rest of algorithm through subroutines.

 use Variables
 implicit none
 
 open(10, file = 'parameters.dat')
 read(10,*) N, t, BC
 print*, 'Initial Parameters: '
 print*, ''
 print*, 'N = ', N
 print*, 't = ', -t, 'eV'
 close(10)
   
 LDA = N
  
 LWORK = 3*LDA-1

  
 allocate(M(LDA,N),W(N),E(N),P(N),k(N),WORK(LWORK))
 
! ----------- Determination of Symmetric TB1D Matrix ------------------ !


 do i = 1,N-1 ! For next neighbour
  M(i+1,i) = -t !For upper triangle use  M(i,i+1)=t
 enddo
 
 if ( BC == 0 ) then
  M(N,1) = -t !Periodic conditions
 elseif (BC == 1) then
  M(N,1) = t !Antiperiodic conditions
 elseif (BC == 2) then
  M(N,1) = 0 !Open boundary condition
 endif

! -------- Analytic Eigenvalues ---------- !

 if (MOD(N,2) == 0) then !N is even
  do i = -N/2,N/2-1
   E(i+N/2+1) = -2*t*cos(i*2*pi/N)
   k(i+N/2+1) = i*2*Pi/N
  enddo
 elseif (MOD(N,2) > 0) then !N is odd
  do i = -(N-1)/2,(N-1)/2
    E(i+(N-1)/2+1) = -2*t*cos(i*2*pi/N)
    k(i+(N-1)/2+1) = i*2*Pi/N
  enddo
 endif
 open(60, file = 'dispanalyt.dat')
 do i = 1,N
  write(60,*) E(i)
 end do
 close(60)
!Note: If Pi is real the division of Pi/N will be real

 call order(E,N) ! to order E for comparison

! ----------- Determination of Eigenvalues and Eigenvectors ---------- !

 call DSYEV('v','L', N, M, LDA, W, WORK, LWORK, INFO)

! ----------- Reordering K with the eigenvalues, respectively -------- !

 call matching

! -------------- Calculation of z-Spin-Spin correlation -------------- !

 call spincor
 call spincoranalyt
 call ascendingorder

! ----------- Ground State energy  and dispersion data --------------- !

 call Dispenergy

 deallocate(W,E,P,K)
 

end program tb1d

subroutine Order(matrix,dims)

 ! B.D.: Order an Array in Ascending order by the bubble sort method.
 
 use Variables
 implicit none
 integer :: INC, dims
 real(real64) :: VY, matrix(dims)

 INC = 1
 do while(INC .LE. dims)
   INC = 3*INC+1
 enddo
 do while(INC .GT. 1)
  INC = INC/3
  do I = INC+1,dims
   VY = matrix(I)
   J = I
   do while(matrix(J-INC) .GT. VY)
    matrix(J) = matrix(J-INC)
    J = J-INC
    if (J .LE. INC) then
     exit
    endif
   enddo
   matrix(J)=VY
  enddo 
 enddo
 
 return
end subroutine Order

subroutine spincor

! B.D.: Gives the spin spin correlation between site i=0 and i=i+1. It's
! necessary to remark thatwe find eigenvalues and eigenvector for both 
! the site index positive and negative. We need to be consistent and
! choose only one. In this case, the positive ones.
! 

 use Variables
 implicit none
 real(real64) :: G
 integer :: l

 open(40, file = 'spincor.dat')
 
 j = 0

 do i = 2,N/2+1 !to avoid self correlation we begin with i = 1
  G = 0
  do l = 1, N/2
   G = G + M(1,l)*M(i,l)
  enddo
   G = -4*G*G
  write(40,*) i-1, G
 enddo
 close(40)

end subroutine spincor

subroutine Matching

! B.D.: Matches the momenta k with how the eigenvalues are printed by the 
! DSYEV function.
 
 use Variables
 implicit none

 if (MOD(N,2) == 0) then !N is even
  P(N) = k(1)
  do i = 2,N,2
     if (i == N) then
      cycle ! next loop
     endif
     P(N-i) = k(1+i/2)
     P(i+1) = k(N/2+i/2+1)
  enddo
  P(1) = k(N/2+1)
  elseif (MOD(N,2) > 0) then !N is odd, we cannot put it in one loop
   P(N) = k(1)
  do i = 2,(N-1),2
     if (i == (N-1)) then
      cycle
     endif
     P(N-i) = k(1+i/2)
  enddo
  do i = 2,(N-1),2
     P(i) = k((N-1)/2+i/2+1)
  enddo
  P(1) = k((N-1)/2+1)
 endif

end subroutine Matching

subroutine Dispenergy

! B.D.: Gives dispersion curve and groundstate of simple TB model, 
! after the ascending order Subroutine is called.

 use Variables
 implicit none

 GS = 0
 do i = N/4+1,3*N/4+1
  GS = GS+W(i)
 enddo

 print*, 'Energy/t = ', 2*GS/N  ! Ground state

 open(30, file = 'disp.dat')
 ! W(i) is the initial one and P(i) ordered with increansing value.
 do i = 1,N
  write(30,*) P(i), W(i)
 enddo
 close(30)

end subroutine Dispenergy

subroutine ascendingorder

! B.D.: Reorder Energy eigenvalue matrix M with increasing   
! momenta P(i) in the 1st Brillouin zone

 use Variables
 implicit none

 do i = 1, UBOUND(M,1)-1
  do j = i+1, UBOUND (M,1)
   if (P(i) > P(j)) then
            temp1 = P(j)
            temp2 = W(j)
            P(j) = P(i)
            W(j) = W(i)
            P(i) = temp1 
            W(i) = temp2
   endif
  enddo
 enddo
end subroutine ascendingorder

subroutine spincoranalyt

 ! B.D.: Calculates the analytic values of z component of the spin-spin 
 ! correlation 
 
 use Variables
 use ISO_FORTRAN_ENV
 implicit none

    !Variables to calculate the analytical value (sum(exp(ikr)))
 complex(real64) :: Im
 complex(real64) :: expo1,expo2,expoa,expob,expo3,expo4
 complex(real64) :: expoanti1,expoanti2,expoanti3,expoanti4,expoantia,expoantib
 real(real64), ALLOCATABLE :: Ga(:)
 Im = (0.0d0,1.0d0) ! definition of i

ALLOCATE(Ga(N/2))

!------------------Analytical Result--------------------
 open(50, file = 'spincoranalyt.dat')

 do i = 1, N/2 ! here we don't need to start from 2 due to analytical expression
  expo1 = EXP(Im*(i*(N+2)*pi)/N)
  expo2 = EXP(Im*(i*2.0*pi)/N)
  expo3 = EXP(-Im*(i*(N+2)*pi)/N)
  expo4 = EXP(-Im*(i*2.0*pi)/N)
  expoa = (1-expo1)/(1-expo2)
  expob = (1-expo3)/(1-expo4)

  expoanti1 = EXP(Im*pi*i)
  expoanti2 = EXP(-Im*(pi*(N-2)*i)/(2*N))
  expoanti3 = EXP(-Im*pi*i)
  expoanti4 = EXP(Im*(pi*(N-2)*i)/(2*N))
  expoantia = (1-expoanti1)*expoanti2/(-1+expo2)
  expoantib = (1-expoanti3)*expoanti4/(-1+expo4)
  Ga(i) = real(-4*expoa*expob/(N*N),8) !periodic
 write(50,*) i,Ga(i) !Store the correlations to be ploted
 enddo
 close(50)
 
 deallocate(Ga)

END SUBROUTINE
