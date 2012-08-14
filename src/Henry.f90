!     Last change:  SW   31 Oct 2005   12:16 pm
!********************************************************************************************
! This program solves Henry's problem
! Psi = sum{sum{A_{m,n}*sin(m*pi*y)*cos(n*pi*x/xi), n=0..infinity}, m=1..infinity}
! Psi = sum{sum{B_{r,s}*sin(s*pi*x/xi)*cos(r*pi*y), s=1..infinity}, r=0..infinity}
! A_{g,h} = frac{sum{B_{r,h}*h*N(g,r),r=0..infinity} + 4/pi*W(g,h)}
!      {epsilon_2*a*pi**2*[g**2+h**2/xi**2]*xi}
! B_{g,h} = frac{pi/4*sum{sum{sum{sum{A_{m,n}*B_{r,s}*(m*s*L*R - n*r*F*G),
!      s=1..infinity},r=0..infinity},n=0..infinity},m=0..infinity}
!      + sum{A_{g,n}*g*N(h,n), n=0..infinity}
!      + epsilon_1*sum{B_{g,s}*s*N(h,s),s=1..infinity} + 4/pi*W(h,g)}
!      {epsilon_1*b*pi**2*[g**2+h**2/xi**2]*xi}
! Where A and B are both two-dimensional matrices.
! A      : Fourier coefficients for Psi
! B      : Fourier coefficients for C
! a=Q/(k_1*d): dimensionless
! b=D/Q    : dimensionless
! d      : thickness of aquifer, m
! D      : dispersion coefficient, m**2/sec
! E      : dimensionless coefficient relating the concentration to
!      : the density of the solution
! k      : permeability of sand, m**2
! kbar=k*rho_0*g/mu  : transmission coefficient of aquifer, m/sec
! k_1=kbar*(rho_s-rho_o)/rho_o : transmission coefficient times the density-difference
!      : ratio
! l      : length of aquifer, m
! Q      : net freshwater discharge per unit length of beach,
!      : m**2/sec
! rho    : density of solution, kg/m**3
! rho_o    : density of freshwater, kg/m**3
! rho_s    : density of saltwater, kg/m**3
! mu    : viscosity of water, N-sec/m**2
! xi=l/d  : aspect ratio
!********************************************************************************************
PROGRAM Henrys_Problem
  
  USE Henry_funcs
  USE Henry_sums
  USE LU
  IMPLICIT NONE

  REAL, DIMENSION(:, :), ALLOCATABLE :: MATRIX_A, MATRIX_B, Matrix !Matrices of Fourier Coefficients
  REAL, DIMENSION(:, :), ALLOCATABLE :: quad, linear !Quadratic and linear terms
  REAL, DIMENSION(:, :), ALLOCATABLE :: Psi, C !Stramlines and Isochlors
  INTEGER, DIMENSION(:), ALLOCATABLE :: indx
  REAL, DIMENSION(:), ALLOCATABLE :: x, y !Grid
  REAL, PARAMETER :: a = 0.263, b = 0.1, d = 1.0, l = 2.0, dx = 0.05, dy= 0.05 !Problem parameters
  REAL :: xi, pid4, fourdpi, xi2, api2, bpi2 !Constants
  REAL :: Psi_err, C_err, ERROR
  INTEGER :: h_a(5), h_b(0:4) !Size of h depends on both i and which coefficient A or B
  INTEGER :: i_a, i_b, j_a, j_b, i_y, j_x !Size of the Matrices
  INTEGER :: g, h, m, n, p, q, r, s !Loop counters
  INTEGER :: AllocateStatus !Status variable for ALLOCATE
  REAL :: EPSILON = 5E-4

  !WRITE (*,*), "Enter in the size of aquifer (d = thickness l = length) in meters."
  !READ *, d, l
  !WRITE (*,*), "Enter in the partition size (dx, dy)."
  !READ *, dx, dy
  !WRITE (*,*), "Now enter in your values for a and b."
  !READ *, a, b

  OPEN(30, FILE='Psi0_1.txt')
  OPEN(31,FILE='C0_1.txt')
  !OPEN (32, "B2h.txt")

  h_a = (/ 15, 10, 5, 3, 2 /)
  h_b = (/ 20, 10, 5, 3, 2 /)

  i_a = 5
  i_b = 4
  j_a = MAXVAL(h_a)
  j_b = MAXVAL(h_b)
  i_y = NINT(d/dy)
  j_x = NINT(l/dx)
  xi = l/d

  pid4 = pi/4.
  fourdpi = 4./pi
  xi2 = xi**2
  api2 = a*pi**2
  bpi2 = b*pi**2



  ALLOCATE(MATRIX_A(i_a, 0:j_a), MATRIX_B(0:i_b, j_b), quad(0:i_b, j_b), linear(i_b, j_b), Psi(0:i_y, 0:j_x), &
    C(0:i_y, 0:j_x), x(0:j_x), y(0:i_y), Matrix(i_a*(j_a+1)+(i_b+1)*j_b,i_a*(j_a+1)+(i_b+1)*j_b), STAT = AllocateStatus)
  IF(AllocateStatus /= 0) STOP "*** NOT ENOUGH MEMORY ***"

  MATRIX_A = 0. !Initializes array so that A = B = 0
  MATRIX_B = 0. 

  x_loop0: DO h = 0, j_x     !Initializing the values for x
    x(h) = dx*h
  END DO x_loop0

  y_loop0: DO g = 0, i_y     !Initializing the values for y
    y(g) = dy*g
  END DO y_loop0
  
  error = 1. !Initialize error terms
  WRITE (*,*) 'b=', b

  CALL build_system() 

  CALL WRITETOFILE(Psi, C, x, y, i_y, j_x)

  CLOSE (30)
  CLOSE (31)
  !CLOSE (32)

  DEALLOCATE(MATRIX_A, MATRIX_B, quad, linear, Psi, C, x, y)

!**************************************** Begin supporting FUNCTIONS ************************
        
  CONTAINS

  SUBROUTINE build_system()
    INTEGER :: i, j
    
    DO i=1, SIZE(Matrix,1)
      DO j=1, SIZE(Matrix,2)
        IF(i .LE. i_a*(j_b+1) .AND. j .LE. i_a*(j_b+1)) THEN
          Matrix(i,j) = eps(j-1)*api2*(i**2 + (j-1)**2/xi2)*xi
        ELSE IF(i .LE. i_a*(j_b+1)) THEN

        ELSE IF(j .LE. i_a*(j_b+1) THEN

        ELSE

        END IF
      END DO
    END DO

  END SUBROUTINE

  SUBROUTINE WRITETOFILE(Psi, C, x, y, i_y, j_x)

    INTEGER, INTENT(IN) :: i_y, j_x
    REAL, DIMENSION(0:i_y, 0:j_x), INTENT(IN) :: Psi, C
    REAL, DIMENSION(0:j_x), INTENT(IN) :: x
    REAL, DIMENSION(0:i_y), INTENT(IN) :: y
    INTEGER :: g, h

101     FORMAT(1x, 3(F7.3,3x))

!    WRITE (30,*) Psi
!    WRITE (31,*) C
!    WRITE (90,*) x
!    WRITE (91,*) y
    DO h = 0, j_x !Write Psi and C out to file
      DO g = 0, i_y
        WRITE (30, 101) x(h), y(g), Psi(g, h)
        WRITE (31, 101) x(h), y(g), C(g, h)
      END DO
      WRITE (30,*)
      WRITE (31,*)
    END DO

  END SUBROUTINE

END PROGRAM Henrys_Problem
