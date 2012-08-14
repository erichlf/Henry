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

  REAL, DIMENSION(:, :), ALLOCATABLE :: MATRIX_A, MATRIX_B !Matrices of Fourier Coefficients
  REAL, DIMENSION(:, :), ALLOCATABLE :: B0h, Bl1h, Bl2h !Matrices
  REAL, DIMENSION(:, :), ALLOCATABLE :: quad, linear !Quadratic and linear terms
  REAL, DIMENSION(:, :), ALLOCATABLE :: Psi, C, Psi_old, C_old !Stramlines and Isochlors
  REAL, DIMENSION(:), ALLOCATABLE :: B0h_W, Bl1h_W, Bl2h_W !Constant vectors containing lin, quad and W_FUNC terms
  INTEGER, DIMENSION(:), ALLOCATABLE :: indx0, indx1, indx2
  REAL, DIMENSION(:), ALLOCATABLE :: x, y !Grid
  REAL, PARAMETER :: a = 0.263, b = 0.1, d = 1.0, l = 2.0, dx = 0.05, dy= 0.05 !Problem parameters
  REAL :: xi, pid4, fourdpi, xi2, api2, bpi2 !Constants
  REAL :: Psi_err, C_err, ERROR
  INTEGER :: h_a(5), h_b(0:4) !Size of h depends on both i and which coefficient A or B
  INTEGER :: i_a, i_b, j_a, j_b, i_y, j_x, Bl1h_size, Bl2h_size !Size of the Matrices
  INTEGER :: g, h, m, n, p, q, r, s, count, start !Loop counters
  INTEGER :: AllocateStatus !Status variable for ALLOCATE
  REAL :: EPSILON = 5E-4
  INTEGER :: loop = 4
  REAL :: d0, d1, d2

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

  Bl1h_size = h_b(1) + h_b(2)
  Bl2h_size = h_b(3) + h_b(4)

  ALLOCATE(MATRIX_A(i_a, 0:j_a), MATRIX_B(0:i_b, j_b), quad(0:i_b, j_b), linear(i_b, j_b), Psi(0:i_y, 0:j_x), &
    Psi_old(0:i_y, 0:j_x), C(0:i_y, 0:j_x), C_old(0:i_y, 0:j_x), x(0:j_x), y(0:i_y), STAT = AllocateStatus)
  ALLOCATE(B0h(h_b(0), h_b(0)), B0h_W(h_b(0)), Bl1h(Bl1h_size, Bl1h_size), Bl1h_W(Bl1h_size), Bl2h(Bl2h_size, Bl2h_size), &
    Bl2h_W(Bl2h_size), indx0(h_b(0)), indx1(Bl1h_size), indx2(Bl2h_size), STAT = AllocateStatus)
  IF(AllocateStatus /= 0) STOP "*** NOT ENOUGH MEMORY ***"

  MATRIX_A = 0. !Initializes array so that A = B = 0
  MATRIX_B = 0. 
  Psi_old = 1.
  C_old = 1.

  x_loop0: DO h = 0, j_x     !Initializing the values for x
    x(h) = dx*h
  END DO x_loop0

  y_loop0: DO g = 0, i_y     !Initializing the values for y
    y(g) = dy*g
  END DO y_loop0
  
  DO count = 0, 0!5!6
    error = 1. !Initialize error terms
    bpi2 = (b - count*0.01)*pi**2
    WRITE (*,*) 'b=', b-count*0.01
    DO g = 0, 4 !Initialize matrices and take their inverses
      IF (g == 0 .OR. (g/2)*2 /= g) THEN 
        DO m = 1, h_b(g)
          h = m
          DO n = 1, h_b(g) !Non-diagonal terms
            IF (m /= n) THEN
              s = n
              IF (g == 0) THEN
                B0h(m, n) = -eps(g)*s*N_FUNC(h, s)
              ELSE IF (g == 1) THEN
                Bl1h(m, n) = -eps(g)*s*N_FUNC(h, s)
              ELSE
                Bl2h(m, n) = -eps(g)*s*N_FUNC(h, s)
              END IF
            END IF
          END DO

          IF (g == 0) THEN !Diagonal terms
            B0h(m, m) = eps(g)*bpi2*(g**2 + h**2/xi2)*xi
          ELSE IF (g == 1) THEN
            Bl1h(m, m) = eps(g)*bpi2*(g**2 + h**2/xi2)*xi
          ELSE
            Bl2h(m, m) = eps(g)*bpi2*(g**2 + h**2/xi2)*xi
          END IF
        END DO
      ELSE
        DO m = h_b(g - 1) + 1, h_b(g - 1) + h_b(g)
          h = m - h_b(g - 1)
          DO n = h_b(g - 1) + 1, h_b(g - 1) + h_b(g)
            IF (m /= n) THEN !Diagonal terms
              s = n - h_b(g - 1)
              IF (g == 2) THEN
                Bl1h(m, n) = -eps(g)*s*N_FUNC(h, s)
              ELSE
                Bl2h(m, n) = -eps(g)*s*N_FUNC(h, s)
              END IF
            END IF
          END DO
          IF (g == 2) THEN !Diagonal terms
            Bl1h(m, m) = eps(g)*bpi2*(g**2 + h**2/xi2)*xi
          ELSE
            Bl2h(m, m) = eps(g)*bpi2*(g**2 + h**2/xi2)*xi
          END IF
        END DO
      END IF
    END DO

    CALL ludcmp(B0h, h_b(0), h_b(0), indx0, d0)
    CALL ludcmp(Bl1h, Bl1h_size, Bl1h_size, indx1, d1)
    CALL ludcmp(Bl2h, Bl2h_size, Bl2h_size, indx2, d2)
    !WRITE (*,*) 'b=', b
        
    Error_Loop: DO
      WRITE (*,*) "Error=", error
      IF (error <= Epsilon) THEN !Terminate loop when ABS(Error) <= Epsilon
        EXIT Error_Loop
      ELSE
        DO q = 1, loop
          DO g = 0, i_b !Update linear and quadratic terms in B(g, h)
             DO h = 1, h_b(g)
               quad(g, h) = pid4*SUM_AB(MATRIX_A, MATRIX_B, i_a, j_a, i_b, j_b, g, h) !Quadratic terms
             END DO
          END DO

          g = 0
          DO m = 1, h_b(g) !Setup the linear system B0h
            h = m
            !Constant terms for B0h including Quadratic and Linear terms
            B0h_W(m) = fourdpi*W_FUNC(h, g) + quad(g, h)
          END DO
          CALL lubksb(B0h, h_b(g), h_b(g), indx0, B0h_W)

          DO h = 1, h_b(g)
            MATRIX_B(g, h) = B0h_W(h) !Place Boh_solution into corresponding B(g, h)
          END DO

          DO g = 1, i_a !Update A(g, h) with new B(0, h)
            IF ((g*2)/2 /= g) THEN
              start = 1
            ELSE 
              start = 0
            END IF
            DO h = start , h_a(g)
              MATRIX_A(g, h) = 1./(eps(h)*api2*(g**2 + h**2/xi2)*xi)* &
                (SUM_BN(MATRIX_B, i_b, j_b, g, h) + fourdpi*W_FUNC(g, h))
            END DO
          END DO
        END DO

        DO q = 1, loop
          DO g = 0, i_b !Update linear and quadratic terms in B(g, h)
            DO h = 1, h_b(g)
              IF (g > 0) THEN
                Linear(g, h) = SUM_AN(MATRIX_A, i_a, j_a, g, h)
              END IF
              quad(g, h) = pid4*SUM_AB(MATRIX_A, MATRIX_B, i_a, j_a, i_b, j_b, g, h) !Quadratic terms
            END DO
          END DO

          DO m = 1, Bl1h_size
            IF (m <= h_b(1)) THEN
              g = 1
              h = m
            ELSE
              g = 2
              h = m - h_b(1)
            END IF
            !Constant terms for Bl1h including Quadratic and Linear terms
            Bl1h_W(m) = fourdpi*W_FUNC(h, g) + quad(g, h) + linear(g, h)
          END DO

          CALL lubksb(Bl1h, Bl1h_size, Bl1h_size, indx1, Bl1h_W)

          DO m = 1, Bl1h_size
            IF (m <= h_b(1)) THEN
              g = 1
              h = m
            ELSE
              g = 2
              h = m - h_b(1)
            END IF
            MATRIX_B(g, h) = Bl1h_W(m) !Place Boh_solution into corresponding B(g, h)
          END DO

          DO g = 1, i_a !Update A(g, h) with new B(0, h)
            IF ((g*2)/2 /= g) THEN
              start = 1
            ELSE 
              start = 0
            END IF
            DO h = 1 , h_a(g)
              MATRIX_A(g, h) = 1./(eps(h)*api2*(g**2 + h**2/xi2)*xi)* &
                (SUM_BN(MATRIX_B, i_b, j_b, g, h) + fourdpi*W_FUNC(g, h))
            END DO
          END DO
        END DO

        DO q = 1, loop
          DO g = 0, i_b !Update linear and quadratic terms in B(g, h)
            DO h = 1, h_b(g)
              IF (g > 0) THEN
                Linear(g, h) = SUM_AN(MATRIX_A, i_a, j_a, g, h)
              END IF
              quad(g, h) = pid4*SUM_AB(MATRIX_A, MATRIX_B, i_a, j_a, i_b, j_b, g, h) !Quadratic terms
            END DO
          END DO

          DO m = 1, Bl2h_size
            IF (m <= h_b(3)) THEN
              g = 3
              h = m
            ELSE
              g = 4
              h = m - h_b(3)
            END IF
            !Constant terms for Bl2h including Quadratic and Linear terms
            Bl2h_W(m) = fourdpi*W_FUNC(h, g) + quad(g, h) + linear(g, h)
          END DO

          CALL lubksb(Bl2h, Bl2h_size, Bl2h_size, indx2, Bl2h_W)

          DO m = 1, Bl2h_size
            IF (m <= h_b(3)) THEN
              g = 3
              h = m
            ELSE
              g = 4
              h = m - h_b(3)
            END IF
            MATRIX_B(g, h) = Bl2h_W(m) !Place Boh_solution into corresponding B(g, h)
          END DO

          DO g = 1, i_a !Update A(g, h) with new B(0, h)
            IF ((g*2)/2 /= g) THEN
              start = 1
            ELSE 
              start = 0
            END IF
            DO h = start , h_a(g)
              MATRIX_A(g, h) = 1./(eps(h)*api2*(g**2 + h**2/xi2)*xi)* &
                      (SUM_BN(MATRIX_B, i_b, j_b, g, h) + fourdpi*W_FUNC(g, h))
            END DO
          END DO
        END DO
                        
        x_loop: DO g = 0, i_y !Find Psi(g, h) and C(g, h)
          y_loop: DO h = 0, j_x
            IF (g == 0 .AND. h == 0) THEN
              Psi(g, h) = 0.
              C(g, h) = 0.
            ELSE
              Psi(g, h) = SUM_A(MATRIX_A, x, y, xi, g, h, i_a, j_a, i_y, j_x) + y(g)/d
              C(g, h) = SUM_B(MATRIX_B, x, y, xi, g, h, i_b, j_b, i_y, j_x) + x(h)/(d*xi)
            END IF
          END DO y_loop
        END DO x_loop

        !Find the error
        C_err = MAXVAL(ABS(C - C_old))
        Psi_err = MAXVAL(ABS(Psi - Psi_old))

        error = MAX(Psi_err, C_err)
        C_old = C
        Psi_old = Psi
      END IF
    END DO Error_Loop
  END DO

  CALL WRITETOFILE(Psi, C, x, y, i_y, j_x)

  CLOSE (30)
  CLOSE (31)
  !CLOSE (32)

  DEALLOCATE(MATRIX_A, MATRIX_B, quad, linear, Psi, C, Psi_old, C_old, x, y)
  DEALLOCATE(B0h, B0h_W, indx0, Bl1h, Bl1h_W, indx1, Bl2h, Bl2h_W, indx2)

!**************************************** Begin supporting FUNCTIONS ************************
        
  CONTAINS

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
