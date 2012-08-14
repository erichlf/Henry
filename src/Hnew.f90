!     Last change:  SW    6 Jun 2005    1:39 pm
!********************************************************************************************
! This program solves Henry's problem
!       Psi = sum{sum{A_{m,n}*sin(m*pi*y)*cos(n*pi*x/xi), n=0..infinity}, m=1..infinity}
!       Psi = sum{sum{B_{r,s}*sin(s*pi*x/xi)*cos(r*pi*y), s=1..infinity}, r=0..infinity}
!	A_{g,h} = frac{sum{B_{r,h}*h*N(g,r),r=0..infinity} + 4/pi*W(g,h)}
!			{epsilon_2*a*pi**2*[g**2+h**2/xi**2]*xi}
!       B_{g,h} = frac{pi/4*sum{sum{sum{sum{A_{m,n}*B_{r,s}*(m*s*L*R - n*r*F*G),
!			s=1..infinity},r=0..infinity},n=0..infinity},m=0..infinity}
!			+ sum{A_{g,n}*g*N(h,n), n=0..infinity}
!			+ epsilon_1*sum{B_{g,s}*s*N(h,s),s=1..infinity} + 4/pi*W(h,g)}
!			{epsilon_1*b*pi**2*[g**2+h**2/xi**2]*xi}
! Where A and B are both two-dimensional matrices.
!       A					: Fourier coefficients for Psi
!       B					: Fourier coefficients for C
!       a=Q/(k_1*d)				: dimensionless
!       b=D/Q					: dimensionless
!       d					: thickness of aquifer, m
!       D					: dispersion coefficient, m**2/sec
!       E					: dimensionless coefficient relating the concentration to
!						: the density of the solution
!       k					: permeability of sand, m**2
!       kbar=k*rho_0*g/mu			: transmission coefficient of aquifer, m/sec
!       k_1=kbar*(rho_s-rho_o)/rho_o		: transmission coefficient times the density-difference
!						: ratio
!       l					: length of aquifer, m
!       Q					: net freshwater discharge per unit length of beach,
!						: m**2/sec
!       rho					: density of solution, kg/m**3
!       rho_o					: density of freshwater, kg/m**3
!       rho_s					: density of saltwater, kg/m**3
!       mu					: viscosity of water, N-sec/m**2
!       xi=l/d					: aspect ratio
!********************************************************************************************
PROGRAM Henrys_Problem !Use Newtown-Raphson to solve Henry's problem.

  USE Henry_funcs
  USE Henry_sums
!  USE linpacks
  USE LU
  IMPLICIT NONE

  REAL, DIMENSION(:, :), ALLOCATABLE :: MATRIX_A, MATRIX_B !Matrices of Fourier Coefficients
  REAL, DIMENSION(:, :), ALLOCATABLE :: Psi, C !Stramlines and Isochlors
  REAL :: a = 0.263, b = 0.04, xi, d = 1.0, l = 2.0, dx = 0.05, dy= 0.05 !Problem parameters
  REAL :: pid4, fourdpi, xi2, api2, bpi2 !Constants
  INTEGER*4 :: h_a(5), h_b(0:4) !Size of h depends on both i and which coefficient A or B
  INTEGER*4 :: i_a, i_b, j_a, j_b, i_y, j_x, i !Size of the Matrices
  INTEGER*4 :: g, h, m, n, o, p, q, r, s !Loop counters
  INTEGER*4 :: AllocateStatus !Status variable for ALLOCATE
  INTEGER*4 :: linearsize
  REAL :: EPSILON = 5E-4

!  COMMON :: a, b, xi, d, l, dx, dy, pid4, fourdpi, xi2, api2, bpi2
!  COMMON :: i_a, i_b, j_a, j_b, i_y, j_x, linearsize
!  COMMON :: EPSILON

  OPEN (30, FILE="Psi0_05.txt")
  OPEN (31, FILE="C0_05.txt")
  OPEN (32, FILE="A.txt")
  OPEN (33, FILE="B.txt")
  OPEN (34, FILE="A0_05.txt")
  OPEN (35, FILE="B0_05.txt")
  OPEN (40, FILE="alpha0_05.txt")
  OPEN (90, FILE="x_05.txt")
  OPEN (91, FILE="y_05.txt")

  !h_a = (/ 15, 10, 10, 10, 10 /)! !Size of each row of A and B
  !h_b = 25!(/ 20, 10, 5, 3, 2 /)!
  h_a = (/ 15, 10, 5, 3, 2 /)
  h_b = (/ 20, 10, 5, 3, 2 /)

  i_a = 5
  i_b = 4
  j_a = MAXVAL(h_a) !Size of largest row in A
  j_b = MAXVAL(h_b) !Size of largest row in B
  i_y = NINT(d/dy) !Number of partitions in y for grid
  j_x = NINT(l/dx) !Number of partitions in x for grid
  xi = l/d

  pid4 = pi/4.
  fourdpi = 4./pi
  xi2 = xi**2
  api2 = a*pi**2

  linearsize = i_a*(j_a + 1) + (i_b + 1)*j_b

  ALLOCATE(MATRIX_A(i_a, 0:j_a), MATRIX_B(0:i_b, j_b), Psi(0:i_y, 0:j_x), C(0:i_y, 0:j_x), STAT = AllocateStatus)
  IF(AllocateStatus /= 0) THEN
    STOP "*** NOT ENOUGH MEMORY ***"
  END IF

  MATRIX_A = 0.
  MATRIX_B = 0.
!  READ (32, *) MATRIX_A !Reads in A and B from file for initial guess
!  READ (33, *) MATRIX_B

!	DO i = 0, 4
!    IF (i /= 0) THEN
!      b = b - 0.01
!    END IF
    bpi2 = b*pi**2
    WRITE (*, *) "b=", b
    CALL NEWTON() !Performs Newton-Raphson
!	END DO

  CALL PSIANDC(Psi, C) !Assigns Psi and C based on A and B
  CALL Dispersivity()

  WRITE (34, *) MATRIX_A !Writes A and B to file
  WRITE (35, *) MATRIX_B

  CLOSE (30)
  CLOSE (31)
  CLOSE (32)
  CLOSE (33)
  CLOSE (34)
  CLOSE (35)

  DEALLOCATE(MATRIX_A, MATRIX_B, Psi, C)

!**************************************** Begin supporting FUNCTIONS ************************
        
  CONTAINS

  FUNCTION dFdA(g, h, m, n) !A derivative of quadratic portion of F
    INTEGER*4, INTENT(IN) :: g, h, m, n
    INTEGER*4 :: r, s
    REAL :: dFdA

    dFdA = 0. !Initialize derivative to 0

    DO r = 0, i_b
      DO s = 1, j_b
        dFdA = dFdA + MATRIX_B(r,s)*(REAL(m*s*L_FUNC(m, r, g)*R_FUNC(h, n, s)) &
          - REAL(n*r*F_FUNC(m, r, g)*G_FUNC(h, n, s)))
      END DO
    END DO
  END FUNCTION dFdA

  FUNCTION dFdB(g, h, r, s) !B derivative of quadratic portion of F
    INTEGER*4, INTENT(IN) :: g, h, r, s
    INTEGER*4 :: m, n
    REAL :: dFdB

    dFdB = 0. !Initialize derivative to 0

    DO m = 1, i_a
      DO n = 0, j_a
        dFdB = dFdB + MATRIX_A(m,n)*(REAL(m*s*L_FUNC(m, r, g)*R_FUNC(h, n, s)) &
          - REAL(n*r*F_FUNC(m, r, g)*G_FUNC(h, n, s)))
      END DO
    END DO
  END FUNCTION dFdB

  FUNCTION F() RESULT(F_VECTOR)
    REAL, DIMENSION (:), ALLOCATABLE :: F_VECTOR
    INTEGER*4 :: o, q

    ALLOCATE(F_VECTOR(linearsize), STAT = AllocateStatus)
    IF(AllocateStatus /= 0) THEN
      STOP "*** NOT ENOUGH MEMORY ***"
    END IF

    F_VECTOR = 0.

    DO o = 1, linearsize !Place values in both F and DF
      IF (o <= i_a*(j_a + 1)) THEN !A portion of F and DF
        CALL AGANDH(o, g, h) !Determine what g and h are for this row
        F_vector(o) = eps(h)*api2*MATRIX_A(g,h)*(g**2 + h**2/xi2)*xi &
          - SUM_BN(MATRIX_B, i_b, j_b, g, h) - fourdpi*W_FUNC(g, h) !A portion of F
      ELSE !B portion of F and DF
        CALL BGANDH(o, g, h) !Determine what g and h are for this row
        F_vector(o) = eps(g)*bpi2*MATRIX_B(g, h)*(g**2 + h**2/xi2)*xi &
          - SUM_AN(MATRIX_A, i_a, j_a, g, h) - eps(g)*SUM_BsN(MATRIX_B, i_b, j_b, g, h) &
          - pid4*SUM_AB(MATRIX_A, MATRIX_B, i_a, j_a, i_b, j_b, g, h) - fourdpi*W_FUNC(h, g) !B portion of F
      END IF
    END DO
  END FUNCTION F

  FUNCTION DF() RESULT(DF_MATRIX)
    REAL, DIMENSION (:, :), ALLOCATABLE :: DF_MATRIX
    INTEGER*4 :: o, p, q, r, s

    ALLOCATE(DF_MATRIX(linearsize,linearsize), STAT = AllocateStatus)
    IF(AllocateStatus /= 0) THEN
      STOP "*** NOT ENOUGH MEMORY ***"
    END IF

    DF_MATRIX = 0. !Initializes DF

    DO o = 1, linearsize !Place values in both F and DF
      IF (o <= i_a*(j_a + 1)) THEN !A portion of F and DF
        CALL AGANDH(o, g, h) !Determine what g and h are for this row
        DO p = i_a*(j_a + 1) + 1, linearsize
          CALL BGANDH(p, r, s) !Determine what r and s are for this row and column
          IF (s == h) THEN
            IF (r /= g) THEN
              DF_MATRIX(o, p) = -h*N_FUNC(g, r) !B derivatives of A portion of F placed into DF
            END IF
          END IF
        END DO
        DF_MATRIX(o, o) = eps(h)*api2*(g**2 + h**2/xi2)*xi !A derivative of A portion of F placed into DF
      ELSE !B portion of F and DF
        CALL BGANDH(o, g, h) !Determine what g and h are for this row
        DO p = 1, linearsize
          IF (p <= i_a*(j_a + 1)) THEN !A derivatives of the B portion of F placed into DF
            CALL AGANDH(p, m, n) !Determine what m and n are for this row and column
            IF (m == g) THEN
              IF (n /= h) THEN
                DF_MATRIX(o, p) = -g*N_FUNC(h, n) &
                  - pid4*dFdA(g, h, m, n) !A derivatives of B portion of F placed into DF
              ELSE
                DF_MATRIX(o, p) = -pid4*dFdA(g, h, m, n) !A derivatives of B portion of F placed into DF
              END IF
            ELSE
              DF_MATRIX(o, p) = -pid4*dFdA(g, h, m, n) !A derivatives of B portion of F placed into DF
            END IF
          ELSE !B derivatives of the B portion F placed into DF
            CALL BGANDH(p, r, s) !Determine what r and s are for this row and column
            IF (r == g) THEN
              IF (s /= h) THEN
                DF_MATRIX(o, p) = -eps(g)*s*N_FUNC(h, s) &
                  - pid4*dFdB(g, h, r, s) !B derivatives of B portion of F placed into DF
              ELSE
                DF_MATRIX(o, p) = eps(g)*bpi2*(g**2 + h**2/xi2)*xi &
                  - pid4*dFdB(g, h, r, s)
              END IF
            ELSE
              DF_MATRIX(o, p) = -pid4*dFdB(g, h, r, s) !B derivatives of B portion of F placed into DF
            END IF
          END IF
        END DO
        DF_MATRIX(o, o) = eps(g)*bpi2*(g**2 + h**2/xi2)*xi - pid4*dFdB(g, h, g, h) !Diagonals of the B portion of DF
      END IF
    END DO
  END FUNCTION DF

  SUBROUTINE AGANDH(o, g, h) !Determines g and h given o inside A section of matrix or vector
    INTEGER*4, INTENT(IN) :: o
    INTEGER*4, INTENT(OUT) :: g, h
    INTEGER*4 :: q

    DO q = 1, i_a !Determine g and h
      IF (o >= (q - 1)*j_a + q .AND. o <= q*j_a + q) THEN
        g = q
      END IF
    END DO
    h = MOD(o, j_a + 1) - 1
    IF (h < 0) THEN
      h = j_a
    END IF
  END SUBROUTINE AGANDH

  SUBROUTINE BGANDH(o, g, h) !Determines g and h given o inside B section of matrix or vector
    INTEGER*4, INTENT(IN) :: o
    INTEGER*4, INTENT(OUT) :: g, h
    INTEGER*4 :: q, r

    r = o - i_a*(j_a + 1)
    DO q = 0, i_b
      IF (r >= q*j_b + 1 .AND. r <= (q + 1)*j_b) THEN
        g = q
      END IF
    END DO
    h = MOD(r, j_b)
    IF (h == 0) THEN
      h = j_b
    END IF
  END SUBROUTINE BGANDH

  SUBROUTINE NEWTON()
    REAL, DIMENSION (linearsize) :: F_VECTOR, Solution_new, Solution_old
    REAL, DIMENSION (linearsize, linearsize) :: DF_MATRIX
    INTEGER*4, DIMENSION(linearsize) :: indx
    INTEGER :: info, job
    REAL :: d0, ERROR

    Solution_new = 0. !Initializes Solution_new
    CALL ASSIGNSOLOLD(Solution_old) !Initial guess for solution

    ERROR = 1.

    DO
      F_VECTOR = F() !Assigns F
      DF_MATRIX = DF() !Assigns DF

!      CALL sgefa(DF_MATRIX,linearsize,linearsize,indx,info)
!      IF( info /= 0 ) THEN
!        WRITE ( *, '(a)' ) ' '
!        WRITE ( *, '(a,i6)' ) '  SGEFA returned an error flag INFO = ', info
!        RETURN
!      END IF
!      job = 0
!      CALL sgesl (DF_MATRIX,linearsize,linearsize,indx,F_VECTOR,job)
      CALL ludcmp(DF_MATRIX, linearsize, linearsize, indx, d0) !Replaces DF with its LU decomposition
      CALL lubksb(DF_MATRIX, linearsize, linearsize, indx, F_VECTOR) !Solves DF*F=F using ludcmp
!      CALL REGULARIZATION(DF_MATRIX, linearsize, F_VECTOR)

      Solution_new = Solution_old - F_VECTOR !Finds new solution

      CALL ASSIGNAANDB(Solution_new) !Translates solution into A and B

      !Find the L2 error for this step of Newton's method
      ERROR = L2NORM(F_VECTOR, linearsize) 

      Solution_old = Solution_new
      F_VECTOR = F() !Reevaluates F since it was replaced by solution to DF*F=F
      
      WRITE (*, *) "error=", ERROR !Writes error to screen
      IF (ERROR < EPSILON) THEN !Checks for convergence
        !Check second Criterion
        IF (L2NORM(F_vector,linearsize) < EPSILON*100) THEN
          EXIT
        ELSE !States why it still has not converged
          WRITE (*,*) "Second Criterion not met L2NORM(F)=", L2NORM(F_VECTOR,linearsize)
        END IF
      END IF
    END DO
  END SUBROUTINE NEWTON

  FUNCTION L2NORM(F, n) RESULT(Norm)
    INTEGER*4, INTENT(IN) :: n
    REAL, DIMENSION(n), INTENT(IN) :: F
    REAL :: Norm
    INTEGER*4 :: i

    Norm = 0.

    DO i = 1, n
      Norm = Norm + (F(i))**2
    END DO

    Norm = Norm**0.5
  END FUNCTION L2NORM

  SUBROUTINE ASSIGNAANDB(Solution_new)
    REAL, DIMENSION(linearsize), INTENT(IN) :: Solution_new
    INTEGER*4 :: g, h, o, q

    DO o = 1, linearsize
      IF (o <= i_a*(j_a + 1)) THEN !A portion of Solution_new
        CALL AGANDH(o, g, h) !Finds g and h corresponding to counter o
        MATRIX_A(g, h) = Solution_new(o)
      ELSE !B portion of Solution_new
        CALL BGANDH(o, g, h) !Finds g and h corresponding to counter o
        MATRIX_B(g, h) = Solution_new(o)
      END IF
    END DO
  END SUBROUTINE ASSIGNAANDB

  SUBROUTINE ASSIGNSOLOLD(Solution_old)
    REAL, DIMENSION(linearsize), INTENT(OUT) :: Solution_old
    INTEGER*4 :: g, h, o, q

    DO o = 1, linearsize
      IF (o <= i_a*(j_a + 1)) THEN !A portion of Solution_new
        CALL AGANDH(o, g, h) !Finds g and h corresponding to counter o
        Solution_old(o) = MATRIX_A(g, h)
      ELSE !B portion of Solution_new
        CALL BGANDH(o, g, h) !Finds g and h corresponding to counter o
        Solution_old(o) = MATRIX_B(g, h)
      END IF
    END DO
  END SUBROUTINE ASSIGNSOLOLD

  SUBROUTINE PSIANDC(Psi, C)
    REAL, DIMENSION(0:i_y, 0:j_x) :: Psi, C
    REAL, DIMENSION(0:j_x) :: x
    REAL, DIMENSION(0:i_y) :: y
    INTEGER*4 :: g, h

    DO h = 0, j_x     !Initializing the values for x
      x(h) = dx*h
    END DO

    DO g = 0, i_y     !Initializing the values for y
            y(g) = dy*g
    END DO

    DO g = 0, i_y !Find Psi(g, h) and C(g, h)
      DO h = 0, j_x
        IF (g == 0 .AND. h == 0) THEN
          Psi(g, h) = 0.
          C(g, h) = 0.
        ELSE
          Psi(g, h) = SUM_A(MATRIX_A, x, y, xi, g, h, i_a, j_a, i_y, j_x) + y(g)/d
          C(g, h) = SUM_B(MATRIX_B, x, y, xi, g, h, i_b, j_b, i_y, j_x) + x(h)/(d*xi)
        END IF
      END DO
    END DO

    CALL WRITETOFILE(x, y) !Writes Psi and C out to file
  END SUBROUTINE PSIANDC

  SUBROUTINE WRITETOFILE(x, y)
    REAL, DIMENSION(0:j_x), INTENT(IN) :: x
    REAL, DIMENSION(0:i_y), INTENT(IN) :: y
    INTEGER*4 :: g, h
      
    WRITE(30,*) "#  X  Y  Z"
    WRITE(31,*) "#  X  Y  Z"
    DO h =0, j_x
      DO g = 0, i_y
        WRITE(30,101) x(h), y(g), Psi(g,h)
        WRITE(31,101) x(h), y(g), C(g,h)
      END DO
      WRITE(30,*) ""
      WRITE(30,*) ""
    END DO

    !DO h = 0, j_x !Write Psi and C out to file
!      WRITE (91,*) y
    !  DO g = 0, i_y
!        WRITE (90,*) x
    !  END DO
    !END DO

    101     FORMAT(1x, 3(F7.3,3x))
  END SUBROUTINE

  !Regularizes the Matrix DF, so it is less singular and solve ((D^-1/2)*AT*A*(D^-1/2) + alpha*I)*(D^1/2)x=(D^-1/2)*AT*b
  SUBROUTINE REGULARIZATION(DF_MATRIX, F_VECTOR) 
    REAL, DIMENSION (linearsize) :: F_VECTOR
    REAL, DIMENSION (linearsize, linearsize) :: DF_MATRIX, DF_T, D, alphaI
    INTEGER*4, DIMENSION(linearsize) :: indx
    REAL :: d0, alpha = 0.
    INTEGER*4 :: i

    DO i = 1, linearsize !Creates the diagonal matrix of DF and takes the square root then inverts it and creates alpha*I
      D(i,i) = (DF_MATRIX(i,i))**(-0.5)
      alphaI(i,i) = alpha
    END DO

    !Regularize the matrix DF
    DF_T = TRANSPOSE(DF_MATRIX)
    DF_MATRIX = MATMUL(DF_MATRIX, D)
    DF_MATRIX = MATMUL(DF_T, DF_MATRIX)
    DF_MATRIX = MATMUL(D, DF_MATRIX)
    DF_MATRIX = DF_MATRIX + alphaI

    !Regularize the vector F
    F_VECTOR =  MATMUL(DF_T, F_VECTOR)
    F_VECTOR =  MATMUL(F_VECTOR, D)

    CALL ludcmp(DF_MATRIX, linearsize, linearsize, indx, d0) !Replaces DF with its LU decomposition
    CALL lubksb(DF_MATRIX, linearsize, linearsize, indx, F_VECTOR) !Solves DF*F=F using ludcmp

    F_VECTOR = MATMUL(D, F_VECTOR) !obtain solution from the regularized solution
  END SUBROUTINE

  SUBROUTINE Dispersivity()
    REAL :: Q = 10., depth=50., D
    REAL, DIMENSION(0:i_y, 1:j_x) :: alpha, dPsi, v
    INTEGER*4 :: i, j
    REAL, DIMENSION(1:j_x) :: x
    REAL, DIMENSION(0:i_y) :: y

    DO j = 1, j_x     !Initializing the values for x
      x(j) = dx*j
    END DO

    DO i = 0, i_y     !Initializing the values for y
      y(i) = dy*i
    END DO

    D = b*Q

    DO i = 0, i_y
      DO j = 1, j_x
        dPsi(i,j) = Psi(i,j) - Psi(i, j-1)
        v(i,j) = -(dPsi(i,j)/dx)*(Q/depth)
        alpha(i,j) = D/v(i,j)
        WRITE (40, 101) x(j), y(i), alpha(i, j)
      END DO
    END DO

    101     FORMAT(1x, 3(F7.3,3x))

  END SUBROUTINE

END PROGRAM Henrys_Problem
