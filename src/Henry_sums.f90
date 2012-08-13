MODULE Henry_sums

  USE Henry_funcs
  IMPLICIT none

  CONTAINS

  FUNCTION SUM_AB(MATRIX_A, MATRIX_B, i_a, j_a, i_b, j_b, g, h)

    INTEGER, INTENT(IN) :: i_a, j_a, i_b, j_b, g, h
    INTEGER :: m, n, r, s
    REAL, DIMENSION (i_a, 0:j_a), INTENT(IN) :: MATRIX_A
    REAL, DIMENSION (0:i_b, j_b), INTENT(IN) :: MATRIX_B
    REAL :: SUM_AB

    SUM_AB = 0.                !Initializes SUM_AB to 0

    m_loop: DO m = 1, i_a
      n_loop: DO n = 0, j_a
        r_loop: DO r = 0, i_b
          s_loop: DO s = 1, j_b
            SUM_AB = SUM_AB + MATRIX_A(m,n)*MATRIX_B(r,s) &
                *(REAL(m*s*L_FUNC(m, r, g)*R_FUNC(h, n, s)) &
                  - REAL(n*r*F_FUNC(m, r, g)*G_FUNC(h, n, s)))
          END DO s_loop
        END DO r_loop
      END DO n_loop
    END DO m_loop

  END FUNCTION SUM_AB

  FUNCTION SUM_AN(MATRIX_A, i, j, g, h)
            
    INTEGER, INTENT(IN) :: i, j, g, h
    INTEGER :: n
    REAL, DIMENSION (i, 0:j), INTENT(IN) :: MATRIX_A
    REAL :: SUM_AN

    SUM_AN = 0.                !Initializes SUM_AN to 0

    IF (g /= 0) THEN
      DO n = 0, j
        SUM_AN = SUM_AN + MATRIX_A(g, n)*g*N_FUNC(h, n)
      END DO
    END IF

  END FUNCTION SUM_AN

  FUNCTION SUM_BN(MATRIX_B, i, j, g, h)

    INTEGER, INTENT(IN) :: i, j, g, h
    INTEGER :: r
    REAL, DIMENSION (0:i, j), INTENT(IN) :: MATRIX_B
    REAL :: SUM_BN

    SUM_BN = 0.                !Initializes SUM_BN to 0

    IF (h /= 0) THEN
      DO r = 0, i
        SUM_BN = SUM_BN + MATRIX_B(r, h)*h*N_FUNC(g, r)
      END DO
    END IF

  END FUNCTION SUM_BN

  FUNCTION SUM_A(MATRIX_A, x, y, xi, g, h, i, j, i_y, j_x)
    
    INTEGER, INTENT(IN) :: g, h, i, j, i_y, j_x
    REAL, INTENT(IN) :: xi
    INTEGER :: m, n
    REAL, DIMENSION (i, 0:j), INTENT(IN) :: MATRIX_A
    REAL, DIMENSION (0:i_y), INTENT(IN) :: y
    REAL, DIMENSION (0:j_x), INTENT(IN) :: x
    REAL :: SUM_A

    SUM_A = 0.       !Initializes SUM_A to zero

    m_loop: DO m = 1, i
      n_loop: DO n = 0 , j
        SUM_A = SUM_A + MATRIX_A(m, n)*sin(m*pi*y(g))*cos(n*pi*x(h)/xi)
      END DO n_loop
    END DO m_loop

  END FUNCTION SUM_A

  FUNCTION SUM_B(MATRIX_B, x, y, xi, g, h, i, j, i_y, j_x)

    INTEGER, INTENT(IN) :: g, h, i, j, i_y, j_x
    REAL, INTENT(IN) :: xi
    INTEGER :: r, s
    REAL, DIMENSION (0:i, j), INTENT(IN) :: MATRIX_B
    REAL, DIMENSION (0:i_y), INTENT(IN) :: y
    REAL, DIMENSION (0:j_x), INTENT(IN) :: x
    REAL :: SUM_B

    SUM_B = 0.       !Initializes SUM_A to zero

    r_loop: DO r = 0, i
      s_loop: DO s = 1 , j
        SUM_B = SUM_B + MATRIX_B(r, s)*cos(r*pi*y(g))*sin(s*pi*x(h)/xi)
      END DO s_loop
    END DO r_loop

  END FUNCTION SUM_B

  FUNCTION SUM_BsN(MATRIX_B, i, j, g, h)
    INTEGER*4, INTENT(IN) :: i, j, g, h
    INTEGER*4 :: s
    REAL, DIMENSION (0:i, j), INTENT(IN) :: MATRIX_B
    REAL :: SUM_BsN

    SUM_BsN = 0. !Initializes SUM_BsN to 0

    IF (h /= 0) THEN
      DO s = 1, i
        SUM_BsN = SUM_BsN + MATRIX_B(g, s)*s*N_FUNC(h, s)
      END DO
    END IF
  END FUNCTION SUM_BsN

END MODULE Henry_sums
