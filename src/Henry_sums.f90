MODULE Henry_sums

  USE Henry_funcs
  IMPLICIT none

  CONTAINS

  FUNCTION Sum_AB(A_Matrix, B_Matrix, i_a, j_a, i_b, j_b, g, h)
    INTEGER, INTENT(IN) :: i_a, j_a, i_b, j_b, g, h
    INTEGER :: m, n, r, s
    REAL, DIMENSION (i_a, 0:j_a), INTENT(IN) :: A_Matrix
    REAL, DIMENSION (0:i_b, j_b), INTENT(IN) :: B_Matrix
    REAL :: Sum_AB

    Sum_AB = 0.                !Initializes Sum_AB to 0

    m_loop: DO m = 1, i_a
      n_loop: DO n = 0, j_a
        r_loop: DO r = 0, i_b
          s_loop: DO s = 1, j_b
            Sum_AB = Sum_AB + A_Matrix(m,n)*B_Matrix(r,s) &
                *(REAL(m*s*L_FUNC(m, r, g)*R_FUNC(h, n, s)) &
                  - REAL(n*r*F_FUNC(m, r, g)*G_FUNC(h, n, s)))
          END DO s_loop
        END DO r_loop
      END DO n_loop
    END DO m_loop
  END FUNCTION Sum_AB

  FUNCTION Sum_AN(A_Matrix, i, j, g, h)
    INTEGER, INTENT(IN) :: i, j, g, h
    INTEGER :: n
    REAL, DIMENSION (i, 0:j), INTENT(IN) :: A_Matrix
    REAL :: Sum_AN

    Sum_AN = 0.                !Initializes Sum_AN to 0

    IF (g /= 0) THEN
      DO n = 0, j
        Sum_AN = Sum_AN + A_Matrix(g, n)*g*N_FUNC(h, n)
      END DO
    END IF
  END FUNCTION Sum_AN

  FUNCTION Sum_BN(B_Matrix, i, j, g, h)
    INTEGER, INTENT(IN) :: i, j, g, h
    INTEGER :: r
    REAL, DIMENSION (0:i, j), INTENT(IN) :: B_Matrix
    REAL :: Sum_BN

    Sum_BN = 0.                !Initializes Sum_BN to 0

    IF (h /= 0) THEN
      DO r = 0, i
        Sum_BN = Sum_BN + B_Matrix(r, h)*h*N_FUNC(g, r)
      END DO
    END IF
  END FUNCTION Sum_BN

  FUNCTION Sum_BsN(B_Matrix, i, j, g, h)
    INTEGER*4, INTENT(IN) :: i, j, g, h
    INTEGER*4 :: s
    REAL, DIMENSION (0:i, j), INTENT(IN) :: B_Matrix
    REAL :: Sum_BsN

    Sum_BsN = 0. !Initializes Sum_BsN to 0

    IF (h /= 0) THEN
      DO s = 1, i
        Sum_BsN = Sum_BsN + B_Matrix(g, s)*s*N_FUNC(h, s)
      END DO
    END IF
  END FUNCTION Sum_BsN

  FUNCTION Sum_A(A_Matrix, x, y, xi, g, h, i, j, i_y, j_x)
    INTEGER, INTENT(IN) :: g, h, i, j, i_y, j_x
    REAL, INTENT(IN) :: xi
    INTEGER :: m, n
    REAL, DIMENSION (i, 0:j), INTENT(IN) :: A_Matrix
    REAL, DIMENSION (0:i_y), INTENT(IN) :: y
    REAL, DIMENSION (0:j_x), INTENT(IN) :: x
    REAL :: Sum_A

    Sum_A = 0.       !Initializes Sum_A to zero

    m_loop: DO m = 1, i
      n_loop: DO n = 0 , j
        Sum_A = Sum_A + A_Matrix(m, n)*SIN(m*pi*y(g))*COS(n*pi*x(h)/xi)
      END DO n_loop
    END DO m_loop
  END FUNCTION Sum_A

  FUNCTION Sum_B(B_Matrix, x, y, xi, g, h, i, j, i_y, j_x)
    INTEGER, INTENT(IN) :: g, h, i, j, i_y, j_x
    REAL, INTENT(IN) :: xi
    INTEGER :: r, s
    REAL, DIMENSION (0:i, j), INTENT(IN) :: B_Matrix
    REAL, DIMENSION (0:i_y), INTENT(IN) :: y
    REAL, DIMENSION (0:j_x), INTENT(IN) :: x
    REAL :: Sum_B

    Sum_B = 0.       !Initializes Sum_A to zero

    r_loop: DO r = 0, i
      s_loop: DO s = 1 , j
        Sum_B = Sum_B + B_Matrix(r, s)*COS(r*pi*y(g))*SIN(s*pi*x(h)/xi)
      END DO s_loop
    END DO r_loop
  END FUNCTION Sum_B

END MODULE Henry_sums
