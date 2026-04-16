module duchon_tps_module
  use, intrinsic :: iso_fortran_env, only : wp => real64
  implicit none
  private

  public :: run_duchon_tps_demo

  interface
     subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
       import :: wp
       integer, intent(in) :: n, nrhs, lda, ldb
       integer, intent(out) :: ipiv(*)
       integer, intent(out) :: info
       real(wp), intent(inout) :: a(lda,*), b(ldb,*)
     end subroutine dgesv
  end interface

contains

  subroutine run_duchon_tps_demo()
    implicit none

    ! ------------------------------------------------------------
    ! DATA DICTIONARY
    !
    ! n                 number of data sites
    ! mpoly             number of polynomial basis terms
    ! neq               total saddle-point system size
    !
    ! x(i), y(i)        coordinates of data site i in R^2
    ! f(i)              observed scalar datum at site i
    !
    ! K(i,j)            TPS kernel matrix entry phi(||xi-xj||)
    ! M(i,j)            regularized kernel matrix K + n*lambda*I
    ! T(i,k)            polynomial/design matrix at data site i
    !                   columns = [1, x, y]
    !
    ! A                 full block saddle-point matrix
    ! rhs               right-hand side [f; 0]
    ! sol               solution vector [c; d]
    !
    ! c(i)              radial basis coefficients
    ! d(k)              polynomial coefficients
    !
    ! lambda            smoothing parameter
    !
    ! fit(i)            fitted value at training site i
    ! resid(i)          residual f(i) - fit(i)
    !
    ! quad_penalty      c^T K c
    ! rss               residual sum of squares
    ! objective         0.5 * rss + lambda * c^T K c
    !
    ! constraint_norm2  squared Euclidean norm of T^T c
    ! ------------------------------------------------------------

    integer, parameter :: n = 9
    integer, parameter :: mpoly = 3
    integer, parameter :: neq = n + mpoly

    real(wp) :: x(n), y(n), f(n)
    real(wp) :: lambda
    real(wp) :: K(n,n), M(n,n), T(n,mpoly)
    real(wp) :: A(neq,neq), rhs(neq), sol(neq)
    real(wp) :: c(n), d(mpoly)
    real(wp) :: fit(n), resid(n), ttc(mpoly)
    real(wp) :: quad_penalty, rss, objective, constraint_norm2
    integer :: ipiv(neq), info
    integer :: i

    call init_sample_data(n, x, y, f)
    lambda = 1.0e-4_wp

    call build_tps_kernel_matrix(n, x, y, K)
    call build_polynomial_matrix(n, x, y, T)
    call build_regularized_matrix(n, lambda, K, M)
    call build_block_system(n, mpoly, M, T, f, A, rhs)

    sol = rhs
    call dgesv(neq, 1, A, neq, ipiv, sol, neq, info)
    if (info /= 0) then
       write(*,'(A,I0)') 'DGESV failed, INFO = ', info
       stop 1
    end if

    c = sol(1:n)
    d = sol(n+1:neq)

    call evaluate_fit_at_sites(n, x, y, x, y, c, d, fit)

    resid = f - fit
    ttc = matmul(transpose(T), c)
    constraint_norm2 = dot_product(ttc, ttc)
    quad_penalty = quadratic_form(n, c, K)
    rss = dot_product(resid, resid)
    objective = 0.5_wp * rss + lambda * quad_penalty

    call print_header(lambda, n, mpoly)
    call print_input_data(n, x, y, f)
    call print_solution(n, mpoly, c, d)
    call print_site_fits(n, x, y, f, fit, resid)
    call print_diagnostics(mpoly, ttc, constraint_norm2, quad_penalty, rss, objective)

  end subroutine run_duchon_tps_demo


  subroutine init_sample_data(n, x, y, f)
    implicit none
    integer, intent(in) :: n
    real(wp), intent(out) :: x(n), y(n), f(n)

    if (n /= 9) then
       write(*,*) 'This demo expects n = 9.'
       stop 1
    end if

    ! 3x3 grid in [0,1]x[0,1]
    x = [ 0.0_wp, 0.5_wp, 1.0_wp, &
          0.0_wp, 0.5_wp, 1.0_wp, &
          0.0_wp, 0.5_wp, 1.0_wp ]

    y = [ 0.0_wp, 0.0_wp, 0.0_wp, &
          0.5_wp, 0.5_wp, 0.5_wp, &
          1.0_wp, 1.0_wp, 1.0_wp ]

    ! Smooth target field with slight deterministic perturbation
    f(1) = field_function(x(1), y(1)) + 0.00_wp
    f(2) = field_function(x(2), y(2)) + 0.02_wp
    f(3) = field_function(x(3), y(3)) - 0.01_wp
    f(4) = field_function(x(4), y(4)) + 0.01_wp
    f(5) = field_function(x(5), y(5)) - 0.02_wp
    f(6) = field_function(x(6), y(6)) + 0.00_wp
    f(7) = field_function(x(7), y(7)) - 0.01_wp
    f(8) = field_function(x(8), y(8)) + 0.01_wp
    f(9) = field_function(x(9), y(9)) + 0.00_wp
  end subroutine init_sample_data


  pure real(wp) function field_function(x, y) result(val)
    implicit none
    real(wp), intent(in) :: x, y

    val = sin(3.141592653589793_wp * x) * cos(3.141592653589793_wp * y) + &
          0.30_wp * x - 0.20_wp * y
  end function field_function


  subroutine build_tps_kernel_matrix(n, x, y, K)
    implicit none
    integer, intent(in) :: n
    real(wp), intent(in) :: x(n), y(n)
    real(wp), intent(out) :: K(n,n)

    integer :: i, j
    real(wp) :: dx, dy, r

    do i = 1, n
       do j = 1, n
          dx = x(i) - x(j)
          dy = y(i) - y(j)
          r = sqrt(dx*dx + dy*dy)
          K(i,j) = thin_plate_spline_phi(r)
       end do
    end do
  end subroutine build_tps_kernel_matrix


  pure real(wp) function thin_plate_spline_phi(r) result(phi)
    implicit none
    real(wp), intent(in) :: r
    real(wp), parameter :: tiny = 1.0e-14_wp

    ! Thin-plate spline basis in 2D:
    ! phi(r) = r^2 log(r), with phi(0) = 0 by continuity
    if (r <= tiny) then
       phi = 0.0_wp
    else
       phi = (r*r) * log(r)
    end if
  end function thin_plate_spline_phi


  subroutine build_polynomial_matrix(n, x, y, T)
    implicit none
    integer, intent(in) :: n
    real(wp), intent(in) :: x(n), y(n)
    real(wp), intent(out) :: T(n,3)

    integer :: i

    do i = 1, n
       T(i,1) = 1.0_wp
       T(i,2) = x(i)
       T(i,3) = y(i)
    end do
  end subroutine build_polynomial_matrix


  subroutine build_regularized_matrix(n, lambda, K, M)
    implicit none
    integer, intent(in) :: n
    real(wp), intent(in) :: lambda
    real(wp), intent(in) :: K(n,n)
    real(wp), intent(out) :: M(n,n)

    integer :: i

    M = K
    do i = 1, n
       M(i,i) = M(i,i) + real(n, wp) * lambda
    end do
  end subroutine build_regularized_matrix


  subroutine build_block_system(n, mpoly, M, T, f, A, rhs)
    implicit none
    integer, intent(in) :: n, mpoly
    real(wp), intent(in) :: M(n,n), T(n,mpoly), f(n)
    real(wp), intent(out) :: A(n+mpoly, n+mpoly), rhs(n+mpoly)

    integer :: i, j

    A = 0.0_wp
    rhs = 0.0_wp

    ! Top-left: M
    A(1:n, 1:n) = M

    ! Top-right: T
    A(1:n, n+1:n+mpoly) = T

    ! Bottom-left: T^T
    A(n+1:n+mpoly, 1:n) = transpose(T)

    ! Bottom-right: 0
    rhs(1:n) = f

    do i = 1, n + mpoly
       do j = 1, n + mpoly
          if (.not. isfinite(A(i,j))) then
             write(*,*) 'Non-finite entry detected in system matrix.'
             stop 1
          end if
       end do
       if (.not. isfinite(rhs(i))) then
          write(*,*) 'Non-finite entry detected in RHS.'
          stop 1
       end if
    end do
  end subroutine build_block_system


  subroutine evaluate_fit_at_sites(nsrc, xs, ys, xq, yq, c, d, fq)
    implicit none
    integer, intent(in) :: nsrc
    real(wp), intent(in) :: xs(nsrc), ys(nsrc)
    real(wp), intent(in) :: xq(:), yq(:)
    real(wp), intent(in) :: c(nsrc), d(3)
    real(wp), intent(out) :: fq(size(xq))

    integer :: i, j, nq
    real(wp) :: dx, dy, r

    nq = size(xq)
    if (size(yq) /= nq .or. size(fq) /= nq) then
       write(*,*) 'Query array size mismatch.'
       stop 1
    end if

    do i = 1, nq
       fq(i) = d(1) + d(2) * xq(i) + d(3) * yq(i)
       do j = 1, nsrc
          dx = xq(i) - xs(j)
          dy = yq(i) - ys(j)
          r = sqrt(dx*dx + dy*dy)
          fq(i) = fq(i) + c(j) * thin_plate_spline_phi(r)
       end do
    end do
  end subroutine evaluate_fit_at_sites


  pure real(wp) function quadratic_form(n, c, K) result(val)
    implicit none
    integer, intent(in) :: n
    real(wp), intent(in) :: c(n), K(n,n)

    val = dot_product(c, matmul(K, c))
  end function quadratic_form


  elemental logical function isfinite(x)
    implicit none
    real(wp), intent(in) :: x
    isfinite = (x == x) .and. (abs(x) < huge(x))
  end function isfinite


  subroutine print_header(lambda, n, mpoly)
    implicit none
    real(wp), intent(in) :: lambda
    integer, intent(in) :: n, mpoly

    write(*,'(/,A)') '=== Research-grade Duchon / Thin-Plate Spline Demo ==='
    write(*,'(A,I0)') 'Number of sites n              = ', n
    write(*,'(A,I0)') 'Polynomial basis dimension     = ', mpoly
    write(*,'(A,ES14.6)') 'Smoothing parameter lambda     = ', lambda
    write(*,'(A)') 'Kernel                          = phi(r) = r^2 log(r), phi(0)=0'
    write(*,'(A)') 'Polynomial basis                = [1, x, y]'
  end subroutine print_header


  subroutine print_input_data(n, x, y, f)
    implicit none
    integer, intent(in) :: n
    real(wp), intent(in) :: x(n), y(n), f(n)
    integer :: i

    write(*,'(/,A)') 'Input data:'
    write(*,'(A)') '   i            x            y            f'
    do i = 1, n
       write(*,'(I4,3F13.6)') i, x(i), y(i), f(i)
    end do
  end subroutine print_input_data


  subroutine print_solution(n, mpoly, c, d)
    implicit none
    integer, intent(in) :: n, mpoly
    real(wp), intent(in) :: c(n), d(mpoly)
    integer :: i

    write(*,'(/,A)') 'Radial coefficients c:'
    do i = 1, n
       write(*,'(I4,ES22.10)') i, c(i)
    end do

    write(*,'(/,A)') 'Polynomial coefficients d:'
    do i = 1, mpoly
       write(*,'(I4,ES22.10)') i, d(i)
    end do
  end subroutine print_solution


  subroutine print_site_fits(n, x, y, f, fit, resid)
    implicit none
    integer, intent(in) :: n
    real(wp), intent(in) :: x(n), y(n), f(n), fit(n), resid(n)
    integer :: i

    write(*,'(/,A)') 'Fit at training sites:'
    write(*,'(A)') '   i            x            y        observed         fitted       residual'
    do i = 1, n
       write(*,'(I4,5F13.6)') i, x(i), y(i), f(i), fit(i), resid(i)
    end do
  end subroutine print_site_fits


  subroutine print_diagnostics(mpoly, ttc, constraint_norm2, quad_penalty, rss, objective)
    implicit none
    integer, intent(in) :: mpoly
    real(wp), intent(in) :: ttc(mpoly), constraint_norm2, quad_penalty, rss, objective
    integer :: i

    write(*,'(/,A)') 'Constraint check T^T c:'
    do i = 1, mpoly
       write(*,'(A,I0,A,ES18.8)') '  component ', i, ' = ', ttc(i)
    end do

    write(*,'(/,A,ES18.8)') '||T^T c||^2                 = ', constraint_norm2
    write(*,'(A,ES18.8)')   'c^T K c                    = ', quad_penalty
    write(*,'(A,ES18.8)')   'Residual sum of squares    = ', rss
    write(*,'(A,ES18.8)')   'Objective                  = ', objective
  end subroutine print_diagnostics

end module duchon_tps_module


program research_grade_duchon_tps
  use duchon_tps_module
  implicit none

  call run_duchon_tps_demo()
end program research_grade_duchon_tps
