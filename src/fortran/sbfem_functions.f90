module sbfem_functions
    !! This module contains the functions used in the sbfem code.
    !! The functions are the main functions used for
    !use iso_fortran_env
    use, intrinsic :: iso_c_binding, only : c_double, c_int, c_ptr, c_loc, c_bool, c_f_pointer
    use functional, only : arange
    use helper_functions, only : point_table
    use math, only : lagrange, lagrange_diff
    implicit none
    private

    public :: n_vec_f, n_mat_f, r_hat, r, j_mat, j_hat_mat, g_xi, g_eta, n_xi, &
    n_eta, n_vec_f_c

contains


    function n_vec_f(eta, poly_ord, deriv) result(n_vec)
        !! Computes the shape function vector or its derivative for a
        !!given polynomial order and evaluation point.

        ! Arguments
        real(c_double), intent(in) :: eta
        !! The natural coordinate (real number) at which the shape function
        !! or its derivative is evaluated.
        integer(c_int), intent(in) :: poly_ord
        !! The natural coordinate (real number) at which the shape function
        logical, intent(in) :: deriv
        !! A logical type which determines whether the shape function or its
        !! derivative is computed. If deriv is true, the derivative is computed.
        real(c_double) :: n_vec(poly_ord + 1)
        !! The shape function vector or its derivative.
        real(c_double) :: val
        integer(c_int) :: i
        real(c_double), allocatable :: xm(:)
        xm = point_table(poly_ord)
        do i = 1, poly_ord + 1
            if (deriv .eqv. .false.) then
                val = lagrange(eta, i, xm)
            else
                val = lagrange_diff(eta, i, xm)
            end if
            n_vec(i) = val
        end do

    end function n_vec_f

!    function n_vec_f_c(eta, poly_ord, deriv, nn) result(ptr) bind(C)
!        real(c_double), intent(in) :: eta
!        integer(c_int), intent(in) :: poly_ord
!        integer(c_bool), intent(in) :: deriv
!        integer(c_int), intent(out) :: nn
!        type(c_ptr) :: ptr
!
!        real(c_double), allocatable, target :: n_vec(:)
!        integer(c_int) :: i
!        real(c_double), allocatable :: xm(:)
!        print*, 'poly_ord: ', poly_ord
!        print*, 'eta:     ', eta
!        print*, 'deriv:   ', deriv
!
!
!
!
!        nn = poly_ord + 1
!        allocate(n_vec(nn))
!        if (.not. allocated(n_vec)) then
!            print *, "Allocation failed for n_vec"
!            return
!        end if
!        xm = point_table(poly_ord)
!        do i = 1, nn
!            if (deriv == 0_c_bool) then
!                ! Assuming lagrange and lagrange_diff are defined elsewhere
!                n_vec(i) = lagrange(eta, i, xm)
!            else
!                n_vec(i) = lagrange_diff(eta, i, xm)
!            end if
!        end do
!        ptr = c_loc(n_vec(1))
!    end function n_vec_f_c

    subroutine n_vec_f_c(eta, poly_ord, deriv, n_vec) bind(C, name="n_vec_f_c")
        use iso_c_binding, only: c_double, c_int
        real(c_double), intent(in) :: eta
        integer(c_int), intent(in) :: poly_ord
        integer(c_bool), intent(in) :: deriv
        real(c_double), intent(out) :: n_vec(*)

        integer :: i
        real(c_double), allocatable :: xm(:)

        xm = point_table(poly_ord)  ! Assuming point_table is defined elsewhere and works as expected.
        do i = 1, poly_ord + 1
            if (deriv == 0) then
                n_vec(i) = lagrange_diff(eta, i, xm)  ! Assuming lagrange_diff is defined elsewhere.
            else
                n_vec(i) = lagrange(eta, i, xm)  ! Assuming lagrange is defined elsewhere.
            end if
        end do
    end subroutine n_vec_f_c

!    function n_vec_f_c(eta, poly_ord, deriv, nn) result(ptr) bind(C)
!        real(c_double), intent(in) :: eta
!        integer(c_int), intent(in) :: poly_ord
!        integer(c_bool), intent(in) :: deriv
!        integer(c_int), intent(out) :: nn  ! Now nn will be an input and output variable
!        type(c_ptr) :: ptr
!
!        real(c_double), target :: n_vec(:)
!        ! NOTE
!        ! .... The target attribute in Fortran is used to indicate that a variable
!        ! or array may be the target of a pointer assignment.
!        ! This is necessary when you want to obtain a C pointer to a Fortran variable or
!        ! array using c_loc, as c_loc requires its argument to have a defined memory address
!        ! that does not change—something that is not guaranteed for allocatable arrays
!        ! or variables without the target attribute.
!
!        ! Allocate n_vec based on poly_ord, fill it, then return its size in 'n'
!        ! Example allocation and assignment
!
!        integer(c_int) :: i
!        real(c_double) :: val
!        real(c_double), allocatable :: xm(:)
!        print*, 'deriv:   ', deriv
!        print*, 'eta:     ', eta
!        print*, 'poly_ord: ', poly_ord
!
!
!        nn = poly_ord + 1;
!        allocate(n_vec(nn))
!        xm = point_table(poly_ord)
!        do i = 1, nn
!            if (deriv == 0) then
!                val = lagrange(eta, i, xm)
!            else
!                val = lagrange_diff(eta, i, xm)
!            end if
!            n_vec(i) = val
!        end do
!        ptr = c_loc(n_vec(1))
!    end function n_vec_f_c




    function n_mat_f(eta, poly_ord, deriv) result(n_mat)
        !! create the shape function Vector and Matrix N
        real(c_double), intent(in) :: eta
        integer(c_int), intent(in) :: poly_ord
        logical, intent(in) :: deriv
        real(c_double) :: n_mat(2, (poly_ord + 1) * 2)
        real(c_double) :: val
        integer(c_int) :: i, j ! loop variables
        real(c_double), allocatable :: xm(:)
        integer(c_int), allocatable :: range_list(:)
        xm = point_table(poly_ord)
        range_list = arange(0, poly_ord + 1)
        do i = 1, (poly_ord + 1)
            j = i + range_list(i)
            if (deriv .eqv. .false.) then
                val = lagrange(eta, i, xm)
            else
                val = lagrange_diff(eta, i, xm)
            end if
            n_mat(1, j) = val
            n_mat(1, j + 1) = 0._c_double
            n_mat(2, j) = 0._c_double
            n_mat(2, j + 1) = val
        end do
    end function n_mat_f

    function r_hat(xi, eta, poly_ord, coord_vec, centre) result(r_vec)
        !! create the shape function Vector and Matrix N
        real(c_double), intent(in) :: xi
        real(c_double), intent(in) :: eta
        integer(c_int), intent(in) :: poly_ord
        real(c_double), intent(in) :: coord_vec(:)
        real(c_double), optional, intent(in) :: centre(2)
        real(c_double) :: r_vec(2)
        real(c_double) :: cen(2)
        if(present(centre)) then
            cen = centre
        else
            cen = [0._c_double, 0._c_double]
        endif

        r_vec = xi * matmul(n_mat_f(eta, poly_ord, .false.), coord_vec) + cen

    end function r_hat

    function r(eta, poly_ord, coord_vec, centre) result(r_vec)
        !! create the shape function Vector and Matrix N
        real(c_double), intent(in) :: eta
        integer(c_int), intent(in) :: poly_ord
        real(c_double), intent(in) :: coord_vec(:)
        real(c_double), optional, intent(in) :: centre(2)
        real(c_double) :: r_vec(2)
        real(c_double) :: cen(2)
        if(present(centre)) then
            cen = centre
        else
            cen = [0._c_double, 0._c_double]
        endif

        r_vec = matmul(n_mat_f(eta, poly_ord, .false.), coord_vec) + cen

    end function r

    function j_mat(eta, coord_vec, poly_ord, shape_function_type, centre) result (j_matrix)
        ! Arguments
        real(c_double), intent(in) :: eta
        real(c_double), intent(in), dimension(2) :: coord_vec
        integer(c_int), intent(in) :: poly_ord
        character(len = :), allocatable, intent(in) :: shape_function_type
        real(c_double), dimension(2, 2) :: j_matrix
        real(c_double), intent(in), dimension(2) :: centre
        real(c_double), dimension(2) :: r_result
        real(c_double), dimension(2) :: shape_dN_result
        ! Call the r function
        r_result = r(eta, poly_ord, coord_vec, centre)

        ! Call the shape_dN function
        shape_dN_result = n_vec_f(eta, poly_ord, .true.)

        ! Populate the matrix
        j_matrix(1, 1) = r_result(1)
        j_matrix(1, 2) = r_result(2)
        j_matrix(2, 1) = dot_product(shape_dN_result, coord_vec)
        j_matrix(2, 2) = dot_product(shape_dN_result, coord_vec)

    end function j_mat

    ! J_Hat Matrix function
    function j_hat_mat(xi, eta, coord_vec, poly_ord, shape_function_type, centre) result (j_hat_matrix)
        ! Arguments
        real(c_double), intent(in) :: xi, eta
        real(c_double), intent(in), dimension(2) :: coord_vec
        integer(c_int), intent(in) :: poly_ord
        character(len = :), allocatable, intent(in) :: shape_function_type
        real(c_double), dimension(2, 2) :: j_hat_matrix
        real(c_double), intent(in), dimension(2) :: centre
        real(c_double), dimension(2) :: r_result
        real(c_double), dimension(2) :: shape_dN_result

        ! Call the r function
        r_result = r(eta, poly_ord, coord_vec, centre)

        ! Call the shape_dN function
        shape_dN_result = n_vec_f(eta, poly_ord, .true.)

        ! Populate the matrix
        j_hat_matrix(1, 1) = r_result(1)
        j_hat_matrix(1, 2) = r_result(2)
        j_hat_matrix(2, 1) = xi * dot_product(shape_dN_result, coord_vec)
        j_hat_matrix(2, 2) = xi * dot_product(shape_dN_result, coord_vec)

    end function j_hat_mat


    function g_xi(eta, coord_vec, poly_ord, centre) result(g_vec)
        real(c_double), intent(in) :: eta
        real(c_double), intent(in) :: coord_vec(:)
        integer(c_int), intent(in) :: poly_ord
        real(c_double), intent(in), optional :: centre(2)
        real(c_double) :: g_vec(2)
        real(c_double), allocatable :: dN(:, :)

        dN = n_mat_f(eta, poly_ord, .true.)
        g_vec = [-dN(2, :) * coord_vec(2), dN(2, :) * coord_vec(1)]
    end function g_xi

    function g_eta(eta, coord_vec, poly_ord, centre) result(g_vec)
        real(c_double), intent(in) :: eta
        real(c_double), intent(in) :: coord_vec(:)
        integer(c_int), intent(in) :: poly_ord
        real(c_double), intent(in), optional :: centre(2)
        real(c_double) :: g_vec(2)
        real(c_double) :: r_vec(2)

        r_vec = r(eta, poly_ord, coord_vec, centre)
        g_vec = [-r_vec(2), r_vec(1)]
    end function g_eta

    function n_xi(eta, coord_vec, poly_ord, centre) result(n_vec)
        real(c_double), intent(in) :: eta
        real(c_double), intent(in) :: coord_vec(:)
        integer(c_int), intent(in) :: poly_ord
        real(c_double), intent(in), optional :: centre(2)
        real(c_double) :: n_vec(2)
        real(c_double) :: g_vec(2)

        g_vec = g_xi(eta, coord_vec, poly_ord, centre)
        n_vec = g_vec / sqrt(dot_product(g_vec, g_vec))
    end function n_xi

    function n_eta(eta, coord_vec, poly_ord, centre) result(n_vec)
        real(c_double), intent(in) :: eta
        real(c_double), intent(in) :: coord_vec(:)
        integer(c_int), intent(in) :: poly_ord
        real(c_double), intent(in), optional :: centre(2)
        real(c_double) :: n_vec(2)
        real(c_double) :: g_vec(2)

        g_vec = g_eta(eta, coord_vec, poly_ord, centre)
        n_vec = g_vec / sqrt(dot_product(g_vec, g_vec))
    end function n_eta

    !     def r_hat_c(xi, eta, poly_ord, coord_vec, centre=np.array([0, 0])):
    !     return xi * np.dot(shape_N(eta, poly_ord)[1], coord_vec) + centre


    ! def r_hat(xi, eta, coord_vec, poly_ord):
    !     return xi * np.dot(shape_N(eta, poly_ord)[1], coord_vec)


    ! def r_c(eta, poly_ord, coord_vec, shape_function_type, centre=np.array([0, 0])):
    !     if shape_function_type == "standard shape functions":
    !         return (r_hat_c(1, eta, coord_vec, poly_ord, centre) - centre)
    !     elif shape_function_type == "hierarchical shape functions":
    !         return (r_hat_c(1, eta, coord_vec, poly_ord, centre) - centre)  # JET TO IMPLEMENT !!!!!
    !     else:
    !         print("ERROR: No valid shape function type in FUNCTION  rc(eta,coord_vec,centre)")
    !         return np.array([0, 0])


    ! def r(eta, poly_ord, coord_vec, shape_function_type, centre=np.array([0, 0])):
    !     if shape_function_type == "standard shape functions":
    !         return r_hat(1, eta, coord_vec, poly_ord)
    !     elif shape_function_type == "hierarchical shape functions":
    !         return r_hat(1, eta, coord_vec, poly_ord)  # JET TO IMPLEMENT !!!!!
    !     else:
    !         print("ERROR: No valid shape function type in FUNCTION -r(eta,coord_vec)")
    !         return np.array([0, 0])


end module sbfem_functions
