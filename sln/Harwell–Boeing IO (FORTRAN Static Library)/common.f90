module CommonMod
    use, intrinsic :: iso_c_binding, only : c_char, c_size_t, c_null_char, c_double
    implicit none
    integer(c_size_t) :: MAX_NAME_LENGTH
    parameter (MAX_NAME_LENGTH = 250) ! max string length for file path
    type, bind(c) :: HBHeader
        ! line 1
        character(c_char) :: title(73), key(9) ! +1 element for null–character of C–strings
        ! line 2
        integer(c_size_t) :: totcrd, ptrcrd, indcrd, valcrd, rhscrd
        ! line 3
        character(c_char) :: mxtype(4)
        integer(c_size_t) :: nrow, ncol, nnzero, neltvl;
        ! line 4
        character(c_char) :: ptrfmt(17), indfmt(17), valfmt(21), rhsfmt(21)
        ! line 5 (optional)
        character(c_char) :: rhstyp(4)
        integer(c_size_t) :: nrhs, nrhsix
    end type HBHeader
    ! FORTRAN–strings
    character(kind = c_char, len = MAX_NAME_LENGTH) :: fname
    character(c_char) :: title*72, key*8, mxtype*3, rhstyp*3, ptrfmt*16, indfmt*16, valfmt*20, rhsfmt*20
    ! formats
    character(kind = c_char, len = *), parameter :: headerFormat1to4 = '(a72, a8 / 5i14 / a3, 11x, 4i14 / 2a16, 2a20)', headerFormat5 = '(a3, 11x, 2i14)'
    contains 
    !
    ! convert FORTRAN string @s to array (C–string) @a
    !
    subroutine string2array(s, a)
        character(kind = c_char, len = *) :: s    ! input string
        character(c_char)                 :: a(*) ! output array
        integer  (c_size_t)               :: i    ! dummy index
        do i = 1, len(s)
            a(i) = s(i : i)
        end do
        a(i) = c_null_char
    end subroutine string2array
    !
    ! convert array (C–string) @a to FORTRAN string @s
    !
    subroutine array2string(a, s)
        character(c_char)                 :: a(*)  ! input array
        character(kind = c_char, len = *) :: s     ! output string
        integer(c_size_t)                 :: i     ! dummy index
        ! make string empty
        do i = 1, len(s)
            s(i : i) = ' '
        end do
        ! fill string
        i = 1
        do while (a(i) .ne. c_null_char)
            s(i : i) = a(i)
            i = i + 1
        end do
    end subroutine array2string
end module CommonMod
    