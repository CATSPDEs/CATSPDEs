    module Input
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
        character(c_char) :: name_f*MAX_NAME_LENGTH, title*72, key*8, mxtype*3, rhstyp*3, ptrfmt*16, indfmt*16, valfmt*20, rhsfmt*20
    contains
        !
        ! read header
        !
        subroutine readHeader(name_c, header) bind(c, name = 'readHarwellBoeingHeader') 
            ! i/o data
            character(c_char) :: name_c(MAX_NAME_LENGTH)
            type(HBHeader)    :: header
            ! convert C–string (input file name) to FORTRAN ″
            name_f = transfer(name_c, name_f)
            ! read data
            open(56, file = name_f(1 : index(name_f, c_null_char) - 1))
            read(56, 1000) title, key, header.totcrd, header.ptrcrd, header.indcrd, header.valcrd, header.rhscrd, mxtype, header.nrow, header.ncol, header.nnzero, header.neltvl, ptrfmt, indfmt, valfmt, rhsfmt
            if (header.rhscrd .gt. 0) read(56, 1001) rhstyp, header.nrhs, header.nrhsix
            ! convert FORTRAN strings to С ″
            call string2array(title, header.title, 72)
            call string2array(key, header.key, 8)
            call string2array(mxtype, header.mxtype, 3)
            call string2array(rhstyp, header.rhstyp, 3)
            call string2array(ptrfmt, header.ptrfmt, 16)
            call string2array(indfmt, header.indfmt, 16)
            call string2array(valfmt, header.valfmt, 20)
            call string2array(rhsfmt, header.rhsfmt, 20)
            ! formats
            1000 format(a72, a8 / 5i14 / a3, 11x, 4i14 / 2a16, 2a20)
            1001 format(a3, 11x, 2i14)
        end subroutine readHeader
        !
        ! read structure 
        !
        subroutine readStruct(name_c, header, colptr, rowind, values) bind(c, name = 'readHarwellBoeingStruct') 
            ! i/o data
            character(c_char) :: name_c(MAX_NAME_LENGTH)
            type(HBHeader)    :: header
            integer(c_size_t) :: colptr(*), rowind(*), i
            real(c_double)    :: values(*)
            ! convert C–string (input file name etc.) to FORTRAN ″
            name_f = transfer(name_c, name_f)
            ptrfmt = transfer(header.ptrfmt, ptrfmt)
            indfmt = transfer(header.indfmt, indfmt)
            valfmt = transfer(header.valfmt, valfmt)
            ! read data
            open(56, file = name_f(1 : index(name_f, c_null_char) - 1))
            ! read matrix structure
            read(56, ptrfmt(1 : 16)) (colptr(i), i = 1, header.ncol + 1)
            read(56, indfmt(1 : 16)) (rowind(i), i = 1, header.nnzero)
            if (header.valcrd .gt. 0 ) then
                ! read matrix values
                if (header.mxtype(3) .eq. 'A')  then
                    read(56, valfmt(1 : 20)) (values(i), i = 1, header.nnzero)
                else ! actually, we are interested in working w/ assembled matrices only, so…
                    read(56, valfmt(1 : 20)) (values(i), i = 1, header.neltvl) 
                endif
            ! …
            endif
            close(56)
        end subroutine readStruct
        ! helper
        subroutine string2array(s, a, l)
            integer(c_size_t) :: i, l
            character :: s*l
            character :: a(l + 1)
            do i = 1, l
                a(i) = s(i : i)
            end do
            a(l + 1) = c_null_char
        end subroutine string2array
    end module Input