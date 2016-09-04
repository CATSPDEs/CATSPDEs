module LoadMod
    use CommonMod
    implicit none
contains
    !
    ! read header
    !
    subroutine loadHeader(fname_array, header) bind(c, name = 'loadHarwellBoeingHeader_f90') 
        ! i/o data
        character(c_char) :: fname_array(*)
        type(HBHeader)    :: header
        ! convert C–string (input file name) to FORTRAN ″
        call array2string(fname_array, fname)
        ! read data
        open(56, file = fname, action = 'read')
        read(56, headerFormat1to4) title, key, header.totcrd, header.ptrcrd, header.indcrd, header.valcrd, header.rhscrd, mxtype, header.nrow, header.ncol, header.nnzero, header.neltvl, ptrfmt, indfmt, valfmt, rhsfmt
        if (header.rhscrd .gt. 0) read(56, headerFormat5) rhstyp, header.nrhs, header.nrhsix
        ! convert FORTRAN strings to С ″
        call string2array(title,  header.title)
        call string2array(key,    header.key)
        call string2array(mxtype, header.mxtype)
        call string2array(rhstyp, header.rhstyp)
        call string2array(ptrfmt, header.ptrfmt)
        call string2array(indfmt, header.indfmt)
        call string2array(valfmt, header.valfmt)
        call string2array(rhsfmt, header.rhsfmt)
    end subroutine loadHeader
    !
    ! read structure 
    !
    subroutine loadStruct(fname_array, header, colptr, rowind, values) bind(c, name = 'loadHarwellBoeingStruct_f90') 
        ! i/o data
        character(kind = c_char) :: fname_array(*)
        type(HBHeader)    :: header
        integer(c_size_t) :: colptr(*), rowind(*), i
        real(c_double)    :: values(*)
        ! convert C–string (input file name etc.) to FORTRAN ″
        call array2string(fname_array,   fname)
        call array2string(header.ptrfmt, ptrfmt)
        call array2string(header.indfmt, indfmt)
        call array2string(header.valfmt, valfmt)
        ! read data
        open(56, file = fname, action = 'read')
        ! read matrix structure
        read(56, ptrfmt) (colptr(i), i = 1, header.ncol + 1)
        read(56, indfmt) (rowind(i), i = 1, header.nnzero)
        if (header.valcrd .gt. 0 ) then
            ! read matrix values
            if (header.mxtype(1) .eq. 'R') then ! real
                read(56, valfmt) (values(i), i = 1, header.nnzero)
            else if (header.mxtype(1) .eq. 'C') then ! complex
                read(56, valfmt) (values(i), i = 1, 2 * header.nnzero)
            endif
        endif
        close(56)
    end subroutine loadStruct
    !
    ! read structure (symmetric case) 
    !
    !subroutine loadStructSym(name_c, header, colptr, rowind, diag, lval) bind(c, name = 'loadHarwellBoeingStructSym_f90') 
    !    ! i/o data
    !    character(c_char) :: name_c(MAX_NAME_LENGTH)
    !    type(HBHeader)    :: header
    !    integer(c_size_t) :: colptr(*), rowind(*), i, j, k, m
    !    real(c_double)    :: diag(*), lval(*)
    !    ! convert C–string (input file name etc.) to FORTRAN ″
    !    name_f = transfer(name_c, name_f)
    !    ptrfmt = transfer(header.ptrfmt, ptrfmt)
    !    indfmt = transfer(header.indfmt, indfmt)
    !    valfmt = transfer(header.valfmt, valfmt)
    !    ! read data
    !    open(56, file = name_f(1 : index(name_f, c_null_char) - 1), action = 'read')
    !    ! read matrix structure
    !    ! read col pointers
    !    read(56, ptrfmt(1 : 16)) (colptr(i), i = 1, header.ncol + 1)
    !    ! read row indicies, fix col pointers
    !    k = 0 ! numb of stored diag elements
    !    do i = 1, header.ncol
    !        diag(i) = 0.
    !        read(56, indfmt(1 : 16)) m
    !        if (m .eq. i) then
    !            diag(i) = 1.
    !            k = k + 1
    !            colptr(i + 1) = colptr(i + 1) - k
    !            read(56, indfmt(1 : 16)) (rowind(j), j = colptr(i), colptr(i + 1) - 1)
    !        else
    !            colptr(i + 1) = colptr(i + 1) - k
    !            rowind(colptr(i)) = m
    !            read(56, indfmt(1 : 16)) (rowind(j), j = colptr(i) + 1, colptr(i + 1) - 1)
    !        end if 
    !    end do
    !    if (header.valcrd .gt. 0 ) then
    !        ! read matrix values
    !        do i = 1, header.ncol
    !            if (diag(i) .eq. 1.) read(56, valfmt(1 : 20)) diag(i) ! diagonal element
    !            read(56, valfmt(1 : 16)) (lval(j), j = colptr(i), colptr(i + 1) - 1)
    !        end do
    !        ! … complex case
    !    endif
    !    close(56)
    !end subroutine loadStructSym
end module LoadMod