!
! Žilyakov Alexander, Sep 2016
!
module LoadMod
    use CommonMod
    implicit none
contains
    !
    ! read header
    !
    subroutine loadHeader(fname_array, header) bind(c, name = 'loadHarwellBoeingHeader_f90') 
        ! input data
        character(c_char) :: fname_array(*)
        ! output data
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
        ! input data data
        character(c_char) :: fname_array(*)
        type(HBHeader)    :: header
        ! output data
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
        if (header.valcrd .gt. 0) then
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
    subroutine loadStructSym(fname_array, header, colptr, rowind, lval, diag, diagSize) bind(c, name = 'loadHarwellBoeingStructSym_f90') 
        ! input data
        character(c_char), intent(in) :: fname_array(*)
        type(HBHeader),    intent(in)    :: header
        ! output data
        integer(c_size_t), intent(out) :: colptr(*), rowind(*), diagSize
        real(c_double),    intent(out) :: diag(*), lval(*)
        ! helpers
        integer(c_size_t) :: i, j, rowIndex
		integer           :: statusIO
        ! convert C–string (input file name etc.) to FORTRAN ″
        call array2string(fname_array,   fname)
        call array2string(header.ptrfmt, ptrfmt)
        call array2string(header.indfmt, indfmt)
        call array2string(header.valfmt, valfmt)
        ! read data
        open(56, file = fname, action = 'read')
        ! read matrix structure
        ! read column pointers
        read(56, ptrfmt) (colptr(i), i = 1, header.ncol + 1)        
        ! read row indicies, fix column pointers
        diagSize = 0 ! numb of stored diag elements
        do i = 1, header.ncol ! loop over columns
            diag(i) = 0.
            read(56, indfmt, advance = 'no') rowIndex
            if (rowIndex .eq. i) then ! if the 1st element of ith column is in ith row (i.e. diagonal)
                diag(i) = 1.
                diagSize = diagSize + 1
                colptr(i + 1) = colptr(i + 1) - diagSize
                read(56, indfmt, advance = 'no') (rowind(j), j = colptr(i), colptr(i + 1) - 1)
            else
                colptr(i + 1) = colptr(i + 1) - diagSize
                rowind(colptr(i)) = rowIndex
                read(56, indfmt, advance = 'no') (rowind(j), j = colptr(i) + 1, colptr(i + 1) - 1)
            end if 
        end do
        if (header.valcrd .gt. 0 ) then
            ! read matrix values
            do i = 1, header.ncol
				read(56, *) ! to line
                if (diag(i) .eq. 1.) read(56, valfmt, advance = 'no') diag(i) ! diagonal element
                read(56, valfmt, advance = 'no') (lval(j), j = colptr(i), colptr(i + 1) - 1)
            end do
            ! … complex case
        endif
        close(56)
    end subroutine loadStructSym
end module LoadMod