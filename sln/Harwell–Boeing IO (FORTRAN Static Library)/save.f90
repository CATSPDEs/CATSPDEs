module SaveMod
    use CommonMod
    implicit none
contains
    !
    ! write header
    !
    subroutine saveHeader(fname_array, header) bind(c, name = 'saveHarwellBoeingHeader_f90') 
        ! input data
        character(c_char) :: fname_array(*)
        type(HBHeader)    :: header
        ! convert C–string (input file name and header data) to FORTRAN ″
        call array2string(fname_array,   fname)
        call array2string(header.title,  title)
        call array2string(header.key,    key)
        call array2string(header.mxtype, mxtype)
        call array2string(header.rhstyp, rhstyp)
        call array2string(header.ptrfmt, ptrfmt)
        call array2string(header.indfmt, indfmt)
        call array2string(header.valfmt, valfmt)
        call array2string(header.rhsfmt, rhsfmt)    
        ! write data
        open(57, file = fname, action = 'write')
        write(57, headerFormat1to4) title, key, &
                                    header.totcrd, header.ptrcrd, header.indcrd, header.valcrd, header.rhscrd, &
                                    mxtype, header.nrow, header.ncol, header.nnzero, header.neltvl, &
                                    ptrfmt, indfmt, valfmt, rhsfmt
        if (header.rhscrd .gt. 0) write(57, headerFormat5) rhstyp, header.nrhs, header.nrhsix
    end subroutine saveHeader
    !
    ! write structure 
    !
    subroutine saveStruct(fname_array, header, colptr, rowind, values) bind(c, name = 'saveHarwellBoeingStruct_f90') 
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
        ! write data
        open(57, file = fname, action = 'write')
        ! write matrix structure
        write(57, ptrfmt) (colptr(i) + 1, i = 1, header.ncol + 1)
        write(57, indfmt) (rowind(i) + 1, i = 1, header.nnzero)
        if (header.valcrd .gt. 0 ) then
            ! write matrix values
            if (header.mxtype(1) .eq. 'R') then ! real
                write(57, valfmt) (values(i), i = 1, header.nnzero)
            else if (header.mxtype(1) .eq. 'C') then ! complex
                write(57, valfmt) (values(i), i = 1, 2 * header.nnzero) 
            endif
        endif
        close(57)
    end subroutine saveStruct
end module SaveMod