subroutine readHeader(input, nrow, ncol, nnzero) bind(c, name = 'readHarwellBoeingHeader') 
    
    use, intrinsic :: iso_c_binding
    implicit none
    
    ! input data (see https://software.intel.com/node/525336)
    character(len=1), dimension(100), intent(in) :: input ! an array of single characters (C-like string) to be converted into…
    integer length
    character(100), pointer :: fname ! …FORTRAN-like string
    
    ! return types
    integer(c_size_t) :: nrow, ncol, nnzero
    
    ! other info
    character title*72, key*8, mxtype*3, rhstyp*3, ptrfmt*16, indfmt*16, valfmt*20, rhsfmt*20
    integer   totcrd, ptrcrd, indcrd, valcrd, rhscrd, neltvl, nrhs, nrhsix
    
    ! get input file name
    ! (convert the array to a character variable)
    call c_f_pointer(c_loc(input), fname)
    length = index(fname, c_null_char) - 1
    open(unit = 56, file = fname(1:length));
    
    ! read data
    print *, 'reading HB header from "', fname(1:length), '" . . .'
    read(56, 100) title, key, totcrd, ptrcrd, indcrd, valcrd, rhscrd, mxtype, nrow, ncol, nnzero, neltvl, ptrfmt, indfmt, valfmt, rhsfmt
    if (rhscrd .gt. 0) read(56, 101) rhstyp, nrhs, nrhsix
        
    ! formats
100 format(a72, a8 / 5i14 / a3, 11x, 4i14 / 2a16, 2a20)
101 format(a3, 11x, 2i14)
    
end subroutine readHeader