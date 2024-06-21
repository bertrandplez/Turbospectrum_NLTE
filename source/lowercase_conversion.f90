subroutine lowercase_conversion(input_string,output_string)

    implicit none
    character(len=100) :: input_string
    character(len=100) :: output_string
    integer :: i

    ! Initialize output string
    output_string = input_string

    ! Convert each character
    do i = 1, len_trim(input_string)
        if (input_string(i:i) >= 'A' .and. input_string(i:i) <= 'Z') then
            ! Convert uppercase to lowercase
            output_string(i:i) = achar(iachar(input_string(i:i)) + 32)
        else 
            output_string(i:i) = input_string(i:i)
        endif
    end do

end subroutine lowercase_conversion
