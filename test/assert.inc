pure subroutine assert(test, msg)
    logical, intent(in) :: test
    character(*), intent(in) :: msg
    if (.not. test) then
        error stop msg
    end if
end subroutine
