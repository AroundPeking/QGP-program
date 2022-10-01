module test
    integer :: a(6)!=(/(1,2,3,4,5,6)/)
    a(5)=1
    write(*,*) a(5)
end