program testpg

implicit none
integer i,j
real x,y

i=1
do j=1,10
x=real(i/j)
y=real(i)/real(j)
print*,i,j,x,y
enddo

end
