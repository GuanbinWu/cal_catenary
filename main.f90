module tools
    implicit none
contains


function f1(x1,x2,y1,y2,k,c1) result(retval)
    real, intent(in) :: x1,x2,y1,y2,k,c1
    real :: retval
    retval=k*cosh((x2+c1)/k)-k*cosh((x1+c1)/k)-(y2-y1)
end function f1

function f2(x1,x2,l,k,c1) result(retval)
    real, intent(in) :: x1,x2,l,k,c1
    real :: retval
    retval=k*(sinh((x2+c1)/k)-sinh((x1+c1)/k))-l
end function f2

function pf1pk(x1,x2,k,c1) result(retval)
    real, intent(in) :: x1,x2,k,c1
    real :: retval
    retval=cosh((x2+c1)/k)+k*sinh((x2+c1)/k)*((-x2-c1)/(k**2)) - &
            (cosh((x1+c1)/k)+k*sinh((x1+c1)/k)*((-x1-c1)/(k**2)))
end function pf1pk

function pf2pk(x1,x2,k,c1) result(retval)
    real, intent(in) :: x1,x2,k,c1
    real :: retval
    retval= sinh((x2+c1)/k)+k*cosh((x2+c1)/k)*((-x2-c1)/(k**2)) - &
            (sinh((x1+c1)/k)+k*cosh((x1+c1)/k)*((-x1-c1)/(k**2)))
end function pf2pk

function pf1pc1(x1,x2,k,c1) result(retval)
    real, intent(in) :: x1,x2,k,c1
    real :: retval
    retval=sinh((x2+c1)/k)-sinh((x1+c1)/k)
end function pf1pc1

function pf2pc1(x1,x2,k,c1) result(retval)
    real, intent(in) :: x1,x2,k,c1
    real :: retval
    retval=cosh((x2+c1)/k)-cosh((x1+c1)/k)
end function pf2pc1

function distance(x1,x2,y1,y2,l,k,c1) result(retval)
    real, intent(in) :: x1,x2,y1,y2,l,k,c1
    real :: retval
    retval=sqrt(f1(x1,x2,y1,y2,k,c1)**2+f2(x1,x2,l,k,c1)**2)
end function distance

function cal_para(x1,x2,y1,y2,l,max_iter) result(para)
    real, intent(in) :: x1,x2,l,y1,y2
    integer, intent(in)::max_iter
    real ::percision,para(2),delta_c1,delta_k,k,c1,p1,p2,p3,p4,p5,p6
    integer ::counter,unit
    unit = 11

    open(unit=unit, file="./log.txt", status="replace")
    
    print *,"Itering..."
    c1=-(x2-x1)/2-x1
    k=((y2-y1)/(x2-x1))/2
    percision=0.0000001
    counter=0
    write(unit, "(a,f12.6,a,f12.6,a,f12.6)") &
    "Init guess: k = ",k,",c1 = ",c1,",error=",distance(x1,x2,y1,y2,l,k,c1)
    ! print *
    do while(distance(x1,x2,y1,y2,l,k,c1) > percision)
        ! print *,"Step: ",counter
        write(unit, "(a,i5)") "Iter step: ",counter
        p1=f1(x1,x2,y1,y2,k,c1)
        p2=f2(x1,x2,l,k,c1)
        p3=pf1pk(x1,x2,k,c1)
        p4=pf2pk(x1,x2,k,c1)
        p5=pf1pc1(x1,x2,k,c1)
        p6=pf2pc1(x1,x2,k,c1)

        ! print *,"f1=",p1
        ! print *,"f2=",p2
        ! print *,"pf1pk=",p3
        ! print *,"pf2pk=",p4
        ! print *,"pf1pc1=",p5
        ! print *,"pf2pc1=",p6

        delta_k=(p1*p6-p2*p5)/(p3*p6-p4*p5)
        delta_c1=(p1*p4-p2*p3)/(p5*p4-p6*p3)
        k=k-delta_k
        c1=c1-delta_c1
        counter = counter+1
        
        ! print *,
        ! print *,
        write(unit, "(a,f12.6,a,f12.6)")  "Delta k = ",delta_k,"Delta c1 = ",delta_c1
        write(unit, "(a,f12.6,a,f12.6,a,f12.6)")  "Update: k = ",k,",c1 = ",c1,",error=",distance(x1,x2,y1,y2,l,k,c1)
        if (counter>max_iter) then
            print *,"Convergence failed."
            exit
        endif
    end do
    print *,"Convergence sucessed."
    para(1)=k
    para(2)=c1
    close(unit)
end function cal_para

function cal_c2(x1,y1,k,c1) result(c2)
    real, intent(in) :: x1,y1,k,c1
    real :: c2
    c2=y1-k*cosh((x1+c1)/k)
end function cal_c2




subroutine sampling(x1, x2,k,c1,c2, steps, filename)
        real, intent(in) :: x1, x2,k,c1,c2
        integer, intent(in) :: steps
        character(len=*), intent(in) :: filename
        integer :: i, unit
        real :: dx, x
        dx = (x2 - x1) / (steps - 1)
        unit = 10
        open(unit=unit, file=filename, status="replace")
        do i = 0, steps - 1
            x = x1 + i * dx
            write(unit, '(F12.6, 1x, F12.6)') x, k*cosh((x+c1)/k)+c2
        end do
        close(unit)

    end subroutine sampling
end module tools




program main
use tools
implicit none


real::x1,y1,x2,y2,l,k_straight,b_straight,para(2),c2
integer :: max_iter
! print *,"Enter initial parameters as format : x1,y1,x2,y2,l,max_iter"
print *,"Enter initial parameters as format : x1,y1,x2,y2,l"
print *,"Notes: x1,y1,x2,y2 are coordinates of the two hanging points,"
print *,"       l is the length of rope,"
! print *,"       max_iter(optional) controls convergence conditions:"
read (*,*)x1,y1,x2,y2,l

max_iter=500

if ((l**2)<((y1-y2)**2+(x1-x2)**2)) then
    print *,"The rope is shorter than distance, impossible"
elseif ((l**2)==((y1-y2)**2+(x1-x2)**2)) then
    k_straight=(y2-y1)/(x2-x1)
    b_straight=y1-k_straight*x1
    print *,"The rope is just straight, and the function is:"
    print *,"y = ",k_straight,"* x + ",b_straight,", x in [ ",x1,",",x2," ] "
else
    para(:)=cal_para(x1,x2,y1,y2,l,max_iter)
    c2=cal_c2(x1,y1,para(1),para(2))
    print *,"The rope forms a catenary, and the function is:"
    print *,"y = ",para(1)," * cosh((x+",para(2),")/",para(1),")+",c2
    call sampling(x1,x2,para(1),para(2),c2,500,"./curve.txt")
endif



end program main