      subroutine cross(n1,n2,n3)
        implicit none
        real*8::n1(3),n2(3),n3(3)
        n3(1)=n1(2)*n2(3)-n1(3)*n2(2)
        n3(2)=n1(3)*n2(1)-n1(1)*n2(3)
        n3(3)=n1(1)*n2(2)-n1(2)*n2(1)
        return
      end subroutine cross
