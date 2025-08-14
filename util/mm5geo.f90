      subroutine nhgeosig(psc,sigma,ptop,ps0,ts0,tlp,i1,j1,k1,h)
      implicit none
      integer i1,j1,k1
      real psc(i1,j1),sigma(k1),h(i1,j1,k1)
      real ptop,p0s,ts0,tlp

      integer i,j,k
      real r,g,ps0,p0,term1,term2

      R=287.04
      G=9.81

      p0s=ps0/100.

      do k=1,k1
      do j=1,j1-1
      do i=1,i1-1
        p0=(sigma(k)*psc(i,j)+ptop)/100.
        term1=r*tlp/(2.*g)*(alog(p0/p0s))**2
        term2=r*ts0/g*alog(p0/p0s)
        h(i,j,k)=-1.*(term1+term2)
      enddo
      enddo
      enddo

      return
      end
