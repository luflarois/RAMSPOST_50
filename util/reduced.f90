!-----------------------------------------------------------------------

      subroutine RAMS_reduced_temp(n1,n2,n3,n4,tempnew,speed,ustar  &
              ,tstar,znew,zold,zrough,patfrac,cantemp,theta,topo,ztop)
      implicit none
      integer :: n1,n2,n3,n4,i,j,np
      real :: tempnew(n1,n2),speed(n1,n2,n3),ustar(n1,n2,n4),znew,zold  &
               ,zrough(n1,n2,n4),patfrac(n1,n2,n4),cantemp(n1,n2,n4)  &
               ,theta(n1,n2,n3),topo(n1,n2),ztop,tstar(n1,n2,n4)
      include 'rconstants.h'
      
      real:: richno,rtgt,zagl,rtemp,rtempw,z0,a2,spd
      
      

      do j=1,n2
         do i=1,n1
            
            rtgt=1.-topo(i,j)/ztop
            zagl=zold*rtgt
            
            rtempw=0.
            
            do np=1,n4
            
               z0=zrough(i,j,np)
               if(np==1) z0=.001
               spd=max(speed(i,j,2),.25)

               richno=g*zagl*(theta(i,j,2)-cantemp(i,j,np))  &
                           /(theta(i,j,2)*spd**2)
               a2 = (vonk / log(znew / z0)) ** 2
      
               if(richno.gt.0.) then
                  rtemp=cantemp(i,j,np)  &
                   +(ustar(i,j,np)*tstar(i,j,np))/(a2*spd)  &
!srf-                 *(1.+15.*richno/sqrt(1+5*richno))  
                      *(1.+15.*richno*sqrt(1+5*richno))  
                  rtemp=min(max(rtemp, cantemp(i,j,np)),theta(i,j,2))
               else
                  rtemp=cantemp(i,j,np)  &
                   +((ustar(i,j,np)*tstar(i,j,np))/(a2*spd))  &
                     / (1.- 15.*richno/(1.+75.*a2   &
                                   * sqrt(-znew*richno/z0)))
                  rtemp=max(min(rtemp, cantemp(i,j,np)),theta(i,j,2))
               endif
               
               !if((i==50.and.j==25)) then
               !   print*,'====tempf2m:',i,j
               !   print*,np,patfrac(i,j,np),cantemp(i,j,np)
               !   print*,np,ustar(i,j,np),zrough(i,j,np),tstar(i,j,np)
               !   print*,np,theta(i,j,2),speed(i,j,2),rtemp
               !endif
               
               rtempw=rtempw+rtemp*patfrac(i,j,np)
            
            enddo
            
            tempnew(i,j)=rtempw
            
  
         enddo
      enddo

      return
      end

!-----------------------------------------------------------------------

subroutine RAMS_reduced_wind(n1,n2,n3,n4,velnew,speed,ustar &
         ,znew,zold,zrough,patfrac,cantemp,theta,pi,topo,ztop)
implicit none
integer :: n1,n2,n3,n4,i,j,np
real :: velnew(n1,n2),speed(n1,n2,n3),ustar(n1,n2,n4),znew,zold  &
          ,zrough(n1,n2,n4),patfrac(n1,n2,n4),cantemp(n1,n2,n4)  &
          ,theta(n1,n2,n3),pi(n1,n2,n3),topo(n1,n2),ztop
include 'rconstants.h'

real:: richno,rtgt,zagl,rwind,rwindw,z0,a2,spd,cantheta,sfcpi



do j=1,n2
   do i=1,n1
      
      rtgt=1.-topo(i,j)/ztop
      zagl=zold*rtgt
      sfcpi=.5*(pi(i,j,1)+pi(i,j,2))
      
      rwindw=0.
      
      do np=1,n4
      
         z0=zrough(i,j,np)
         if(np==1) z0=.001
         spd=max(speed(i,j,2),.25)
         cantheta=cantemp(i,j,np)*cp/sfcpi

         richno=g*zagl*(theta(i,j,2)-cantheta)  &
                      /(theta(i,j,2)*spd**2)
         a2 = (vonk / log(znew / z0)) ** 2

         if(richno.gt.0.) then
            rwind=sqrt(ustar(i,j,np)**2/a2   &
                     *(1.+10.*richno/sqrt(1+5*richno)) )
         else
            rwind=sqrt( ustar(i,j,np)**2/a2  &
                / (1.- 10.*richno/(1.+75.*a2  &
                              * sqrt(-znew*richno/z0))))
         endif
         
         rwind=max(min(rwind,speed(i,j,2)),0.)
         
         !if(i==50.and.j==25) then
         !   print*,'====speed10m'
         !   print*,np,patfrac(i,j,np),cantemp(i,j,np)
         !   print*,np,ustar(i,j,np),zrough(i,j,np)
         !   print*,np,theta(i,j,2),speed(i,j,2),rwind
         !endif
         
         rwindw=rwindw+rwind*patfrac(i,j,np)
      
      enddo
      
      velnew(i,j)=rwindw
      

   enddo
enddo

return
end



----------------------- rconstants.h---------------------------------------

!---------------------------------------------------------------------------
real, parameter ::                    &
        rgas     = 287.               &
    ,   cp       = 1004.              &
    ,   cv       = 717.               &
    ,   rm       = 461.               &
    ,   p00      = 1.e5               &
    ,   t00      = 273.16             &
    ,   g        = 9.80               &
    ,   pi180    = 3.1415927 / 180.   &
    ,   pi4      = 3.1415927 * 4.     &
    ,   spcon    = 111120.            &
    ,   erad     = 6367000.           &
    ,   vonk     = 0.40               &
    ,   tkmin    = 5.e-4              &
    ,   alvl     = 2.50e6             &
    ,   alvi     = 2.834e6            &
    ,   alli     = 0.334e6            &
    ,   alvl2    = 6.25e12            &
    ,   alvi2    = 8.032e12           &
    ,   solar    = 1.3533e3           &
    ,   stefan   = 5.6696e-8          &
    ,   cww      = 4218.              &
    ,   c0       = 752.55 * 4.18684e4 &
    ,   viscos   = .15e-4             &
    ,   rowt     = 1.e3               &
    ,   dlat     = 111120.            &
    ,   omega    = 7.292e-5           &
    ,   rocp     = rgas / cp          &
    ,   p00i     = 1. / p00           &
    ,   cpor     = cp / rgas          &
    ,   rocv     = rgas / cv          &
    ,   cpi      = 1. / cp            &
    ,   cpi4     = 4. * cpi           &
    ,   cp253i   = cpi / 253.         & 
    ,   allii    = 1. / alli          &
    ,   aklv     = alvl / cp          &
    ,   akiv     = alvi / cp          &
    ,   gama     = cp / cv            &
    ,   gg       = .5 * g             &
    ,   ep       = rgas / rm          & 
    ,   p00k     = 26.870941          &  !  = p00 ** rocp  
    ,   p00ki    = 1. / p00k
!---------------------------------------------------------------------------





