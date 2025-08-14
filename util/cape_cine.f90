!******************************** FIM *********************************
      subroutine cape_cine(nx,ny,nz,press,TEMPK,UR,dummy)
!      subroutine cape_cine(nz,nt,nomeIN,nomeOUT,&
!                     nSzP,nSzpP,nSzT,nSzpT,nSzUR,nSzpUR,lit,indef)
!      include 'grade_cape.h'
      real PRESS(nx,ny,nz), TEMPK(nx,ny,nz),UR(nx,ny,nz),indef
      real dummy(nx,ny,nz)
      real cine(nx,ny),cape(nx,ny)
      real press2(nz),temp2(nz),ur2(nz)
      character nomeIN*100,nomeOUT*100
      integer t,z,x,y
      erro0=5.e-5

!      OPEN(34,FILE=nomeIN,STATUS='UNKNOWN'&
!            ,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=nx*ny*4)
!      OPEN(35,FILE=nomeOUT,STATUS='UNKNOWN'&
!            ,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=nx*ny*4)
      
!      irec=1
!      irecP=-nSzpP
!      irecT=-nSzpT
!      irecUR=-nSzpUR
!      do t=1,nt
!        print*,t,' de ',nt
!        irecP=irecP + nSzP + nSzpP
!        irecT=irecT + nSzT + nSzpT
!        irecUR=irecUR + nSzUR + nSzpUR
!        do z=1,nz
!          irecP=irecP+1
!          irecT=irecT+1
!          irecUR=irecUR+1
!          read(34,rec=irecP) ((PRESS(i,j,z),i=1,nx),j=1,ny)
!          read(34,rec=irecT) ((TEMPK(i,j,z),i=1,nx),j=1,ny)
!          read(34,rec=irecUR) ((UR(i,j,z),i=1,nx),j=1,ny)
!          
!        enddo

        do i=1,nx
          do j=1,ny
            do z=1,nz
              press2(z)=PRESS(i,j,z)
              temp2(z)=TEMPK(i,j,z)-273.16
              ur2(z)=UR(i,j,z)
            enddo
            cape(i,j)=calccape(nz,1,press2,temp2,ur2,erro0,indef)
            cine(i,j)=calccine(nz,1,press2,temp2,ur2,erro0,indef)
	    print*,'cape-cine',cape(i,j),cine(i,j)
          enddo
        enddo
!        write(35,rec=irec) cape
!        irec=irec+1
!        write(35,rec=irec) cine
!        irec=irec+1
!      enddo

!      close (35)
!      close (34)

      return
      end
      
      
!************************************************************************
!*  Esta função calcula a CINE de uma determinada sondagem              *
!************************************************************************
      real function calccine(num,i0,pres,temp,urel,erro0,indef)
       implicit none !Tens de declarar tudo...
       real pres(*),temp(*),urel(*)
       real indef,pncl0,rmis0,tpot0,tpeq0,tamb,ramb,tvamb         !Ambiente
       real pres1,pres2,inte1,inte2                               !Comuns
       real tpar1,tpar2,rpar1,rpar2,tvpar                         !Parcela
       real presdoncl,tempvirtual,vartvarp,potencial,potencialeq
       real razaodemistura,varrvarp,integrando,epsi,cine,erro0
       real tparcela,rparcela
       integer num,i0,i
       logical fim
       parameter (epsi=0.62198)
       cine=0.
       fim=.false.
       if (pres(i0).eq.indef.or.temp(i0).eq.indef.or.urel(i0).eq.indef) then
           calccine=indef
           return
       endif
!* Tomo os primeiros valores, para depois jogá-los aos valores velhos...
       pncl0=presdoncl(pres(i0),temp(i0),urel(i0),indef)
       rmis0=razaodemistura(pres(i0),temp(i0),urel(i0),indef)
       tpot0=potencial(pres(i0),temp(i0)+273.16,rmis0,indef)
       tpeq0=potencialeq(pres(i0),temp(i0)+273.16,rmis0,indef)
       pres2=pres(i0)
       tamb=temp(i0)+273.16
       ramb=razaodemistura(pres2,tamb-273.16,urel(i0),indef)
       tvamb=tempvirtual(pres2,tamb,ramb,indef)
       tpar2=temp(i0)+273.16 !Começa com mesma temperatura do ambiente
       rpar2=rmis0           !Começa com mesmo rmis do ambiente
       tvpar=tempvirtual(pres2,tpar2,rpar2,indef)
       inte2=0.
       i=i0+1
       do while (.not.fim.and.i.le.num)
         if (pres(i).ne.indef.and.temp(i).ne.indef.and.urel(i).ne.indef)then
! Passo os valores de algumas variáveis para o valor "velho"
             pres1=pres2
             tpar1=tpar2
             rpar1=rpar2
             inte1=inte2
! Recalculo estas variáveis e calculo a contribuição para o CINE
             pres2=pres(i)
             tamb=temp(i)+273.16
             ramb=razaodemistura(pres2,tamb-273.16,urel(i),indef)
             tvamb=tempvirtual(pres2,tamb,ramb,indef)
             tpar2=tparcela(pres2,pncl0,tpot0,tpeq0,rmis0,erro0,indef)
             rpar2=rparcela(pres2,pncl0,rmis0,tpar2,indef)
             tvpar=tempvirtual(pres2,tpar2,rpar2,indef)
             inte2=integrando(tvamb,tvpar,indef)
             if (inte2.lt.0.) then 
                 fim=.true.
               else             
                 cine=cine-0.5*(inte1+inte2)*log(pres2/pres1)
             endif             
         endif 
         i=i+1
       enddo
!   Caso tenha acabado até aqui, indefini-lo-ei, pois na realidade ele
! vale infinito.....
       if (.not.fim) cine=indef 
       calccine=cine
       return
      end





!***********************************************************************
!  Esta função calcula o NCE de uma determinada sondagem               *
!***********************************************************************
       real function calcnce(num,i0,pres,temp,urel,erro0,indef)
       implicit none !Tens de declarar tudo...
       real pres(*),temp(*),urel(*)
       real indef,pncl0,rmis0,tpot0,tpeq0,tamb,ramb,tvamb         !Ambiente
       real pres1,pres2                                           !Comuns
       real tpar1,tpar2,rpar1,rpar2,tvpar                         !Parcela
       real presdoncl,tempvirtual,vartvarp,potencial,potencialeq
       real razaodemistura,varrvarp,integrando,epsi,nce,erro0
       real tparcela,rparcela
       integer num,i0,i
       logical fim
       parameter (epsi=0.62198)
       nce=0.
       fim=.false.
       if (pres(i0).eq.indef.or.temp(i0).eq.indef.or.urel(i0).eq.indef)then
           nce=indef
           return
       endif
! Tomo os primeiros valores, para depois jogá-los aos valores velhos...
       pncl0=presdoncl(pres(i0),temp(i0),urel(i0),indef)
       rmis0=razaodemistura(pres(i0),temp(i0),urel(i0),indef)
       tpot0=potencial(pres(i0),temp(i0)+273.16,rmis0,indef)
       tpeq0=potencialeq(pres(i0),temp(i0)+273.16,rmis0,indef)
       pres2=pres(i0)
       tamb=temp(i0)+273.16
       ramb=razaodemistura(pres2,tamb-273.16,urel(i0),indef)
       tvamb=tempvirtual(pres2,tamb,ramb,indef)
       tpar2=temp(i0)+273.16 !Começa com mesma temperatura do ambiente
       rpar2=rmis0           !Começa com mesmo rmis do ambiente
       tvpar=tempvirtual(pres2,tpar2,rpar2,indef)
       i=i0+1
       do while (.not.fim.and.i.le.num)
         if (pres(i).ne.indef.and.temp(i).ne.indef.and.urel(i).ne.indef)then
! Passo os valores de algumas variáveis para o valor "velho"
             pres1=pres2
             tpar1=tpar2
             rpar1=rpar2
! Recalculo estas variáveis e calculo a contribuição para o CINE
             pres2=pres(i)
             tamb=temp(i)+273.16
             ramb=razaodemistura(pres2,tamb-273.16,urel(i),indef)
             tvamb=tempvirtual(pres2,tamb,ramb,indef)
             tpar2=tparcela(pres2,pncl0,tpot0,tpeq0,rmis0,erro0,indef)
             rpar2=rparcela(pres2,pncl0,rmis0,tpar2,indef)
             tvpar=tempvirtual(pres2,tpar2,rpar2,indef)
             if (tvpar.gt.tvamb) then 
                 fim=.true.
                 nce=0.5*(pres1+pres2)
             endif             
         endif 
         i=i+1
       enddo
!*   Caso tenha acabado até aqui, zero-o, numa forma de dizer que é inatingível
       if (.not.fim) nce=0.
       calcnce=nce
       return
      end





!************************************************************************
!*  Esta função calcula a CAPE de uma determinada sondagem              *
!************************************************************************
      real function calccape(num,i0,pres,temp,urel,erro0,indef)
       implicit none !Tens de declarar tudo...
       real pres(*),temp(*),urel(*)
       real indef,pncl0,rmis0,tpot0,tpeq0,tamb,ramb,tvamb         !Ambiente
       real pres1,pres2,inte1,inte2                               !Comuns
       real tpar1,tpar2,rpar1,rpar2,tvpar                         !Parcela
       real presdoncl,tempvirtual,vartvarp,potencial,potencialeq
       real razaodemistura,varrvarp,integrando,epsi,cape,erro0
       real tparcela,rparcela
       integer num,i0,i
       logical fim,embaixo
       parameter (epsi=0.62198)
       cape=0.
       fim=.false.
       embaixo=.true.
       if (pres(i0).eq.indef.or.temp(i0).eq.indef.or.urel(i0).eq.indef)then
           calccape=indef
           return
       endif
!* Tomo os primeiros valores, para depois jogá-los aos valores velhos...
       pncl0=presdoncl(pres(i0),temp(i0),urel(i0),indef)
       rmis0=razaodemistura(pres(i0),temp(i0),urel(i0),indef)
       tpot0=potencial(pres(i0),temp(i0)+273.16,rmis0,indef)
       tpeq0=potencialeq(pres(i0),temp(i0)+273.16,rmis0,indef)
       pres2=pres(i0)
       tamb=temp(i0)+273.16
       ramb=razaodemistura(pres2,tamb-273.16,urel(i0),indef)
       tvamb=tempvirtual(pres2,tamb,ramb,indef)
       tpar2=temp(i0)+273.16 !Começa com mesma temperatura do ambiente
       rpar2=rmis0           !Começa com mesmo rmis do ambiente
       tvpar=tempvirtual(pres2,tpar2,rpar2,indef)
       inte2=0.
       i=i0+1
       do while (.not.fim.and.i.le.num)
         if (pres(i).ne.indef.and.temp(i).ne.indef.and.urel(i).ne.indef)then
!* Passo os valores de algumas variáveis para o valor "velho"
             pres1=pres2
             tpar1=tpar2
             rpar1=rpar2
             inte1=inte2
!* Recalculo estas variáveis e calculo a contribuição para o CINE
             pres2=pres(i)
             tamb=temp(i)+273.16
             ramb=razaodemistura(pres2,tamb-273.16,urel(i),indef)
             tvamb=tempvirtual(pres2,tamb,ramb,indef)
             tpar2=tparcela(pres2,pncl0,tpot0,tpeq0,rmis0,erro0,indef)
             rpar2=rparcela(pres2,pncl0,rmis0,tpar2,indef)
             tvpar=tempvirtual(pres2,tpar2,rpar2,indef)
             inte2=integrando(tvamb,tvpar,indef)
             if (.not.embaixo.and.inte2.gt.0.) then 
                 fim=.true.
               elseif (inte2.le.0) then
                 embaixo=.false.         
                 cape=cape+0.5*(inte1+inte2)*log(pres2/pres1)
             endif             
         endif 
         i=i+1
       enddo
!*   Caso tenha acabado até aqui, indefini-lo-ei, pois na realidade ele
!* vale infinito.....
       if (.not.fim) cape=0.
       calccape=cape
       return
      end





!************************************************************************
!*  Esta função calcula o NPE de uma determinada sondagem		*
!************************************************************************
      real function calcnpe(num,i0,pres,temp,urel,erro0,indef)
       implicit none !Tens de declarar tudo...
       real pres(*),temp(*),urel(*)
       real indef,pncl0,rmis0,tpot0,tpeq0,tamb,ramb,tvamb         !Ambiente
       real pres1,pres2,inte1,inte2                               !Comuns
       real tpar1,tpar2,rpar1,rpar2,tvpar                         !Parcela
       real presdoncl,tempvirtual,vartvarp,potencial,potencialeq
       real razaodemistura,varrvarp,integrando,epsi,npe,erro0
       real tparcela,rparcela
       integer num,i0,i
       logical fim,embaixo
       parameter (epsi=0.62198)
       npe=0.
       fim=.false.
       embaixo=.true.
       if (pres(i0).eq.indef.or.temp(i0).eq.indef.or.urel(i0).eq.indef)then
           calcnpe=indef
           return
       endif
!* Tomo os primeiros valores, para depois jogá-los aos valores velhos...
       pncl0=presdoncl(pres(i0),temp(i0),urel(i0),indef)
       rmis0=razaodemistura(pres(i0),temp(i0),urel(i0),indef)
       tpot0=potencial(pres(i0),temp(i0)+273.16,rmis0,indef)
       tpeq0=potencialeq(pres(i0),temp(i0)+273.16,rmis0,indef)
       pres2=pres(i0)
       tamb=temp(i0)+273.16
       ramb=razaodemistura(pres2,tamb-273.16,urel(i0),indef)
       tvamb=tempvirtual(pres2,tamb,ramb,indef)
       tpar2=temp(i0)+273.16 !Começa com mesma temperatura do ambiente
       rpar2=rmis0           !Começa com mesmo rmis do ambiente
       tvpar=tempvirtual(pres2,tpar2,rpar2,indef)
       inte2=0.
       i=i0+1
       do while (.not.fim.and.i.le.num)
         if (pres(i).ne.indef.and.temp(i).ne.indef.and.urel(i).ne.indef)then
!* Passo os valores de algumas variáveis para o valor "velho"
             pres1=pres2
             tpar1=tpar2
             rpar1=rpar2
!* Recalculo estas variáveis e calculo a contribuição para o CINE
             pres2=pres(i)
             tamb=temp(i)+273.16
             ramb=razaodemistura(pres2,tamb-273.16,urel(i),indef)
             tvamb=tempvirtual(pres2,tamb,ramb,indef)
             tpar2=tparcela(pres2,pncl0,tpot0,tpeq0,rmis0,erro0,indef)
             rpar2=rparcela(pres2,pncl0,rmis0,tpar2,indef)
             tvpar=tempvirtual(pres2,tpar2,rpar2,indef)
             if (.not.embaixo.and.tvpar.lt.tvamb) then 
                 fim=.true.
                 npe=0.5*(pres1+pres2)
               elseif (tvpar.ge.tvamb) then
                 embaixo=.false.         
             endif             
         endif 
         i=i+1
       enddo
!*   Caso tenha acabado até aqui, indefini-lo-ei, pois na realidade ele
!* vale infinito.....
       if (.not.fim) npe=indef
       calcnpe=npe
       return
      end







!************************************************************************
!* Função que calcula a pressão do NCL a partir de p,T,Urel		*
!************************************************************************  
      real function presdoncl(pres0,temp0,urel0,indef)
       implicit none !Tens de declarar tudo....
       real pres0,temp0,urel0,indef,tempk,tpot,tncl
       if (pres0.eq.indef.or.temp0.eq.indef.or.urel0.eq.indef) then
           presdoncl=indef
         else
           tempk=temp0+273.16
           tpot=tempk*((1000./pres0)**0.286)
           tncl=1/(1/(tempk-55.)-log(urel0/100.)/2840.)+55.
           presdoncl=1000.*((tncl/tpot)**3.4965035)
       endif
       return
      end
      

!
!
!************************************************************************
!* Função que calcula a temperatura virtual do ar                       *
!************************************************************************
      real function tempvirtual(pres0,temp0,rmis0,indef)
       implicit none !Tens de declarar tudo....
       real pres0,temp0,rmis0,umes,epsi,pvap,indef
       parameter (epsi=0.62198)
       if (pres0.eq.indef.or.temp0.eq.indef.or.rmis0.eq.indef) then
           tempvirtual=indef
         else
          umes=rmis0/(rmis0+1)
          tempvirtual=temp0*(1+0.61*umes)
       endif
       return
      end

!
!
!************************************************************************
!* Função que calcula a razão de mistura em kg/kg                       *
!************************************************************************
      real function razaodemistura(pres0,temp0,urel0,indef)
       implicit none !Tens de declarar tudo....
       real pres0,temp0,urel0,pvap,indef,epsi
       parameter (epsi=0.62198)
       if (pres0.eq.indef.or.temp0.eq.indef.or.urel0.eq.indef) then
           razaodemistura=indef
         else
           pvap=0.01*urel0*6.112*exp(17.67*temp0/(temp0+243.5))
           razaodemistura=epsi*pvap/(pres0-pvap)
       endif
       return
      end

!
!
!************************************************************************
!* Função que calcula a temperatura potencial da parcela                *
!************************************************************************
      real function potencial(pres0,temp0,rmis0,indef)
       implicit none !Tens de declarar tudo
       real pres0,temp0,rmis0,indef,epsi
       parameter (epsi=0.62198)
       if (pres0.eq.indef.or.temp0.eq.indef) then
           potencial=indef
         elseif (rmis0.eq.indef) then
           potencial=temp0*((1000./pres0)**0.2854)
         else
           potencial=temp0*((1000./pres0)**(0.2854*(1-0.28*rmis0)))
       endif
       return
      end



!************************************************************************
!* Função que calcula a temperatura potencial equivalente da parcela    *
!************************************************************************
      real function potencialeq(pres0,temp0,rmis0,indef)
       implicit none !Tens de declarar tudo...
       real pres0,temp0,rmis0,indef,pvap,tncl,epsi,tpot
       parameter (epsi=0.62198)
       if (pres0.eq.indef.or.temp0.eq.indef.or.rmis0.eq.indef) then
           potencialeq=indef
         else
           pvap=pres0*rmis0/(epsi+rmis0)
           tncl=2840./(3.5*log(temp0)-log(pvap)-4.805)+55.
           tpot=temp0*((1000./pres0)**(0.2854*(1-0.28*rmis0)))
           potencialeq=tpot*exp((3.376/tncl-0.00254)*&
                       1000.*rmis0*(1+0.81*rmis0))
       endif
       return
      end

!************************************************************************
!* Função iterativa que calcula a temperatura da parcela                *
!************************************************************************
      real function tparcela(pres0,pncl0,tpot0,tpeq0,rmis0,erro0,indef)
       implicit none !Tens de declarar tudo...
       real pres0,pncl0,tpot0,tpeq0,rmis0,erro0,indef
       real erro,epsi,tparnovo,tparm0,tparm1,esatm0,esatm1,rsatm0
       real rsatm1,tpotm0,tpotm1,tpeqm0,tpeqm1
       parameter (epsi=0.62198)
       tparnovo=273.16
       erro=2.*erro0
       if (pres0.eq.indef.or.pncl0.eq.indef) then
           tparnovo=indef
         elseif(pncl0.eq.0.) then !Pressão do NCL inatingível
           tparnovo=indef
         elseif(pres0.gt.pncl0) then !Iterage com temperatura potencial
           do while (erro.gt.erro0)
             tparm0=tparnovo
             tparm1=tparm0+1.
             rsatm0=1000.*rmis0 !Só por via das dúvidas
             tpotm0=tparm0*((1000./pres0)**(0.2854*(1-0.28e-3*rsatm0)))
             tpotm1=tparm1*((1000./pres0)**(0.2854*(1-0.28e-3*rsatm0)))
             tparnovo=tparm0+(tpot0-tpotm0)/(tpotm1-tpotm0)
             erro=200.*(tparnovo-tparm0)/(tparnovo+tparm0) !O erro está em %
           enddo           
         else !Iterage com temperatura potencial equivalente
           do while (erro.gt.erro0)
             tparm0=tparnovo
             tparm1=tparm0+1.
             esatm0=6.112*exp(17.67*(tparm0-273.16)/(tparm0-29.66))
             esatm1=6.112*exp(17.67*(tparm1-273.16)/(tparm1-29.66))
             rsatm0=1000.*epsi*esatm0/(pres0-esatm0)
             rsatm1=1000.*epsi*esatm1/(pres0-esatm1)
             tpotm0=tparm0*((1000./pres0)**(0.2854*(1-0.28e-3*rsatm0)))
             tpotm1=tparm1*((1000./pres0)**(0.2854*(1-0.28e-3*rsatm1)))
             tpeqm0=tpotm0*exp((3.376/tparm0-0.00254)*rsatm0*(1+0.81e-3*rsatm0))
             tpeqm1=tpotm1*exp((3.376/tparm1-0.00254)*rsatm1*(1+0.81e-3*rsatm1))
             tparnovo=tparm0+(tpeq0-tpeqm0)/(tpeqm1-tpeqm0)
             erro=abs(200.*(tparnovo-tparm0)/(tparnovo+tparm0)) !O erro está em %
           enddo
       endif
       tparcela=tparnovo
       return
      end

!************************************************************************
!* Função que calcula a razão de mistura da parcela                     *
!************************************************************************
      real function rparcela(pres0,pncl0,rmis0,tpar,indef)
       implicit none !Tens de declarar tudo...
       real pres0,pncl0,rmis0,tpar,indef,razaodemistura
       if (pres0.eq.indef.or.pncl0.eq.0) then
           rparcela=indef
         elseif (pres0.gt.pncl0) then
           rparcela=rmis0
         else
           rparcela=razaodemistura(pres0,tpar-273.16,0.*tpar+100.,indef)
       endif
       return
      end


!************************************************************************
!* Função que calcula o termo integrando da CINE e CAPE                 *
!************************************************************************
      real function integrando(tvamb,tvpar,indef)
       implicit none !Tens de declarar tudo....
       real tvamb,tvpar,ra,indef
       parameter (ra=287.04)
       if (tvamb.eq.indef.or.tvpar.eq.indef) then
           integrando=indef
         else
           integrando=-ra*(tvpar-tvamb)
       endif
       return
      end

!***********************************************************************


        
