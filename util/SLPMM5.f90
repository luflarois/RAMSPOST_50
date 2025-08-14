
elseif(cvar(1:lv).eq.'slp') then
   ivar_type=2
   ierr= RAMS_getvar('TOPT',idim_type,ngrd,c,b,flnm)
   ierr= RAMS_getvar('PI',idim_type,ngrd,d,b,flnm)
   ierr= RAMS_getvar('THETA',idim_type,ngrd,e,b,flnm)

   call RAMS_comp_slpmm5(n1,n2,n3,e,d,c,a)
   cdname='sea level pressure;'
   cdunits='mb;'







!***************************************************************************
!-------------------------------------------------------------------------
!*rmc Will Cheng's code for calculating slp with mm5's GRAPH method
! ------- added for calculating SLP from MM5 algorithm ------

      subroutine RAMS_comp_slpmm5(n1,n2,n3,theta,pp,z,slp)

!    The subroutine calculates SLP from an algorithm taken from
!    GRAPH, a post-processing packing of MM5 V3.3
!
!    Input: theta - potential temperature (K)         3D
!           pp    - Exner function        (J/kg K)    3D
!           z     - terrain               (m)         2D
!
!    Ouput: SLP   - sea-level pressure    (hPa)       2D

! ------ define dimension of arrays ----------

      dimension theta(n1,n2,n3), pp(n1,n2,n3), z(n1,n2),slp(n1,n2)

! ------ input variables to GRAPH subroutine ----------

      dimension sfp(n1,n2), ts(n1,n2), t_mm5(n1,n2,n3-1), p_mm5(n1,n2,n3-1)

! -----------------------------------------------------

      cp = 1004
      rgas = 287
      cpor = cp / rgas
      p00 = 1.e5

      do j = 1,n2
        do i = 1,n1
!! calculate surface pressure
        sfp(i,j) = (0.5*(pp(i,j,1)+pp(i,j,2))/cp)**cpor*p00*.01
!! calculate surface temp
        ts(i,j) = (0.5/cp)*(theta(i,j,1)*pp(1,j,1)+&
                            theta(i,j,2)*pp(1,j,2))
        enddo
      enddo

      do k = 2,n3
        kk = n3-k+1
        do j = 1,n2
!! flip array upside down for input to GRAPH subroutine
          do i = 1,n1
            t_mm5(i,j,kk) = theta(i,j,k)*pp(i,j,k)/cp
            p_mm5(i,j,kk) = (pp(i,j,k)/cp)**cpor*p00*.01
          enddo
        enddo
      enddo

      call SEAPRS_0(t_mm5,p_mm5,z,sfp,ts,n1,n2,n3-1,slp)

      return
      end


!------------------------------------------------------------------------
      SUBROUTINE SEAPRS_0(T,PP,TER,SFP,TS,IMX,JMX,KX,SLP)
!
!     SECTION  DIAGNOSTIC
!     PURPOSE  COMPUTES SEA LEVEL PRESSURE FROM THE RULE
!              T1/T2=(P1/P2)**(GAMMA*R/G).
!
!     *** LEVELS GO FROM TOP-DOWN ***
!
!     INPUT       T        TEMPERATURE (Kelvin)                3D
!                 TER      TERRAIN     (m)                     2D
!                 SFP      SURFACE PRESSURE (hPa)              2D
!                 IMX      DOT POINT DIMENSION N-S
!                 JMX      DOT POINT DIMENSION E-W
!                 KX       NUMBER OF VERTICAL LEVELS
!
!     OUTPUT      SLP      SEA LEVEL PRESSURE (hPa)            2D
!
      DIMENSION T(IMX,JMX,KX), PP(IMX,JMX,KX),&
                PS(IMX,JMX)  ,SFP(IMX,JMX) , &
                TER(IMX,JMX)
      DIMENSION PL(IMX,JMX),T0(IMX,JMX),TS(IMX,JMX),&
                XKLEV(IMX,JMX)
      DIMENSION SLP(IMX,JMX)
      PARAMETER (R=287.04,G=9.8,GAMMA=6.5E-3)
      PARAMETER (TC=273.16+17.5) ! T CRITICAL IN PSFC/PSLV
      PARAMETER (PCONST=100.)
!
      LOGICAL L1,L2,L3,L4
!
!
!
!
!     ... SEA LEVEL PRESSURE
!
      XTERM=GAMMA*R/G
!
!     ... COMPUTE PRESSURE AT PCONST MB ABOVE SURFACE (PL)
!
      KUPTO=KX/2
99    CONTINUE
      DO 100 J=1,JMX
      DO 100 I=1,IMX
         PL(I,J)=SFP(I,J)-PCONST
         XKLEV(I,J)=0.
100   CONTINUE
!
!     ... FIND 2 LEVELS ON SIGMA SURFACES SURROUNDING PL AT EACH I,J
!
      DO 150 J=1,JMX
      DO 150 I=1,IMX
         DO 125 K=KX-1,KUPTO,-1
            XK=FLOAT(K)
            XKHOLD=XKLEV(I,J)
!srf            XKLEV(I,J)=CVMGT(XK,XKHOLD,   &
            XKLEV(I,J)=merge(XK,XKHOLD,   &
              (((PP(I,J,K)).LT.PL(I,J)) .AND.  &
               ((PP(I,J,K+1)).GE.PL(I,J))))
125      CONTINUE
         IF(XKLEV(I,J).LT.1.) THEN
            PRINT *,'ERROR FINDING PRESSURE LEVEL ',PCONST,' MB ',&
                   'ABOVE THE SURFACE'
            PRINT *,'LAST K LEVEL =',KUPTO
            IF(KUPTO.NE.1) THEN
               PRINT *,'TRYING AGAIN WITH KUPTO=1'
               KUPTO=1
               GOTO 99
            ELSE
               PRINT *,'I,J=',I,J
               PRINT *,'PL=',PL(I,J)
               PRINT *,'PSFC=',SFP(I,J)
               STOP
            END IF
         END IF
150   CONTINUE
!
!     ... GET TEMPERATURE AT PL (TL), EXTRAPOLATE T AT SURFACE (TS)
!         AND T AT SEA LEVEL (T0) WITH 6.5 K/KM LAPSE RATE
!
      DO 200 J=1,JMX
      DO 200 I=1,IMX
         KLO=NINT(XKLEV(I,J))+1
         KHI=NINT(XKLEV(I,J))
         PLO=PP(I,J,KLO)
         PHI=PP(I,J,KHI)
         TLO=T(I,J,KLO)
         THI=T(I,J,KHI)
         TL=THI-(THI-TLO)*ALOG(PL(I,J)/PHI)/ALOG(PLO/PHI)
         TS(I,J)=TL*(SFP(I,J)/PL(I,J))**XTERM
         TBAR=(TS(I,J)+TL)*0.5
         HL=TER(I,J)-R/G*ALOG(PL(I,J)/SFP(I,J))*TBAR
         T0(I,J)=TL+GAMMA*HL
200   CONTINUE
!
!     ... CORRECT SEA LEVEL TEMPERATURE IF TOO HOT
!
      DO 400 J=1,JMX
      DO 400 I=1,IMX
         L1=T0(I,J).LT.TC
         L2=TS(I,J).LE.TC
         L3=.NOT.L1
         T0HOLD=T0(I,J)
!srf         T0(I,J)=CVMGT(T0HOLD,&
!srf           CVMGT(TC,TC-0.005*(TS(I,J)-TC)**2,L2.AND.L3),L1.AND.L2)
         T0(I,J)=merge(T0HOLD,&
           merge(TC,TC-0.005*(TS(I,J)-TC)**2,L2.AND.L3),L1.AND.L2)
400   CONTINUE
!
!     ... COMPUTE SEA LEVEL PRESSURE
!
      DO 600 J=1,JMX
      DO 600 I=1,IMX
         SLP(I,J)=SFP(I,J)*EXP(2.*G*TER(I,J)/(R*(TS(I,J)+T0(I,J))))
600   CONTINUE
      RETURN
      END

!-------------------------------------------------------------------------

