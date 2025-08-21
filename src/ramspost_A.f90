! -------------------------------------------------------------
! -                                                           -
! - RAMSPOST - RAMS Post Processor for GrADS.                 -
! -                                                           -
! -------------------------------------------------------------
! - Adapted for RAMS  4.3 by Saulo R Freitas (SP/14/05/1998)   -
! - Adapted for RAMS  6.0 by Saulo R Freitas (SP/06/03/2005)   -
! - Adapted for BRAMS 5.1 by Saulo R Freitas (SJC/09/09/2015)  -
! -------------------------------------------------------------
program ramspost

    ! -----------------------
    ! -   BASIC PARAMETERS  -
    ! -----------------------

    include 'rcommons.h'
    parameter(maxfiles=2000, nplmax=50, maxpatch=5)

    character*240 fln(maxfiles), fprefix, cfln, gprefix
    character*40 vpln(200), cdum1
    character*16 vp(200), vpun(200), cdum2, proj, anl2gra, mean_type, ascii_data
    character*1 cgrid, patchnumber
    character*2 ccgrid
    character*3 cmo(12)
    character*15 chdate, chstep, xchstep
    integer nvp, nfiles, nzvp(200), nrec, ipresslev, iplevs(nplmax)
    integer inplevs, zlevmax(maxgrds), ndim(200), iproj, ianl2gra
    logical there

    real a(nxpmax, nypmax, nzpmax), b(nxpmax, nypmax, nzpmax), &
        rout(nxpmax, nypmax, nzpmax), &
        zplev(nxpmax, nypmax, nplmax), &
        mytopo(nxpmax, nypmax), mypi(nxpmax, nypmax, nzpmax), &
        rlat(nxpmax, nypmax), rlon(nxpmax, nypmax), &
        lati(maxgrds), latf(maxgrds), loni(maxgrds), lonf(maxgrds), &
        a2(nxpmax, nypmax, nzgmax, maxpatch), &
        rout2(nxpmax, nypmax, nzgmax, maxpatch)
    real site_lat, site_lon

    !  ---  Grid of GRADS
    ! Define maxima dimensao da grade do GRADS =  120% da grade do RAMS
    !      parameter(maxgx=int(1.2*float(nxpmax)),maxgy=int(1.2*float(nypmax)))
    parameter(maxgx=900, maxgy=900)
    integer nxgrads(maxgrds), nygrads(maxgrds), &
        iinf(maxgx, maxgy), jinf(maxgx, maxgy), &
        nxa(maxgrds), nxb(maxgrds), nya(maxgrds), nyb(maxgrds)
    real rmi(maxgx, maxgy, 4), routgrads(maxgx, maxgy, nzpmax)
    real glatg(maxgy), glong(maxgx)
    !  ---

    namelist /rp_input/ fprefix, nvp, vp, gprefix, anl2gra, proj, mean_type, lati, latf, &
        loni, lonf, zlevmax, ipresslev, inplevs, iplevs, ascii_data, site_lat, site_lon

    data cmo/'jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', &
        'oct', 'nov', 'dec'/

    dimension dep_zlev(nzpmax, maxgrds), iep_nx(maxgrds), &
        iep_ny(maxgrds), iep_nz(maxgrds), &
        iep_stdate(6), iep_step(6), &
        dep_glat(2, maxgrds), dep_glon(2, maxgrds)

    common/rams_data/ftimes(maxfiles), nfgpnts(4, maxgrds, maxfiles) &
        , nfgrids(maxfiles), ifdates(maxfiles), iftimes(maxfiles) &
        , flevels(nzpmax, maxgrds, maxfiles), startutc, httop &
        , fdelx(maxgrds, maxfiles), fdely(maxgrds, maxfiles)
    character*240 wfln(maxfiles)

    ! -----------------------------
    ! -   INITIALIZING ROUTINES   -
    ! -----------------------------
    print *, '############################################'
    print *, ' RamsPost - GrADS Visualization for RAMS    '
    print *, '############################################'

    print *, ' '
    print *, 'Opening ramspost.inp file'

    open (5, file='ramspost.inp', status='old')
    read (5, rp_input)

    cgrid = '0'
    iproj = lastchar(proj)
    ianl2gra = lastchar(anl2gra)
    if (anl2gra(1:ianl2gra) == 'one') anl2gra = 'ONE'
    if (ascii_data(1:lastchar(ascii_data)) == 'only') ascii_data = 'ONLY'
    if (ascii_data(1:lastchar(ascii_data)) == 'yes') ascii_data = 'YES'
    ! --- frequencia com as analises serao escrita
    nstep = 1

    nrec = 0
    ic = lastchar(gprefix)
    call RAMS_anal_init(nfiles, fln, fprefix, &
                        dep_zlev, iep_nx, iep_ny, iep_nz, iep_ng, iep_np, &
                        iep_stdate, iep_step, iep_ngrids)
    chdate = '00:00z00mmm1900'
    call RAMS_get_time_init(1, iyear, imonth, idate, ihour, imin)
    call RAMS_get_time_step(iistep, hunit, nfiles)

!  print*,iyear,imonth,idate,ihour,imin
    write (chdate(1:2), '(i2.2)') ihour
    write (chdate(4:5), '(i2.2)') imin
    write (chdate(7:8), '(i2.2)') idate
    write (chdate(12:15), '(i4.2)') iyear
    chdate(9:11) = cmo(imonth) (1:3)

    if (hunit == 1) chstep = '          sec'
    if (hunit == 2) chstep = '          mn'
    if (hunit == 3) chstep = '          hr'
    write (chstep(8:10), '(i3)') iistep

    ! -----------------
    ! -   GRID LOOP   -
    ! -----------------

    ! do ng=1,iep_ngrids
    do ng = iep_ngrids, 1, -1

        print *, '=========================================='
        print *, '           Writing Grid ', ng
        !.................
        iv = 1
        nfn = 1
        cfln = fln(nfn)
        ip = lastchar(cfln) - 9
        !.................
        !   rlat and rlon = lat and lon of "thermodynamic points" of RAMS model.
        !
        call ep_getvar('lat', &
                       rlat, a, b, iep_nx(ng), iep_ny(ng), &
                       1, ng, cfln(1:ip), vpln(iv), &
                       vpun(iv), n, iep_np, iep_ng, a2, rout2)
        call ep_getvar('lon', &
                       rlon, a, b, iep_nx(ng), iep_ny(ng), &
                       1, ng, cfln(1:ip), vpln(iv), &
                       vpun(iv), n, iep_np, iep_ng, a2, rout2)
        !.................
        call geo_grid(iep_nx(ng), iep_ny(ng), rlat, rlon, &
                      dep_glon(1, ng), dep_glon(2, ng), &
                      dep_glat(1, ng), dep_glat(2, ng), &
                      rlatmin, rlatmax, rlonmin, rlonmax, &
                      nxgrads(ng), nygrads(ng), &
                      proj(1:iproj), glatg, glong, maxgx, maxgy)
        !print*,'XXX',nxgrads(ng),nygrads(ng),iep_nx(ng),iep_ny(ng)
        !print*, dep_glon(1,ng),dep_glon(2,ng),dep_glat(1,ng),dep_glat(2,ng)
        !pause
        !.................
        !
        call Matriz_interp(ng, nxgrads(ng), nygrads(ng), &
                           iep_nx(ng), iep_ny(ng), &
                           !dep_glat(1,ng),dep_glat(2,ng),     &
                           !dep_glon(1,ng),dep_glon(2,ng),     &
                           iinf, jinf, rmi, &
                           proj(1:iproj), mean_type(1:lastchar(mean_type)) &
                           , glatg, glong, maxgx, maxgy)

        !_................
        call define_lim(ng, nxgrads(ng), nygrads(ng), &
                        dep_glat(1, ng), dep_glat(2, ng), &
                        dep_glon(1, ng), dep_glon(2, ng), &
                        lati(ng), latf(ng), loni(ng), lonf(ng), &
                        nxa(ng), nxb(ng), nya(ng), nyb(ng), proj(1:iproj), &
                        iep_nx(ng), iep_ny(ng), rlat, rlon, glatg, glong, maxgx, maxgy)
        !            print*,nxgrads(ng),nygrads(ng),iep_nx(ng),iep_ny(ng) &
        !                   ,iep_nz(ng),iv,vpln(iv), nxa(ng),nxb(ng),nya(ng),nyb(ng),&
        !                         vpun(iv),n
        !     ----------------------
        !     - WRITE GRADS BINARY -
        !     ----------------------
        !     ----------------------
        write (cgrid, '(i1)') ng

        if (anl2gra(1:ianl2gra) /= 'ONE' .and. &
            ascii_data(1:lastchar(ascii_data)) /= 'ONLY') then
            !print *,'Pos 174: open ',gprefix(1:ic)//'_g'//cgrid//'.gra'
            open (19, file=gprefix(1:ic)//'_g'//cgrid//'.gra', &
                  form='unformatted', access='direct', status='replace', &
                  !use this for intel compiler
                  !           recl=1*(nxb(ng)-nxa(ng)+1)*(nyb(ng)-nya(ng)+1))
                  !use this for pgi compiler
                  recl=4*(nxb(ng) - nxa(ng) + 1)*(nyb(ng) - nya(ng) + 1))
            nrec = 0
        end if

        nfn_start = 0

        do nfn = 1, nfiles, nstep
            print *, ' --- > timestep ', nfn

            cfln = fln(nfn)
            ip = lastchar(cfln) - 9
            if (nfn_start == 0) nfn_start = 1
            call date1(ifdates(nfn), iyear, imon, idate)
            !.................
            !.................
            !     ----------------------
            !     - WRITE GRADS BINARY -
            !     ----------------------
            if (anl2gra(1:ianl2gra) == 'ONE' .and. &
                ascii_data(1:lastchar(ascii_data)) /= 'ONLY') then
                call makefnam(wfln(nfn), gprefix(1:ic)//' ', 0, iyear, imon, idate, &
                              iftimes(nfn), 'A', 'g'//cgrid, 'gra')

                !print*,iyear,imon,idate,iftimes(nfn),wfln(nfn),ifdates(nfn)

                !inquire (file=wfln(nfn), exist=there)
                !if (.not. there) then
                    iunit = 19
                    !print *,'Pos 208: open ',wfln(nfn)
                    open (iunit, file=wfln(nfn), form='unformatted', access='direct' &
                          , status='replace', recl=4*(nxb(ng) - nxa(ng) + 1)*(nyb(ng) - nya(ng) + 1))
                    nrec = 0
                    print *, ' --- > doing timestep ', nfn
                    print *, ' --- >'
                !else
                !    print *, 'file=', wfln(nfn) (1:lastchar(wfln(nfn))), ' already exist'
                !    nfn_start = nfn + 1
                !    !print*,'nfn_start=',nfn_start,nfiles
                !    cycle
                !end if

            end if
!.................

            if (ipresslev > 0) then

                call ep_getvar('topo', mytopo, a, b, iep_nx(ng), iep_ny(ng), &
                               1, ng, cfln(1:ip), vpln(iv), &
                               vpun(iv), n, iep_np, iep_ng, a2, rout2)

                call ep_getvar('pi', mypi, a, b, iep_nx(ng), iep_ny(ng), &
                               iep_nz(ng), ng, cfln(1:ip), vpln(iv), &
                               vpun(iv), n, iep_np, iep_ng, a2, rout2)

            end if

            nnvp = 0
            do iv = 1, nvp
                !.................
                print *, 'Get variable:  ', iv, vp(iv)
                call ep_getvar(vp(iv), &
                               rout, a, b, iep_nx(ng), iep_ny(ng), &
                               iep_nz(ng), ng, cfln(1:ip), vpln(iv), &
                               vpun(iv), ndim(iv), iep_np, iep_ng, a2, rout2)

                if (ndim(iv) == 0) cycle
                nnvp = nnvp + 1
                !if(ndim(iv)==0) then
                !   print*,'Variable "',trim(vp(iv)),'" not available in this run/grid'
                !else
                !   print*,'Variable "',trim(vp(iv)),'" ok, dimension=',ndim(iv)
                !endif

                if (ascii_data == 'YES' .or. ascii_data == 'ONLY') then

                    call printout(vp(iv), iep_nx(ng), iep_ny(ng), iep_nz(ng), ng, rlat, rlon, rout &
                                  , cgrid, gprefix, iyear, imon, idate, iftimes(nfn), ifdates(nfn), nfn, n &
                                  , site_lat, site_lon)

                    if (ascii_data == 'ONLY') cycle
                end if
                !................
                !.....
                if (ndim(iv) == 8) then
                    nzvp(iv) = iep_ng
                    do ipatch = 1, iep_np
                        call S4d_to_3d(iep_nx(ng), iep_ny(ng), nzvp(iv), iep_np, &
                                       ipatch, rout, rout2)
                        call proj_rams_to_grads(vp(iv), ndim(iv), &
                                                iep_nx(ng), iep_ny(ng), nzvp(iv), &
                                                nxgrads(ng), nygrads(ng), &
                                                rmi, iinf, jinf, &
                                                rout, routgrads, rlat, rlon, proj(1:iproj))
                        call ep_putvar(routgrads, a, nxgrads(ng), nygrads(ng), &
                                       nxa(ng), nxb(ng), nya(ng), nyb(ng), &
                                       nzvp(iv), nrec, 1, iep_ng)
                        nnvp = nnvp + 1
                    end do
                    nnvp = nnvp - 1
                    !.....
                    !
                elseif (ndim(iv) == 2) then
                    nzvp(iv) = 1
                    call proj_rams_to_grads(vp(iv), ndim(iv), &
                                            iep_nx(ng), iep_ny(ng), nzvp(iv), &
                                            nxgrads(ng), nygrads(ng), &
                                            rmi, iinf, jinf, &
                                            rout, routgrads, rlat, rlon, proj(1:iproj))
                    call ep_putvar(routgrads, a, nxgrads(ng), nygrads(ng), &
                                   nxa(ng), nxb(ng), nya(ng), nyb(ng), &
                                   nzvp(iv), nrec, 1, 1)
              !!.....
                    !
                elseif (ndim(iv) == 7) then
                    nzvp(iv) = iep_np
                    call proj_rams_to_grads(vp(iv), ndim(iv), &
                                            iep_nx(ng), iep_ny(ng), nzvp(iv), &
                                            nxgrads(ng), nygrads(ng), &
                                            rmi, iinf, jinf, &
                                            rout, routgrads, rlat, rlon, proj(1:iproj))
                    do ipatch = 1, iep_np
                        call ep_putvar(routgrads, a, nxgrads(ng), nygrads(ng), &
                                       nxa(ng), nxb(ng), nya(ng), nyb(ng), &
                                       nzvp(iv), nrec, ipatch, ipatch)
                        nnvp = nnvp + 1
                    end do
                    !if(nfn.eq.1) nnvp=nnvp-1
                    !if(nfn == nfn_start)
                    nnvp = nnvp - 1
              !!.....
                    !
                elseif (ndim(iv) == 5) then
                    !             nzvp(iv)=iep_ng
                    !               Call proj_rams_to_grads(vp(iv),ndim(iv),       &
                    !                             iep_nx(ng),iep_ny(ng),nzvp(iv),  &
                    !                             nxgrads(ng),nygrads(ng),              &
                    !                             rmi,iinf,jinf,                      &
                    !                             rout,routgrads,rlat,rlon,proj(1:iproj))
                    !
                    !                Call ep_putvar(routgrads,a,nxgrads(ng),nygrads(ng),  &
                    !                             nxa(ng),nxb(ng),nya(ng),nyb(ng),            &
                    !                              nzvp(iv),nrec,1,iep_ng)
              !!.....
                    !
                elseif (ndim(iv) == 3) then
                    nzvp(iv) = iep_nz(ng)
                    !..
                    if (ipresslev == 1) then
                        !print *, 'before ptransvar, rout', maxval(rout), minval(rout)
                        inplevsef = inplevs
                        !do i=1,iep_nx(ng)
                        !   do j=1,iep_ny(ng)
                        !      write(86,*) i,j,rlon(i,j),rlat(i,j),rout(i,j,1)
                        !   end do
                        !end do
                        call ptransvar(rout, iep_nx(ng), iep_ny(ng), nzvp(iv), &
                                       inplevs, iplevs, mypi, dep_zlev(1, ng), zplev, mytopo)
                        !print *, 'before proj_rams_to_grads, rout', maxval(rout), minval(rout)
                        !do i=1,iep_nx(ng)
                        !   do j=1,iep_ny(ng)
                        !      write(87,*) i,j,rlon(i,j),rlat(i,j),rout(i,j,1)
                        !   end do
                        !end do
                        call proj_rams_to_grads(vp(iv), ndim(iv), &
                                                iep_nx(ng), iep_ny(ng), nzvp(iv), &
                                                nxgrads(ng), nygrads(ng), &
                                                rmi, iinf, jinf, &
                                                rout, routgrads, rlat, rlon, proj(1:iproj))
                        !print *, 'before eputvar, routgrads', maxval(routgrads), minval(routgrads)
                        !do i=1,iep_nx(ng)
                        !   do j=1,iep_ny(ng)
                        !      write(88,*) i,j,rlon(i,j),rlat(i,j),rout(i,j,1)
                        !   end do
                        !end do
                        call ep_putvar(routgrads, a, nxgrads(ng), nygrads(ng), &
                                       nxa(ng), nxb(ng), nya(ng), nyb(ng), &
                                       inplevsef, nrec, 1, inplevsef)

                    elseif (ipresslev == 2) then

                        inplevsef = inplevs
                        call ctransvar(iep_nx(ng), iep_ny(ng), iep_nz(ng), rout &
                                       , mytopo, inplevs, iplevs, ztn(1, ng), zmn(nnzp(1) - 1, 1))

                        call proj_rams_to_grads(vp(iv), ndim(iv), &
                                                iep_nx(ng), iep_ny(ng), nzvp(iv), &
                                                nxgrads(ng), nygrads(ng), &
                                                rmi, iinf, jinf, &
                                                rout, routgrads, rlat, rlon, proj(1:iproj))

                        call ep_putvar(routgrads, a, nxgrads(ng), nygrads(ng), &
                                       nxa(ng), nxb(ng), nya(ng), nyb(ng), &
                                       inplevsef, nrec, 1, inplevsef)

                    elseif (ipresslev == 3) then

                        inplevsef = inplevs
                        call select_sigmaz(iep_nx(ng), iep_ny(ng), iep_nz(ng), rout &
                                           , inplevs, iplevs)

                        call proj_rams_to_grads(vp(iv), ndim(iv), &
                                                iep_nx(ng), iep_ny(ng), nzvp(iv), &
                                                nxgrads(ng), nygrads(ng), &
                                                rmi, iinf, jinf, &
                                                rout, routgrads, rlat, rlon, proj(1:iproj))

                        call ep_putvar(routgrads, a, nxgrads(ng), nygrads(ng), &
                                       nxa(ng), nxb(ng), nya(ng), nyb(ng), &
                                       inplevsef, nrec, 1, inplevsef)
                 !!..
                    else
                 !!..
                        call proj_rams_to_grads(vp(iv), ndim(iv), &
                                                iep_nx(ng), iep_ny(ng), nzvp(iv), &
                                                nxgrads(ng), nygrads(ng), &
                                                rmi, iinf, jinf, &
                                                rout, routgrads, rlat, rlon, proj(1:iproj))
                        if (zlevmax(ng) >= iep_nz(ng)) zlevmax(ng) = iep_nz(ng) - 1

                        call ep_putvar(routgrads, a, nxgrads(ng), nygrads(ng), &
                                       nxa(ng), nxb(ng), nya(ng), nyb(ng), &
                                       nzvp(iv), nrec, 2, zlevmax(ng) + 1)
                 !!..
                    end if
              !!
                end if
           !!................
                !
                !
                !

                !print*,nfn,nfn_start,nnvp

            end do
            print *, '-------------------------------------------------------------'
            print *, 'Grid=', ng, ' Bytes=', 4*nrec* &
                (nxb(ng) - nxa(ng) + 1)*(nyb(ng) - nya(ng) + 1)
            print *, '-------------------------------------------------------------'
            !.................

            if (anl2gra(1:ianl2gra) == 'ONE' .and. ascii_data(1:lastchar(ascii_data)) /= 'ONLY') close (19)
        end do ! enddo do NFILES

        if (anl2gra(1:ianl2gra) /= 'ONE' .and. ascii_data(1:lastchar(ascii_data)) /= 'ONLY') close (19)

        !.................
        !     -----------------------
        !     - WRITE GRADS CONTROL -
        !     -----------------------

! if only ascii data => no grads control files needed => cycle
        if (ascii_data == 'ONLY' .or. nnvp == 0) cycle

        write (cgrid, '(i1)') ng
        iunit = 20
! if all files already done => cycle
        if (nfn_start - 1 == nfiles) cycle

        do nfn = 1, nfiles, nstep

            iuniti = iunit
            iunitf = iunit

! case 1 : all analysis at only one grads file
!----
            if (anl2gra(1:ianl2gra) /= 'ONE' .and. anl2gra(1:ianl2gra) /= 'one' .and. &
                nfn == 1) then
                !print *,'Pos 447: open ',gprefix(1:ic)//'_g'//cgrid//'.ctl'
                open (iunit, file=gprefix(1:ic)//'_g'//cgrid//'.ctl', &
                      status='replace')
                write (iunit, 2001) '^'//gprefix(1:ic)//'_g'//cgrid//'.gra'
! case 2 : one analysis at one grads file
            elseif (anl2gra(1:ianl2gra) == 'ONE' .or. anl2gra(1:ianl2gra) == 'one') then

! case 2 with template
                if (nfn == 1 .and. nfiles > 1) then
                    call date1(ifdates(nfn), iyear, imon, idate)
                    call makefnam(wfln(nfn), gprefix(1:ic)//'-template'//' ', 0, iyear, imon, idate, &
                                  iftimes(nfn), 'A', 'g'//cgrid, 'ctl')
                    open (iunit + 1, file=wfln(nfn), status='unknown')
! valido somente para hora cheia  --------------------------------------vvvv-
                    write (iunit + 1, 2001) '^'//gprefix(1:ic)//'-A-'//'%y4-%m2-%d2-%h20000-'//'g'//cgrid//'.gra'
                    write (iunit + 1, 2002) 'options template'
                    iunitf = iunit + 1
                end if
!
! individual files
                call date1(ifdates(nfn), iyear, imon, idate)
                call makefnam(wfln(nfn), gprefix(1:ic)//' ', 0, iyear, imon, idate, &
                              iftimes(nfn), 'A', 'g'//cgrid, 'ctl')
                inquire (file=wfln(nfn), exist=there)
                if (there .and. nfn > 1) cycle

                open (iunit, file=wfln(nfn), status='unknown')
                write (iunit, 2001) '^'//wfln(nfn) (1:lastchar(wfln(nfn)) - 3)//'gra'

            end if
!----

            do iunit = iuniti, iunitf

                write (iunit, 2002) 'undef -9.99e33'
                write (iunit, 2022) gprefix(1:ic)
!srf   write(iunit,2002) 'title RAMS 4.2 Output'
                write (iunit, 2003) nxb(ng) - nxa(ng) + 1, (dep_glon(i, ng), i=1, 2)
                write (iunit, 2004) nyb(ng) - nya(ng) + 1, (dep_glat(i, ng), i=1, 2)

                !print*,' lon-lat=',dep_glon(:,ng),dep_glat(:,ng)

                if (ipresslev > 0 .and. ipresslev <= 2) then
                    write (iunit, 2005) inplevs, (iplevs(i)*1.0, i=1, inplevs)
                elseif (ipresslev == 3) then
                    write (iunit, 2005) inplevs, (dep_zlev(iplevs(i), ng), i=1, inplevs)
                else
                    if (zlevmax(ng) + 1 < 15) then
                        write (iunit, 2005) zlevmax(ng), &
                            (dep_zlev(i, ng), i=2, zlevmax(ng) + 1)
                    else
                        write (iunit, 2005) zlevmax(ng), (dep_zlev(i, ng), i=2, 15)
                        write (iunit, 2055) (dep_zlev(i, ng), i=16, zlevmax(ng) + 1)
                    end if
                end if
! case 1
                if (anl2gra(1:ianl2gra) /= 'ONE' .and. anl2gra(1:ianl2gra) /= 'one' .and. &
                    nfn == 1) then
                    write (iunit, 2006) nfiles, chdate, chstep
! case 2
                elseif (anl2gra(1:ianl2gra) == 'ONE' .or. anl2gra(1:ianl2gra) == 'one') then

                    !call RAMS_get_time_init(nfn,iyear,imonth,idate,ihour,imin)
                    call date1(ifdates(nfn), iyear, imon, idate)
                    write (chdate(1:2), '(i2.2)') iftimes(nfn)/10000!ihour
                    write (chdate(4:5), '(i2.2)') imin
                    write (chdate(7:8), '(i2.2)') idate
                    write (chdate(12:15), '(i4.2)') iyear
                    chdate(9:11) = cmo(imon) (1:3)

                    ! print*,'chdate',nfn,chdate,ifdates(nfn),iftimes(nfn)

                    if (iunitf == iuniti) write (iunit, 2006) 1, chdate, chstep
                    if (iunitf == iuniti + 1) write (iunit, 2006) nfiles, chdate, chstep ! para template

                end if
!----

                write (iunit, 2007) nnvp
                do i = 1, nvp
                    !print*,'ndim(i)=',i,vp(i),ndim(i)
                    if (ndim(i) == 0) cycle

                    if (ipresslev > 0 .and. nzvp(i) == iep_nz(ng)) then
                        write (iunit, 2008) vp(i), inplevs, vpln(i), vpun(i)
                    else
                        if (ndim(i) == 3) &
                            write (iunit, 2008) vp(i), zlevmax(ng), vpln(i), vpun(i)

                        if (ndim(i) == 2) &
                            write (iunit, 2008) vp(i), 0, vpln(i), vpun(i)

                        if (ndim(i) == 7) then
                            il = lastchar(vp(i))
                            il2 = lastchar(vpln(i))
                            do ipatch = 1, iep_np
                                write (patchnumber, '(i1)') ipatch
                                cdum2 = vp(i) (1:il)//patchnumber
                                cdum1 = vpln(i) (1:il2)//': patch # '//patchnumber
                                write (iunit, 2008) cdum2, 0, cdum1, vpun(i)
                            end do
                        end if
                        if (ndim(i) == 8) then
                            il = lastchar(vp(i))
                            il2 = lastchar(vpln(i))
                            do ipatch = 1, iep_np
                                write (patchnumber, '(i1)') ipatch
                                cdum2 = vp(i) (1:il)//patchnumber
                                cdum1 = vpln(i) (1:il2)//': patch # '//patchnumber
                                write (iunit, 2008) cdum2, nzvp(i), cdum1, vpun(i)
                            end do
                        end if

                        if (ndim(i) == 5) &
                            write (iunit, 2008) vp(i), nzvp(i), vpln(i), vpun(i)

                    end if
                end do
                write (iunit, 2002) 'endvars'
                close (iunit)

            end do ! enddo nas unidades de escrita

!     endif

            if (anl2gra(1:ianl2gra) /= 'ONE' .and. &
                nfn == 1) exit
        end do ! enddo do NFILES

    end do  !enddo do NGRIDS

2001 format('dset ', a)
2002 format(a)
2022 format("title ", a)
2003 format('xdef ', i4, ' linear ', 2f15.4)
2004 format('ydef ', i4, ' linear ', 2f15.4)
2005 format('zdef ', i4, ' levels ', 60f10.1)
2006 format('tdef ', i4, ' linear ', 2a15)
2007 format('vars ', i4)
2008 format(a16, i4, ' 99    - RAMS : ', a40, '[', a8, ']')
2055 format(60f7.0)
    close (iunit + 1)
    stop
end program ramspost

! ---------------------------------------------------------------
! -   SUBROUTINE EP_GETVAR : LOAD RAMS VARIABLE FROM ANALYSIS   -
! ---------------------------------------------------------------

subroutine ep_getvar(cvar, rout, a, b, &
                     nx, ny, nz, ng, fn, cdname, cdunits, itype, &
                     npatch, nzg, a2, rout2)
    dimension a(nx, ny, nz), b(nx, ny, nz), rout(nx, ny, nz), &
        a2(nx, ny, nzg, npatch), rout2(nx, ny, nzg, npatch)
    character*(*) cvar, fn, cdname, cdunits

    call RAMS_varlib(cvar, nx, ny, nz, nzg, npatch, ng, fn, &
                     cdname, cdunits, itype, a, b, a2)

    if (itype == 1) return

    if (itype == 8) then
        do i = 1, nx
            do j = 1, ny
                do izg = 1, nzg
                    do ip = 1, npatch
                        rout2(i, j, izg, ip) = a2(i, j, izg, ip)
                        !             print*,i,j,izg,ip,a2(i,j,izg,ip)
                    end do
                end do
            end do
        end do
    else
        do k = 1, nz
            do i = 1, nx
                do j = 1, ny
                    rout(i, j, k) = a(i, j, k)
                    !cc
                    !         if(a(i,j,k).gt.1.) &
                    !print*,'Var= ',cvar,i,j,k,a(i,j,k)
                    !cc
                end do
            end do
        end do
    end if
    !print *, 'LFR-DBG rout:', maxval(rout), minval(rout)
    !xxxxxxxxxxxxxxxxxxxxxxxxx
    !      k=1
    !        write(17,'(52f10.3)')((rout(ii,jj,k),ii=1,nx),jj=1,ny)
    !
    !xxxxxxxxxxxxxxxxxxxxxxxxx
    return
end subroutine ep_getvar
!
!--------------------------------------------------------------------------
subroutine printout(cvar, nxr, nyr, nzr, ng, rlat, rlon, rout, cgrid, gprefix, &
                    iyear, imon, idate, iftimes, ifdates, nfn, ndim, xlat, xlon)
    integer nxr, nyr, nzr, ng, nfn, ndim, ix, jx
    dimension rlat(nxr, nyr), rlon(nxr, nyr), rout(nxr, nyr, nzr)
    real xlat, xlon
    character*(*) cvar, cgrid, gprefix
    integer iyear, imon, idate, iftimes, ifdates, ihour, imin, isec, klevel
    include 'rcommons.h'
    character*256 filename
    klevel = 1
    if (ndim == 3) klevel = 2

    filename = cvar(1:lastchar(cvar))//'_'//gprefix(1:lastchar(gprefix))//'_' &
               //'g'//cgrid//'.dat'
    ihour = int(iftimes/10000)
    imin = int((iftimes - ihour*10000)/100)

    open (2, file=filename, status='unknown', position='append')

    if (nfn == 1) then
!       xlat=-10.76
!       xlon=-62.36
        call get_latlon(nxr, nyr, ng, xlat, xlon, ix, jx)
!       write(2,'(A6,A20)') 'var=',cvar(1:lastchar(cvar))
        write (2, '(A25,2f10.3)') 'Site : lat,lon           =', xlat, xlon
        write (2, '(A25,2f10.3,3i4)') 'Model: lat,lon,ix,jx,klev=', rlat(ix, jx), rlon(ix, jx), ix, jx &
            , klevel
        write (2, '(A21,A20)') 'year,mon,day,hour,min', cvar(1:lastchar(cvar))
    end if

!     print*,'ponto mais proximo de (lat,lon):',xlat,xlon
!     print*,'� (i,j,lat,lon):', ix,jx,rlat(ix,jx),rlon(ix,jx)
!     print*,'valor da var=',rout(ix,jx,1)
!     print*,'======================================================'
!     print*,'======================================================'

    write (2, '(5i4,e18.6)') iyear, imon, idate, ihour, imin, rout(ix, jx, klevel)
    close (2)

    return
    return
end subroutine printout
!***************************************************************
subroutine get_latlon(nxr, nyr, ng, xlat, xlon, ix, jx)
    include 'rcommons.h'
    ix = 1
    jx = 1

    call ge_to_xy(polelat, polelon, xlon, xlat, x, y)
!     print*,xlat,xlon,x,y,ng,nxr,nyr
!     print*,xtn(1,ng),xtn(nxr,ng),ytn(1,ng),ytn(nyr,ng)

    if (x < xtn(1, ng) .or. x > xtn(nxr, ng)) goto 777
    if (y < ytn(1, ng) .or. y > ytn(nyr, ng)) goto 777

    do i = 1, nxr
        if (x <= xtn(i, ng)) go to 555
    end do
555 continue
    i1 = i - 1
    i2 = i

    do j = 1, nyr
        if (y <= ytn(j, ng)) go to 666
    end do
666 continue
    j1 = j - 1
    j2 = j

    if (x - xtn(i1, ng) <= 0.5*deltaxn(ng)) then
        ix = i1
    else
        ix = i2
    end if
    if (y - ytn(j1, ng) <= 0.5*deltayn(ng)) then
        jx = j1
    else
        jx = j2
    end if
    return
777 continue
    print *, 'ponto fora'
    stop 3344

end subroutine get_latlon

! -------------------------------------------------
! -   SUBROUTINE EP_SETDATE : PARSE DATE RECORD   -
! -------------------------------------------------

subroutine ep_setdate(iyear1, imonth1, idate1, strtim, itrec)
    real time

    integer itrec(6)
    itrec(1) = iyear1
    itrec(2) = imonth1
    itrec(3) = idate1
    itrec(4) = int(mod(strtim, 24.))
    itrec(5) = int(mod(strtim, 1.)*60)
    itrec(6) = int(mod((strtim)*3600., 60.))

    !      print*,'---------------------------------------'
    !      print*,itrec(1),itrec(2),itrec(3),itrec(4),itrec(5),
    !     +       itrec(6)
    !      print*,'---------------------------------------'

    return
end subroutine ep_setdate

!***************************************************************
!     New version that just uses ls and the /tmp directory

subroutine RAMS_filelist(fnames, file_prefix, nfile)
    character fnames(*)*(*), file_prefix*(*)
    character file*240, command*240

    ! This version uses nfile as flag for whether to stop if no files exist.
    !    If nfile.ge.0, then stop

    iflag = nfile

    nfile = 0
    print *, ' '
    print *, 'RAMS file analysis: Checking directory - ', file_prefix

    iprelen = index(file_prefix, ' ')
    if (iprelen == 0) iprelen = len(file_prefix)

    command = &
        !-srf temporary filelist is now written in the local directory
        !      '/bin/ls -1 '//file_prefix(1:iprelen)//' >/tmp/RAMS_filelist'
        '/bin/ls -1 '//file_prefix(1:iprelen)//' >./RAMS_filelist'
    call system(command)
    command = 'chmod 777 ./RAMS_filelist'
! command= 'chmod 777 /tmp/RAMS_filelist'
    call system(command)

    !     Open the directory list and read through the files
    iun = 98
! open(unit=iun,file='/tmp/RAMS_filelist',status='old',err=15)
    open (unit=iun, file='./RAMS_filelist', status='old', err=15)
    rewind iun

    do nf = 1, 1000000
        read (iun, '(a128)', end=30, err=30) file
        fnames(nf) = file
    end do

30  continue

    close (iun)
! command= '/bin/rm -f /tmp/RAMS_filelist'
    command = '/bin/rm -f ./RAMS_filelist'
    call system(command)

    nfile = nf - 1

    if (nfile == 0) then
        print *, 'No RAMS files for prefix:', file_prefix
        if (iflag >= 0) stop 'RAMS_filelist-no_files'
    end if

    return

15  print *, 'RAMS_filelist: Error opening ./RAMS_filelist'
!15 print *, 'RAMS_filelist: Error opening /tmp/RAMS_filelist'
    stop 'RAMS_filelist-/tmp file error : run again'
    return

100 continue

    return
end subroutine RAMS_filelist

!*****************************************************************************

! --------------------------------------------------------
! -   SUBROUTINE EP_PUTVAR : WRITE ARRAY TO GRADS FILE   -
! --------------------------------------------------------

subroutine ep_putvar(rout, a, nx, ny, nxa, nxb, nya, nyb, &
                     nz, nrec, istartz, iendz)
    dimension a(nx, ny), rout(nx, ny, nz)
    integer istartz, iendz
    !
    do k = istartz, iendz
        do j = 1, ny
            do i = 1, nx
                a(i, j) = rout(i, j, k)
                !cc
                !            print*,'PUT VAR=',i,j,k,a(i,j)
                !cc
            end do
        end do
        nrec = nrec + 1
        write (19, rec=nrec) ((a(i, j), i=nxa, nxb), j=nya, nyb)
        !       write(19,rec=nrec) a
    end do

    !      k=1
    !        write(18,'(59f10.3)')((rout(ii,jj,k),ii=1,nx),jj=1,ny)
    !
    return
end subroutine ep_putvar

!-------------------------------------------------------------------
!
subroutine Matriz_interp(ng, nxg, nyg, nxr, nyr, &!rlat1,dlat, &
                         !rlon1,dlon,
                         iinf, jinf, rmi, proj, mean_type, glatg, glong, maxgx, maxgy)
    implicit none
    include 'rcommons.h'
    character*(*) proj, mean_type
    integer itype_proj ! = 0 normal, =1 usa o ponto mais proximo
    integer :: maxgx, maxgy, ng, nxg, nyg, nxr, nyr
    real rm
    real :: glatg(maxgy), glong(maxgx)

!  common/grid2/ glatg(nypmax),glong(nxpmax)
    real:: rmi(nxg, nyg, 4)
    integer ::iinf(nxg, nyg), jinf(nxg, nyg)
    !real rlat1,dlat,rlon1,dlon
    integer i, j, l, ix, i1, i2, iy, j1, j2
    real x, y, undef
    !
    if (proj /= 'YES' .and. proj /= 'yes') return

!!  itype_proj = 1

    !       Construcao da matriz de interpolacao.
    !       Flag para pontos do grads fora do dominio do modelo

!!!  if(itype_proj==0) then

    undef = -9.99e+15
    do i = 1, nxg
        do j = 1, nyg
            iinf(i, j) = 1
            jinf(i, j) = 1
            do l = 1, 4
                rmi(i, j, l) = undef
            end do
        end do
    end do

!!  else
!!      rmi(:,:,:)=-9.99e+15
!!  endif

    do i = 1, nxg
        do j = 1, nyg
            !       Encontra posicao do ponto de grade do GRADS na grade do RAMS
            !
            !        glatg(i)=-37.113
            !        glong(j)=-79.128
            !        xlat=glatg(i)
            !        xlon= glong(j)
            !        Call ll_xy(glatg(i),glong(j),polelat,polelon,x,y)
            !       call getops(pla,plo,xlat,xlon,polelat,polelon)
            !       call pstoxy(x,y,pla,plo,6376000.)
            !
            call ge_to_xy(polelat, polelon, glong(i), glatg(j), x, y)
            !
            !        print*,x,y,glong(i),glatg(j),polelat,polelon

            !       Elimina pontos fora:
            if (x < xtn(1, ng) .or. x > xtn(nxr, ng)) go to 777
            if (y < ytn(1, ng) .or. y > ytn(nyr, ng)) go to 777
            !
            do ix = 1, nxr
                if (x <= xtn(ix, ng)) go to 555
            end do
555         continue
            i1 = ix - 1
            i2 = ix
            iinf(i, j) = i1

            do iy = 1, nyr
                if (y <= ytn(iy, ng)) go to 666
            end do
666         continue
            j1 = iy - 1
            j2 = iy
            jinf(i, j) = j1
            !

!!!!       if(itype_proj==0) then
            if (mean_type == 'BAV' .or. mean_type == 'bav') then
!projecao pela media ponderada
                rmi(i, j, 1) = (x - xtn(i1, ng))/deltaxn(ng)
                rmi(i, j, 2) = 1.-rmi(i, j, 1)
                !
                rmi(i, j, 3) = (y - ytn(j1, ng))/deltayn(ng)
                rmi(i, j, 4) = 1.-rmi(i, j, 3)

            elseif (mean_type == 'VMP' .or. mean_type == 'vmp') then

!-srf-set2005 - projecao pelo ponto mais proximo
                if (x - xtn(i1, ng) <= 0.5*deltaxn(ng)) then
                    rmi(i, j, 1) = 0.
                else
                    rmi(i, j, 1) = 1.
                end if

                if (y - ytn(j1, ng) <= 0.5*deltayn(ng)) then
                    rmi(i, j, 3) = 0.
                else
                    rmi(i, j, 3) = 1.
                end if
                rmi(i, j, 2) = 1.-rmi(i, j, 1)
                rmi(i, j, 4) = 1.-rmi(i, j, 3)

            else
                print *, 'MEAN_TYPE not available:'
                print *, 'Use VMP or BAV'
                print *, 'stop at Matriz_interp routine'
                stop 1234

            end if

777         continue
            !
!                 if(rmi(i,j,1)+rmi(i,j,2)+rmi(i,j,3)+rmi(i,j,4)> 0)then
!                 if(abs(rmi(i,j,1)+rmi(i,j,2)+rmi(i,j,3)+rmi(i,j,4)).ne.2)then
!                 print*,'XXXXXXXXXXXXXXXXXXXXXXX'
!                 stop 33333
!                 print*,i,j,rmi(i,j,1)+rmi(i,j,2)+rmi(i,j,3)+rmi(i,j,4)
!                 print*,rmi(i,j,1),rmi(i,j,2),rmi(i,j,3),rmi(i,j,4)
!                 endif
!                 endif
            !
            !
        end do
    end do
    return
end subroutine Matriz_interp

!*************************************************************************

subroutine proj_rams_to_grads(vp, n, nxr, nyr, nzz, nxg, nyg, &
                              rmi, iinf, jinf, &
                              rout, routgrads, rlat, rlon, proj)

    character*(*) proj
    character*10 vp
    dimension rlat(nxr, nyr), rlon(nxr, nyr)
    dimension rout(nxr, nyr, nzz), routgrads(nxg, nyg, nzz)
    dimension rmi(nxg, nyg, 4), iinf(nxg, nyg), jinf(nxg, nyg)

    if (proj /= 'YES' .and. proj /= 'yes') then
        if (nxg /= nxr .and. nyg /= nyr) then
            print *, 'Projection with problems nxr nxg ...'
            stop
        end if
        call rout_to_routgrads(nxr*nyr*nzz, rout, routgrads)
        return
    end if

    do i = 1, nxg
        do j = 1, nyg

            !

            !if(rmi(i,j,1)+rmi(i,j,2)+rmi(i,j,3)+rmi(i,j,4) < 0.) then
            if (rmi(i, j, 1) < 0. .or. &
                rmi(i, j, 2) < 0. .or. &
                rmi(i, j, 3) < 0. .or. &
                rmi(i, j, 4) < 0.) then

                routgrads(i, j, :) = -9.99e+33
                cycle
            end if
            !print*,'-----',i,j
            !print*,i,j,rmi(i,j,1),rmi(i,j,2),rmi(i,j,3),rmi(i,j,4)
            !
            !
            r1 = rmi(i, j, 1)
            r2 = rmi(i, j, 2)
            r3 = rmi(i, j, 3)
            r4 = rmi(i, j, 4)
            i1 = iinf(i, j)
            i2 = i1 + 1
            j1 = jinf(i, j)
            j2 = j1 + 1

            do k = 1, nzz

                if (rout(i1, j1, k) == -9.99e+33 .or. &
                    rout(i2, j1, k) == -9.99e+33 .or. &
                    rout(i1, j2, k) == -9.99e+33 .or. &
                    rout(i2, j2, k) == -9.99e+33) then
                    routgrads(i, j, k) = -9.99e+33; cycle; end if

                rr1 = rout(i1, j1, k)*(1.-r1) + rout(i2, j1, k)*(1.-r2)
                rr2 = rout(i1, j2, k)*(1.-r1) + rout(i2, j2, k)*(1.-r2)
                routgrads(i, j, k) = rr1*(1.-r3) + rr2*(1.-r4)
                ! if(k==1)print*,i,j,rout(i,j,k),routgrads(i,j,k)

                !  if(abs(routgrads(i,j,k)).gt.1.E+06)  then
                !   if(abs(routgrads(i,j,k)).lt. 1.e+10) &
                !      print*,rr1,rr2,rout(i,j,k),routgrads(i,j,k)
                !      routgrads(i,j,k)=-9.99E+33
                ! endif
                !
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                !   write(2,0998)
                !   write(2,0999) i,j,i1,j1,glatg(j),glong(i)
                !   write(2,1000) rlat(i1,j1),rlat(i2,j1),rlat(i1,j2),rlat(i2,j2)
                !   write(2,1001) rlon(i1,j1),rlon(i2,j1),rlon(i1,j2),rlon(i2,j2)
                !   write(2,1002) rout(i1,j1,k),rout(i2,j1,k),rout(i1,j2,k),&
                !                  rout(i2,j2,k), routgrads(i,j,k)
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end do
        end do
    end do

0998 format(1x, '---------------------------------------------')
0999 format(1x, 4i3, 2f10.2)
1000 format(1x, 4f10.2)
1001 format(1x, 4f10.2)
1002 format(1x, 4f10.2, f16.3)

    !xxxxxxxxxxxxxxxxxxxxxxxxx
    !      k=1
    !      do jj=1,nyr
    !      do ii=1,nxr
    !         write(10,'(2i3,3f8.1)')ii,jj,rlat(ii,jj),rlon(ii,jj)
    !     +             ,rout(ii,jj,k)
    !      enddo
    !      enddo
    !      do jj=1,nyg
    !      do ii=1,nxg
    !         write(11,'(2i3,3f8.1)')ii,jj,glatg(jj),glong(ii)
    !     +             ,routgrads(ii,jj,k)
    !      enddo
    !      enddo
    !xxxxxxxxxxxxxxxxxxxxxxxxx
    return
end subroutine proj_rams_to_grads

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine ge_to_xy(polelat, polelon, xlon, xlat, x, y)

    parameter(rt=6367000.00)
    p = 3.14159265360/180.00

    !       transformacao horizontal:
    b = 1.0 + sin(p*xlat)*sin(p*polelat) + &
        cos(p*xlat)*cos(p*polelat)*cos(p*(xlon - polelon))

    f = 2.00*rt/b

    y = f*(cos(p*polelat)*sin(p*xlat) - &
           sin(p*polelat)*cos(p*xlat)*cos(p*(xlon - polelon)))

    x = f*(cos(p*xlat)*sin(p*(xlon - polelon)))

    return
end subroutine ge_to_xy

!---------------------------------------------------------------------

subroutine geo_grid(nx, ny, rlat, rlon, dep_glon1, dep_glon2, &
                    dep_glat1, dep_glat2, &
                    rlatmin, rlatmax, rlonmin, rlonmax, &
                    nxg, nyg, proj, glatg, glong, maxgx, maxgy)
    implicit none
    include 'rconfig.h'
    integer :: maxgx, maxgy
    real rlat(nx, ny), rlon(nx, ny)
    character*(*) proj
    real glatg(maxgy), glong(maxgx)
    integer nx, ny, x, y, n, i, j, nxg, nyg
    real dep_glon1, dep_glon2, dep_glat1, dep_glat2, rlon1, rlat1, dlon_min, dlat_min
    real rlatmin, rlatmax, rlonmin, rlonmax, xx

    dep_glon1 = rlon(1, 1)
    dep_glon2 = rlon(nx, 1)
    ! RMF TEST
    if (dep_glon1 > dep_glon2) then
        do y = 1, ny
            do x = 1, nx
                rlon(x, y) = mod((rlon(x, y) + 360.), 360.)
            end do
        end do
    end if

    if (dep_glat1 > dep_glat2) then
        do y = 1, ny
            do x = 1, nx
                rlat(x, y) = mod((rlat(x, y) + 180.), 180.)
            end do
        end do
    end if
    ! RMF TEST

    do n = 1, ny
        if (rlon(1, n) > dep_glon1) dep_glon1 = rlon(1, n)
        if (rlon(nx, n) < dep_glon2) dep_glon2 = rlon(nx, n)
    end do
    dep_glon2 = (dep_glon2 - dep_glon1)/(nx - 1)

    dep_glat1 = rlat(1, 1)
    dep_glat2 = rlat(1, ny)
    do n = 1, nx
        if (rlat(n, 1) > dep_glat1) dep_glat1 = rlat(n, 1)
        if (rlat(n, ny) < dep_glat2) dep_glat2 = rlat(n, ny)
    end do
    dep_glat2 = (dep_glat2 - dep_glat1)/(ny - 1)

!---versao 1
    !10/08/98
    x = 0.
    xx = 0.
    do n = 1, ny
        x = x + rlon(1, n)
        xx = xx + (rlon(nx, n) - rlon(1, n))/(nx - 1)
    end do
    dep_glon1 = x/ny
    dep_glon2 = xx/ny

    x = 0.
    xx = 0.
    do n = 1, nx
        x = x + rlat(n, 1)
        xx = xx + (rlat(n, ny) - rlat(n, 1))/(ny - 1)
    end do
    dep_glat1 = x/nx
    dep_glat2 = xx/nx
    go to 222
!- versao 2
!-2007  grade telescopica------------------------
    dlon_min = 1.e10
    dlat_min = 1.e10
!print*,nx,ny
    do i = 2, nx; do j = 2, ny
            dlon_min = min(dlon_min, abs(rlon(i, j) - rlon(i - 1, j)))
            dlat_min = min(dlat_min, abs(rlat(i, j) - rlat(i, j - 1)))
!print*,i,j,rlat(i,j),rlat(i,j-1)
        end do; end do
    dep_glat2 = dlat_min
    dep_glon2 = dlon_min
!print*,dlon_min,dlat_min
!pause
!stop 333
!-2007------------------------------------------
222 continue
    if (proj /= 'YES' .and. proj /= 'yes') then
        nxg = nx
        nyg = ny

    else

        !...... Grade para o GRADS:

        rlatmin = rlat(1, 1)
        rlatmax = rlat(1, 1)
        rlonmin = rlon(1, 1)
        rlonmax = rlon(1, 1)
        do i = 1, nx
            do j = 1, ny
                rlatmin = min(rlatmin, rlat(i, j))
                rlatmax = max(rlatmax, rlat(i, j))
                rlonmin = min(rlonmin, rlon(i, j))
                rlonmax = max(rlonmax, rlon(i, j))
            end do
        end do

        !...... Definicao da grade do GRADS
        !
        ! Para testar dependencia com a resolucao da grade
        !
        !            dep_glon2=0.5*dep_glon2
        !            dep_glat2=0.5*dep_glat2

        nxg = int((rlonmax - rlonmin)/dep_glon2 + 0.5) - 1
        nyg = int((rlatmax - rlatmin)/dep_glat2 + 0.5) - 1
        rlon1 = rlonmin - (nxg - nx - 1)*dep_glon2
        rlat1 = rlatmin - (nyg - ny - 1)*dep_glat2

        !            rlon1=rlonmin-2*dep_glon2
        !            rlat1=rlatmin-2*dep_glat2

        rlon1 = rlonmin - dep_glon2
        rlat1 = rlatmin - dep_glat2
        dep_glat1 = rlat1
        dep_glon1 = rlon1

        !print*,rlonmin,rlonmax,rlatmin,rlatmax
        !print*,nxg,nyg,dep_glon2,dep_glat2,dep_glat1,dep_glon1
        !pause     !            stop
    end if

    !        Define        grade do GRADS

    do i = 1, nxg
        glong(i) = dep_glon1 + float(i - 1)*dep_glon2
        !        print*,' i lon=',i,glong(i)
    end do

    do j = 1, nyg
        glatg(j) = dep_glat1 + float(j - 1)*dep_glat2
        !        print*,' j lat=',j,glatg(j)
    end do

    return
end subroutine geo_grid

!---------------------------------------------------------------------
subroutine rout_to_routgrads(nxyz, rinp, rout)
    dimension rinp(nxyz), rout(nxyz)
    do i = 1, nxyz
        rout(i) = rinp(i)
    end do
    return
end subroutine rout_to_routgrads

!---------------------------------------------------------------------

