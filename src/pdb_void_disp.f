      program pdb_ASA_void 

      implicit		real*4 (a-h,o-z), integer (i-n)
      integer		idloop, iyn15v, iyn15h, numatom, maxatom,
     &			numsteps, i, j, iatm
      parameter		(maxatom=100000)
      double precision	sitime, sec, et, ek, ep, temp, rmsf, rmsd
      real*4            sitime4, sec4, et4, ek4, ep4,
     +                  temp4, rmsf4, rmsd4
      real*8            cord8(3)
      real*4            cord4(maxatom,3)
      real*4            cord(maxatom,3)
      character*80      inppdb, inpfil, outasc
      character*80      line
      character*32      OFILE
      character         sord*1
      character         verbose*1
      character*3       type(maxatom)
      integer  no_center_atom(100)
      real*4   ave_contact4(100)   
      real*4   ave_contact6(100)   
      real*4   ave_contact8(100)   
      integer no_contact4(100)
      integer no_contact6(100)
      integer no_contact8(100)
      real*4  RPROBE
      dimension sph_area(maxatom),radius(maxatom)
      dimension vdw_expose_1(maxatom),vdw_expose_2(maxatom)

      character*4 atomtype(maxatom)
      character*3 sol
      dimension charge(maxatom)
      dimension qq(maxatom)
      dimension co(maxatom,3)
      dimension nbond(maxatom,2)
      dimension orig_asa(maxatom),after_asa(maxatom),diff_asa(maxatom)
      real*4 EU
      character*80 file_name_list(10000)
      integer ident(maxatom)
      dimension void(maxatom,3),void_new(maxatom,3)
      logical simple,switch_vol,volume_calc

c Read options
c
      write(*,*) 'input presto-PDB file= '
      read(5,'(a80)') inppdb
      write(*,*) 'input presto-binary file= '
      read(5,'(a80)') inpfil
      write(*,*) 'output pdb-ascii file= '
      read(5,'(a80)') outasc
      write(*,*) 'structure number= :"-1" means pdb only'
      read(5,*) initstr
      write(*,*) 'Single (s) or Double (d) precision= '
      read(5,'(a1)') sord
      if( sord.eq.'S') then
       sord='s'
      else if( sord.eq.'D') then
       sord='d'
      end if
      write(*,*) 'Verbose mode ? Y/[N] '
      read(5,'(a1)') verbose
      if( verbose .eq. 'Y') then
         verbose='y'
      else if (verbose .eq. 'y') then
         verbose='y'
      endif

      write(*,*) 'No of centers'
      read(*,*) no_center
      do i = 1, no_center
         read(*,*) no_center_atom(i)
      enddo

c     write(*,*) 'structure sampling step= '
c     read(5,*) no_step

      OFILE = 'void.pdb'

      no_step = 100
      eps     = 1.0e-4
      dx      = 1.0 ! shinji's note: mesh grid width [Angstrom] 
      simple  = .true.
      void_display_scale = 0.35 ! void size? [Angstrom]
      volume_calc = .true.

      intot = 0

      no_contact4(:)=0
      no_contact6(:)=0
      no_contact8(:)=0
      do i = 1, 100       
        ave_contact4(i) = 0.0                                  
        ave_contact6(i) = 0.0                                  
        ave_contact8(i) = 0.0                                  
      enddo
      do i = 1, maxatom
        type(i) = '   '
      enddo

      call make_name(file_name_list)
      do i= 1, 10
       write(*,*) file_name_list(i)
      enddo

c---------------------------------------------
c Read template PDB
c---------------------------------------------
      !(written by shinji): 
      !    Args:
      !        inppdb
      !    Outputs: 
      !        numatom: The number of atoms
      !        type : (specifying pro, lig, na, ... which is specific in this program)
      !        atomtype : atom type e.g. H1, CA, ...
      !        co : coordinates
      call  read_pdb_file(inppdb,numatom,type,atomtype,co)

      !(written by shinji): 
      !    Args:
      !        numatom, type, atomtype, co (details are shown in the above comment)
      !    Outputs: 
      !        indent :  atom nunmbers from 1 to the end written in a pdb
      call  identical_atom(numatom,type,atomtype,co,ident)

      write(*,*) 'number of atoms= ', numatom


c DUM = 0.2 in setasa_parm
       sol = 'DUM'
       x = 1.0
       y = 1.0
       charge(:) = 0.0
       call setasa_parm(numatom,type,atomtype,charge,co,x,y,sol
     &  ,numbond,nbond ,qq,radius,sph_area)


       RPROBE=0.1
       call asa_15tab(numatom,co,RPROBE)
       call FENASA(co, RPROBE, EU, orig_asa)

c---------------------------------------------
c Read trajectory file
c---------------------------------------------

      if(initstr .le. 0 ) then

       intot=1
       cord(:,:)=co(:,:)

c DUM = 0.1 in setasa_parm
       sol = 'DUM'
       x = 1.0
       y = 1.0
       charge(:) = 0.0
       call setasa_parm(numatom,type,atomtype,charge,co,x,y,sol
     &  ,numbond,nbond ,qq,radius,sph_area)
       RPROBE = 0.1
       call asa_15tab(numatom,co,RPROBE)
      write(*,*) 'step1c = ', EU  ,initstr      
       call FENASA(co, RPROBE, EU, orig_asa)

       do  i = 1, numatom
           vdw_expose_1(i) =  orig_asa(ident(i))/sph_area(ident(i))
       enddo

c WAT = 1.4 in setasa_parm
       sol = 'WAT'
       x = 1.0
       y = 1.0
       charge(:) = 0.0
       print*, "type, atomtype", type(1:3), atomtype(1:3)
       call setasa_parm(numatom,type,atomtype,charge,co,x,y,sol
     &  ,numbond,nbond ,qq,radius,sph_area)
       RPROBE = 1.4
       call asa_15tab(numatom,co,RPROBE)
       call FENASA(co, RPROBE, EU, after_asa)


       do  i = 1, numatom
           vdw_expose_2(i) =  after_asa(ident(i))/sph_area(ident(i))
       enddo

       do i = 1, numatom
        diff_asa(i) = 4.0*3.14159*radius(i)**2
     &    * (vdw_expose_1(ident(i)) - vdw_expose_2(ident(i)))
       enddo

      call cal_contactNo(numatom, co,no_center,no_center_atom,
     &  type,no_contact4, no_contact6, no_contact8)

       ave_contact4(:)=real(no_contact4(:))
       ave_contact6(:)=real(no_contact6(:))
       ave_contact8(:)=real(no_contact8(:))

      write(*,*) 'before'
      switch_vol=.false.

      !@@@shinji checked
      print*, "--Radius--"
      do i = 1, size(radius)
        if (radius(i) /= real(0.0, kind=4)) print*, i, radius(i) 
      enddo

      call make_mesh(numatom,co,radius,type,
     &          dx, switch_vol,no_void,void)
      write(*,*) 'after '

      else !------------------

      no_out_pdb = 0

      open(unit=9, file=inpfil, status='old', form='unformatted')
100   continue

      if( sord.eq.'d' ) then
        read(9,end=200) idloop, sitime, sec, et, ek, temp,
     +                  ep, rmsf, iyn15v, iyn15h, rmsd
      intot = intot + 1
         if(verbose.eq.'y') then
              write(*,*) idloop, sitime, sec, et, ek, temp,
     +                  ep, rmsf, iyn15v, iyn15h, rmsd
         endif

        read(9,end=200) ((cord8(j), j=1,3), iatm=1,numatom)
	do iatm=1, numatom
          cord(iatm,1) = real(cord8(1))
          cord(iatm,2) = real(cord8(2))
          cord(iatm,3) = real(cord8(3))
	enddo

      else
        read(9,end=200) idloop, sitime4, sec4, et4, ek4, temp4,
     +                  ep4, rmsf4, iyn15v, iyn15h, rmsd4
      intot = intot + 1
          if(verbose.eq.'y') then
             write(*,*) idloop, sitime4, sec4, et4, ek4, temp4,
     +                  ep4, rmsf4, iyn15v, iyn15h, rmsd4
          endif

      read(9,end=200) ((cord4(iatm,j), j=1,3), iatm=1,numatom)

      do iatm=1, numatom
         cord(iatm,1) = cord4(iatm,1)
         cord(iatm,2) = cord4(iatm,2)
         cord(iatm,3) = cord4(iatm,3)
      enddo

c---------------------------------------------
c calculate distance
c---------------------------------------------
      call cal_contactNo(numatom, cord,no_center,no_center_atom,
     &  type,no_contact4, no_contact6, no_contact8)

       write(*,* ) 'CN4=',idloop,real(sitime4),(no_contact4(i),i=1, 4)        
       write(*,* ) 'CN6=',idloop,real(sitime4),(no_contact6(i),i=1, 4)        
       write(*,* ) 'CN8=',idloop,real(sitime4),(no_contact8(i),i=1, 4)        
 97    format(a4,1x,i9,2x,f8.3,'4(2x,i6)')

       do i = 1, no_center
        ave_contact4(i) = ave_contact4(i)+real(no_contact4(i))
        ave_contact6(i) = ave_contact6(i)+real(no_contact6(i))
        ave_contact8(i) = ave_contact8(i)+real(no_contact8(i))
       enddo


c---------------------------------------------
c calculate ASA      
c---------------------------------------------

c DUM = 0.1 in setasa_parm
       sol = 'DUM'
       x = 1.0
       y = 1.0
       charge(:) = 0.0
       call setasa_parm(numatom,type,atomtype,charge,co,x,y,sol
     &  ,numbond,nbond ,qq,radius,sph_area)
       RPROBE = 0.1
       call asa_15tab(numatom,cord,RPROBE)
       call FENASA(cord, RPROBE, EU, orig_asa)

       do  i = 1, numatom
           vdw_expose_1(i) =  orig_asa(i)/sph_area(i)
       enddo

c WAT = 1.4 in setasa_parm
       sol = 'WAT'
       x = 1.0
       y = 1.0
       charge(:) = 0.0
       call setasa_parm(numatom,type,atomtype,charge,co,x,y,sol
     &  ,numbond,nbond ,qq,radius,sph_area)
       RPROBE = 1.4
       call asa_15tab(numatom,cord,RPROBE)
       call FENASA(cord, RPROBE, EU, after_asa)
        

       do  i = 1, numatom
           vdw_expose_2(i) =  after_asa(i)/sph_area(i)
       enddo

       do i = 1, numatom
        diff_asa(i) = 4.0*3.14159*radius(i)**2 
     &    * (vdw_expose_1(i) - vdw_expose_2(i))
       enddo

      end if

c---------------------------------------------
c Write PDB
c---------------------------------------------
      
      write(*,*) 'aaa', intot,initstr,no_step,no_out_pdb
      if( mod(intot,no_step).eq.1) then
      no_out_pdb = no_out_pdb+1
      write(*,*) ' write pdb ',no_out_pdb

      open(unit=1, file=inppdb, status='old')
      open(unit=3,file=file_name_list(no_out_pdb),status='unknown')
      write(3,'("REMARK ",i7,"TH STRUCTURE")') no_out_pdb

      ia = 0
      icng = 0
   11 continue
        read(1,'(a80)',end=19) line
        if(line(1:5).eq.'ATOM ' .or. line(1:5).eq.'HETAT') then
          ia= ia + 1
          write(line(31:54),'(3f8.3)')  (cord(ia,j),j=1,3)
	  write(line(62:66),'(f5.2)') diff_asa(ident(ia))/10.0 
          write(3,'(a80)') line
        end if
        go to 11
   19 continue
      write(3,'("END")')
      close(1)
      close(3)
      endif


c---------------------------------------------
c Write PDB end
c---------------------------------------------

      if( intot.lt.initstr ) then
         go to 100
      else
         go to 200 
      end if

 200  continue
      close(9)
c---------------------------------------------
c Loop end 
c---------------------------------------------

      endif


       write(*,*) '-----------------------------------'
       write(*,*) '- Contact number   center 1, 2,.. -'
       write(*,*) '-----------------------------------'
       do i = 1, no_center
        ave_contact4(i) = ave_contact4(i)/real(intot)          
        ave_contact6(i) = ave_contact6(i)/real(intot)          
        ave_contact8(i) = ave_contact8(i)/real(intot)          
       enddo

       write(*,* ) 'CN4f=',idloop,sitime4,(ave_contact4(i),i=1, 4)        
       write(*,* ) 'CN6f=',idloop,sitime4,(ave_contact6(i),i=1, 4)        
       write(*,* ) 'CN8f=',idloop,sitime4,(ave_contact8(i),i=1, 4)        
 98    format(a5,1x,i8,1x,f8.3,'4(1x,f8.3)')

c---------------------------------------------
c Write PDB           
c---------------------------------------------
      open(unit=1, file=inppdb, status='old')
      open(unit=3, file=outasc, status='unknown')
      write(3,'("REMARK ",i7,"TH STRUCTURE")') initstr

      open(unit=1, file=inppdb, status='old')
      ia = 0
      icng = 0
   21 continue
        read(1,'(a80)',end=29) line
        if(line(1:5).eq.'ATOM ') then
          ia= ia + 1
	  write(line(31:54),'(3f8.3)')  (cord(ia,j),j=1,3)
ccc       write(line(56:60),'(f5.2)') diff_asa(ia)/10.0
	  write(line(62:66),'(f5.2)') diff_asa(ident(ia))/1.0 
          write(3,'(a80)') line
        end if
        go to 21
   29 continue
      write(3,'("END")')
      close(1)
      close(3)


      if (simple) then

       Rth  = 2.1
       print*, "Rth =", Rth
       Rth2 = Rth*Rth
       n_th = int(void_display_scale*(2.0*Rth/dx)**3) !n threshold
       no_void_new = 0
       
       !call writepdb2(no_void,void,OFILE)
       !, which generate all lattice points saved in void.pdb.
       !Shinji just used it to get my head around how it works.
       !stop

       do i = 1, no_void
        n_kosuu = 0

        do j = 1, no_void
          if(abs(void(j,1)-void(i,1)).gt.Rth) cycle
          if(abs(void(j,2)-void(i,2)).gt.Rth) cycle
          if(abs(void(j,3)-void(i,3)).gt.Rth) cycle
          dd=(void(i,1)-void(j,1))**2 +
     &       (void(i,2)-void(j,2))**2 + (void(i,3)-void(j,3))**2

          if(dd <= Rth2) then
            n_kosuu = n_kosuu+1 !the number of mesh points (voids) within Rth**2 
          endif
        enddo
        
        if(n_kosuu > n_th) then
           !print*, "N mesh points within a point = ", n_kosuu
           no_void_new = no_void_new + 1
           void_new(no_void_new,1:3) = void(i,1:3)
           !print*, "Counting No. of voids...", no_void_new
        endif

       enddo

       no_void = no_void_new
       void(1:no_void,1:3) = void_new(1:no_void,1:3)

      endif 

      call writepdb2(no_void,void,OFILE)
    
      if( volume_calc) then
        dx = 1.0 !0.2
        print*, "dx is set to", dx, "for volume calculation."
        switch_vol = .true. 
        call make_mesh(numatom,co,radius,type,
     &                 dx,switch_vol, no_void, void)
      endif

      stop
      end

c---------------------------------------------
c calculate distance
c---------------------------------------------
      subroutine cal_contactNo(numatom, cord,no_center,no_center_atom,
     &  type,no_contact4, no_contact6, no_contact8)
      implicit real*4(a-h,o-z)
      parameter		(maxatom=100000)
      real*4            cord(maxatom,3)
      real*4  R_thre4, R_thre6, R_thre8, x, y, z, rr
      character*3       type(maxatom)
      integer i, numatom, no_center
      integer no_center_atom(100)
      integer no_contact4(100)
      integer no_contact6(100)
      integer no_contact8(100)

      R_thre4 = 4.0
      R_thre6 = 6.0
      R_thre8 = 8.0

      do i = 1, no_center
        no_contact4(i) = 0    
        no_contact6(i) = 0    
        no_contact8(i) = 0    
      enddo

      do i = 1, no_center
       x = cord(no_center_atom(i),1)      
       y = cord(no_center_atom(i),2)      
       z = cord(no_center_atom(i),3)      
      do j = 1, numatom
       if( type(j) .eq. 'wao'.or.type(j).eq.'nat'.or.
     &     type(j) .eq. 'clo' ) then
       rr=(cord(j,1)-x)**2+(cord(j,2)-y)**2+(cord(j,3)-z)**2
       rr = sqrt(rr)

       if ( rr .lt. R_thre4 ) then
         no_contact4(i) = no_contact4(i) + 1
       endif
       if ( rr .lt. R_thre6 ) then
         no_contact6(i) = no_contact6(i) + 1
       endif
       if ( rr .lt. R_thre8 ) then
         no_contact8(i) = no_contact8(i) + 1
       endif
       endif
      enddo
      enddo

      return
      end
c---------------------------------------------
c Read template PDB
c---------------------------------------------
      subroutine identical_atom(numatom,type,atomtype,co,ident)
      implicit		real*4 (a-h,o-z), integer (i-n)
      integer		i, numatom 
      parameter		(maxatom=100000)
      real*4            co(maxatom,3)
      character*80      inppdb 
      character*80      line
      character*1       at1            
      character*2       at2             
      character*3       type(maxatom)
      character*4       atomtype(maxatom),dum
      integer           ident(maxatom)

      do i =1, numatom

       if(type(i) .ne. 'pro') then
        ident(i) = i ! indent stores non-proten atom numbers.
        cycle
       endif

       dum = atomtype(i)
       at2 = dum(1:2)
       at1 = dum(1:1)

       if(at2 .eq. ' H' .or. at2 .eq. 'H ') at1 = 'H' 
       !This line would be rewritten by just using trim()

       !indent stores non-hydrogen atom numbers in a protein here unlike the above. 
       if(at1 .ne. 'H' .or. i .le. 1 ) then
         ident(i) = i
       else 
         ident(i) = ident(i-1) !@@@why this written?
       endif

      enddo 

      return
      end
c---------------------------------------------
c---------------------------------------------
c Read template PDB
c---------------------------------------------
      subroutine read_pdb_file(inppdb,numatom,type,atomtype,co)
      implicit		real*4 (a-h,o-z), integer (i-n)
      integer		ia, numatom 
      parameter		(maxatom=100000)
      real*4            co(maxatom,3)
      character*80      inppdb 
      character*80      line
      character*3       type(maxatom)
      character*4       atomtype(maxatom),dum

c Read options
c
c
c
      open(unit=1, file=inppdb, status='old')
      ia = 0
    1 continue
        read(1,'(a80)',end=9) line


        if(line(1:6).eq.'ATOM  ' .or. line(1:6).eq.'HETATM') then
          ia= ia + 1
             type(ia) = 'pro'

          if( line(1:6).eq.'HETATM') then

            if(line(14:15).eq.'NA'.or.line(15:16).eq.'NA') then
              type(ia) = 'nat'
            endif

            if(line(14:15).eq.'Na'.or.line(15:16).eq.'Na') then
              type(ia) = 'nat'
            endif

            if(line(14:15).eq.'CL'.or.line(15:16).eq.'CL') then
              type(ia) = 'clo'
            endif

            if(line(14:15).eq.'Cl'.or.line(15:16).eq.'Cl') then
              type(ia) = 'clo'
            endif

          endif

          if(line(1:6).eq.'HETATM') then

            if(line(18:18).eq.'L') then
              type(ia) = 'lig'

            else
              type(ia) = 'het'
            endif

          endif

          if(line(18:20).eq.'WAT'.or.line(18:20).eq.'TIP'
     &       .or. line(18:20).eq.'HOH') then

            if(line(14:14) .eq. 'O') then
              type(ia) = 'wao'

            else
              type(ia) = 'wah'

            endif

          endif

	read(line(31:54),'(3f8.3)')  (co(ia,j),j=1,3)
	read(line(13:16),'(a)')  dum             
         atomtype(ia) = trim(dum)

        end if
        go to 1
    9 continue
      close(1)
      numatom = ia

      return
      end
c-----------------------------------
       subroutine make_name(file_name_list)
       integer i
       character*3 file_suf
       character*80 line,file_name_list(10000)
 
       file_suf='pdb'

       max_conf=  9000
      do i=1, max_conf
        write(line,'(a4)') 'asad'

      if(i.le.9)
     &  write(line(5:9),'(i1,a1,a3)')i,'.',file_suf

      if(i.ge.10.and.i.le.99)
     &  write(line(5:10),'(i2,a1,a3)')i,'.',file_suf

      if(i.ge.100.and.i.le.999)
     &  write(line(5:11),'(i3,a1,a3)')i,'.',file_suf

      if(i.ge.1000.and.i.le.9999)
     &  write(line(5:12),'(i4,a1,a3)')i,'.',file_suf

        file_name_list(i)=line
      enddo

      return
      end
c-----------------------------------

      subroutine make_mesh(numatom_in,cord_in,radius_in,type,
     &        dx, switch_vol,no_void,void)
      implicit          real*4 (a-h,o-z), integer (i-n)
      parameter         (maxatom=100000)
      character*3       type(maxatom)
      real*4  RPROBE
      dimension sph_area(maxatom),radius(maxatom)
      character*3 sol
      dimension cord_in(maxatom,3),radius_in(maxatom)
      dimension cord(maxatom,3)

      dimension co(1000,3),rad(1000),dist(1000),no(1000)
      dimension coo(1000,3),rado(1000),rr(1000)
      dimension xyz_min(3),xyz_max(3),xx(3),yy(3)
      logical point_is_void
      logical on_surface   
      logical probe_exist  
      logical switch_vol
      dimension void(maxatom,3),vec(3),vec2(3),vec3(3)
      dimension vec4(3),cp(3)

      !Save no. of atoms, coordinates, and radius only for protein atoms.
      numatom=0
      do i = 1, numatom_in
        if(type(i) == 'pro') then
          numatom           = numatom+1
          cord(numatom,1:3) = cord_in(i,1:3)
          radius(numatom)   = radius_in(i)
        endif
      enddo

      Rth1    = 7.0 !cut-off for neighbouring grid points from a grid point. ??
      r_probe = 1.4

      xyz_min(:) = 999.9
      xyz_max(:) = -999.9

      if(switch_vol) then
        write(*,*) "#=================================="
        write(*,*) "# calculate void_volume            "
        write(*,*) "#=================================="
      endif

      do i = 1, numatom
       do j = 1, 3
         if(cord(i,j) < xyz_min(j)) xyz_min(j) = cord(i,j)
         if(cord(i,j) > xyz_max(j)) xyz_max(j) = cord(i,j)
       enddo
      enddo

      !making min and max side of a grid points (generated later) with probe-based margin 
      xyz_min(1:3) = xyz_min(1:3) - 1.0*r_probe
      xyz_max(1:3) = xyz_max(1:3) + 1.0*r_probe

      !number of grid points for each axis
      nx_max = (xyz_max(1)-xyz_min(1))/dx
      ny_max = (xyz_max(2)-xyz_min(2))/dx
      nz_max = (xyz_max(3)-xyz_min(3))/dx

      write(*,*) '## ',numatom,nx_max,ny_max,nz_max
      write(*,*) '###min ',xyz_min(1),xyz_min(2),xyz_min(3)
      write(*,*) '###max ',xyz_max(1),xyz_max(2),xyz_max(3)

      no_void=0
      n1_tag_old= 0

      !start to three loop for grid points 
      do n1 = 1, nx_max
        xx(1) = xyz_min(1)+dx*real(n1)
      do n2 = 1, ny_max
        xx(2) = xyz_min(2)+dx*real(n2)
      do n3 = 1, nz_max
        xx(3) = xyz_min(3)+dx*real(n3)

      point_is_void=.false.

      !save coordinates of grid points within Rth1
      noatom=0 !the number of grid points within Rth1
      do i = 1, numatom
        if( abs(cord(i,1) - xx(1)) > Rth1) cycle
        if( abs(cord(i,2) - xx(2)) > Rth1) cycle
        if( abs(cord(i,3) - xx(3)) > Rth1) cycle
        noatom = noatom + 1
        coo(noatom,1) = cord(i,1)
        coo(noatom,2) = cord(i,2)
        coo(noatom,3) = cord(i,3)
        rado(noatom)  = radius(i) !which comes from ASA.f: RADASA(i) = 1.8
      enddo

      if(noatom <= 1) cycle

      !calculate squred distance  of the saved coordinates above
      dist(:) = 0.0 !initialised by shinji
      do i = 1, noatom
        no(i) = i
        !note xx(1:3) is a grid point along each axis
        dist(i) = (coo(i,1)-xx(1))**2 + (coo(i,2)-xx(2))**2
     &          + (coo(i,3)-xx(3))**2
      enddo

      !sorting distance as ascending order like
      !-
      !--
      !-------
      do j = 1, noatom
        do i = 1, noatom - 1
          if(dist(i) > dist(i+1)) then
            x         = dist(i)
            k         = no(i)
            dist(i)   = dist(i+1) !The bigger one replaced w. smaller dist
            no(i)     = no(i+1)   !index as well
            dist(i+1) = x !the smaller one is replaced w. the bigger one.
            no(i+1)   = k 
          endif
        enddo
      enddo
      !check if it is really ascending order
      !do i = 1, noatom
      !   print*, n1, n2, n3, dist(i)
      !enddo

      !sorting the coodinates of grid points according to the sorting 
      !of distance done above. And radisu (rado) as well
      do i = 1, noatom
        co(i,1:3) = coo(no(i),1:3)
        rad(i)    = rado(no(i))
      enddo

      !calculate distance btwn the first coordinates of grid and those of a grid point (xx) 
      !Note that this distance is the SMALLEST one from dist saved above
      d1=sqrt((co(1,1)-xx(1))**2+(co(1,2)-xx(2))**2
     &      + (co(1,3)-xx(3))**2 )

      !btwn the 2nd coordinates of grid and those of a bin (xx) 
      d2=sqrt((co(2,1)-xx(1))**2+(co(2,2)-xx(2))**2
     &      + (co(2,3)-xx(3))**2 )

      !a grid point is not a void if within rad(ius) (1.8 A). 
      !1.8 is, perhaps, the lower limit of a vdw sphere.
      if(d1 < rad(1) .or. d2 < rad(2)) then
        point_is_void=.false.
        cycle
      endif

      !r_probe = 1.4 A for water. 2*r_probe is diameter
      !a point is not a void if the 1st&2nd smallest distance from surroundings > (vdw lower lim) + (wat. diameter)
      !In other words, as contraposition,
      !if a point is a void, the distance <= (vdw lower lim) + (wat.diameter).
      !Meaning 
      if(d1 > rad(1)+2.0*r_probe .and. d2 > rad(2)+2.0*r_probe) then
        point_is_void=.false.
        cycle
      endif
      !Note that the false volds are not processed in the subsequent lines 

      !true voids are processed by the susequent lines
      !---??
      vec(1:3) = xx(1:3) - co(1,1:3) !vector from a grid point to the 1st surrounding grid point within Rth1
      xnorm    = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
      vec(1:3) = vec(1:3)/xnorm !unitfied
      yy(1:3)  = co(1,1:3) + vec(1:3)*(rad(1) + r_probe) !it directs co 
      do i = 1, noatom
        !rr(i) is a radius (vdw + wat) of a sphere 
        rr(i) = sqrt( (co(i,1)-yy(1))**2 + (co(i,2)-yy(2))**2
     &              + (co(i,3)-yy(3))**2 )
      enddo

      !---??
      on_surface = .true.
      do i = 2, noatom !why from 2?
        if(rr(i) <= rad(i) + 1.0*r_probe) on_surface = .false. 
      enddo
      if(on_surface) cycle !because grids on surface should be omitted from the def of void. 

      !---??
      R =sqrt((co(2,1)-co(1,1))**2+(co(2,2)-co(1,2))**2
     &      + (co(2,3)-co(1,3))**2 )

      Ar = rad(1)+r_probe
      Br = rad(2)+r_probe

      ct1  = acos((d1*d1+R*R-d2*d2)/(2.0*d1*R))
      ct12 = acos((Ar*Ar+R*R-Br*Br)/(2.0*Ar*R))
      ct2  = ct12 - ct1
      d3   = sqrt(Ar*Ar+d1*d1-2.0*Ar*d1*cos(ct2))


      if(d3 > r_probe) then
        point_is_void=.true. 
      else
        vec(1:3)  = xx(1:3)-co(1,1:3)
        vec2(1:3) = co(2,1:3)-co(1,1:3)
        x_norm    = sqrt(vec2(1)**2+vec2(2)**2+vec2(3)**2)
        vec2(1:3) = vec2(1:3)/x_norm
        vec3(1:3) = d1*cos(ct1)*vec2(1:3)
        vec4(1:3) = vec(1:3)-vec3(1:3)
        x_norm    = sqrt(vec4(1)**2+vec4(2)**2+vec4(3)**2)
        vec4(1:3) = vec4(1:3)/x_norm

        cp(1:3)   = co(1,1:3)+Ar*cos(ct12)*vec2(1:3)
     &   + Ar*sin(ct12)*vec4(1:3)

        do i = 1, noatom
        rr(i) = sqrt( (co(i,1)-cp(1))**2+(co(i,2)-cp(2))**2
     &          +(co(i,3)-cp(3))**2 )
        enddo
        probe_exist=.true.

        if(noatom >= 3) then
         do i = 3, noatom
          if(rr(i) <= rad(i)+1.0*r_probe) then
            probe_exist=.false. 
            cycle
          endif
         enddo
        endif

        if(.not.probe_exist) point_is_void=.true.

      endif

      if(point_is_void) then
        no_void = no_void + 1

        if(.not.switch_vol) then
          void(no_void,1:3) = xx(1:3)
        endif

        n1_tag = int(real(n1)/real(nx_max)*100.0)
        if(n1_tag_old .ne. n1_tag) then
          write(*,*) '## ', n1_tag,
     &               ' %-done No_of_void_points=',no_void
        endif
        n1_tag_old=n1_tag
      endif
     
      enddo
      enddo
      enddo
! --------------

      if(switch_vol) then
        vol = real(no_void)*dx**3
        write(*,*) '================================='
        write(*,*) '@vol void_volume= ',vol,' A^3/mol'
        write(*,*) '================================='
      endif

      return
      end

c------------------------------------------------------
      subroutine writepdb2(numatom,co,OFILE)
      implicit real*4(a-h,o-z)
      character*7 dum
      character*2 ato2
      character*32 OFILE
      character*80 line,line2
      parameter         (maxatom=100000)
      dimension co(maxatom, 3)
      dimension ato2(maxatom)

      open(11,file=OFILE,status='unknown',form='formatted')

      dum =  'ATOM   '
      do i=1,numatom
       ato2(i)='CH'
       j = mod(i,10000)
       write(11,95)dum,j,ato2(i),'VID',
     &       ' Z ', 999,co(i,1),co(i,2),co(i,3)
      end do
 95   format(a7,i4,2x,a2,2x,a3,a3,i3,4x,3f8.3)
        write(11,96) 'TER'
 96       format(a3)

      close(11)

      write(*,*)'#============================='
      write(*,*)'# void_points were written already'
      write(*,*)'# output_file name= ',OFILE        
      write(*,*)'#============================='

      return
      end subroutine writepdb2

