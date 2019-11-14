c----------------------------------------
      subroutine setasa_parm(numatom,type,atomtype,charge,co,x,y,sol
     &  ,numbond,nbond ,qn,radius,sph_area)
      implicit real*4(a-h,o-z)
      parameter(maxat=100000 )
      character*3 type(maxat)
      character*4 atomtype(maxat), atom_name
      character*1 at1
      character*2 at2
      character*3 sol               
      dimension charge(maxat)
      dimension qq(maxat)
      dimension qn(maxat)
      dimension co(maxat,3)
      dimension qg(maxat)
      dimension nbond(maxat,2)
      dimension sph_area(maxat),radius(maxat)

C---------------------------------------------------------------------C
C     DIMENSION PARAMETER FOR COMASA
C---------------------------------------------------------------------C

      INTEGER MAXNOH,MAXASA,NONONH

      PARAMETER (MAXNOH=50000,MAXASA=350)

C---------------------------------------------------------------------C
C     COMMON AREA
C---------------------------------------------------------------------C

      COMMON/ASAPO1/ NONONH
      COMMON/ASAPO3/NONHYD(MAXNOH),INASA(MAXNOH),IASATB(MAXNOH,MAXASA)
      COMMON/ASAPO4/ ASPCOF(MAXNOH),RADASA(MAXNOH)

C---------------------------------------------------------------------C
      INTEGER NONHYD,IASATB,INASA
      REAL*4    ASPCOF,RADASA
C---------------------------------------------------------------------C
      
c  x : hydrophobic atom
c  y : hydrophilic atom

      NONONH = 0

      do i = 1, numatom
        qq(i) = charge(i)
      enddo
      radius(:)=0.0
      sph_area(:) = 0.0
c
      do i = 1, numatom

        if( type(i) .ne. 'pro') cycle 

        atom_name = atomtype(i)
        at1 = atom_name(1:1)
        at2 = atom_name(1:2)

        !@@@Shinji commented out:
        !older version (if you wanna push back, remove commented out):
        !if( at2 .eq. ' H') at1 = 'H'           
        !if( at2 .eq. 'H ') at1 = 'H'           

        !@@@Shinji version
        if( at2 .eq. ' H') then
          at1 = 'H'           
        elseif( at2 .eq. 'H ') then 
          at1 = 'H'           
        else
          at1 = at2(2:2)
        endif 
        if (at1 == " ") then 
            print*, "The variable at1 in ASA.f was blank. 
     &       i.e. at1 = ", at1 
            stop "STOPPED" 
        endif 
        !@@@

        print*, "i, at1, at2, atmtyp: ", i, at1, at2, atomtype(i)

        !ignore hydrogen atoms in the calc.
        if(at1 .ne. 'H') then
          NONONH = NONONH + 1
          NONHYD(NONONH) = i
          ASPCOF(i) = x
          RADASA(i) = 1.8

          if(at1 .eq. 'C') then              
            ASPCOF(i) = x 
            RADASA(i) = 1.9
          endif 
          if(at1 .eq. 'N') then           
            ASPCOF(i) = y 
            RADASA(i) = 1.8
          endif 
          if(at1 .eq. 'O') then              
            ASPCOF(i) = y 
            RADASA(i) = 1.7
          endif 
          if(at1 .eq. 'S') then            
            ASPCOF(i) = x 
            RADASA(i) = 2.1
          endif 
          if(at1 .eq. 'P') then           
            ASPCOF(i) = y 
            RADASA(i) = 2.1
          endif 
          if(at1 .eq. 'F') then             
            ASPCOF(i) = y 
            RADASA(i) = 1.2 
          endif 
 
 
          radius(i) = RADASA(i)
 
          if(sol .eq. 'WAT') then
            offset_w = 1.4
            RADASA(i) = RADASA(i) + offset_w
            ASPCOF(i) = x  ! 10.0   
          endif
 
          if(sol .eq. 'OCT') then
            offset_o = 0.1 
            RADASA(i) = RADASA(i) + offset_o
            ASPCOF(i) = y  !  0.0  
          endif
       
          if(sol.eq.'DUM'.or.sol.eq.'DES') then
            ASPCOF(i) =  1.0
          endif
 
          sph_area(i) = 4.0*3.14159*RADASA(i)**2
        endif
      enddo


      return
      end subroutine setasa_parm
c----------------------------------------
      subroutine setGB_parm(numatom,atomtype,
     &  numbond,nbond, vdW, scale)
      implicit real*4(a-h,o-z)
      parameter(maxat=100000 )
      character*6 atomtype(maxat),atom1, atom2
      character*1 at1, at2
      dimension nbond(maxat,2)
      dimension vdW(maxat), scale(maxat)

C---------------------------------------------------------------------C

      do i = 1, numatom
           scale(i) = 0.90
           vdW(i) = 2.0

         if(atomtype(i) .eq. 'H') then              
           scale(i) = 0.85
           vdW(i) = 1.25
         endif 

         if(atomtype(i) .eq. 'C.1') then              
           scale(i) = 0.72
           vdW(i) = 1.9
         endif 
         if(atomtype(i) .eq. 'C.2') then           
           scale(i) = 0.72
           vdW(i) = 1.825
         endif 
         if(atomtype(i) .eq. 'C.3') then           
           scale(i) = 0.72
           vdW(i) = 1.875
         endif 
         if(atomtype(i) .eq. 'C.ar') then           
           scale(i) = 0.72
           vdW(i) = 1.875
         endif 
         if(atomtype(i) .eq. 'C.cat') then           
           scale(i) = 0.72
           vdW(i) = 1.875
         endif 
         if(atomtype(i) .eq. 'N.1') then           
           scale(i) = 0.79
           vdW(i) = 1.6
         endif 
         if(atomtype(i) .eq. 'N.2') then           
           scale(i) = 0.79
           vdW(i) = 1.7063
         endif 
         if(atomtype(i) .eq. 'N.3') then           
           scale(i) = 0.79
           vdW(i) = 1.7063
         endif 
         if(atomtype(i) .eq. 'N.4') then           
           scale(i) = 0.79
           vdW(i) = 1.625
         endif 
         if(atomtype(i) .eq. 'N.am') then            
           scale(i) = y 
           vdW(i) = 1.7063
         endif 
         if(atomtype(i) .eq. 'N.ar') then            
           scale(i) = 0.79
           vdW(i) = 1.7063
         endif 
         if(atomtype(i) .eq. 'N.pl3') then           
           scale(i) = y 
           vdW(i) = 1.7063
         endif 
         if(atomtype(i) .eq. 'O.2') then            
           scale(i) = 0.85
           vdW(i) = 1.535
         endif 
         if(atomtype(i) .eq. 'O.3') then           
           scale(i) = 0.85
           vdW(i) = 1.535
         endif 
         if(atomtype(i) .eq. 'O.co2') then           
           scale(i) = 0.85
           vdW(i) = 1.535
         endif 
         if(atomtype(i) .eq. 'S.2') then            
           scale(i) = 0.96
           vdW(i) = 1.775
         endif 
         if(atomtype(i) .eq. 'S.3') then           
           scale(i) = 0.96
           vdW(i) = 1.775
         endif 
         if(atomtype(i) .eq. 'S.o') then           
           scale(i) = 0.96
           vdW(i) = 1.775
         endif 
         if(atomtype(i).eq.'S.o2'.or.atomtype(i).eq.'S.O2')then           
           scale(i) = 0.96
           vdW(i) = 1.775
         endif 
         if(atomtype(i) .eq. 'P.2') then           
           scale(i) = 0.86
           vdW(i) = 1.87
         endif 
         if(atomtype(i) .eq. 'P.3') then           
           scale(i) = 0.86
           vdW(i) = 1.87
         endif 
         if(atomtype(i) .eq. 'F') then             
           vdW(i) = 1.47
         endif 
         if(atomtype(i) .eq. 'Cl') then             
           vdW(i) = 1.735
         endif 
         if(atomtype(i) .eq. 'Br') then            
           vdW(i) = 1.9
         endif 
         if(atomtype(i) .eq. 'I') then             
           vdW(i) = 2.1
         endif 

      enddo

      do i = 1, numbond
         n1=nbond(i,1)
         n2=nbond(i,2)
         atom2=atomtype(n2)
        write(atom2(1:1),'(a1)')  at2    

         if(atomtype(n1) .eq. 'H') then
           if(at2 .eq. 'N') then
             vdW(n1) = 1.15
           endif
           if(at2 .eq. 'O') then
             vdW(n1) = 1.05
           endif
         endif

        atom1=atomtype(n1)
        write(atom1(1:1),'(a1)')  at1    

         if(atomtype(n2) .eq. 'H') then
           if(at1 .eq. 'N') then
             vdW(n2) = 1.15
           endif
           if(at1 .eq. 'O') then
             vdW(n2) = 1.05
           endif
         endif
      enddo

      do i=1,numatom 
       vdW(i)=vdW(i)+1.4
      enddo

      return
      end subroutine setGB_parm
c----------------------------------------
c----------------------------------------
      subroutine shuffle(numatom,charge)
      implicit real*4(a-h,o-z)
      parameter(maxat=100000 )
      PARAMETER (MAXNOH=50000,MAXASA=350)
      dimension charge(maxat)
C---------------------------------------------------------------------C
      COMMON/ASAPO1/ NONONH
      COMMON/ASAPO3/NONHYD(MAXNOH),INASA(MAXNOH),IASATB(MAXNOH,MAXASA)
C     COMMON AREA
C---------------------------------------------------------------------C
      COMMON/ASAPO4/ ASPCOF(MAXNOH),RADASA(MAXNOH)

      do ncyc = 1, NONONH *100
        do i = 1, NONONH
        n1 = NONHYD(i)
        call random_number(rnd)
        j2 = int(rnd * real(NONONH))
        if (j2. le. 0 ) j2 = 1
        if (j2. ge. NONONH  ) j2 = NONONH 
           n2 = NONHYD(j2)
           x = ASPCOF(n1)  
           ASPCOF(n1) = ASPCOF(n2)
           ASPCOF(n2) = x  
        enddo
      enddo

      return
      end subroutine shuffle
c----------------------------------------
      subroutine shuffle_atom(numatom,atomtype,charge)
      implicit real*4(a-h,o-z)
      parameter(maxat=100000 )
      character*6 atomtype(maxat)
      character*6 x
      dimension charge(maxat)
      real*4 ch

      PARAMETER (MAXNOH=50000,MAXASA=350)
C---------------------------------------------------------------------C
      COMMON/ASAPO1/ NONONH
      COMMON/ASAPO3/NONHYD(MAXNOH),INASA(MAXNOH),IASATB(MAXNOH,MAXASA)
C---------------------------------------------------------------------C
      COMMON/ASAPO4/ ASPCOF(MAXNOH),RADASA(MAXNOH)

      do ncyc = 1, NONONH *100
        do i = 1, NONONH
        n1 = NONHYD(i)
        call random_number(rnd)
        j2 = int(rnd * real(NONONH))
        if (j2. le. 0 ) j2 = 1
        if (j2. ge. NONONH  ) j2 = NONONH
           n2 = NONHYD(j2)
           x = atomtype(n1)
           atomtype(n1) = atomtype(n2)
           atomtype(n2) = x
           ch = charge(n1)
           charge(n1) = charge(n2)
           charge(n2) = ch
        enddo
      enddo


      return
      end subroutine shuffle_atom
c---------------------------------------
      subroutine asa_15tab(numatom,co,RPROBE)
      implicit real*4(a-h,o-z)
      parameter(maxat=100000 )

C---------------------------------------------------------------------C
C     DIMENSION PARAMETER FOR COMASA
C---------------------------------------------------------------------C

      INTEGER MAXNOH,MAXASA,NONONH
      REAL*4    RPROBE, CTLASA

      PARAMETER (MAXNOH=50000,MAXASA=350)
c     PARAMETER (RPROBE=1.4,CTLASA=4.5)
      PARAMETER (           CTLASA=4.5)
      dimension co(maxat,3)

C---------------------------------------------------------------------C
C     COMMON AREA
C---------------------------------------------------------------------C
      COMMON/ASAPO3/NONHYD(MAXNOH),INASA(MAXNOH),IASATB(MAXNOH,MAXASA)
C---------------------------------------------------------------------C
      NONONH = 0

      do i = 1, numatom
       INASA(i) = 0
       x1=co(i,1)
       y1=co(i,2)
       z1=co(i,3)
       do j = 1, numatom

       x2=co(j,1)
       y2=co(j,2)
       z2=co(j,3)
       r=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
       if(r .lt. 8.0 .and. i.ne.j ) then
       INASA(i) = INASA(i) + 1
       IASATB(i,INASA(i))=j
       endif
      
       enddo
c      write(6,*) 'asa15',INASA(i)
      enddo

      return
      end subroutine asa_15tab
c----------------------------------------



c     SUBROUTINE FENASA(CORD,EU)
      subroutine  FENASA(CORD, RPROBE, EU,part_asa )
C
C  ANALYTICAL ACCESSIBLE SURFACE AREA CALCULATION 
C  BASED ON THE PROGRAM WRITTEN BY
C   T.J.RICHMOND
C   J. Mol. Biol. (1984) 178:63-89 
C
C   Modified 9/86 by Morgan Wesson
C
C   Modified by Akira Kinjo May, 1999
C   This version is only for EMBOSS.
C
C*******************************************************************
C REGARDING EMBOSS VERSION!!!
C      Calculation of ASA energy is done when 
C      IYEFLG(23) = 1(See COMBAS & COMERG)
C      
C      Weight of ASA energy is 
C      FUWASA with default value of 1.0
C 
C*******************************************************************
C
C**  SIG=.01 is used in the test for significant 
C    sphere overlaps.
C**
c   Arguments passed in:
C     CORD - Coordinates for energy evaluation
C     EU   - User energy to be returned
C** Indices:
C**   MARC = max. no. of partial arcs for IR sphere (ith sphere)
C**   MAXA = max. no. of atoms
C**   MOV  = max. no. of IN overlapping spheres (j & k spheres for the ith)
C**   MPT  = max. no. of overlap end pts. on a circle of intersection
C
      implicit real*4(a-h,o-z)
c     include 'COMBAS'
c     include 'COMERG'
c     include 'COMMIS'
c     include 'COMASA'

      PARAMETER (MAXATM = 100000)
C   definitions of passed arguments...
      real*4 EU
      real*4 CORD(MAXATM,3)
      real*4 part_asa(MAXATM)

      real*4 AREA, L_WEIGHT

C** Variable dimension arrays:
      PARAMETER (MARC = 201)
      PARAMETER (MOV = 200)
      PARAMETER (MPT = 300) 

      LOGICAL ISKIP(200)
      INTEGER INTAG1(200),INTAG(200),ITAG(200),IDER(200),
     1   SIGN_YDER(200)
      REAL*4  XC1(MOV),YC1(MOV),ZC1(MOV),
     1       BG(MOV),THER(MOV),CSTHER(MOV),SNTHER(MOV),
     1       RI(MOV),RISQ(MOV),
     1       B1(MOV),DSQ1(MOV),BSQ1(MOV),GR(MOV),
     1       XC(MOV),YC(MOV),ZC(MOV),
     1       UX(MOV),UY(MOV),UZ(MOV),
     1       DSQ(MOV),BSQ(MOV),B(MOV)
      REAL*4  KENT(MARC),KOUT(MARC)
      REAL*4  ARCI(MPT),ARCF(MPT),EX(MPT),LT(MPT)

C**   local logical variables...
      logical LONE,LTOP,ISKIPS   


C   local integer variables...
      integer IB_LOCAL, JB_LOCAL, I, IR, IO, IN, K, L, IO1, K1, 
     1    KARC, MI, N, J, M, II, IFAIL, JBURIE

c   local real variables...
      real*4 SIG, SIGSQ, PI, PIX2, PIX4, PID2, ARCLEN, EXANG,
     1  XR, YR, ZR, RR, RRX2, RRSQ, RPLUS, TX, TY, TZ, XYSQ,
     2   CCSQ, CC, RMINUS, TXK, TYK, TZK, BSQK, BK, GL, THERK, TD,
     3   DK, GK, RISQK, RIK, T1, AXX, AXY, AXZ, AYX, AYY, AZX, AZY,
     4   AZZ, TXL, TYL, TZL, UXL, UYL, UZL, DSQL, TB, TXB, TYB, TR,
     5   TXR, TYR, TK1, TK2, THE, TF, ARCSUM, T, TT, RCN, BGL, BSQL,
     6   RISQL, WXLSQ, WXL, P, V, DEAL, DECL, DTKAL, DTKCL, S, 
     7   T2, DTLAL, DTLCL, GACA, GACB, FACA, FACB, FACC, DAX, DAY, DAZ, 
     8   TI, ASP_IR, ACOS_INPUT, SQRT_INPUT, DGK, CSTHK,SNTHK,TC,TS,
     9   CSSIG
      PARAMETER (PI=3.14159265)
      PARAMETER (PIX2=2.*PI,PIX4=4.*PI,PID2=PI/2.)
      PARAMETER (SIG=.01,SIGSQ=SIG*SIG,CSSIG=.99998)


C---------------------------------------------------------------------C
C     DIMENSION PARAMETER FOR COMASA
C---------------------------------------------------------------------C

      INTEGER MAXNOH,MAXASA,NONONH
      REAL*4    RPROBE, CTLASA

      PARAMETER (MAXNOH=50000,MAXASA=350)
c     PARAMETER (RPROBE=1.4,CTLASA=4.5)
      PARAMETER (           CTLASA=4.5)

C---------------------------------------------------------------------C
C     COMMON AREA
C---------------------------------------------------------------------C

      COMMON/ASAPO1/ NONONH
      COMMON/ASAPO3/NONHYD(MAXNOH),INASA(MAXNOH),IASATB(MAXNOH,MAXASA)
      COMMON/ASAPO4/ ASPCOF(MAXNOH),RADASA(MAXNOH)

C---------------------------------------------------------------------C
      INTEGER NONHYD,IASATB,INASA
      REAL*4    ASPCOF,RADASA
C---------------------------------------------------------------------C

      EU=0.0
C
c     L_WEIGHT = FUWASA
      L_WEIGHT = 1.0   

      DO 3 I=1,MOV
          IDER(I)=0
          SIGN_YDER(I)=0
3     CONTINUE  

C** Process each atom
C** Find the IN spheres which overlap the IR sphere

      DO 4 JR=1,NONONH
          IR = NONHYD(JR)
          ASP_IR = ASPCOF(IR) 
          ASP_IR = 1.0

          LONE=.FALSE.
          IO=1
          JB_LOCAL=0
          IB_LOCAL=0
          ARCLEN=0.
          EXANG=0.
          AREA = 0.
          XR=CORD(IR,1)
          YR=CORD(IR,2)
          ZR=CORD(IR,3)
          RR=RADASA(IR)
          RRX2=RR*2.
          RRSQ=RR*RR 

          JBURIE=0
          DO 12 JN = 1, INASA(JR)
C            IN = NONHYD(JN)
             IN = IASATB(JR,JN) 
            if (IN.eq.IR) go to 12        !exit IN loop
C
C**       Is the IN sphere next to the IR sphere
C
            RPLUS=RR+RADASA(IN)
            TX=CORD(IN,1)-XR
            IF(ABS(TX).GE.RPLUS)GO TO 12 !exit IN loop

            TY=CORD(IN,2)-YR
            IF(ABS(TY).GE.RPLUS)GO TO 12 !exit IN loop

            TZ=CORD(IN,3)-ZR
            IF(ABS(TZ).GE.RPLUS)GO TO 12 !exit IN loop
C
C**       Check for overlap of spheres by testing center to center distance
C**       Against sum and difference of radii
C
            XYSQ=TX*TX+TY*TY 
            IF (XYSQ.lt.SIGSQ) then
              TX=SIG
              TY=0.
              XYSQ=SIGSQ
            END IF

            CCSQ=XYSQ+TZ*TZ 
            CC=SQRT(CCSQ)
            IF(RPLUS-CC.LE.SIG)GO TO 12  !exit IN loop
            RMINUS=RR-RADASA(IN)
            if (CC-ABS(RMINUS).le.SIG) then
              IF(RMINUS.LE.0.) JBURIE = 1
              go to 12  !IN atom is buried, exit IN loop
            end if
C
C**       Calc. overlap parameters
C
            XC1(IO)=TX
            YC1(IO)=TY
            ZC1(IO)=TZ
   
            DSQ1(IO)=XYSQ
            BSQ1(IO)=CCSQ
            B1(IO)=CC
            GR(IO)=(CCSQ+RPLUS*RMINUS)/(RRX2*CC) !this is the 
         !distance from the IR sphere to the plane containing its 
         !intersection with the IN sphere, divided by the radius of the 
         !IR sphere.
            INTAG1(IO)=IN
            IO=IO+1
12        CONTINUE
          IF(JBURIE.EQ.1) GOTO 4 

          IO=IO-1

C
C** No sphere overlaps the sphere <IR>
C
          if (IO.eq.0) then
            AREA =PIX4*RRSQ 
            EU = EU + L_WEIGHT*ASP_IR*AREA
         part_asa(IR)=L_WEIGHT*ASP_IR*AREA
            GO TO 4 
          end if

C
C** Only one sphere overlaps the sphere <IR>
C
          if (IO.eq.1) then
            AREA = PIX2 + GR(1)*PIX2
            AREA = MOD(AREA,PIX4)
            AREA = AREA*RRSQ
            IN = INTAG1(1)
            EU = EU + L_WEIGHT*ASP_IR*AREA
         part_asa(IR)=L_WEIGHT*ASP_IR*AREA
            DGK = CCSQ - RRSQ + RADASA(IN)*RADASA(IN)
            DGK = DGK/(2*BSQ1(1)*B1(1))
            DAX = DGK*XC1(1)
            DAY = DGK*YC1(1)
            DAZ = DGK*ZC1(1)

C            K=1
C            LONE=.TRUE.
C            TXK=XC1(1)
C            TYK=YC1(1)
C            TZK=ZC1(1)
C            BSQK=BSQ1(1)
C            BK=B1(1)
C            INTAG(1)=INTAG1(1)
C            ARCSUM=PIX2
C            IB_LOCAL=IB_LOCAL+1

            GO TO 4
          end if
C
C**  Two or more sphere overlap the sphere <IR> 
C 
C
C**        Sort IN spheres by degree of overlap with IR sphere
C
          IFAIL=0
          CALL SORTAG(GR,IO,ITAG)
          DO 110 L=1,IO
            K=ITAG(L)
            IN=INTAG1(K)
            INTAG(L)=IN

            XC(L)=XC1(K)
            YC(L)=YC1(K)
            ZC(L)=ZC1(K)

            DSQ(L)=DSQ1(K)
            B(L)=B1(K)
            BSQ(L)=BSQ1(K)
            ISKIP(L)=.FALSE.
110       CONTINUE      

          DO 120 L=1,IO
            GL=GR(L)*RR
            BG(L)=B(L)*GL
            RISQ(L)=RRSQ-GL*GL 
            RI(L)=SQRT(RISQ(L))
C
C**       Radius of the IN circle on the surface of the sphere
C
            THER(L)=PID2-ASIN(GR(L))
            CSTHER(L) = COS(THER(L))
            SNTHER(L) = SIN(THER(L))
120       CONTINUE 

C
C** Find boundary of inaccessible area on IR sphere
C
          IO1=IO-1
          JBURIE=0
          DO 42 K=1,IO1
          if (ISKIP(K)) GO TO 42 
            TXK=XC(K)
            TYK=YC(K)
            TZK=ZC(K)
            BK=B(K)
            THERK=THER(K)
            CSTHK=CSTHER(K)
            SNTHK=SNTHER(K)
            K1=K+1
            DO 31 L=K1,IO
            if (ISKIP(L)) GO TO 31 
C
C** Is L circle intersecting K circle?
C** Distance between circle centers and sum of radii
C
               ACOS_INPUT = (TXK*XC(L)+TYK*YC(L)+TZK*ZC(L))/(BK*B(L))
               IF(ABS(ACOS_INPUT) .GT.1) 
     +            ACOS_INPUT = SIGN(1.,ACOS_INPUT)
C
C** CC = acos( ACOS_INPUT )
C
               TC = CSTHK*CSTHER(L)
               TS = SNTHK*SNTHER(L)
               TD=THERK+THER(L)
               IF(TD.GT.PI) GO TO 23 
C
C** Circles enclose separate regions?
C               IF(CC.GE.TD)GO TO 31
C
               IF(ACOS_INPUT.LE.TC-TS) GO TO 31
C
C** Circle L completely inside circle K?
C               if (CC+THER(L).lt.THER(K)) then
C
                IF(ACOS_INPUT.GT.TC+TS) GOTO 10
C
C** IF(CC.gt.SIG) then, Circles <L> and <K> are essentially parallel.
C   If so, <L> is completely inside circle <K> by virtue of 
C   sorting carried out above.
C

23          IF(ACOS_INPUT.LT.CSSIG) GOTO 25
10          ISKIP(L)=.TRUE.

C
C** See if the sphere <IR> is completely buried
C   IF(PIX2-CC.LE.TD) JBURIE=1
C

25          IF(ACOS_INPUT.LE.TC-TS) JBURIE=1

C
31          CONTINUE
            IF(JBURIE .EQ.1) GOTO 4 
42        CONTINUE 
C
C** Find T value of circle intersections
C
          DO 140 K=1,IO
          if (ISKIP(K)) GOTO 140 
            ISKIPS=ISKIP(K)
            ISKIP(K)=.TRUE.
            KARC=0
            LTOP=.FALSE.
            TXK=XC(K)
            TYK=YC(K)
            TZK=ZC(K)
            DK=SQRT(DSQ(K))
            BSQK=BSQ(K)
            BK=B(K)
            GK=GR(K)*RR
            RISQK=RISQ(K)
            RIK=RI(K)
            THERK=THER(K)
C
C** Rotation matrix elements
C
            T1=TZK/(BK*DK)
            AXX=TXK*T1
            AXY=TYK*T1
            AXZ=DK/BK
            AYX=TYK/DK
            AYY=TXK/DK
            AZX=TXK/BK
            AZY=TYK/BK
            AZZ=TZK/BK
C
*VOCL LOOP,NOVREC(ARCF,ARCI,LT,EX)
C
            DO 150 L=1,IO
            if (ISKIP(L)) GO TO 150 
              TXL=XC(L)
              TYL=YC(L)
              TZL=ZC(L)
C
C** Rotate spheres so K vector colinear with z-axis
C
              UXL=TXL*AXX+TYL*AXY-TZL*AXZ
              UYL=TYL*AYY-TXL*AYX
              UZL=TXL*AZX+TYL*AZY+TZL*AZZ
C
C              ACOS_INPUT = UZL/B(L)
C
C If (acos( ACOS_INPUT ).lt.THERK+THER(L)) GOTO 150
C
              TD=THERK+THER(L)
              IF(TD.LT.PI.AND.UZL/B(L).LT.COS(TD)) GOTO 150
C
              GL=GR(L)*RR
              DSQL=UXL*UXL+UYL*UYL 
              TB=UZL*GK-BG(L)
              TXB=UXL*TB
              TYB=UYL*TB
              TD=RIK*DSQL
              SQRT_INPUT = RISQK * DSQL - TB*TB 
              IF(SQRT_INPUT .lt. 0.0) SQRT_INPUT = 0.0 
              TR = SQRT( SQRT_INPUT )
              TXR=UXL*TR
              TYR=UYL*TR
C
C** T values of intersection for K circle
C
              TB=(TXB+TYR)/TD
              if (ABS(TB).gt.1.) TB=SIGN(1.,TB)
              TK1=ACOS(TB)
              if (TYB-TXR.lt.0.) TK1=PIX2-TK1

              TB=(TXB-TYR)/TD
              if (ABS(TB).gt.1.) TB=SIGN(1.,TB)
              TK2=ACOS(TB)
              if (TYB+TXR.lt.0.) TK2=PIX2-TK2

              ACOS_INPUT = (RRSQ*UZL-GK*BG(L))/(RIK*RI(L)*B(L))
              IF(ABS(ACOS_INPUT).GT.1.) 
     +           ACOS_INPUT = SIGN(1.,ACOS_INPUT)
C              THE = -acos( ACOS_INPUT )
               THE = ACOS_INPUT 
C
C** Is TK1 entry or exit point?  check T=0 point.
C** TI IS EXIT PT., TF IS ENTRY PT.
C
              ACOS_INPUT = (UZL*GK-UXL*RIK)/(B(L)*RR)
              IF(ABS(ACOS_INPUT).GT.1.) 
     +           ACOS_INPUT = SIGN(1.,ACOS_INPUT)

              if ((ACOS(ACOS_INPUT)-THER(L))*(TK2-TK1).le.0) then
                TI=TK2
                TF=TK1
              else
                TI=TK1
                TF=TK2
              end if

              KARC=KARC+1

              if (TF.le.TI) then
                ARCF(KARC)=TF
                ARCI(KARC)=0.
                TF=PIX2
                LT(KARC)=L
                EX(KARC)=THE
                LTOP=.TRUE.
                KARC=KARC+1
              end if

              ARCF(KARC)=TF
              ARCI(KARC)=TI
              LT(KARC)=L
              EX(KARC)=THE
              UX(L)=UXL             
              UY(L)=UYL             
              UZ(L)=UZL             
150         CONTINUE    
C           end of <L> loop
            ISKIP(K)=ISKIPS
C
C** Special case: K circle without intersections?
C
            if (KARC.LE.0) then
              ARCSUM=PIX2
              IB_LOCAL=IB_LOCAL+1
              go to 19
            end if
C
C** General case: sum up arclength and set connectivity code
C
            IFAIL=0
            CALL SORTAG(ARCI,KARC,ITAG)
            ARCSUM=ARCI(1)
            MI=ITAG(1)
            T=ARCF(MI)
            N=MI
            if (KARC.ne.1) then
              DO 152 J=2,KARC
                 M=ITAG(J)
                 if (T.lt.ARCI(J)) then
                   ARCSUM=ARCSUM+ARCI(J)-T
                   EXANG=EXANG-ACOS(EX(N))
                   JB_LOCAL=JB_LOCAL+1
                   L=LT(N)
                   IDER(L)=IDER(L)+1
                   SIGN_YDER(L) = SIGN_YDER(L) + 1
                   KENT(JB_LOCAL)=L*1024+K
                   L=LT(M)
                   IDER(L)=IDER(L)+1
                   SIGN_YDER(L) = SIGN_YDER(L) - 1
                   KOUT(JB_LOCAL)=K*1024+L
                 end if
                 TT=ARCF(M)
                 if (TT.ge.T) then
                   T=TT
                   N=M
                 end if
152           CONTINUE
            end if
            ARCSUM=ARCSUM+PIX2-T
            if (.not.LTOP) then
              EXANG=EXANG-ACOS(EX(N))
              JB_LOCAL=JB_LOCAL+1
              L=LT(N)
              IDER(L)=IDER(L)+1
              SIGN_YDER(L) = SIGN_YDER(L) + 1
              KENT(JB_LOCAL)=L*1024+K
              L=LT(MI)
              IDER(L)=IDER(L)+1
              SIGN_YDER(L) = SIGN_YDER(L) - 1
              KOUT(JB_LOCAL)=K*1024+L
            end if
C
C** Calculate derivatives.  Comments refer to the equivalent notation in
C Tim Richmond's paper.
C
            DO 153 L=1,IO   !K ; K index is J
            if (IDER(L).eq.0) GOTO 153 
              RCN=IDER(L)*RRSQ      !IDER(L)* rho**2
              IDER(L)=0
              UZL=UZ(L)         !C(K)
              GL=GR(L)*RR         !G(K)
              BGL=BG(L)         !D(K)*G(K)
              BSQL=BSQ(L)         !D(K)**2
              RISQL=RISQ(L)      !R(K)**2
              WXLSQ=BSQL-UZL*UZL      !A(K)**2
              WXL=SQRT(WXLSQ)      !A(K)
              P=BGL-GK*UZL      !D(K)*G(K) - G(J)*C(K)
              V=RISQK*WXLSQ-P*P      !E ** 2
              IF(V .lt. 0.000001) V = 0.000001
              V=SQRT(V)         !E
              T1=RR*(GK*(BGL-BSQL)+UZL*(BGL-RRSQ))/(V*RISQL*BSQL)
c           rho*( G(J)* (D(K)*G(K) - D(K)**2) + C(K)( D(K)*G(K)-rho**2))
c            / (E*R(K)**2 * D(K) **2)
              DEAL=-WXL*T1      !d(omega)/dA(K)
              DECL=-UZL*T1-RR/V      !d(omega)/dC(K)
              DTKAL=(WXLSQ-P)/(WXL*V)   !d( T(lambda+1))/dA(K)
              DTKCL=(UZL-GK)/V      !d( T(lambda+1))/dC(K)
              S=GK*B(L)-GL*UZL      !G(J)*D(K) - G(K)*C(K)
              T1=2.*GK-UZL      !2*(G(J)-C(K)
              T2=RRSQ-BGL         !rho **2 - D(K)*G(K)
              DTLAL=-(RISQL*WXLSQ*B(L)*T1-S*(WXLSQ*T2+RISQL*BSQL))
     1        /(RISQL*WXL*BSQL*V)      !d( T(lambda))/dA(K)
              DTLCL=-(RISQL*B(L)*(UZL*T1-BGL)-UZL*T2*S)/(RISQL*BSQL*V)
c                   !d( T(lambda))/dC(K)
              GACA=RCN*(DEAL-(GK*DTKAL-GL*DTLAL)/RR)/WXL
                     !x component of derivative in
                     !mu coordinate system
              GACB = GK - UZL*GL/B(L)   !y component of 
              GACB = GACB * SIGN_YDER(L) * RR / WXLSQ !derivative in
                         !mu coordinate system
              SIGN_YDER(L) = 0
              FACA=UX(L)*GACA - UY(L)*GACB !rotate back around z axis.
              FACB=UY(L)*GACA + UX(L)*GACB
              FACC=RCN*(DECL-(GK*DTKCL-GL*DTLCL)/RR) !z component of 
                     !derivative in
                     !mu coordinate system
              DAX=AXX*FACA-AYX*FACB+AZX*FACC !transform back to the
              DAY=AXY*FACA+AYY*FACB+AZY*FACC !real coordinates
              DAZ=AZZ*FACC-AXZ*FACA
              IN=INTAG(L)
153         CONTINUE 
C
19      ARCLEN=ARCLEN+GR(K)*ARCSUM
            IN=INTAG(K)
            T1=ARCSUM*RRSQ*(BSQK-RRSQ+RADASA(IN)*RADASA(IN))
     +          /(RRX2*BSQK*BK)

            if (LONE) go to 56
140       CONTINUE
C         end of <K> loop
C
          if (ARCLEN.eq.0.) go to 4
          if (JB_LOCAL.eq.0) go to 56
C
C** Find number of independent boundaries
C
          J=0
          DO 60 K=1,JB_LOCAL
          if (KOUT(K).eq.0) GOTO 60
            I=K
62          N=KOUT(I)
            KOUT(I)=0
            J=J+1
            DO 61 II=1,JB_LOCAL
            if (N.ne.KENT(II)) GOTO 61 
               if (II.eq.K) then
                 IB_LOCAL=IB_LOCAL+1
                 if (J.eq.JB_LOCAL) go to 56
                 go to 60
               end if
             I=II
             go to 62
61          CONTINUE 
60       CONTINUE
C
         IB_LOCAL=IB_LOCAL+1
56       AREA=IB_LOCAL*PIX2+EXANG+ARCLEN
         AREA=MOD(AREA,PIX4)
16       AREA=AREA*RRSQ

         EU = EU + L_WEIGHT*ASP_IR*AREA
         
         part_asa(IR)=L_WEIGHT*ASP_IR*AREA
c        write(6,*) 'asp',L_WEIGHT,ASP_IR,AREA
 
4     CONTINUE

      return
      end subroutine  FENASA

c----------------------------------


C
C************SORTAG************
C called from DFEASA and FENASA
C
      SUBROUTINE SORTAG(A,N,TAG)
      INTEGER TAG,TG
      REAL*4 A 
      DIMENSION A(N),IU(16),IL(16),TAG(N)
      DO 1  I=1,N
      TAG(I)=I
    1 CONTINUE
      M=1
      I=1
      J=N
  5   IF(I.GE.J) GO TO 70
 10   K=I
      IJ=(J+I)/2
      T=A(IJ)
      IF(A(I).LE.T) GO TO 20
      A(IJ)= A(I)
      A(I)=T
      T=A(IJ)
      TG=TAG(IJ)
      TAG(IJ)=TAG(I)
      TAG(I)=TG
 20   L=J
      IF(A(J).GE.T) GO TO 40
      A(IJ)=A(J)
      A(J)=T
      T=A(IJ)
      TG=TAG(IJ)
      TAG(IJ)=TAG(J)
      TAG(J)=TG
      IF(A(I).LE.T) GO TO 40
      A(IJ)=A(I)
      A(I)=T
      T=A(IJ)
      TG=TAG(IJ)
      TAG(IJ)=TAG(I)
      TAG(I)=TG
      GO TO 40
 30   A(L)=A(K)
      A(K)=TT
      TG=TAG(L)
      TAG(L)=TAG(K)
      TAG(K)=TG
 40   L=L-1
      IF(A(L).GT.T) GO TO 40
      TT=A(L)
 50   K=K+1
      IF(A(K).LT.T) GO TO 50
      IF(K.LE.L) GO TO 30
      IF(L-I.LE.J-K) GO TO 60
      IL(M)=I
      IU(M)=L
      I=K
      M=M+1
      GO TO 80
 60   IL(M)=K
      IU(M)=J
      J=L
      M=M+1
      GO TO 80
 70   M=M-1
      IF(M.EQ.0) RETURN
      I=IL(M)
      J=IU(M)
 80   IF(J-I.GE.1) GO TO 10
      IF(I.EQ.1) GO TO 5
      I=I-1
 90   I=I+1
      IF(I.EQ.J) GO TO 70
      T=A(I+1)
      IF(A(I).LE.T) GO TO 90
      TG=TAG(I+1)
      K=I
 100  A(K+1)=A(K)
      TAG(K+1)=TAG(K)
      K=K-1
      IF(T.LT.A(K)) GO TO 100
      A(K+1)=T
      TAG(K+1)=TG
      GO TO 90
      END SUBROUTINE SORTAG
c------------------------
