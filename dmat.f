C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C:::::::::::::::: DMATA  Create harmonic matrix  ::::::::::::::::::
      SUBROUTINE  DMATA
      PARAMETER (NA1=1995,NA2=249,NA3=755,NA4=755,NA5=400,NA6=121)
      COMMON  AAA, BBB, NONC, CCC, NON2, NCOIL, NHOWA, 
     &        NPO1,NOM, NPO2, NPO3, NPOR, NELEM, NB, NDEG, NDE, DF,
     &        DH(NA1,NA2), DK(NA1), DAA(NA1), DA(NA1),
     &        DS(NA3,18), DCS(NA3), DB(NA3,22), DD(NA3,NA6), NOD(0:NA4,4),           
     &        XY(NA5,2), DENRYU(1,11), DC(NA3,11), DCPRE(NA3,11),
     &        DMUO, DOMEG, DPI, ITR, TOTAL, 
     &        DN(11,11), DBH(11,4)
      DIMENSION  AB(0:72), VSC(4), VSS(4)
C
      DO 10  NE=1, NELEM
        IF (INT(NOD(NE,4) / 100) .EQ. 2)  GO TO 1000
        DO 12  K=1, 25
          DD(NE,K) = 0.0
   12   CONTINUE
        DD(NE,1) = 1.0 / DMUO
        DD(NE,7) = DD(NE,1)
        DD(NE,13) = DD(NE,1)
        DD(NE,19) = DD(NE,1)
        DD(NE,25) = DD(NE,1)
        GO TO  9999
 1000   CONTINUE
        HD = 2.0 * DPI / 72.0
        CO = 1.0 / DPI
        DO 30  NT=0, 72
          DT = FLOAT(NT)
          SIN1 = SIN(HD*DT)
          COS1 = COS(HD*DT)
          SIN2 = SIN(2.0*HD*DT)
          COS2 = COS(2.0*HD*DT)
          BX = DB(NE,1)+DB(NE,2)*SIN1+DB(NE,3)*COS1+DB(NE,4)*SIN2 
     &          +DB(NE,5)*COS2
          BY = DB(NE,6)+DB(NE,7)*SIN1+DB(NE,8)*COS1+DB(NE,9)*SIN2
     &          +DB(NE,10)*COS2
          BT = SQRT(BX*BX + BY*BY)
          AB(NT) = DBH(1,1)
   30   CONTINUE
        SO = AB(0)
        DO 40  NT=1, 71, 2
          SO = SO + 4.0 * AB(NT) + 2.0 * AB(NT+1)
   40   CONTINUE
        VO=(SO - AB(72)) * HD/3.0 * CO * 0.5
C------------------------------------------------------------------
C       IF (NE .EQ. 48)  THEN
C         PRINT  '(1X, " VO=", E14.6)', VO
C       END IF
C------------------------------------------------------------------
        DO 50  N=1, 4
          DRN = FLOAT(N)
          SC = AB(0)
          DO 60  NT=1, 71, 2
            DT = FLOAT(NT)
            YC = AB(NT) * COS(DRN*HD*DT)
            SC = SC+4. *YC
            YC = AB(NT+1)*COS(DRN*HD*(DT+1.0))
            SC = SC + 2.0 * YC
   60     CONTINUE
          VSC(N) = (SC - YC) * HD/3.0 * CO
C         SS = AB(0)
          SS = 0.0
          DO 70  NT=1, 71, 2
            DT = FLOAT(NT)
            YS = AB(NT) * SIN(DRN * HD * DT)
            SS = SS + 4.0 * YS
            YS = AB(NT+1) * SIN(DRN*HD*(DT+1.0))
            SS = SS + 2.0 * YS
   70     CONTINUE
          VSS(N) = (SS-YS) * HD /3.0 * CO
C--------------------------------------------------------------------------
C         IF(NE .EQ. 48)  THEN
C           PRINT '(1X, "VSC(", I1, ")=", E14.6)', N, VSC(N)
C           PRINT '(1X, "VSS(", I1, ")=", E14.6)', N, VSS(N)
C         END IF 
C------------------------------------------------------------------------------
  50    CONTINUE
        DD(NE,1) = VO
        DD(NE,2) = VSS(1) * 0.5
        DD(NE,3) = VSC(1) * 0.5
        DD(NE,4) = VSS(2) * 0.5
        DD(NE,5) = VSC(2) * 0.5
C
        DD(NE,6) = DD(NE,2) * 2.0 
        DD(NE,7) = VO - VSC(2) * 0.5
        DD(NE,8) = VSS(2) * 0.5
        DD(NE,9) = (VSC(1) - VSC(3)) * 0.5
        DD(NE,10) =(-VSS(1) + VSS(3)) * 0.5
C
        DD(NE,11) = DD(NE,3) * 2.0
        DD(NE,12) = DD(NE,8)
        DD(NE,13) = VO + VSC(2) * 0.5
        DD(NE,14) = (VSS(1) + VSS(3)) * 0.5
        DD(NE,15) = (VSC(1) + VSC(3)) * 0.5
C      
        DD(NE,16) = DD(NE,4) * 2.0
        DD(NE,17) = DD(NE,9)
        DD(NE,18) = DD(NE,14)
        DD(NE,19) = VO - VSC(4) * 0.5
        DD(NE,20) = VSS(4) * 0.5
C
        DD(NE,21) = DD(NE,5) * 2.0
        DD(NE,22) = DD(NE,10)
        DD(NE,23) = DD(NE,15)
        DD(NE,24) = DD(NE,20)
        DD(NE,25) = VO + VSC(4) * 0.5
 9999   CONTINUE
   10 CONTINUE
      PRINT  '(10X, " NDEG=5")'
      PRINT  '(10X, " DMATA")'
      RETURN
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C::::::::::::::::: DMATB  Create harmonic matrix ::::::::::::::::::
      SUBROUTINE  DMATB
      PARAMETER (NA1=1995,NA2=249,NA3=755,NA4=744,NA5=400,NA6=121)
      COMMON  AAA, BBB, NONC, CCC, NONC2, NCOIL, NHOWA,
     &        NPO1, NOM, NPO2, NPOR, NELEM, NB, NDEG, NDE, DF, 
     &        DH(NA1,NA2), DK(NA1), DAA(NA1), DA(NA1), 
     &        DS(NA3,18), DCS(NA3), DB(NA3,22), DD(NA3,NA6), NOD(0:NA4,4),
     &        XY(NA5,2), DENRYU(1,11), DC(NA3,11), DCPRE(NA3,11),
     &        DMUO, DOMEG, DPI, ITR, TOTAL, 
     &        DN(11,11), DBH(11,4)
      DIMENSION  AB(0:72), VSC(6), VS(6)
C
      DO 10 NE=1, NELEM
        IF (INT(NOD(NE,4)/100) .EQ. 2)  GO TO 1000
        DO 12  K=1, 49
          DO (NE,K) = 0.0
   12   CONTINUE
        DD(NE,1) = 1.0 / DMUO
        DD(NE,9) = DD(NE,1)
        DD(NE,17) = DD(NE,1)
        DD(NE,25) = DD(NE,1)
        DD(NE,33) = DD(NE,1)
        DD(NE,41) = DD(NE,1)
        DD(NE,49) = DD(NE,1)
        GO TO  9999
 1000   CONTINUE
        HD = 2.0 * DPI / 72.0
        CO = 1.0 / DPI
        DO 30  T=0, 72
          SIN1 = SIN(HD*T)
          COS1 = COS(HD*T)
          SIN2 = SIN(2.0*HD*T)
          COS2 = COS(2.0*HD*T)
          SIN3 = SIN(3.0*HD*T)
          COS3 = COS(3.0*HD*T)
          BX = DB(NE,1)+DB(NE,2)*SIN1+DB(NE,3)*COS1+DB(NE,4)*SIN2
     &         +DB(NE,5)*COS2+DB(NE,6)*SIN3+DB(NE,7)*COS3
          BY = DB(NE,8)+DB(NE,9)*SIN1+DB(NE,10)*COS1+DB(NE,11)*SIN2
     &         +DB(NE,12)*COS2+DB(NE,13)*SIN3+DB(NE,14)*COS3
          BT = SQRT(BX*BX+BY*BY)
          AB(T) = DBH(1,1)
C-------------------------------------------------FERRITE  H7C1
C         IF (ABS(BT) .GT. 0.4)  GOTO 510
C         AB(T) = 90.5 + 1582.4*BT**4
C         GO TO 590
C 510     IF (ABD(BT) .GT. 0.435)  GOTO 520
C         AB(T) = 50.0 + 1428.0 * (ABS(BT) - 0.4)
C         GOTO 580
C 520     IF (ABS(BT) .GT. 0.46)  GOTO  530
C         AB(T) = 100.0 + 4.0E3 * (ABS(BT) - 0.435)
C         GOTO 580
C 530     IF (ABS(BT) .GT. 0.49)  GOTO 540
C         AB(T) = 200.0 + 1.0E4*(ABS(BT) - 0.460)
C         GOTO  580 
C 540     IF(ABS(BT) .GT. 0.51)  GOTO 550
C         AB(T) = 500. +2.0E4 * (ABS(BT) - 0.490)
C         GOTO 580
C 550     IF (ABS(BT) .GT. 0.52)  GOTO 560
C         AB(T) = 900.0 + 6.0E4 * (ABS(BT) - 0.510)
C         GOTO 580
C 560     AB(T) = 1500.0 +1.0E5 * (ABS(BT) - 0.520)
C 580     AB(T) = AB(T) / ABS(BT)
C 590     CONTINUE
C------------------------------------------------------------------------------------
   30   CONTINUE
        SO = AB(0)
        DO 40  T=1, 71, 2
          SO = SO + 4.0 * AB(T) + 2.0 * AB(T+1)
   40   CONTINUE
        VO = (SO-AB(72))* HD / 3.0 * CO * 0.5
        DO 50  N=1, 6
          SC = AB(0)
          DO 60  T=1, 71, 2
            YC = AB(T) * COS(N*HD*T)
            SC = SC + 4.0 * YC
            YC = AB(T+1) * COS(N*HD*(T+1))
            SC = SC +2.0 * YC
   60     CONTINUE
          VSC(N) = (SC-YC) * HD / 3. *CO
          SS=0.0
          DO 70  T=1, 71, 2
            YS = AB(T) * SIN(N*HD*T)
            SS = SS + 4.0 *YS
            YS = AB(T+1) * SIN(N*HD*(T+1))
            SS = SS + 2.0 * YS
   70     CONTINUE
          VSS(N) = (SS-YS) * HD / 3.0 * CO
   50   CONTINUE
        DD(NE,1) = VO
        DD(NE,2) = VSS(1) * 0.5
        DD(NE,3) = VSC(1) * 0.5
        DD(NE,4) = VSS(2) * 0.5
        DD(NE,5) = VSC(2) * 0.5
        DD(NE,6) = VSS(3) * 0.5
        DD(NE,7) = VSC(3) * 0.5
        DD(NE,9) = VO - VSC(2) * 0.5
        DD(NE,10) = VSS(2) * 0.5
        DD(NE,11) = (VSC(1) - VSC(3)) * 0.5
        DD(NE,12) = (-VSS(1) + VSS(3)) * 0.5
        DD(NE,13) = (VSC(2) - VSC(4)) * 0.5
        DD(NE,14) =(-VSS(2) +VSS(4)) * 0.5
        DD(NE,17) = VO + VSC(2) * 0.5
        DD(NE,18) = (VSS(1) + VSS(3)) * 0.5
        DD(NE,19) = (VSC(1) + VSC(3)) * 0.5
        DD(NE,20) = (VSS(2) + VSS(4)) * 0.5
        DD(NE,21) = (VSC(2) + VSC(4)) * 0.5
        DD(NE,25) = VO - VSC(4) * 0.5
        DD(NE,26) = VSS(4) * 0.5
        DD(NE,27) = (VSC(1)- VSC(5)) * 0.5
        DD(NE,28) = (-VSS(1) + VSS(5)) * 0.5
        DD(NE,33) = VO + VSC(4) * 0.5
        DD(NE,34) = (VSS(1) + VSS(5)) * 0.5
        DD(NE,35) = (VSC(1) + VSC(5)) * 0.5
        DD(NE,41) = VO - VSC(6) * 0.5
        DD(NE,42) = VSS(6) * 0.5
        DD(NE,49) = VO + VSC(6) * 0.5
        DD(NE,8)  = DD(NE,2) * 2.0
        DD(NE,15) = DD(NE,3) * 2.0
        DD(NE,16) = DD(NE,10)
        DD(NE,22) = DD(NE,4) * 2.0
        DD(NE,23) = DD(NE,11)
        DD(NE,24) = DD(NE,18)
        DD(NE,29) = DD(NE,5) * 2.0
        DD(NE,30) = DD(NE,12)
        DD(NE,31) = DD(NE,19)
        DD(NE,32) = DD(NE,26)
        DD(NE,36) = DD(NE,6) * 2.0
        DD(NE,37) = DD(NE,13)
        DD(NE,38) = DD(NE,20)
        DD(NE,39) = DD(NE,27)
        DD(NE,40) = DD(NE,34)
        DD(NE,43) = DD(NE,7) * 2.0
        DD(NE,44) = DD(NE,14)
        DD(NE,45) = DD(NE,21)
        DD(NE,46) = DD(NE,28)
        DD(NE,47) = DD(NE,35)
        DD(NE,48) = DD(NE,42)      
 9999   CONTINUE
   10 CONTINUE
      PRINT  '(10X, "NDEG=7")'
      PRINT  '(10X, "DMATB")'
      RETURN
      END
C:::::::::::::::::::::::::::::DMATC ?????????????????????????::::::::::::::::::::::::::::::::::::::::::
          SUBROUTINE        DMATC
           PARAMETER  (NA1=1995, NA2=249, NA3=755, NA4=744, NA5=400, NA6=121)
           COMMON     AAA, BBB, NONC, CCC, NONC2, NCOIL, NHOWA,
         &                      NPO1, NOM, NPO2, NPOR, NELEM, NB, NDEG, NDE, DF, 
         &                      DH(NA1, NA2), DK(NA1), DAA(NA1), DA(NA1), 
         &                 DS(NA3, 18), DCS(NA3), DB(NA3, 22), DD(NA3, NA6), NOD(0:NA4,4),
         &                 XY(NA5, 2), DENRYU(1, 11), DC(NA3, 11), DCPRE(NA3, 11),
         &                 DMUO, DOMEG, DPI, ITR, TOTAL, 
         &                 DN(11, 11), DBH(11, 4)
            DIMENSION    AB(0: 72), VSC(10), VSS(10)
C***********???????????????????*************
C         IMPLICIT  REAL* 8(D)
C
C <<<<<<<<<<<<<<???????????LOOP>>>>>>>>>>>>>>>
           DO  10  NE=1, NELEM
           IF  ( INT(NOD(NE, 4)  / 100). EQ. 2)              GO    TO   1000
C
C------------------?????????????????-----------------------------
                 DO   12  K=1, NA6
                 DD(NE, K)= 0.0
       12      CONTINUE
                 DD(NE, 1) = 1. 0 / DMUO
                 DD(NE, 13) = DD(NE, 1)
                 DD(NE, 25) = DD(NE,  1)
                 DD(NE, 37) = DD(NE, 1)
                 DD(NE, 49) = DD(NE,  1)
                 DD(NE, 61) = DD(NE, 1)
                 DD(NE, 73) = DD(NE,  1)
                 DD(NE, 85) = DD(NE, 1)
                 DD(NE, 97) = DD(NE,  1)
                 DD(NE, 109) = DD(NE, 1)
                 DD(NE, 121) = DD(NE,  1)
      GO      TO    9999
C
C--------------------??????????????200 < = NOD(I,4) < =299-----------------------
 1000   CONTINUE
            HD = 2. *DPI / 72.
            CO = 1. / DPI
            DO  30  T=0, 72
            SIN1= SIN (HD* T)
            COS1 = COS(HD*T)
            SIN2= SIN (2. *HD* T)
            COS2 = COS(2. *HD*T)
            SIN3= SIN (3. *HD* T)
            COS3 = COS(3. *HD*T)
            SIN4= SIN (4. *HD* T)
            COS4 = COS(4. *HD*T)
            SIN5= SIN (5. *HD* T)
            COS5 = COS(5. *HD*T)
            BX = DB(NE, 1) + DB(NE, 2) * SIN1+ DB(NE, 3) * COS1 + DB (NE, 4) * SIN 2
     &           + DB (NE, 5) * COS 2 + DB (NE, 6) * SIN 3 + DB (NE, 7) * COS 3
     & + DB(NE, 8) * SIN4 + DB(NE,9) * COS4 + DB (NE,10) * SIN5 + DB (NE,11) * COS5
           BY = DB(NE, 12)+ DB(NE, 13) * SIN1 + DB(NE, 14) *COS1 + DB(NE, 15) *SIN2
     &            +DB(NE, 16) * COS2 + DB (NE, 17) * SIN3 + DB(NE, 18) * COS3
     & + DB(NE,19) * SIN4 + DB(NE,20) * COS4 + DB(NE,21) * SIN5 + DB(NE, 22) * COS5 
C        DBX =HD * (DB(NE, 2) *COS1 - DB(NE, 3) * SIN1 + 2 *DB(NE, 4) * COS2
C   &              - 2 * DB(NE,5) *SIN2 + 3* DB(NE,6) *COS3 - 3 *DB(NE,7) * SIN 3)
C        DBY =HD * (DB(NE, 9) *COS1 - DB(NE, 10) * SIN1 + 2 *DB(NE,11) * COS2
C   &              - 2 * DB(NE,12) *SIN2 + 3* DB(NE,13) *COS3 - 3 *DB(NE,14) * SIN 3)
           BT =  SQRT(BX * BX + BY *BY)
C        JIKITEKOURITU            H/B = F(B) /B
          AB( T) = DBH(1, 1)
C       AB( T) = 90.5 + 1582.4 * (BT**4)
C       AB( T) = 15.9 + 166. *(BT**2)
C   ******  FERRITE  H7C1 ******
C
C                IF  (ABS (BT) . GT. 0. 4)                    GO   TO  510
C                       AB(T) = 90.5 + 1582.4 * BT ** 4
C                                   GO   TO    590
C   510      IF (ABS( BT). GT.  0. 435)                  GO   TO  520
C                     AB(T) = 50. + 1428. *(ABS(BT) - 0.4)
C                                    GO   TO   580
C  520      IF  (ABS(BT). GT. 0. 46)                      GO   TO  530
C                    AB(T) =  100. + 4. 0E3 * (ABS(BT) - 0. 435)
C                                   GO   TO  580
C  530      IF  (ABS(BT). GT. 0. 49)                      GO   TO  540
C                    AB(T) =  200. + 1. 0E4 * (ABS(BT) - 0. 460)
C                                   GO   TO  580
C  540      IF  (ABS(BT). GT. 0. 51)                      GO   TO  550
C                    AB(T) =  500. + 2. 0E4 * (ABS(BT) - 0. 490)
C                                   GO   TO  580
C  550      IF  (ABS(BT). GT. 0. 52)                      GO   TO  560
C                    AB(T) =  900. + 6. 0E4 * (ABS(BT) - 0. 510)
C                                   GO   TO  580
C  560            AB(T) = 1500. + 1. 0E5 * (ABS (BT) - 0.520)
C  580            AB(T) = AB(T) / ABS(BT)
C 590      CONTINUE
C*******************************************************
C
C         AB(T) = 20. +100. *BT *BT + . 5 * (2. *BX*DBX + 2. *BY*DBY)
   30    CONTINUE
            SO = AB(0)
           DO    40   T=1, 71, 2
           SO = SO + 4.*AB(T)  +2. *AB(T+1)
   40   CONTINUE
            VO+ (SO - AB(72)) * HD/ 3.0 * CO * 0.5
            DO   50   N=1, 10
           SC = AB(0)
           DO   60  T=1, 71, 2
           YC = AB(T) * COS(N*HD*T)
           SC = SC + 4. * YC
           YC = AB(T +1) * COS(N*HD *(T+1))
           SC = SC + 2.. *YC
   60   CONTINUE
          VSC(N) = (SC - YC) *HD/ 3.*CO
          SS=0.
          DO  70     T=1, 71, 2
          YS = AB (T) * SIN(N * HD* (T)
          SS= SS + 4.*YS
          YS = AB(T +1) * SIN(N*HD*(T+1))
          SS = SS +2. *YS
  70   CONTINUE
         VSS(N) = (SS-YS) * HD/3.0 *CO
  50    CONTINUE
         DD(NE, 1) = VO
         DD(NE, 2) = VSS(1) *. 5
         DD(NE, 3) = VSS(1) *. 5
         DD(NE, 4) = VSS(2) *. 5
         DD(NE, 5) = VSS(2) *. 5
         DD(NE, 6) = VSS(3) *. 5
         DD(NE, 7) = VSS(3) *. 5
         DD(NE, 8) = VSS(4) *. 5
         DD(NE, 9) = VSS(4) *. 5
         DD(NE, 10) = VSS(5) *. 5
         DD(NE, 11) = VSS(5) *. 5
C
        DD(NE, 12) = DD (NE, 2) * 2.
        DD(NE, 13) =  VO - VSC(2) *. 5
        DD(NE, 14) = VSS(2) *. 5
        DD(NE, 15) = (VSC(1) - VSC(3))*. 5
        DD(NE, 16) = ( - VSS(1) + VSS(3))*. 5
        DD(NE, 17) = (VSC(2) - VSC(4))*. 5
        DD(NE, 18) = (- VSS(2) + VSS(4))*. 5
        DD(NE, 19) = (VSC(3) - VSC(5))*. 5
        DD(NE, 20) = ( - VSS(3) + VSS(5))*. 5
        DD(NE, 21) = (VSC(4) - VSC(6))*. 5
        DD(NE, 22) = (- VSS(4) + VSS(6))*. 5
C
        DD(NE, 23) = DD(NE, 3) *2. 
        DD(NE, 24) = DD(NE, 14)
        DD(NE, 25) = VO + VSC( 2)*. 5
        DD(NE, 26) = (VSS(1) + VSS(3)) *. 5
        DD(NE, 27) = (VSC(1) + VSC(3)) *.5
        DD(NE, 28) = (VSS(2) + VSS(4)) *. 5
        DD(NE, 29) = (VSC(2) + VSC(4)) *.5
        DD(NE, 30) = (VSS(3) + VSS(5)) *. 5
        DD(NE, 31) = (VSC(3) + VSC(5)) *.5
        DD(NE, 32) = (VSS(4) + VSS(6)) *. 5
        DD(NE, 33) = (VSC(4) + VSC(6)) *.5
C
        DD(NE, 34) = DD(NE, 4) * 2.
        DD(NE, 35) = DD( NE, 15)
        DD(NE, 36) = DD( NE, 26)
        DD(NE, 37) = VO - VSC(4) *. 5
        DD(NE, 38) = VSS(4) *. 5
        DD(NE, 39) = (VSC(1) - VSC(5)) *. 5
        DD(NE, 40) = ( - VSS(1) + VSS(5)) *. 5
        DD(NE, 41) = (VSC(2) - VSC(6)) *. 5
        DD(NE, 42) = ( - VSS(2) + VSS(6)) *. 5
        DD(NE, 43) = (VSC(3) - VSC(7)) *. 5
        DD(NE, 44) = ( - VSS(3) + VSS(7)) *. 5
C
        DD(NE, 45) = DD(NE,5) * 2.
        DD(NE,  46) = DD(NE, 16)
        DD(NE, 47) = DD(NE,27) 
        DD(NE, 48) = DD(NE,38)
        DD(NE, 49) = VO + VSC(4) *. 5
        DD(NE, 50) = (VSS(1) + VSS(5)) *. 5
        DD( NE, 51) = (VSC(1) + VSC(5)) *. 5 
        DD(NE, 52) = (VSS(2) + VSS(6)) *. 5
        DD( NE, 53) = (VSC(2) + VSC(6)) *. 5 
        DD(NE, 54) = (VSS(3) + VSS(7)) *. 5
        DD( NE, 55) = (VSC(3) + VSC(7)) *. 5 
C
        DD(NE, 56) = DD(NE, 6) *2.
        DD(NE, 57) = DD(NE, 17)
        DD(NE, 58) = DD(NE, 28)
        DD(NE, 59) = DD(NE, 39)
        DD(NE, 60) = DD(NE, 50)
        DD(NE, 61) = VO - VSC(6) *. 5
        DD(NE, 62) = VSS(6) *. 5
        DD(NE, 63) = (VSC(1) - VSC(7)) *. 5
        DD(NE, 64) = ( - VSS( 1) + VSS(7)) *. 5
        DD(NE, 65) = ( VSC(2) - VSC(8)) *. 5
        DD(NE, 66) = (- VSS(2) + VSS(8)) *. 5
C
        DD(NE, 67) = DD(NE, 7) *2.
        DD(NE, 68) = DD(NE, 18)
        DD(NE, 69) = DD( NE, 29)
        DD(NE, 70) = DD(NE, 40)
        DD(NE, 71) = DD(NE, 51)
        DD(NE, 72) = DD(NE, 62)
        DD(NE, 73) = VO + VSC(6) *. 5
        DD(NE, 74) = (VSS(1) + VSS(7))*. 5
        DD(NE, 75) = (VSC(1) + VSC(7))*. 5
        DD(NE, 76) = (VSS(2) + VSS(8))*. 5
        DD(NE, 77) = (VSC(2) + VSC(8))*. 5
C
        DD(NE, 78) = DD(NE, 8) *2. 
        DD(NE, 79) = DD(NE, 19)
        DD(NE, 80) = DD(NE, 30)
        DD(NE, 81) = DD(NE, 41)
        DD(NE, 82) = DD(NE, 52)
        DD(NE, 83) = DD(NE, 63)
        DD(NE, 84) = DD(NE, 74)
        DD(NE, 85) = VO - VSC(8)*.5
        DD(NE, 86) = VSS(8)*.5
        DD(NE, 87) = (VSC(1) - VSC(9))*.5
        DD(NE, 88) = (-VSS(1) + VSS(9))*.5
C
        DD(NE, 89) = DD(NE, 9)*2.
        DD(NE, 90) = DD(NE, 20)
        DD(NE, 91) = DD(NE, 31)
        DD(NE, 92) = DD(NE, 42)
        DD(NE, 93) = DD(NE, 53)
        DD(NE, 94) = DD(NE, 64)
        DD(NE, 95) = DD(NE, 75)
        DD(NE, 96) = DD(NE, 86)
        DD(NE, 97) = VO + VSC(8)*.5
        DD(NE, 98) = (VSS(1) + VSS(9))*.5
        DD(NE, 99) = (VSC(1) + VSC(9))*.5
C
       DD(NE, 100) = DD(NE, 10)*2.
       DD(NE, 101) = DD(NE, 21)
       DD(NE, 102) = DD(NE, 32)
       DD(NE, 103) = DD(NE, 43)
       DD(NE, 104) = DD(NE, 54)
       DD(NE, 105) = DD(NE, 65)
       DD(NE, 106) = DD(NE, 76)
       DD(NE, 107) = DD(NE, 87)
       DD(NE, 108) = DD(NE, 98)
       DD(NE, 109) = VO - VSC(10)*.5
       DD(NE, 110) = VSS(10) *.5
C
       DD(NE, 111) = DD(NE, 11)*2.
       DD(NE, 112) = DD(NE, 22)
       DD(NE, 113) = DD(NE, 33)
       DD(NE, 114) = DD(NE, 44)
       DD(NE, 115) = DD(NE, 55)
       DD(NE, 116) = DD(NE, 66)
       DD(NE, 117) = DD(NE, 77)
       DD(NE, 118) = DD(NE, 88)
       DD(NE, 119) = DD(NE, 99)
       DD(NE, 120) = DD(NE, 110)
       DD(NE, 121) = VO + VSC(10)*.5
C
 9999  CONTINUE
   10  CONTINUE
       PRINT '(10X,"NDEG=11")'
       PRINT '(10X,"DMATC")'
       RETURN
       END
:::::::::::::::::::::::::::::::::::::::::::::::::
C::::::::::::::::: DMAT?????????????::::::::::::::
       SUBROUTINE   DMATH 
           PARAMETER  (NA1=1995, NA2=249, NA3=755, NA4=744, NA5=400, NA6=121)
           COMMON     AAA, BBB, NONC, CCC, NONC2, NCOIL, NHOWA,
         &                      NPO1, NOM, NPO2, NPOR, NELEM, NB, NDEG, NDE, DF, 
         &                      DH(NA1, NA2), DK(NA1), DAA(NA1), DA(NA1), 
         &                 DS(NA3, 18), DCS(NA3), DB(NA3, 22), DD(NA3, NA6), NOD(0:NA4,4),
         &                 XY(NA5, 2), DENRYU(1, 11), DC(NA3, 11), DCPRE(NA3, 11),
         &                 DMUO, DOMEG, DPI, ITR, TOTAL, 
         &                 DN(11, 11), DBH(11, 4)
            DIMENSION    AB(0: 72), VSC(4), VSS(4)
C
            DO  10  NE=1, NELEM
               IF (INT(NOD(NE,4)/100).EQ.2)     GO  TO  1000
               DO 12  K=1,25
                  DD(NE,K) = 0.
      12  CONTINUE
          DD(NE,1) = 1.0 /DMUO
          DD(NE,7) = DD(NE, 1)
          DD(NE,13)= DD(NE,1)
          DD(NE,19)= DD(NE,1)
          DD(NE,25)= DD(NE,1)
          GO  TO  9999
000       CONTINUE
          HD=2.*DPI/72.
          CO=1./DPI
          DO  30  NT = 0, 72
              DT= FLOAT (NT)
              SIN1 = SIN(HD*DT)
              COS1 = COS(HD*DT)
              SIN2 = SIN(2.*HD*DT)
              COS2 = COS(2.*HD*DT)
              BX = DB(NE,1) + DB(NE,2)*SIN1+ DB(NE,3)*COS1+ DB(NE,4)*SIN2
     &            + DB(NE,5)*COS2
              BY = DB(NE,6)+DB(NE,7)*SIN1+DB(NE,8)*COS1+ DB(NE,9)*SIN2
     &            + DB(NE, 10)* COS2
              BT = SQRT(BX*BX+BY*BY)
              DBXT=(DB(NE,2)*COS1 - 1.*DB(NE,3)*SIN1+ 2.*DB(NE,4)*COS2
     &            - 2.*DB(NE,5)*SIN2)*DOMEG
              DBYT=(DB(NE,7)*COS1 - 1.*DB(NE, 8)*SIN1+2.*DB(NE,9)*COS2
     &            - 2.*DB(NE,10)*SIN2)*DOMEG
              IF (BT.NE.0.) THEN
              AB(NT) = DBH(1,1) +DBH(1,2)*(BX*DBXT+BY*DBYT) / (BT*BT)           
              ELSE
              AB(NT) = AB(NT - 1)
              END IF
              IF (AB(NT).GT. 7000.)  AB(NT) =7000.
              IF (AB(NT).LT.-7000.)  AB(NT) =-7000.
C
   30       CONTINUE
            SO = AB(0)
            DO   40  NT=1,71, 2
                 SO =SO+4.*AB(NT) +2.*AB(NT+1)
   40       CONTINUE
            VO=(S0-AB(72))*HD/3.0*CO*0.5
C-------------------------------------------------------
C           IF(NE.EQ.48)  THEN
C           PRINT '(1X,"VO=", E14.6)',VO
C           END IF 
C--------------------------------------------------------
            DO   50   N=1,4
               DRN= FLOAT(N)
               SC=AB(0)
               DO 60 NT=1, 71, 2
                 DT=FLOAT(NT)
                 YC=AB(NT)*COS(DRN*HD*DT)
                 SC= SC+4.*YC
                 YC= AB(NT+1)*COS(DRN*HD*(NT+1.))
                 SC= SC+2.*YC
  60          CONTINUE
              VSC(N)=(SC - YC) * HD/3.*CO
C             SS= AB(0)
              SS=0
              DO  70  NT=1, 71, 2
                  DT=FLOAT(NT)
                  YS = AB(NT)*SIN(DRN*HD*DT)
                  SS=SS+4.*YS
                  YS=AB(NT+1)*SIN(DRN*HD*(DT+1.))
                  SS=SS+2.*YS
  70          CONTINUE
              VSS(N)=(SS-YS)*HD/3.0*CO
C-----------------------------------------------------------
C        IF (NE.EQ.48)  THEN
C        PRINT '(1X, "VSC(", I1, ") = ",E14.6)',N,VSC(N)
C        PRINT '(1X, "VSS(", I1, ") = ",E14.6)',N,VSS(N)
C        END IF
c---------------------------------------------------------
 50         CONTINUE
         DD(NE,1) = VO
         DD(NE,2) = VSS(1)*.5
         DD(NE,3) = VSC(1)*.5
         DD(NE,4) = VSS(2)*.5
         DD(NE,5) = VSC(2)*.5
C
         DD(NE,6) = DD(NE,2)*2.
         DD(NE,7) = VO - VSC(2)*.5
         DD(NE,8) = VSS(2)*.5
         DD(NE,9) = (VSC(1) - VSC(3))*.5
         DD(NE,10) = (-VSS(1) + VSS(3))*.5
C
         DD(NE,11) = DD(NE,3)*2.
         DD(NE,12) = DD(NE,8)
         DD(NE,13) = VO+VSC(2)*.5
         DD(NE,14) = (VSS(1) + VSS(3))*.5
         DD(NE,15) = (VSC(1) + VSC(3))*.5
C
         DD(NE,16) = DD(NE,4)*2.
         DD(NE,17) = DD(NE,9)
         DD(NE,18) = DD(NE,14)
         DD(NE,19) = VO - VSC(4)*.5
         DD(NE,20) = VSS(4)*.5
C
         DD(NE,21) = DD(NE,5)*2.
         DD(NE,22) = DD(NE,10)
         DD(NE,23) = DD(NE,15) 
         DD(NE,24) = DD(NE,20)
         DD(NE,25) = VO + VSC(4)*.5 
 9999    CONTINUE
   10    CONTINUE
         PRINT '(10X,"NDEG=5")'
         PRINT '(10X,"DMATH")'
         RETURN
         END