C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C:::::::::::::::  GAUSS Elimination Method    VER .2  :::::::::::::
C
      SUBROUTINE  GAUSS
      PARAMETER (NA1=1995,NA2=249,NA3=755,NA4=755,NA5=400,NA6=121)
      COMMON  AAA, BBB, NONC, CCC, NONC2, NCOIL, NHOWA,
     &        NPO1, NOM, NPO2, NPO3, NPOR, NELEM, NB, NDEG, NDE, DF,
     &        DH(NA1,NA2), DK(NA1), DAA(NA1), DA(NA1),
     &        DS(NA3,18), DCS(NA3), DB(NA3,22), DD(NA3,NA6), NOD(0:NA4,4),
     &        XY(NA5,2), DENRYU(1,11), DC(NA3,11), DCPRE(NA3,11),
     &        DMUO, DOMEG, DPI, ITR, TOTAL,
     &        DN(11,11), DBH(11,4)
C
      INTEGER  NM
C******************************************************************
C     DH(I,J) :  MATRIX          (NOB, NB)
C     DK(I)   :  FORCED VECTOR   (NOB)
C     DA(I)   :  POTENTIAL       (NOB)
C     NOM     :  NUMBER OF UNKNOW POTENTIAL
C     NB      :  BAND  WIDTH
C******************* ELIMINATE L OF MATRIX DH(I,J) ****************
C
C  ---------------- Double precision calculation ------------------
C     DOUBLE  DT1, DT
C
C  ------- Sequence number of centre point of matrix DH() ---------
C     NDIA = (NB + 1) / 2
C     NDIA = NB
C
      DO 100  NN=1, NOM-1
        DT1=1.0 / DH(NN,NDIA)
C ---------------------------------- (NM, NN)
        DO 200  NM=NN+1, NOM
          J2 = NDIA - NM + NN
          I2 = NM
          IF (J2 .LT. 1)  GOTO  100
          DT = DT1 * DH(I2,J2)
C;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
C
          IF (DT .EQ. 0.0)  GOTO  200
C
          I1 = NM
          DK(I1) = DK(I1) - DT * DK(NN)
C ---------------------------------- (NM, NK)
          DO 300  NK=NN, NOM
            II1 = NM
            JJ1= NDIA - NM + NK
            IF (JJ1 .GT. 2*NB-1)  GOTO  200
            II2 = NN
            JJ2 = NDIA -NN + NK
            IF (JJ2 .GT. 2*NB-1)  GOTO  200
C
            DH(II1,JJ1) = DH(II1,JJ1) - DT * DH(II2,JJ2)
C ----------------------------------- NK 
  300     CONTINUE
C ----------------------------------- NM
  200   CONTINUE]
C ----------------------------------- NN
  100 CONTINUE
C
C
      DA(NOM) = DK(NOM) / DH(NOM,NDIA)
C
      DO 400  NN=2, NOM
        NM = NOM - NN + 1
C
        DO 500  NN=2, NOM
          I1 = NM
          J1 = NDIA - NM + NK
          IF (J1 .GT. 2*NB-1)  GOTO  410
          DK(I1) = DK(I1) - DH(I1,J1) * DA(NK)
C ----------------------------------- NK
  500   CONTINUE
  410   IF (DH(NM,NDIA) .EQ. 0.0)  GOTO 600
        DA(NM) = DK(NM) / DH(NM,NDIA)
C ----------------------------------- NN
  400 CONTINUE
      RETURN
  600 CONTINUE
      WRITE (6,9900)  NM
      PRINT '(5X," MATRIX  ERROR AT GAUSS ",I4)', NM
 9900 FORMAT (1H, 2OH MATRIX ERROR , I4)
      STOP
      END