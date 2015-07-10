C:::   RDATA  Read coefficients and data for for vectors ::::::::::
      SUBROUTINE RDATA
	  
      implicit real*8 (A-h,o-z) integer (i-n)
	  
      PARAMETER (NA1=1995,NA2=249,NA3=755,NA4=755,NA5=400,NA6=121)
      COMMON AAA, BBB, NONC, CCC, NONC2, NCOIL, NHOWA,  
     &       NPO1, NOM, NPO2, NPO3, NPOR, NELEM, NB, NDEG, NDE, NF,
     &       DH(NA1,NA2), DK(NA1), DAA(NA1), DA(NA1),
     &       DS(NA3,18), DCS(NA3), DB(NA3,22), DD(NA3,NA6), NOD(0:NA4,4),
     &       XY(NA5,2), DENTYU(1,11), DC(NA3,11), DCPRE(NA3,11),
     &       DMUO, DONEG, DPI, ITR, TOTAL, 
     &       DN(11,11), DBH(11,4)
      INTEGER  TOTAL
C
      CHARACTER  CHAR*80
      READ(15,9900) CHAR
      READ(15,9900) CHAR
 9900 FORMAT(A80)
      READ(15,*) NPO1, NPO2, NPO3, NELEM, NB, NCOIL
      NB=NDEG*NB
      NOM=NDEG*NPO1
C-----Read the vertex coordinates
      READ(15,9900) CHAR
      DO I= 1, NPO3, 2
        READ(15,*) XY(I,1), XY(I,2), XY(I+1,1), XY(I+1,2)
      END DO
C-----Change units of X-Y coordinate
      DO  I=1, NPO3
        DO   J=1, 2
          XY(I,J) = XY(I,J) * .001
        END DO   
      END DO
C-----Read the nodes 
      READ(15,9900) CHAR
      DO I=1, NELEM, 4
        READ(15,*) (NOD(I+J-1, 1), NOD(I+J-1, 2),
     &              NOD (I+J-1, 3), NOD(I+J-1, 4), J=1,4)
      ENDDO
      DO I=1, NELEM
        NOD(I,1) = NOD(I,1)-1
        NOD(I,2) = NOD(I,2)-1
        NOD(I,3) = NOD(I,3)-1
      ENDDO
C-----Read the DBH?
      READ(15,9900) CHAR
      READ(15,*) (DBH(1,J), J=1,2)
      PRINT *
      PRINT '(1P,(2E15.4))', (DBH (1,J), J=1,2)
      PRINT *
C---  Current density source (J) ---
      READ(15,9900) CHAR
      READ(15,9900) NDE
      DO 100  KK=1, NDE
        READ(15,*)  (DENRYU(KK,J), J=1, NDEG)
        IF (NDEG .EQ. 5) THEN 
          PRINT '(1P(3E15.4))', (DENRYU(KK,J), J=1,3)
          PRINT '(1P(2E15.4))', (DENRYU(KK,J), J=4,5)
          ELSE IF (NDEG .EQ. 7) THEN 
            PRINT '(1P(3E15.4))', (DENRYU(KK,J), J=1,3)
            PRINT '(1P(3E15.4))', (DENRYU(KK,J), J=4,7)
          ELSE IF (NDEG .EQ. 11) THEN
            PRINT '(1P(3E15.4))', (DENRYU(KK,J), J=1,3)
            PRINT '(1P(4E15.4))', (DENRYU(KK,J), J=4,7)
            PRINT '(1P(4E15.4))', (DENRYU(KK,J), J=8,11)
        END IF 
  100 CONTINUE
C----- Read data into matrix DAA -----
      IF (NHOWA .EQ. 1)  THEN
        PRINT ' (10X, "Set initial value of DAA to zero")'
        DO 900 I=1, NOM
          DAA(I)=0.0
  900   CONTINUE
        TOTAL =0
      ELSE
        PRINT '(10X, "Set initial value of DAA from file")'
        READ (19,9900) CHAR
        READ (19,9900) CHAR  
        READ (19,9900) CHAR
        READ (19,9900) CHAR
        READ (19,9900) CHAR
        READ (19,*)  TOTAL 
        PRINT '(20X, "TOTAL=",I6)',TOTAL 
        READ (19,9900) CHAR
        DO 200 K=1, NOM-NDEG+1,NDEG
          IF (NDEG .EQ. 5) THEN
            READ (19,*) NDAMY, DAA(K), DAA(K+1), DAA(K+2)
            READ (19,*) DAA(K+3), DAA(K+4)
          ELSE IF (NDEG .EQ. 7) THEN 
            READ (19,*) NDAMY, DAA(K), DAA(K+1), DAA(K+2)
            READ (19,*) DAA(K+3), DAA(K+4), DAA(K+5), DAA(K+6)
          ELSE IF (NDEG. EQ. 11) THEN
            READ (19,*) NDAMY, DAA(K), DAA(K+1), DAA(K+2)
            READ (19,*) DAA(K+3), DAA(K+4), DAA(K+5), DAA(K+6) 
            READ (19,*) DAA(K+7), DAA(K+8), DAA(K+9), DAA(K+10)
          END IF 
  200   CONTINUE 
      END IF 
      RETURN
      END 
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C:::: Set initial values ::::
      SUBROUTINE  PRI
      PARAMETER (NA1=1995,NA2=249,NA3=755,NA4=755,NA5=400, NA6=121)
      COMMON  AAA, BBB, NONC, CCC, NONC2, NCOILL, NHOWA,
     &        NPO1, NOM, NPO2, NPO3, NPOR, NELEM, NB, NDEG, NDE, DF,
     &        DH(NA1,NA2), DK(NA1), DAA(NA1), DA(NA1),
     &        DS(NA3,18), DCS(NA3), DB(NA3,22), DD(NA3,NA6), NOD(0:NA4,4),
     &        XY(NA5,2), DENRYU(1,11), DC(NA3,11), DCPRE(NA3,11),
     &        DMUO, DOMEG, DPI, ITR, TOTAL,
     &        DN(11,11), DBH(11,4)
C----------------------------------------------Initialization
      DO  I=1, NA1
        DO  J=1, NA2
          DH(I,J)=0.0
        END DO
      ENDDO
   
      DO  I=1, NA3
        DO  J=1, NA6
          DD(I,J)=0
        END DO
        DO  K=1, NDEG
          DCPRE(I,K)=0.0
        END DO
      ENDDO
C.... Convergence points
      NONC=0
      NONC2=0
C--- Set values to matrix DN ---
C--- Harmonic Matrix 'N'
C---  |0 -1             .....|
C---  |1  0                  |
C---  |     0 -3             |
C---  |     3  0        .....|
C---  |          0 -5        |
C---  |          5  0   .....|
C---  |               .......|
      DO  J=1, NDEG
        DO  K=1, NDEG
          DN(J,K)=0.0
        END DO
      ENDDO
      DN(3,2)=1.0
      DN(2,3)= -1.0
      DN(5,4)=2.0
      DN(4,5)= -2.0
      IF (NDEG .EQ. 5) GOTO 90
      DN(7,6)= 3.0
      DN(6,7)= -3.0
      IF (NDEG .EQ. 7) GOTO 90
      DN(9,8)=4.0
      DN(8,9)= -4.0
      DN(11,10)=5.0
      DN(10,11)= -5.0
   90 CONTINUE
      RETURN
      END