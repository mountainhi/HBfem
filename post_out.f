C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::::::::::::::: POST  Calculation values ::::::::::::::::::::::::
      SUBROUTINE  POST
      PARAMETER (NA1=1995,NA2=249,NA3=755,NA4=755,NA5=400,NA6=121)
      COMMON  AAA, BBB, NONC, CCC, NONC2, NCOIL, NHOWA,
     &        NPO1, NOM, NPO2, NPO3, NPOR, NELEM, NB, NDEG, NDE, DF,
     &        DH(NA1,NA2), DK(NA1), DAA(NA1), DA(NA1),
     &        DS(A3,18), DCS(NA3), DB(NA3,22), DD(NA3,NA6), NOD(0:NA4,4),
     &        XY(NA5,2), DENRYU(1,11), DC(NA3,111), DCPRE(NA3,11),
     &        DMUO, NOMEG, NPI, ITR, TOTAL,
     &        DN(11,11), DBH(11,4)
      CALL  FLUX
      RETURN
      END
C 
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::::::::::::: ODATA   Storage of output data ::::::::::::::::::::
C
      SUBROUTINE  ODATA
      PARAMETER (NA1=1995,NA2=249,NA3=755,NA4=755,NA5=400,NA6=121)
      COMMON  AAA, BBB, NONC, CCC, NONC2, NCOIL, NHOWA,
     &        NPO1, NOM, NPO2, NPO3,, NPOR, NELEM, NB, NDEG, NDE, DF,
     &        DH(NA1,NA2), DK(NA1), DAA(NA1), DA(NA1),
     &        DS(NA3,18), DCS(NA3), DB(NA3,22), DD(NA3,NA6), NOD(0:NA4,4), 
     &        XY(NA5,2), DENRYU(1,11), DC(NA3,11), DCPRE(NA3,11), 
     &        DMUO, NOMEG, DPI, ITR, TOTAL,
     &        DN(11,11), DBH(11,4)
C
      INTEGER  TOTAL
C
C *************** Double precision calculation ********************
C     IMPLICIT   REAL * 8(D)
C
C --------- WRITE TOTAL ----------
C
C     WRITE(11,'(A)') 'DATA FEM163 THREE-LEGGED CORE WITH SHADING'
      IF (NDEG .EQ. 5)  THEN
        WRITE(11, '(1P(3E15.4))')  (DENRYU(1,J), J=1,3)
        WRITE(11, '(1P(2E15.4))')  (DENRYU(1,J), J=4,5)
        WRITE(11, '(A)')  'DATA'
      ELSE IF (NDEG .EQ. 7)  THEN
        WRITE(11, '(1P(3E15.4))')  (DENRYU(1,J), J=1, 3)
        WRITE(11, '(1P(4E15.4))')  (DENRYU(1,J), J=4, 7)
        WRITE(11, '(A)')  'DATA'
      ELSE IF (NDEG .EQ. 11)  THEN
        WRITE(11, '(1P(3E15.4))')  (DENRYU(1,J), J=1, 3)
        WRITE(11, '(1P(4E15.4))')  (DENRYU(1,J), J=4, 7)
        WRITE(11, '(1P(4E15.4))')  (DENRYU(1,J), J=8, 11)
      END IF
      WRITE(11, '(1P(E15.4))')  (DF)
      WRITE(11, '(16)')  TOTAL + ITR - 1
      WRITE(11, '(16)')  TOTAL + ITR - 1
C
C--------------------- Write matrix DA() --------------------------
C    
      WRITE(11,'(A)') 'VALUES OF VECTOR POTENTIAL  0 SIN1 COS1 SIN3'
      IF (NDEG .EQ. 5)  THEN
        DO 10  N=1, NOM, NDEG
          WRITE(11, '(I4, 3E14.5)')  N, DAA(N), DAA(N+1), DAA(N+2)
          WRITE(11, '(4X, 2E14.5)')  DAA(N+3), DAA(N+4)
   10   CONTINUE
      ELSE IF (NDEG .EQ. 7)  THEN
        DO 20  N=1, NOM, NDEG
          WRITE(11, '(I4, 3E14.5)') N, DAA(N), DAA(N+1), DAA(N+2)
          WRITE(11,' (4X, 4E14.5)') DAA(N+3), DAA(N+4), DAA(N+5), DAA(N+6)
   20   CONTINUE
      ELSE IF (NDEG .EQ. 11)  THEN
        DO 30  N=1, NOM, NDEG
          WRITE(11, '(I4, 3E14.5)') N, DAA(N), DAA(N+1), DAA(N+2)
          WRITE(11, '(4X, 4E14.5)') DAA(N+3), DAA(N+4), DAA(N+5), DAA(N+6)
          WRITE(11, '(4X, 4E14.5)') DAA(N+7), DAA(N+8), DAA(N+9), DAA(N+10)
   30   CONTINUE
      END IF
C
C-------------------- Write matrix DB() --------------------------
      WRITE(11,'(A)') 'Values of flux density DB()  0 SIN1-X COS1-X SIN3-X...'
      WRITE(11,'(A)') '                             0 SIN1-Y COS1-Y SIN3-Y...'
      IF (NDEG .EQ. 5)  THEN
        DO 100  I=1, NELEM
          WRITE(11, '(I4, 3E14.5)')  I, (DB(I,K), K=1, 3)
          WRITE(11, '(4X, 4E14.5)')  (DB(I,K), K=4, 5)
          WRITE(11, '(4X, 4E14.5)')  (DB(I,K), K=6, 8)
          WRITE(11, '(4X, 3E14.5)')  (DB(I,K), K=9, 10)
  100   CONTINUE
      ELSE IF (NDEG .EQ. 7)  THEN
        DO 110  I=1, NELEM
          WRITE(11, '(I4, 3E14.5)')  I, (DB(I,K), K=1, 3)
          WRITE(11, '(4X, 4E14.5)')  (DB(I,K), K=4, 7)
          WRITE(11, '(4X, 4E14.5)')  (DB(I,K), K=8, 10)
          WRITE(11, '(4X, 4E14.5)')  (DB(I,K), K=11, 14)
  110   CONTINUE
      ELSE IF (NDEG .EQ. 11)  THEN
        DO 120  I=1, NELEM
          WRITE(11, '(I4, 3E14.5)')  I, (DB(I,K), K=1, 3)
          WRITE(11, '(4X, 4E14.5)')  (DB(I,K), K=4, 7)
          WRITE(11, '(4X, 4E14.5)')  (DB(I,K), K=8, 11)
          WRITE(11, '(3X, 4E14.5)')  (DB(I,K), K=12, 14)
          WRITE(11, '(4X, 4E14.5)')  (DB(I,K), K=15, 18)
          WRITE(11, '(4X, 4E14.5)')  (DB(I,K), K=19, 22)
  120   CONTINUE
      END IF
      RETURN
      END