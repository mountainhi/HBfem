C      HBFEM.FORT77(FEMXY99B)   9/4       FOR ELEMENT=
C 
C  1100=11*100, (100=UNKNOWNPOINT), 600=(NB*2-1)*NDEG+4,
C  303=ELEMENT+3,  200=ALLPOINT+1,
C  17=3+2+1+11, 22=11+11, 121=11**2
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C::::::::::::::::::::  MAINPROGRAM  :::::::::::::::::::::::::::::::
             PROGRAM HBFEM_2D
C PARAMETER (NA1=1100,NA2=600,NA3=303,NA4=303,NA5=200,NA6=121)
      implicit real*8 (A-h,o-z) integer (i-n)
	   
      PARAMETER (NA1=1995,NA2=249,NA3=755,NA4=755,NA5=400,NA6=121)
      COMMON AAA,BBB,NONC,CCC,NONC2,NCOIL,NHOWA,
     &       NPO1,NOM,NPO2,NPO3,NPOR,NELEM,NB,NDEG,NDE,DF,
     &       DH(NA1,NA2),DK(NA1),DAA(NA1),DA(NA1), 
     &       DS(NA3,18),DCS(NA3),DB(NA3,22),DD(NA3,NA6),NDD(O:NA4,4),
     &       XY(NA5,2),DENRYU(1,11),DC(NA3,11),DCPRE(NA3,11),
     &       DMUO,DOMEG,DPI,ITR,TOTAL,
     &       DN(11,11),DBH(11,4)
      INTEGER   TOTAL
C
C     MAX=80
      MAX=40
C     DF=10.
      DF=1.0363E6
      DF=1.E6
      NDEG=5
C          5...COMPUTE TO 2ND HARMONIC,  7...TO 3RD,  11...TO 5TH 
C      NHOWA=2
      NHOWA = 1  
C          1... DAA = 0.   2... INPUT FROM FILE
      AAA=0.005
      BBB=0.1
      CCC=0.005
      DPI=3.141593
      DMUO=4.E-7*DPI
      DOMEG=2.*DPI*DF
C
      CALL  RDATA
      CALL PRI
      ITR=1
C     CALL SMAT
      CALL SMATP
  100 CONTINUE
      IF (ITR.GT.MAX)  THEN
        PRINT  '(45X, " STOP AT MAIN ") '
        GO TO  9999
      END IF
      PRINT  '(5X,I3,A)', ITR+TOTAL, 'INTERATION TIMES'
	  
      CALL CMATP
C     CALL CMAT2
C     CALL  FLUX
      CALL FLUXP
      IF (NDEG. EQ. 5) THEN
      CALL DMATA
C     CALL CMATH
      ELSE IF (NDEG. EQ. 7 ) THEN  
            CALL  DMATB
      END IF
	  
      CALL  FEMP
C     CALL  FEM
      CALL  GAUSS
      CALL  CONV
C     CALL  CONV2
      ITR = ITR +1
      IF (ITR .EQ. 2)  GOTO 100
      IF ((NONC .NE. O) .OR. (NONC2 .NE. O)) GOTO 100
      IF ((NONC .EQ. O) .AND. (NONC2 .EQ.O)) THEN
        PRINT '(10X,"*********************************")'
        PRINT '(10X,"***        CONVERGENT         ***")'
        PRINT '(10X,"*********************************")'
      END IF
	  
 9999 CONTINUE
      PRINT '(10X, "POST")'
      CALL  POST
      PRINT '(10X,"ODATA")'
      CALL  ODATA
      STOP
	  
      END