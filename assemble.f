C:::: SMAT Computing of coefficients of elements in matrix DS :::::
      SUBROUTINE SMAT

      use param 
      use matvec
      
!      PARAMETER (NA1=1995,NA2=249,NA3=755,NA4=755,NA5=400,NA6=121)
!      COMMON  AAA, BBB, NONC, CCC, NONC2, NCOIL, NHOWA,
!     &        NPO1, NOM, NPO2, NPO3, NPOR, NELEM, NB,NDEG, NDE, DF,
!     &        DH(NA1,NA2), DK(NA1), DAA(NA1), DA(NA1),
!     &        DS(NA3,18), DCS(NA3), DB(NA3,22), DD(NA3,NA6), NDD(0:NA4,4),
!     &        XY(NA5,2), DENRYU(1,11), DC(NA3,11), DCPRE(NA3,11),
!     &        DMUO, DOMEG, DPI, ITR, TOTAL,
!     &        DN(11,11), DBH(11,4)
      real(8), DIMENSION(3)::  DX, DY, DR , DQ
C
      DO 10  I=1, NELEM
        DO  J=1, 3
          DX(J) = XY(NOD(I,J)+1,1)
          DY(J) = XY(NOD(I,J)+1,2)
        ENDDO
        DCS(I) = (DX(2)*DY(3)+DX(1)*DY(2)+DX(3)*DY(1)
     &            -DX(2)*DY(1)-DX(1)*DY(3)-DX(3)*DY(2))*0.5D0
C-------------------------------------------------------------DEBUG
        IF (DCS(I) .EQ. 0.0) THEN 
          PRINT *, 'area =0   ',I
          PRINT *, DX(1), DX(2), DX(3), DY(1), DY(2), DY(3)
          STOP
        END IF
C---------------------------------------------------------------------
C------    b (DQ) and c(DR) in the interpolating function
        DO J=1, 3
          DQ(J) = DY(MOD(J,3)+1) - DY(MOD(J+1,3)+1)
          DR(J) = DX(MOD(J+1,3)+1) - DX(MOD(J, 3)+1)
        ENDDO
   
        N=1
        DO   J=1, 3
          DO   K=J, 3
            DS(I,N) = 0.25/DCS(I)*(DR(J)*DR(K)+DQ(J)* DQ(K))
            N=N+1
          ENDDO
        ENDDO
c...    Current density source
        IF (INT(NOD(I,4)/100) .EQ. 3) THEN
          DC1 = DENRYU(1,1)
          DC2 = DENRYU(1,2)
          DC3 = DENRYU(1,3)
          DC4 = DENRYU(1,4)
          DC5 = DENRYU(1,5)
          IF (NDEG .EQ. 5) GOTO 60
          DC6 = DENRYU(1,6)
          DC7 = DENRYU(1,7)
          IF (NDEG .EQ. 7) GOTO 60
          DC8 = DENRYU(1,8)
          DC9 = DENRYU(1,9)
          DC10 = DENRYU(1,10)
          DC11 = DENRYU(1,11)
   60     CONTINUE
        ELSE
          DC1=0.0
          DC2=0.0
          DC3=0.0
          DC4=0.0
          DC5=0.0
          IF (NDEG .EQ. 5) GOTO 70
          DC6=0.0
          DC7=0.0
          IF (NDEG .EQ. 7) GOTO 70
          DC8=0.0
          DC9=0.0
          DC10=0.0
          DC11=0.0
   70     CONTINUE
        END IF 
C-------------  Calculation of current density  -------------
C---    K = delta * {J1s J1c .....} /3
        DCON = DCS(I) * 0.3333333333333333
        DS(I,7) = DC1 * DCON
        DS(I,8) = DC2 * DCON
        DS(I,9) = DC3 * DCON
        DS(I,10) = DC4 * DCON
        DS(I,11) = DC5 * DCON
        IF (NDEG .EQ. 5) GOTO 80
        DS(I,12) = DC6 * DCON
        DS(I,13) = DC7 * DCON
        IF (NDEG .EQ. 7) GOTO 80
        DS(I,14) = DC8 * DCON
        DS(I,15) = DC9 * DCON
        DS(I,16) = DC10 * DCON
        DS (I,17) = DC11 * DCON
   80   CONTINUE
C--------    Computing coefficients of eddy current   --------
        IF (INT(NOD(I,4)/100) .EQ. 3) THEN
C         DS(I,18) = DBH(NOD(I,4)-199, 4)*DOMEG*DCS(I)/12.0
          DS(I,18) = 5.92E7*DOMEG*DCS(I)/12.0
        ELSE
          DS(I,18) = 0.0
        END IF 
   10 CONTINUE
      RETURN
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C:: SMATP  Computing coefficients in Matrix DS for axially symmetric
      SUBROUTINE SMATP
      use param
      use matvec
C       PARAMETER (NA1=1995,NA2=249,NA3=755,NA4=755,NA5=400,NA6=121)
C      COMMON  AAA, BBB, NONC, CCC, NONC2, NCOIL, NHOWA,
C     &        NPO1, NOM, NPO2, NPO3, NPOR, NELEM, NB, NDEG, NDE, DF, 
C     &        DH(NA1, NA2), DK(NA1), DAA(NA1), DA(NA1),
C     &        DS(NA3,18), DCS(NA3), DB(NA3,22), DD(NA3,NA6), NOD(0:NA4,4),
C     &        XY(NA5,2), DENTYU(1,11), DC(NA3,11), DCPRE(NA3,11), 
C     &        DNUO, DOMEG, DPI, ITR, TOTAL,
C     &        DN(11,11), DBH(11,4)

      real(8), DIMENSION(3)::  DX, DY, DR, DQ
C
      DO 10  I=1, NELEM
	  
        DO  J=1, 3
           DX(J) = XY(NOD(I,J)+1,1)
           DY(J) = XY(NOD(I,J)+1,2)
        END DO
C.......Compute the area of elements
        DCS(I) = ( DX(2)*DY(3)+DX(1)*DY(2)+DX(3)*DY(1)
     &            - DX(2)*DY(1) - DX(1)*DY(3) - DX(3)*DY(2) )*.5D0
C-------------------------------------------------------------DEBUG
        IF (DCS(I) .EQ. 0.0) THEN
          PRINT *, 'Area = 0  ',I
          PRINT *, DX(1), DX(2), DX(3), DY(1), DY(2), DY(3)
          STOP
        END IF
C-------------------------------------------------------------------------------------------------
C------    b (DQ) and c(DR) in the interpolating function
        DO  J=1, 3
          DQ(J) = DY(MOD(J,3)+1) - DY(MOD(J+1,3)+1)
          DR(J) = DX(MOD(J+1,3)+1) - DX(MOD(J,3)+1)
        ENDDO
C-----------   Symmetric field 1   ------------
        DDX = (DX(1)+DX(2)+DX(3))/3.0
C
C.......  1/4*delta*[S]
        N=1
        DO   J=1, 3
          DO  K=J, 3
            DS(I,N) = 0.25/DCS(I)*(DR(J)*DR(K)+DQ(J)*DQ(K))*DDX
     &                +(DQ(J)+DQ(K))/6.0+DCS(K)/9.0/DDX
            N=N+1
          ENDDO
        ENDDO
C--------------------------------------------------------------------------------------
        IF (INT(NOD(I,4)/100) .EQ. 3) THEN
          DC1 = DENRYU(1,1)
          DC2 = DENRYU(1,2)
          DC3 = DENRYU(1,3)
          DC4 = DENRYU(1,4)
          DC5 = DENRYU(1,5)
          IF (NDEG .EQ. 5)  GOTO 60
          DC6 = DENRYU(1,6)
          DC7 = DENRYU(1,7)
          IF (NDEG .EQ. 7)  GOTO 60
          DC8 = DENRYU(1,8)
          DC9 = DENRYU(1,9)
          DC10 = DENRYU(1,10)
          DC11 = DENRYU(1,11)
   60     CONTINUE
        ELSE
          DC1 = 0.0
          DC2 = 0.0
          DC3 = 0.0
          DC4 = 0.0
          DC5 = 0.0
          IF (NDEG .EQ. 5)  GOTO 70
          DC6 = 0.0
          DC7 = 0.0
          IF (NDEG .EQ. 7)  GOTO 70
          DC8 = 0.0
          DC9 = 0.0
          DC10 = 0.0
          DC11 = 0.0
   70     CONTINUE
        END IF 
C --------   current density calculation   ----------
C---------   Symmetric field 2    ----------
        DCON = DCS(I)/4.0
C------------------------------------------------------------------
        DS(I,7) = DC1 * DCON
        DS(I,8) = DC2 * DCON
        DS(I,9) = DC3 * DCON
        DS(I,10) = DC4 * DCON
        DS(I,11) = DC5 * DCON
        IF (NDEG .EQ. 5)  GOTO 80
        DS(I,12) = DC6 * DCON
        DS(I,13) = DC7 * DCON
        IF (NDEG .EQ. 7)  GOTO 80
        DS(I,14) = DC8 * DCON
        DS(I,15) = DC9 * DCON
        DS(I,17) = DC10 * DCON
        DS(I,18) = DC11 * DCON
   80   CONTINUE
C------  Computing coefficients of eddy current  -------
        IF (INT(NOD(I,4)/100) .EQ. 3)  THEN
C         DS(I,18) = DBH(NOD(I,4)-199, 4)*DOMEG*DCS(I)/12.0
C-------  Symmetric field 3  --------
          DS(I,18) = 5.92E7*DOMEG*DCS(I)/12.0*DDX
C------------------------------------------------------
        ELSE
          DS(I,18) = 0.0
        END IF
   10 CONTINUE
      RETURN
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C::::
      SUBROUTINE  CMAT2
      use param
      use matvec
C      PARAMETER (NA1=1995,NA2=249,NA3=755,NA4=755,NA5=400,NA6=121)                                                        
c      COMMON  AAA, BBB, NONC, CCC, NONC2, NCOIL, NHOWA,
C     &        NPO1, NOM, NPO2, NPO3, NPOR, NELEM, NB, NDEG, NDE, DF,
C     &        DH(NA1,NA2), DK(NA1), DAA(NA1), DA(NA1),
C     &        DS(NA3,18), DCS(NA3), DB(NA3,22), DD(NA3,NA6), NOD(0:NA4,4),
C     &        XY(NA5,2), DENRYU(1,11), DC(NA3,11), DCPRE(NA3,11),
C     &        DMUO, NOMEG, DPI, ITR, TOTAL,
C     &        DN(11,11), DBH(11,4)
C      DIMENSION  DCC(NA3,11), DCCC(NA3,11), S(9), DDA(3,11)
      real(8), DIMENSION(NA3,11)::  DCC, DCCC
      real(8), DIMENSION(9)::  S
      real(8), DIMENSION(3,11)::  DDA
	  
C.....initial
      DO  I=1, NELEM
        DO  J=1, NDEG
          DCCC(I,J) = 0.0
        ENDDO
      ENDDO
	  
      DO I=1, NCOIL
        S(I)=0.0
        DO  J=1, NDEG
          DCC(I,J)=0.0
        ENDDO
      ENDDO
C.....	  
      DO 60  I=1, NELEM
C              .....NOD(I,4) = 100 or 200 or 301 	  
        IF (INT(NOD(I,4)/100) .NE. 3)  GOTO 60
        IF (NOD(I,4) .EQ. 300)  GOTO 60
        J = NOD(I,4)-300
        S(J) = S(J)+DCS(I)
  103   CONTINUE
  
        DO J=1,3
          N = NOD(I,J)
          DO K=2,NDEG
            DDA(J, K) = DAA(NDEG * N + K)
          ENDDO
        ENDDO

        DO J=2,NDEG
          DCCC(I, J) = 0.0
          DO K=1,3
            DCCC(I,J) = DCCC(I, J) + DDA(K, J)
          ENDDO
        ENDDO
   60 CONTINUE
C---------------------------------------------------------DEBUG
C     DO 300  I=1, NELEM
C       DO 310  J=1, NDEG
C         PRINT '(10X,"DCCC(",I2,",",I2,")=",E15.4)',I,J,DCCC(I,J)
C 310   CONTINUE
C 300 CONTINUE
C-------------------------------------------------------------------
      DO I=1, NELEM
        IF (INT(NOD(I,4)/100) .NE. 3) exit
        IF (NOD(I,4) .EQ. 300) exit
        J = NOD(I,4) - 300
        DO  K=2, NDEG
          DCC(J,K) = DCC(J,K) + DCCC(I,K) * DCS(I)
        ENDDO
      ENDDO
C--------------------------------------------------------------DEBUG
C     DO 340  I=1, NCOIL
C       DO 330  J=1, NDEG
C         PRINT '(10X,"DCC( ",I2,",",I2,")=",E15.4)',I,J,DCC(I, J)
C 330   CONTINUE
C 340 CONTINUE
C-------------------------------------------------------------------------------------
      DO 107  I=1, NELEM
        DOMEG2 = 5.92E7 * DCS(I)/9.0 * DOMEG
        IF (INT(NOD(I,4)/100) .NE. 3) GOTO 104
        IF (NOD(I,4) .EQ. 300) GOTO 104
        J = NOD(I,4) - 300
        DC(I,1) = 0.0
        DC(I,2) = DOMEG2 * DCC(J,3) / S(J)
        DC(I,3) = -1.0 * DOMEG2 * DCC(J,2) / S(J)
        DC(I,4) =  2.0 * DOMEG2 * DCC(J,5) / S(J)
        DC(I,5) = -2.0 * DOMEG2 * DCC(J,4) / S(J)
        IF (NDEG .EQ. 5) GOTO 107
        DC(I,6) =  3.0 * DOMEG2 * DCC(J,7) / S(J)
        DC(I,7) = -3.0 * DOMEG2 * DCC(J,6) /S( J)
        IF (NDEG .EQ. 7) GOTO 107
        DC(I,8) =  4.0 * DOMEG2 * DCC(J,9) / S(J)
        DC(I,9) = -4.0 * DOMEG2 * DCC(J,8) / S(J)
        DC(I,10) =  5.0 * DOMEG2 * DCC(J,11) / S(J)
        DC(I,11) = -5.0 * DOMEG2 * DCC(J,10) / S(J)
        GOTO   107
  104   CONTINUE
  
        DO  J=1, NDEG
          DC(I,J) = 0.0
        ENDDO
  107 CONTINUE
C ------------------------------------------------------------DEBUG
C     DO 350  I=1, NELEM
C       DO 360  J=1, NDEG
C         PRINT '(10X,"DC(",I2,",",I2,")=",E15.4)',I,J,DCC(I,J)
C 360   CONTINUE
C 350 CONTINUE
C -----------------------------------------------------------------
      RETURN
      END 
C
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C::::Computing coefficients in Matrix DC
      SUBROUTINE  CMATP
      use param
      use matvec
C      PARAMETER (NA1=1995,NA2 = 249,NA3=755,NA4=755,NA5=400,NA6=121)
C      COMMON  AAA, BBB, NONC, CCC, NONC2, NCOIL, NHOWA,
C     &        NPO1, NOM, NPO2, NPO3, NPOR, NELEM, NB, NDEG, NDE, DF, 
C     &        DH(NA1,NA2), DK(NA1), DAA(NA1), DA(NA1),
C     &        DS(NA3,18), DCS(NA3), DB(NA3,22), DD(NA3,NA6), NOD(0:NA4,4)
C     &        XY(NA5,2), DENRYU(1,11), DC(NA3,11), DCPRE(NA3,11), 
C     &        DMUO , DOMEG, DPI, ITR, TOTAL,
c     &        DN(11,11), DBH(11,4)
C      DIMENSION  DCC (NA3,11), DCCC(NA3,11), S(9), DDA(3,11), DX(3)
      real(8), DIMENSION(NA3,11)::  DCC, DCCC
      real(8), DIMENSION(9)::  S
      real(8), DIMENSION(3,11)::  DDA
      real(8), DIMENSION(3):: DX
C 
      DO I=1, NELEM
        DO J=1, NDEG
          DCCC(I,J) = 0.0
        END DO
      END DO
  
      DO I=1, NCOIL
        S(I) = 0.0
        DO J=1, NDEG
          DCC(I,J) = 0.0
        END DO
      END DO
  
      DO 60  I=1, NELEM
        IF (INT(NOD(I,4)/100) .NE. 3) GOTO 60
        IF (NOD(I,4) .EQ. 300) GOTO 60
        J = NOD(I,4) - 300
        S(J) = S(J) + DCS(I)
  
        DO  J=1, 3
          N = NOD(I,J)
          DO  K=2, NDEG
            DDA(J,K) = DAA(NDEG * N + K)
          END DO
          DX(J) = XY(N+1, 1)
        ENDDO

        DDX = (DX(1) + DX(2) + DX(3)) / 3.0

        DO J=2, NDEG
          DCCC(I,J) = 0.0
          DO  K=1, 3
            DCCC(I,J) = DCCC(I,J) + DDA(K,J) * (DDX + DX(K) / 3.0)
          END DO
        ENDDO
   60 CONTINUE
C-------------------------------------------------------------DEBUG
C     DO 300  J=1, NELEM
C       DO 310  J=1, NDEG
C         PRINT '(10X,"DCCC(",I2,",",I2,")=",E15.4)',I,J,DCCC(I,J)
C 310   CONTINUE
C 300 CONTINUE
C ---------------------------------------------------------------------
      DO I=1, NELEM
        IF (INT(NOD(I,4) / 100) .NE. 3) exit
        IF (NOD(I,4) .EQ. 300) exit
        J = NOD(I,4) - 300
        DO K=2, NDEG
          DCC(J,K) = DCC(J,K) + DCCC(K,K) * DCS(I)
        END DO
      ENDDO
C ------------------------------------------------------------DEBUG
C     DO 340  I=1, NCOIL
C       DO 330  J=1, NDEG
C         PRINT '(10X,"DCC(",I2,",",I2,")=",E15.4)',I,J,DCC(I,J)
C 330   CONTINUE
C 340 CONTINUE
C------------------------------------------------------------------
      DO 107  I=1, NELEM
C       DOMEG2=5.92E7 * DCS(I) /16. * DOMEG
        DOMEG2=5.92E7 * DCS(I) /12. * DOMEG
        IF (INT(NOD(I,4) / 100 ) .NE. 3) GOTO 104
        IF (NOD(I,4) .EQ. 300) GOTO 104
C--------------------------------------------------------------DEBUG
C       DO 200 K=1,3
C         DX(K) = XY(NOD(I,K)+1, 1)
C 200   CONTINUE
C       DDX = (DX(1) + DX(2) + DX(3)) / 3.0
C -----------------------------------------------------------------
        J = NOD(I,4) - 300
        DC(I,1) = 0.0
        DC(I,2) = DOMEG2 * DCC(J,3) / S(J)
        DC(I,3) = -1.0 * DOMEG2 * DCC(J,2) / S(J)
        DC(I,4) =  2.0 * DOMEG2 * DCC(J,5) / S(J)
        DC(I,5) = -2.0 * DOMEG2 * DCC(J,4) / S(J)
        IF (NDEG .EQ. 5)  GOTO 107
        DC(I,6) =  3.0 * DOMEG2 * DCC(J,7) / S(J)
        DC(I,7) = -3.0 * DOMEG2 * DCC(J,6) / S(J)
        IF (NDEG .EQ. 7)  GOTO 107
        DC(I,8) =  4.0 * DOMEG2 * DCC(J,9) / S(J)
        DC(I,9) = -4.0 * DOMEG2 * DCC(J,8) / S(J)
        DC(I,10) =  5.0 * DOMEG2 * DCC(J,11) / S(J)
        DC(I,11) = -5.0 * DOMEG2 * DCC(J,10) / S(J)
        GOTO  107
  104   CONTINUE
        DO  J=1,  NDEG
          DC(I,J) = 0.0
        END DO
  107 CONTINUE
C-------------------------------------------------------------DEBUG
C     DO  350  I=1, NELEM
C       DO  360  J=1, NDEG
C         PRINT '(10X,"DC(",I2,",",I2,")=",E15.4)',I,J,DC(I,J)
C 360   CONTINUE
C 350 CONTINUE
C------------------------------------------------------------------
      RETURN
      END
C
