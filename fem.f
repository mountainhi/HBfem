C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C::::: FEM  Computing of coefficients of element matrix DH() ::::::
      SUBROUTINE  FEM
      PARAMETER (NA1=1995,NA2=249,NA3=755,NA4=755,NA5=400,NA6=121)
      COMMON  AAA, BBB, NONC, CCC, NONC2, NCOIL, NHOWA, 
     &        NPO1, NOM, NPO2, NPO3, NPOR, NELEM, NB, NDEG, NDE, DF,
     &        DH(NA1,NA2), DK(NA1), DAA(NA1), DA(NA1),
     &        DS(NA3,18), DCS(NA3), DB(NA3,22), DD(NA3,NA6), NOD(0:NA4,4),
     &        XY(NA5,2), DENRYU(1,11), DC(NA3,11), DCPRE(NA3,11),
     &        DMUO, DOMEG, DPI, ITR, TOTAL,
     &        DN(11,11), DBH(11,4)
      DIMENSION  DHH(11,11), DSS(3,3)
      REAL NNN
C----------------------------------------------------INITIALIZATION
      DO 99  I=1, NA1
        DK(I) = 0.0
   99 CONTINUE
      DO 21  I=1, NA1
        DO 31  J=1, NA2
          DH(I,J) = 0.0
   31   CONTINUE
   21 CONTINUE
C------------------------------------------------------------------
      DO 500  NE=1, NELEM
        NSS=0
        DO 10  J=1, 3
          DO 20  K=J, 3
            NSS = NSS + 1
            DSS(J,K) = DS(NE, NSS)
            DSS(K,J) = DS(NE, NSS)
   20     CONTINUE
   10   CONTINUE
        NSS=0
        DO 30  N=1, NDEG
          DO 40  M=1, NDEG
            NSS = NSS + 1
            DHH(N,M) = DD(NE, NSS)
   40     CONTINUE
   30   CONTINUE
C----------------------------------------- Substitute into Matrix K
        DO 50  J=1, 3
          II = NOD(NE,J) + 1
          NLOW = NDEG * (II-1) + 1
          IF ((II .GT. NPO1) .AND. (II .LE. NPO2)) THEN
            GOTO 50
          ELSE
            DSIG1 = 1.0
          END IF
          IF (INT(NOD(NE,4) / 100) .NE. 3)   GO TO  200
          DO 60  K=1, NDEG
            DK(NLOW+K-1) = DK(NLOW+K-1)+DS(NE,6+K)*DSIG1-DC(NE,K)    
C           DK(NLOW+K-1) = DK(NLOW+K-1)+DS(NE,6+K)*DSIG1
   60     CONTINUE
  200     CONTINUE
          DO 70  K=1, 3
            NNN = 1.0
            IF (J .EQ. K)  NNN=2.0
            JJ = NOD(NE,K) + 1
            NCOL = NDEG * (JJ-1) + 1
            IF ((JJ .GT. NPO1) .AND. (JJ .LE. NPO2)) THEN
              GOTO   70
            ELSE
              DSIG2 = 1.0
            END IF
C------------------------------- Substitute into band matrix DH()
            DCON = DSIG1 * DSIG2 * DSS(J,K)
            DO 80  M=1, NDEG
              DO 90  N=1, NDEG
                NL = NLOW + M -1
                NC = NCOL + N -1
                NC1 = NC - NL + NB
                CON = DHH(M,N) * DCON + DN(M,N) * DS(NE,18) * NNN
C               CON = DHH(M,N) * DCON
                DH(NL,NC1) = DH(NL,NC1) + CON
   90         CONTINUE
   80       CONTINUE
   70     CONTINUE
   50   CONTINUE 
  550 CONTINUE
      RETURN
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C::::::: FEM  Computing coefficients of element matrix DH() :::::::
      SUBROUTINE  FEMP
      include "dm.inc"
CCC      PARAMETER (NA1=1995,NA2=249,NA3=755,NA4=755,NA5=400,NA6=121)
C      COMMON  AAA, BBB, NONC, CCC, NONC2, NCOIL, NOWA,
C     &        NPO1, NOM, NPO2, NPO3, NPOR, NELEM, NB, NDEG, NDE, DF,
C     &        DH(NA1,NA2), DK(NA1), DAA(NA1), DA(NA1),
C     &        DS(NA3,18), DCS(NA3), DB(NA3, 22), DD(NA3, NA6), NOD(0:NA4,4),
C     &        XY(NA5,2), DENRYU(1,11), DC(NA3,11), DCPRE(NA3,11),
C     &        DMUO, DOMEG, DPI, ITR, TOTAL, 
C     &        DN(11,11), DBH(11,4)
C******************************************************************
C     DH(I,J) :  MATRIX          (NOB, NB)
C     DK(I)   :  FORCED VECTOR   (NOB)
C     DA(I)   :  POTENTIAL       (NOB)
C     NOM     :  NUMBER OF UNKNOW POTENTIAL
C     NB      :  BAND  WIDTH
      DIMENSION  DHH(11,11), DSS(3,3), DX(3)
      REAL  NNN
C----------------------------------------------------INITIALIZATION
      DO  I=1, NA1
        DK(I) = 0.0
      ENDDO
      DO  I=1, NA1
        DO  J=1, NA2
          DH(I,J) = 0.0
        ENDDO
      ENDDO
C------------------------------------------------------------------
      DO  500  NE=1, NELEM
        NSS = 0
        DO  J=1, 3
          DO  K=J, 3
            NSS = NSS + 1
            DSS(J,K) = DS(NE,NSS)
            DSS(K,J) = DS(NE,NSS)
          ENDDO
          DX(J) = XY(NID(NE,J)+1, 1)
        ENDDO
        DDX = (DX(1) + DX(2) + DX(3)) /3.0
        NSS = 0
		
        DO  N=1, NDEG
          DO  M=1, NDEG
            NSS = NSS + 1
            DHH(N,M) = DD(NE,NSS)
          ENDDO
        ENDDO
C------------------------------------------Substitute into matrix K
        DO 50  J=1, 3
          II = NDD(NE,J) + 1
          NLOW = NDEG * (II-1) + 1
          IF ((II .GT. NPO1) .AND. (II .LE. NPO2))  THEN
            GOTO 50
          ELSE
            DSIG1 = 1.0
          END   IF
          IF(INT(NOD(NE,4)/100) .NE. 3)  GO TO  200
          DDX2 = DDX + XY(II,1) / 3.0 
          DO  K=1, NDEG
            DK(NLOW+K-1) = DK(NLOW+K-1)+DS(NE,6+K)*DSIG1*DDX2-DC(NE,K)
            DK(NLOW+K-1) = DK(NLOW+K-1)+DS(NE,6+K)*DSIG1*DDX2
          ENDDO
  200     CONTINUE
          DO 70  K=1, 3
            NNN=1.0
            IF (J .EQ. K)  NNN=2.0                                                           
            JJ = NOD(NE,K) + 1
            NCOL = NDEG * (JJ-1) + 1
            IF ((JJ .GT. NPO1) .AND. (JJ .LE. NPO2))  THEN
              GOTO 70
            ELSE
              DSIG2 = 1.0
            END IF
C-------------------------------- Substitute into band matrix DH()
            DCON= DSIG1 * DSIG2 * DSS(J,K)
            DO  M=1, NDEG                        
              DO  N=1, NDEG
                NL = NLOW + M - 1
                NC = NCOL + N - 1
                NC1 = NC - NL + NB
                CON = DHH(M,N) * DCON + DN(M,N) * DS(NE,18) * NNN       
C               CON = DHH(M,N) * DCON
                DH(NL,NC1) = DH(NL,NC1) + CON
              ENDDO
            ENDDO
			
   70      CONTINUE
   50   CONTINUE
  500 CONTINUE
      RETURN
      END