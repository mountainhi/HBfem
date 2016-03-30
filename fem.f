C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C::::: FEM  Assemble the matrices DH() , DK()                ::::::
      SUBROUTINE  FEM
      use param
      use matvec 
C      PARAMETER (NA1=1995,NA2=249,NA3=755,NA4=755,NA5=400,NA6=121)
C      COMMON  AAA, BBB, NONC, CCC, NONC2, NCOIL, NHOWA, 
C     &        NPO1, NOM, NPO2, NPO3, NPOR, NELEM, NB, NDEG, NDE, DF,
C     &        DH(NA1,NA2), DK(NA1), DAA(NA1), DA(NA1),
C     &        DS(NA3,18), DCS(NA3), DB(NA3,22), DD(NA3,NA6), NOD(0:NA4,4),
C     &        XY(NA5,2), DENRYU(1,11), DC(NA3,11), DCPRE(NA3,11),
C     &        DMUO, DOMEG, DPI, ITR, TOTAL,
C     &        DN(11,11), DBH(11,4)

      real(8), DIMENSION(11,11)::DHH 
      real(8), DIMENSION(3,3)::DSS
      REAL(8) ::NNN
C---------------------
C.... Internal
C     DX ()
C     DHH()
C     DDX()
C 
C******************************************************************
C....    DH ( NA1,NA2) * DA () = DK (NA1)
C     DH(I,J) :  MATRIX          (NOB, NB)
C     DK(I)   :  FORCED VECTOR   (NOB)
C     DA(I)   :  POTENTIAL       (NOB)
C     NOM     :  NUMBER OF UNKNOW POTENTIAL
C     NB      :  BAND  WIDTH
C----------------------------------------------------INITIALIZATION
      DO I=1, NA1
        DK(I) = 0.0
      ENDDO
      DO I=1, NA1
        DO  J=1, NA2
          DH(I,J) = 0.0
        ENDDO
      ENDDO
C------------------------------------------------------------------
      DO 500  NE=1, NELEM
        NSS=0		
        DO  J=1, 3
          DO  K=J, 3
            NSS = NSS + 1
            DSS(J,K) = DS(NE, NSS)
            DSS(K,J) = DS(NE, NSS)
          ENDDO
        ENDDO
		
        NSS=0
        DO  N=1, NDEG
          DO  M=1, NDEG
            NSS = NSS + 1
            DHH(N,M) = DD(NE, NSS)
          ENDDO
        ENDDO
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
		  
          DO  K=1, NDEG
C                           ....DS( , 6+k) = delta * {J1s J1c .....} /3
            DK(NLOW+K-1) = DK(NLOW+K-1)+DS(NE,6+K)*DSIG1-DC(NE,K)    
C           DK(NLOW+K-1) = DK(NLOW+K-1)+DS(NE,6+K)*DSIG1
          ENDDO
   
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
C------------------------------- Assemble band matrix DH()
C                              (Reluctivity + Harmonics)
C
C.......    parameters (bb+cc)/4*delta
            DCON = DSIG1 * DSIG2 * DSS(J,K)
			
            DO M=1, NDEG
              DO N=1, NDEG
                NL = NLOW + M -1
                NC = NCOL + N -1
                NC1 = NC - NL + NB
C                                DS(NE,18) = 5.93e7*omeg*delta/12
                CON = DCON * DHH(M,N) + NNN * DS(NE,18) * DN(M,N)
C               CON = DHH(M,N) * DCON
                DH(NL,NC1) = DH(NL,NC1) + CON
              ENDDO
            ENDDO
   
   70     CONTINUE
   50   CONTINUE 
   
  500 CONTINUE
      RETURN
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C::::::: FEM  Computing coefficients of element matrix DH() :::::::
      SUBROUTINE  FEMP
      use param
      use matvec 
CCC      PARAMETER (NA1=1995,NA2=249,NA3=755,NA4=755,NA5=400,NA6=121)
C      COMMON  AAA, BBB, NONC, CCC, NONC2, NCOIL, NOWA,
C     &        NPO1, NOM, NPO2, NPO3, NPOR, NELEM, NB, NDEG, NDE, DF,
C     &        DH(NA1,NA2), DK(NA1), DAA(NA1), DA(NA1),
C     &        DS(NA3,18), DCS(NA3), DB(NA3, 22), DD(NA3, NA6), NOD(0:NA4,4),
C     &        XY(NA5,2), DENRYU(1,11), DC(NA3,11), DCPRE(NA3,11),
C     &        DMUO, DOMEG, DPI, ITR, TOTAL, 
C     &        DN(11,11), DBH(11,4)
C---------------------
C.... Internal
C     DX ()
C     DHH()
C     DDX()
C 
C******************************************************************
C....  
C     DH(I,J) :  MATRIX          (NOB, NB)
C     DK(I)   :  FORCED VECTOR   (NOB)
C     DA(I)   :  POTENTIAL       (NOB)
C     NOM     :  NUMBER OF UNKNOW POTENTIAL
C     NB      :  BAND  WIDTH
      real(8), DIMENSION(11,11)::DHH 
      real(8), DIMENSION(3,3)::DSS
      real(8), DIMENSION(3)::DX

      REAL(8)::NNN
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
          DX(J) = XY(NOD(NE,J)+1, 1)
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
          II= NOD(NE,J)+1
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
C            DK(NLOW+K-1) = DK(NLOW+K-1)+DS(NE,6+K)*DSIG1*DDX2
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