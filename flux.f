C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C :: FLUX  Computing of coefficient of flux density  from DAA() :::
C
      SUBROUTINE  FLUX
      include 'dm.inc'
C      PARAMETER (NA1=1995,NA2=249,NA3=755,NA4=755,NA5=400,NA6=121)
C      COMMON  AAA, BBB, NONC, CCC, NONC2, NCOIL, NHOWA,
C     &        NPO1, NOM, NPO2, NPO3, NPOR, NELEM, NB, NDEG, NDE, DF,
C     &        DH(NA1,NA2), DK(NA1), DAA(NA3,22), DD(NA3,NA6), NOD(0:NA4,4),
C     &        DS(NA3,18), DCS(NA3), DB(NA3,22), DD(NA3,NA6), NOD(0:NA4,4),
C     &        XY(NA5,2), DENRYU(1,11), DC(NA3,11), DCPRE(NA3,11),
c     &        DMUO, DOMEG, DPI, ITR, TOTAL, 
c     &        DN(11,11), DBH(11,4)
      DIMENSION  DX(3), DY(3), DQ(3), DR(3), DDA(3,11)
C
C...  
      DO 10  I=1, NELEM
	  
        DO J=1,3
          DX(J) = XY(NOD(I,J) + 1, 1)
          DY(J) = XY(NOD(I,J) + 1, 2)
          DO  K=1, NDEG
            N = NOD(I,J)
            IF (N .LT. NPO1)  THEN
              DDA(J,K) = DAA(NDEG * N + K)
            ELSE
              DDA(J,K) = 0.0
            END IF
          ENDDO
        ENDDO
   
        DO J=1,3
          DQ(J) = DY(MOD(J,3) + 1) - DY(MOD(J+1,3) + 1)
          DR(J) = DX(MOD(J+1,3) + 1) - DX(MOD(J,3) + 1)
        ENDDO
		
        DO  J=1, NDEG
          JJJ = J + NDEG
          DB(I,J) = 0.0
          DB(I,JJJ) = 0.0
          DO K=1,3
            DB(I,J) = DB(I,J) + DR(K) * DDA(K,J)
            DB(I,JJJ) = DB(I,JJJ) - DQ(K) * DDA(K,J)
          ENDDO
          DB(I,J) = DB(I,J) * 0.5 / DCS(I)
          DB(I,JJJ) = DB(I,JJJ) * 0.5 / DCS(I)
        ENDDO
   
   10 CONTINUE
      RETURN 
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::  FLUX  Computing of coefficient of flux density from DAA() :::
C::                   axial symmetric
C            Matrix DB
      SUBROUTINE  FLUXP
      include 'dm.inc'
C      PARAMETER (NA1=1995,NA2=249,NA3=755,NA4=755,NA5=400,NA6=121)
C      COMMON  AAA,BBB,NONC,CCC,NONC2,NCOIL,NHOWA,
C     &        NPO1, NOM, NPO2, NPO3, NPOR, NELEM, NB, NDEG, NDE, DF,
C     &        DH(NA1,NA2),DK(NA1),DAA(NA1), DA(NA1),
C     &        DS(NA3,18), DCS(NA3), DB(NA3,22),DD(NA3, NA6),NOD(0:NA4,4),
C     &        XY(NA5M2), DENRYU(1,11), DC(NA3,11), DCPRE(NA3,11)
C     &        DMUO, NOMEG, NPI, ITR, TOTAL,
C     &        DN(11,11),  DBH(11,4)
      DIMENSION  DX(3), DY(3), DQ(3), DR(3), DDA(3,11), AVEA(11)
C
C.....For each elements
      DO 10  I=1, NELEM
	  
        DO  J=1,3
          DX(J) = XY(NOD(I,J) + 1, 1)
          DY(J) = XY(NOD(I,J) + 1, 2)
          DO  K=1, NDEG
            N = NOD(I,J)
            IF (N .LT. NPO1)  THEN
              DDA(J,K) = DAA(NDEG * N + K)
            ELSE 
              DDA(J,K) = 0.0
            END IF
          ENDDO
        ENDDO
   
        DDX = (DX(1) + DX(2) +DX (3)) / 3.0
		
        DO J=1, NDEG
          AVEA(J) = (DDA(1, J) + DDA(2,J) + DDA(3,J)) / 3.0
        ENDDO
        DO J=1,3
          DQ(J) = DY(MOD(J,3) + 1) - DY(MOD(J+1,3) + 1)
          DR(J) = DX(MOD(J+1,3) + 1) - DX(MOD(J,3) + 1)
        ENDDO
		
        DO  J=1, NDEG
          JJJ= J + NDEG
          DB(I,J) = 0.0
          DB(I,JJJ) = 0.0
          DO K=1,3
            DB(I,J) = DB(I,J) - DR(K) * DDA(K,J)
            DB(I,JJJ) = DB(I,JJJ) + DQ(K) * DDA(K,J)
          ENDDO
          DB(I,J) = DB(I,J) * 0.5/DCS(I)
          DB(I,JJJ) = DB(I,JJJ) * 0.5/DCS(I) + AVEA(J) / DDX
        ENDDO
   
   10 CONTINUE
      RETURN
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::