C::::::::::::::::::::::::::::::::::::::::::::::::::::::
C::::::: Convergence test and Relaxation ::::::::::::::::::::
        SUBROUTINE  CONV
        include 'dm.inc'  
C		 PARAMETER  (NA1=1995, NA2=249, NA3=755, NA4=744, NA5=400, NA6=121)
C           COMMON     AAA, BBB, NONC, CCC, NONC2, NCOIL, NHOWA,
C         &                      NPO1, NOM, NPO2, NPOR, NELEM, NB, NDEG, NDE, DF, 
C         &                      DH(NA1, NA2), DK(NA1), DAA(NA1), DA(NA1), 
C         &                 DS(NA3, 18), DCS(NA3), DB(NA3, 22), DD(NA3, NA6), NOD(0:NA4,4),
C         &                 XY(NA5, 2), DENRYU(1, 11), DC(NA3, 11), DCPRE(NA3, 11),
C         &                 DMUO, DOMEG, DPI, ITR, TOTAL, 
C         &                 DN(11, 11), DBH(11, 4)
        INTEGER  TOTAL
C
        IF (ITR.EQ.1)  THEN
         PRINT '(30X,"Initial values were used when calculation(current)")'
         PRINT '(30X,"first results for current calculation")'
         PRINT '(30X,"?????????????????????????")'
           GOTO  70
        END  IF
C
        NN1=0
        NN2=0
        NN3=0
        NN4=0
        NONC=0
        NONC2=0
		
        DRR1=0.
        DRR2=0.
        DRR3=0.
        DRR4=0.

        DCMAX=0.
        NSOSU=0
		
        DAMAX=0.
        DO I=1,NOM
           IF  (ABS(DA(I)).GT.DAMAX)   DAMAX=ABS(DA(I))
        ENDDO
        DO I=1,NOM
C--------------------------------------------------------------DEBUG
              IF  (DAA(I).EQ.0.)  THEN
              PRINT *,'DAA(',i,')=0.'
              STOP
              END IF
C--------------------------------------------------------------
              DERR=(DA(I) - DAA(I))/DAMAX
C--------------------------------------------------------------DEBUG
C             IF (DERR.EQ.0.)   THEN 
C             PRINT  *, 'DERR=0. WHEN I=',I
C             END IF
C----------------------------------------------------------
              IF (ABS(DERR). GT.AAA)  THEN
              DRR1 = ABS(DERR) +DRR1
              NN1=NN1+1
              NONC= NONC+1
              END IF
        ENDDO
C----------------------------------------------------------Test on DC
       DO  30 I=1, NELEM
           IF (INT(NOD(I,4) /100).NE.3)  GOTO 30
           IF (NOD(I,4).EQ.300)  GOTO 30
           NSOSU=NSOSU+1
           DO  40 J=2, NDEG
              DERR2=(DC(I,J) - DCPRE(I,J))/DC(I,J)
              IF(ABS(DERR2).GT.CCC)  THEN
              DRR3=DRR3+ABS(DERR2)
              NN3=NN3+1
              NONC2=NONC2+1
C------------------------------------------------------------DEBUG
C                IF  (ITR.GT.30)
C     &          PRINT *,'DC(I,J) does not converge,I=', I,'J=', J
C-----------------------------------------------------------
              END IF
            ENDDO
   30  CONTINUE
   
   
        PRINT '(5X,"  The number of non-convergence point of A =",
     &	  	I6, 5X, "Total =", I6)'
     &       , NONC, NOM
        PRINT '(5X,"  Average number of non-convergence point of A =",
     &		I6, 5X, "Total =", I6)'
     &       , NONC2, NSOSU*NDEG
        IF ((NONC.NE.0). AND. (NONC2.NE.0))  THEN
        PRINT '(30X," A, Aver of A, Not all converged")'
        GOTO 90
        ELSE IF((NONC.EQ.0).AND.(NONC2.EQ.0)) THEN
        PRINT '(30X," A and Aver of A.  All converged")'           
        GOTO 70
        ELSE IF (NONC.EQ. 0)  THEN
        PRINT '(30X, "  A converged")'
        GOTO 90
        ELSE IF (NONC2.EQ.0)  THEN
        PRINT '(30X, "  Aver of A converged")'           
        GOTO 90
        END IF
   90 CONTINUE
   
        IF (NN1.NE.0)  THEN
        PRINT'(10X," Threshold AAA=", E9.2,10X," Average error ="
     &     'E15.4)',AAA, DRR1/NN1
        END IF
        IF (NN3.NE.0)  THEN
        PRINT'(10X," Threshold CCC=", E9.2,10X," Average error ="
     &     'E15.4)',CCC, DRR3/NN3
        END IF
C---------------------------------------------------------deceleration relaxation
   70 CONTINUE
C       IF (NONC.NE.0)  THEN
C----------------------------------------------------DF=300.E3
C       IF (TOTAL+ITR.GE.500)  THEN
C           BBBB=BBB*1.1
C       ELSE IF (TOTAL+ITR.GE.200)  THEN
C           BBBB=BBB*1.3
C       ELSE IF (TOTAL+ITR.GE.100)  THEN
C           BBBB=BBB*1.5
C       ELSE IF (TOTAL+ITR.GE.50)  THEN
C           BBBB=BBB*1.8
C       ELSE IF (TOTAL+ITR.GE.10)  THEN
C           BBBB=BBB*2.0
C       ELSE IF (TOTAL+ITR.GE.4)  THEN
C           BBBB=BBB*3.0
C       ELSE IF (TOTAL+ITR.GE.1)  THEN
C           BBBB=BBB*4.0
C       END IF
C---------------------------------------------------DF=10.E1
C       IF(TOTAL+ITR.GE.50)  THEN
C           BBBB=BBB*0.1
C       ELSE IF(TOTAL+ITR.GE.40)  THEN
C           BBBB=BBB*0.4
C       ELSE IF (TOTAL+ITR.GE.30)  THEN
C           BBBB=BBB*0.7
C       ELSE IF (TOTAL+ITR.GE.20)  THEN
C           BBBB=BBB*1.0
C       ELSE IF (TOTAL+ITR.GE.10)  THEN
C           BBBB=BBB*2.0
C       ELSE IF (TOTAL+ITR.GE.4)  THEN
C           BBBB=BBB*3.0
C       ELSE IF (TOTAL+ITR.GE.1)  THEN
C           BBBB=BBB*4.0
C       END IF
C---------------------------------------------------DF=1000.E0
C       IF (TOTAL+ITR.GE.44)  THEN
C           BBBB=BBB*0.01
C       ELSE IF (TOTAL+ITR.GE.40)  THEN
C           BBBB=BBB*1.2
C       ELSE IF (TOTAL+ITR.GE.35)  THEN
C           BBBB=BBB*1.5
C       ELSE IF (TOTAL+ITR.GE.20)  THEN
C           BBBB=BBB*1.8
C       ELSE IF (TOTAL+ITR.GE.10)  THEN
C           BBBB=BBB*2.0
C       ELSE IF (TOTAL+ITR.GE.1)  THEN
C           BBBB=BBB*3.0
C       END IF
C--------------------------------------------------
C    ELSE
C       IF (TOTAL+ITR.GE.80)  THEN
C           BBBB=BBB*0.04
C       ELSE IF (TOTAL+ITR.GE.1)  THEN
C           BBBB=BBB*0.01
C       END IF
C     END IF
C-------------------------------------------------------
      BBBB=.5
C      IF(TOTAL +ITR.GE.10)  THEN
C         BBBB=0.4
C      END IF
      IF ((NONC.NE.0).OR.(NONC2.NE.0))  THEN
       PRINT '(10X, "Deceleration relaxation =".E9.2)', BBBB
      END IF
C  100 CONTINUE
C..
       DO I=1,NOM
          DAA(I) = (1. -BBBB)*DAA(I)+BBBB*DA(I)
       ENDDO
       DO  L=1,NELEM
          DO  M=1,NDEG
             DCPRE(L,M)=DC(L,M)
          ENDDO
       ENDDO
      RETURN
      END