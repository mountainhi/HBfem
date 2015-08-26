
      PARAMETER (NA1=1995,NA2=249,NA3=755,NA4=755,NA5=400,NA6=121)
C...  from main.f
	  parameter dpi = 3.1415927
	  parameter df  = 1.0363E6
	  parameter MAX = 40
	  
C*****     AAA is relative error bound
	  COMMON  AAA, BBB, CCC, NCOIL, NOWA,
C*****    convergence test
     &     NONC, NONC2,   
C*****   Size of the problem
C*****     NPO3 : number of nodes
     &     NPO1, NOM, NPO2, NPO3, NPOR, NELEM, NB, NDEG, NDE, DF,
C*****   matrix and vectors
     &     DH(NA1,NA2), DK(NA1), DAA(NA1), DA(NA1),
     &     DS(NA3,18), DCS(NA3), DB(NA3, 22), DD(NA3, NA6), 
C*****   geometric info
     &     NOD(0:NA4,4),XY(NA5,2), 
C*****  Current density
     &     DENRYU(1,11), DC(NA3,11), DCPRE(NA3,11),
C*****  computational parameters 

     &     DMUO, DOMEG, DPI, ITR, TOTAL, 
C*****  harmonic balance N matrix
     &     DN(11,11), DBH(11,4)
	 
C******************************************************************
C     DH(I,J) :  MATRIX          (NOB, NB)
C     DK(I)   :  FORCED VECTOR   (NOB)
C     DA(I)   :  POTENTIAL       (NOB)
C     NOM     :  NUMBER OF UNKNOW POTENTIAL
C     NB      :  BAND  WIDTH
C******************* ELIMINATE L OF MATRIX DH(I,J) ****************

