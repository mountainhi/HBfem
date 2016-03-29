
      PARAMETER (NA1=1995,NA2=249,NA3=755,NA4=755,NA5=400,NA6=121)
C...  from main.f
	  parameter dpi = 3.1415927
C     frequency
	  parameter df  = 1.0363E6
	  
C*****     AAA : criterion of convergence
C          BBB : deceleration factor
C          CCC :
           NCOIL : number of coil
	  COMMON  AAA, BBB, CCC, NCOIL, NOWA,
	  
C*****    status return from convergence test
     &     NONC, NONC2,   
C*****   Size of the problem
C          NPO1 : number of unknown potentials(Dof)  
C          NPO2 : last number of node on dirichlet's condition    
C          NPO3 : number of vertex
C          NOM  : size of system matrix
C          NB   : bandwidth of matrix 
C          NPOR : count of node on periodic boundary 
C                             NPOR = NPO3 - NPO2
C          NELEM: number of elements(nodes)

     &     NPO1, NOM, NPO2, NPO3, NPOR, NELEM, NB, NDEG, NDE, DF,
	 
C*****   matrix and vectors
C*****          H * A = K
C           DH :
C           DK :
C	        DA :
C		    DAA: initial values of vector potential 
C           DCS:  cross section of element 

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

