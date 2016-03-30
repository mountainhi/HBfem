C-------------------------------------------------------------------
       module param
            save

          real(8)::    AAA, BBB, CCC
          integer(4):: NCOIL, NHOWA
          integer(4):: NONC , NONC2
          integer(4):: NPO1, NOM, NPO2,NPO3,NPOR,NELEM,NB,NDEG,NDE
          real(8),parameter:: df=1.0363E6, DPI = 3.1415927
          real(8)::    DMUO, DOMEG
          integer(4):: ITR, TOTAL

        end module param
       
C*****     AAA : criterion of convergence
C          BBB : deceleration factor
C          CCC :
C          NCOIL : number of coil
C	         NHOWA :
C	  
C*****    NONC, NONC2 : status return from convergence test
C
C*****   Size of the problem
C          NPO1 : number of unknown potentials(Dof)  
C          NPO2 : last number of node on dirichlet's condition    
C          NPO3 : number of vertex
C          NOM  : size of system matrix
C          NB   : bandwidth of matrix 
C          NPOR : count of node on periodic boundary 
C                             NPOR = NPO3 - NPO2
C          NELEM: number of elements(nodes)
C*******
C*****  computational parameters 
C
C          DMUO, DOMEG, DPI, ITR, TOTAL, 
C
c.... DMUO is the permeability of air
C----------------------------------------------------------------------
        module matvec
          save
            integer,parameter::NA1=1995,NA2=249,NA3=755
            integer,parameter::NA4=755,NA5=400,NA6=121
            real(8):: DH(NA1,NA2), DK(NA1), DAA(NA1), DA(NA1)
            real(8):: DS(NA3,18), DCS(NA3), DB(NA3, 22), DD(NA3, NA6)
            real(8):: DN(11,11), DBH(11,4)
            real(8):: NOD(0:NA4,4),XY(NA5,2)
            real(8):: DENRYU(1,11), DC(NA3,11), DCPRE(NA3,11)
        end module matvec
C*****   matrices and vectors
C*****          H * A = K
C           DH :
C           DK :
C	          DA :
C		        DAA: initial values of vector potential 
C           DCS:  cross section of element 
      
C*****   geometric info
C        NOD(0:NA4,4),XY(NA5,2),
C 
C*****  Current density
C        DENRYU(1,11), DC(NA3,11), DCPRE(NA3,11),
C
C*****  harmonic balance N matrix
C        DN(11,11), DBH(11,4)
	 
C******************************************************************
C     DH(I,J) :  MATRIX          (NOB, NB)
C     DK(I)   :  FORCED VECTOR   (NOB)
C     DA(I)   :  POTENTIAL       (NOB)
C     NOM     :  NUMBER OF UNKNOW POTENTIAL
C     NB      :  BAND  WIDTH
C******************* ELIMINATE L OF MATRIX DH(I,J) ****************
