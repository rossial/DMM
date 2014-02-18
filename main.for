C --------------------------------------------------------------------------------------
C Program DMM: Bayesian and Classical Inference of Dynamic Mixture Models
C Developed by A.Rossi, C.Planas and G.Fiorentini     
C
C State-space format:   y(t) = c(t)z(t) + H(t)x(t)   + G(t)u(t)
C                       x(t) = a(t)     + F(t)x(t-1) + R(t)u(t)
C
C y(t) (ny x 1)          ny  = # of endogenous series
C z(t) (nz x 1)          nz  = # of exogenous series
C x(t) (nx x 1)          nx  = # of continous states 
C u(t) (nu x 1)          nu  = # of shocks
C c(t) (ny x nz x ns1)   ns1 = # of states for c(t)
C H(t) (ny x nx x ns2)   ns2 = # of states for H(t)  
C G(t) (ny x nu x ns3)   ns3 = # of states for G(t)  
C a(t) (nx x ns4)        ns4 = # of states for a(t)  
C F(t) (nx x nx x ns5)   ns5 = # of states for F(t)  
C R(t) (nx x nu x ns6)   ns6 = # of states for R(t)  
C
C Copyright (C) 2010-2014 European Commission 
C
C This file is part of Program DMM
C
C DMM is free software developed at the Joint Research Centre of the 
C European Commission: you can redistribute it and/or modify it under 
C the terms of the GNU General Public License as published by
C the Free Software Foundation, either version 3 of the License, or
C (at your option) any later version.
C
C DMM is distributed in the hope that it will be useful,
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C GNU General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with DMM.  If not, see <http://www.gnu.org/licenses/>.  
C      
C In addition, as a special exception, the copyright holders give
C permission to link the code of portions of this program with the
C NAG Fortran library under certain conditions as described in each
C individual source file, and distribute linked combinations including
C the two.
C
C You must obey the GNU General Public License in all respects for all
C of the code used other than NAG Fortran library. If you modify file(s)
C with this exception, you may extend this exception to your
C version of the file(s), but you are not obligated to do so. If
C you do not wish to do so, delete this exception statement from
C your version. If you delete this exception statement from all
C source files in the program, then also delete it here.      
C --------------------------------------------------------------------------------------
      PROGRAM DMM
	USE dfwin
C DECLARE an "interface block" to the .DLL that contains DESIGN

	INTERFACE
	 SUBROUTINE DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R) 
	 INTEGER ny,nz,nx,nu,ns(6),nt
	 DOUBLE PRECISION theta(nt)
	 DOUBLE PRECISION c(ny,max(1,nz),ns(1)),H(ny,nx,ns(2)),
	1 G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)),R(nx,nu,ns(6))
	 END SUBROUTINE
	END INTERFACE
	CHARACTER*1 fittizia
	POINTER (pdll,fittizia)  ! ASSOCIATE  pointer pdll alla DLL ad una varibile fittizia
	POINTER (pdesign,DESIGN) 
C DECLARE an "interface block" to the .DLL that contains SETFILEM
	INTERFACE
	 SUBROUTINE SETFILEM(mfile,pathmfile) 
	 CHARACTER*200 mfile,pathmfile
       END SUBROUTINE
      END INTERFACE	
      POINTER (psetfilem,SETFILEM) 
C DECLARE an "interface block" to the .DLL that contains GETERRSTR
	INTERFACE
	 SUBROUTINE GETERRSTR(matlaberror) 
	 CHARACTER*1024 matlaberror
       END SUBROUTINE
      END INTERFACE	
      POINTER (pgeterrstr,GETERRSTR)
      
      LOGICAL status
	CHARACTER*200 DLLNAME    ! name of the DLL (defined by the user)

C NAMELIST DECLARATIONS
	INTEGER ny,nz,nx,nu,d(2),nv,nt,nf,nobs
	INTEGER seed,thin,burnin,GGG,HBL
	CHARACTER*1 MargLik,datasim,check 
	CHARACTER*2 thetasampler,estimation
	CHARACTER*3 Ssampler
	DOUBLE PRECISION, ALLOCATABLE:: hypS(:,:,:),hyptheta(:,:),
	1 obs(:)
      CHARACTER(2), ALLOCATABLE:: pdftheta(:)
C LOCALS
	INTEGER, ALLOCATABLE:: Z(:),ZW(:),S(:,:),NEVAL(:),
	1 CUMN(:),IYK(:,:),INDT(:),ACCRATE(:),gibZ(:,:),
     1 IT1(:),IT2(:),DATE_ITIME(:),np(:),ns(:),INFOS(:,:)
	DOUBLE PRECISION, ALLOCATABLE:: yk(:,:),STATE(:,:),theta(:),
	1 theta0(:),thetaprior(:,:),psi(:),psi0(:),psiprior(:,:),
     2 INN(:,:),FORE(:,:),ykmis(:),PTR(:,:,:),PM(:,:),
     3 gibtheta(:,:),MLHM(:,:),MLMW(:,:),thetase(:),AKMSE(:,:),HESS(:),
     4 psise(:),SSMOOTH(:,:)
      DOUBLE PRECISION,ALLOCATABLE::c(:,:,:),H(:,:,:),
	1 G(:,:,:),a(:,:),F(:,:,:),R(:,:,:)
      CHARACTER(12), ALLOCATABLE:: REAL_CLOCK(:)
	INTEGER ntf,nstot,nmis,indmis,IT,I,J,K,L1,jjj,IND,IFAIL,IMAX(1),
     1 IMIN(1),IMSVAR
	DOUBLE PRECISION AUX,lastl,lasth
	CHARACTER*1 DEB
      CHARACTER*3 DLLEXT
      CHARACTER*200 mfile,pathmfile
	CHARACTER*200 FILEIN,NMLNAME,PATH,FILEOUT,DMMTITLE,CURDIR
      CHARACTER*1024 matlaberror
	EXTERNAL GETARG
C TIME	
      ALLOCATE(np(3),ns(6),INFOS(9,6),IT1(7),IT2(7),DATE_ITIME(8),
     1 REAL_CLOCK(3))      
      CALL DATE_AND_TIME(REAL_CLOCK(1),REAL_CLOCK(2),REAL_CLOCK(3),
     1                   DATE_ITIME)      
      IT1(1:3) = DATE_ITIME(1:3)
      IT1(4:7) = DATE_ITIME(5:8)      
      
C GET the namelist specified by FILEIN
	DEB = 'R' 
	IF (DEB.EQ.'R') THEN  
	 CALL GETARG(1,FILEIN)  ! load name of input file
      ELSE
       FILEIN = 'H:\AROSSI\DMM\NILE\nile.nml'      
	ENDIF

C CHECK FILEIN
	IF (TRIM(FILEIN).EQ.'') THEN 
	 TYPE *, ' '
	 TYPE *, ' No input file provided'  
	 TYPE *, ' Program aborting'
	 PAUSE
	 STOP
	ENDIF
	
C LOAD input from FILEIN
	ALLOCATE(obs(30000),hyptheta(4,200),hypS(50,50,6),pdftheta(200))
      CALL input(FILEIN,NMLNAME,PATH,
	1 ny,nz,nx,nu,d,nv,ns,nstot,np,nf,INFOS,
     2 seed,thin,burnin,GGG,thetasampler,datasim,DLLNAME,check,
     3 estimation,nt,pdftheta,hyptheta,hypS,nobs,obs,
     4 Ssampler,HBL,MargLik)	      
      
C CHECK DLL NAME AND FIND FILE EXTENSION (.dll or .m)
      J = SCAN(DLLNAME,'\', BACK = .TRUE.)      
      I = SCAN(DLLNAME,'.', BACK = .TRUE.)	      
      DLLEXT = DLLNAME(I+1:I+3)               
      IF ((DLLEXT.EQ.'M  ').OR.(DLLEXT.EQ.'m  ')) THEN                 
       mfile     = DLLNAME(J+1:I-1)
       pathmfile = DLLNAME(1:J-1)              
C       DLLNAME   = 'H:\arossi\dmm64\matlabdll\debug\matlabdll.dll' ! provvisorio                             
       IND = GETCWD(CURDIR)  ! current directory
       DLLNAME = TRIM(CURDIR) // '\matlabdll.dll'                 ! definitivo
      ENDIF      
 
C FIND the DLL and LOAD it into the memory
      pdll = loadlibrary(DLLNAME)
	IF (pdll.EQ.0) THEN
	 TYPE *, ' '
	 TYPE *, TRIM(DLLNAME) // ' cannot be found or opened'  
	 TYPE *, ' Program aborting'
	 PAUSE
	 STOP
	ENDIF

C SET UP the pointer to the DLL function
	pdesign = getprocaddress(pdll, "design_"C)
	IF (pdesign.EQ.0) THEN
	 TYPE *, ' '
	 TYPE *, ' Sub DESIGN cannot be found into '// DLLNAME
	 TYPE *, ' Program aborting'
	 PAUSE
	 STOP
      ENDIF

C CHECK the MatLab file if needed
      IF ((DLLEXT.EQ.'M  ').OR.(DLLEXT.EQ.'m  ')) THEN
C SET UP the pointer to the DLL function
	 psetfilem = getprocaddress(pdll, "setfilem_"C)
	 IF (psetfilem.EQ.0) THEN
	  TYPE *, ' '
	  TYPE *, ' Sub SETFILEM cannot be found into '// DLLNAME
	  TYPE *, ' Program aborting'
	  PAUSE
	  STOP
       ENDIF

C SET UP the pointer to the DLL function
	 pgeterrstr = getprocaddress(pdll, "geterrstr_"C)
	 IF (pgeterrstr.EQ.0) THEN
	  TYPE *, ' '
	  TYPE *, ' Sub GETERRSTR cannot be found into '// DLLNAME
	  TYPE *, ' Program aborting'
	  PAUSE
	  STOP
       ENDIF
   
C Assign the name of the matlab file        
       ALLOCATE( c(ny,max(nz,1),ns(1)),H(ny,nx,ns(2)),
	1  G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)),R(nx,nu,ns(6)),
     1  theta(nt))
       CALL SETFILEM(mfile,pathmfile)  ! ONLY THE FIRST TIME 
       theta(:) = 1.D0
       CALL DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R)
       DEALLOCATE(c,H,G,a,F,R,theta)       
       IF (ny.EQ.0) THEN
        TYPE *, ' '
	  TYPE *, ' Can''t start MATLAB engine'
	  TYPE *, ' Program aborting'
	  PAUSE
	  STOP         
       ELSEIF (ny.EQ.-1) THEN
        TYPE *, ' '
	  TYPE *, ' Can''t read ny in the MATLAB file'
	  TYPE *, ' Program aborting'
	  PAUSE
	  STOP
       ELSEIF (ny.EQ.-2) THEN
        TYPE *, ' '
	  TYPE *, ' Can''t read nz in the MATLAB file'
	  TYPE *, ' Program aborting'
	  PAUSE
	  STOP
       ELSEIF (ny.EQ.-3) THEN
        TYPE *, ' '
	  TYPE *, ' Can''t read nx in the MATLAB file'
	  TYPE *, ' Program aborting'
	  PAUSE
	  STOP
       ELSEIF (ny.EQ.-4) THEN
        TYPE *, ' '
	  TYPE *, ' Can''t read nu in the MATLAB file'
	  TYPE *, ' Program aborting'
	  PAUSE
	  STOP
       ELSEIF (ny.EQ.-5) THEN
        TYPE *, ' '
	  TYPE *, ' Can''t read ns in the MATLAB file'
	  TYPE *, ' Program aborting'
	  PAUSE
	  STOP
       ELSEIF (ny.EQ.-6) THEN
        TYPE *, ' '
	  TYPE *, ' Can''t read nt in the MATLAB file'
	  TYPE *, ' Program aborting'
	  PAUSE
	  STOP
       ELSEIF (ny.EQ.-7) THEN
        TYPE *, ' '
	  TYPE *, ' Can''t find or open the MatLab function'
	  TYPE *, ' Program aborting'
	  PAUSE
	  STOP
       ELSEIF (ny.EQ.-8) THEN
        CALL GETERRSTR(matlaberror) 
        TYPE *, ' '
	  TYPE *, ' the MATLAB funtion can not be executed:'
        TYPE *, trim(matlaberror)
	  TYPE *, ' Program aborting'
	  PAUSE
	  STOP 
       ELSEIF (ny.LT.-100) THEN        
        TYPE *, ' '
	  TYPE *, ' One of the output canot be assigned during the call '
        TYPE *, ' ' // trim(DLLNAME)    
        TYPE *, ' Program aborting'
	  PAUSE
	  STOP
       ENDIF 
      ENDIF
      
C SET SHELL title        
	DMMTITLE = 'title DMM input:' // TRIM(PATH) // TRIM(NMLNAME) 
     #	     // '.nml' // ' - ' // TRIM(DLLNAME) 
	CALL system(DMMTITLE)

C INITIALISE THE RANDOM NUMBER GENERATOR
      CALL INITRAND(SEED)
	
C ASSIGN DATA and MISSING VALUES
	ALLOCATE(yk(nobs+nf,ny+nz),IYK(nobs,ny+1))
	K = 0
	DO 10 I = 1,nobs+nf
	DO 10 J = 1,ny+nz
	K = K+1
10	yk(I,J) = obs(K)

	IYK(:,:) = 0
	INDMIS   = 1 
	DO 11 I = 1,nobs
	 K = 0
	 DO 11 J = 1,ny
	 IF(yk(I,J).NE.-99999.D0) THEN
        K = K+1
	  IYK(I,K) = J
       ELSE
	  DO JJJ=1,nz
	   IF (yk(I,ny+JJJ).EQ.-99999.D0)indmis=0
        END DO
	 ENDIF
11	 IYK(I,ny+1) = K
	nmis = ny*nobs-SUM(IYK(1:nobs,ny+1))
	DEALLOCATE(obs)

C Allocate and Assign S 
	ALLOCATE(S(nobs,6),Z(nobs))
	S(1:nobs,1:6) = 1
		
C ASSIGN THETA-PRIORS
	ALLOCATE(thetaprior(nt,4))
	DO 30 I = 1,nt
30	thetaprior(I,1:4) = hyptheta(1:4,I)
	DEALLOCATE(hyptheta)
      
C ASSIGN PSI hyperparameters (# ind. Dirichelet x max # hyp)
	IF (nv.GT.0) THEN
	 ALLOCATE(psiprior(np(2),np(3)))	 
	 K = 0
	 DO I = 1,nv	  
	  IF (INFOS(9,I).EQ.1) THEN  ! S~iid
	   psiprior(K+1,1:INFOS(8,I)) = hypS(1:INFOS(8,I),1,I)
	   K = K+1
	  ELSEIF (INFOS(9,I).EQ.2) THEN  ! S~Markov 
	   DO J = 1,INFOS(8,I)
	    psiprior(K+J,1:INFOS(8,I)) = hypS(1:INFOS(8,I),J,I)
	   ENDDO
         K = K+INFOS(8,I)
	  ENDIF
	 END DO
	ENDIF
	DEALLOCATE(hypS)

C THETA STARTING VALUES & TRACK FREE PARAMETERS 
	ALLOCATE(theta0(nt),theta(nt),INDT(nt+2))
      CALL SIMPRIOR(estimation,nt,thetaprior,pdftheta(1:nt),ntf,INDT,
     1              theta0)
	theta(1:nt) = theta0(1:nt)
      
C PSI STARTING VALUES 
	IF (nv.GT.0) THEN
	 ALLOCATE(psi0(np(1)),psi(np(1)),ZW(2*nobs))
	ENDIF
	K = 0
	DO 80 J=1,nv
	 IF (INFOS(9,J).EQ.1) THEN  ! S-IID
	  CALL G05FEF(1.D0,1.D0,INFOS(8,J)-1,psi0(K+1:K+INFOS(8,J)-1),
	1              IFAIL) ! beta 
	  CALL G05FEF(1.D0,1.D0,1,AUX,IFAIL)
	  psi0(K+1:K+INFOS(8,J)-1) = psi0(K+1:K+INFOS(8,J)-1)/	  
     #  (SUM(psi0(K+1:K+INFOS(8,J)-1))+AUX) 
	  K = K + INFOS(8,J)-1
	 ELSE IF (INFOS(9,J).EQ.2) THEN  ! S-MARKOV
	  DO I = 1,INFOS(8,J)
	   CALL G05FEF(1.D0,1.D0,INFOS(8,J)-1,psi0(K+1:K+INFOS(8,J)-1),
	1               IFAIL)
	   CALL G05FEF(1.D0,1.D0,1,AUX,IFAIL)
	   psi0(K+1:K+INFOS(8,J)-1) = psi0(K+1:K+INFOS(8,J)-1)/	  
     #   (SUM(psi0(K+1:K+INFOS(8,J)-1))+AUX) 
	   K = K + INFOS(8,J)-1
	  ENDDO
80	 ENDIF
      
C WRITE HYPERPARAMTERS for THETA and PSI plus DATA
	FILEOUT = TRIM(PATH)//TRIM(NMLNAME)//'.PRI'
 	OPEN(10,FILE = FILEOUT, ACCESS='SEQUENTIAL')       
	 WRITE(10,'(<11+nv>(I6))') nt,np(1:3),nf,nz,seed,nx,ny,nobs,
	1                           nv,INFOS(8,1:nv)       
       WRITE(10,'(A2)') estimation
       DO I =1,nt
	  WRITE(10,1111) thetaprior(I,1:4),pdftheta(I)
       END DO       
	 K = 0
	 DO I = 1,nv
	  IF (INFOS(9,I).EQ.1) THEN
	    WRITE(10,1112) INFOS(8,I),psiprior(K+1,:),INFOS(9,I)
		K = K+1	    
	  ELSEIF (INFOS(9,I).EQ.2) THEN
	    DO J = 1,INFOS(8,I)
	     WRITE(10,1112) INFOS(8,I),psiprior(K+1,:),INFOS(9,I)
           K = K + 1
	    END DO
	  ENDIF
       END DO	 
       DO I =1,nobs+nf
	  WRITE(10,'(<ny+nz>(F20.10))') yk(I,1:ny+nz)
	 END DO
	CLOSE(10)
      
	ALLOCATE(STATE(nobs,nx),NEVAL(nt),CUMN(nt))
      
C CHECK DESIGN.dll
	IF ((check.EQ.'Y').OR.(check.EQ.'y')) THEN
	 CALL CHECKDESIGN(ny,nz,nx,nu,ns,nt,d,theta0,pdll,PATH,NMLNAME)
	 GOTO 7777 
	ENDIF
      
C SIMULATION of DATA and UNOBSERVABLES 
	IF ((datasim.EQ.'Y').OR.(datasim.EQ.'y')) THEN
	 CALL OPENFILES(ESTIMATION,SEED,NV,0,0,datasim,MARGLIK,
	1                PATH,NMLNAME)
	 CALL SIMDATA(nobs,d,ny,nz,nx,nu,ns,nstot,nt,nv,np,INFOS,pdll,
     2              theta0,psi0,Z,STATE,yk)
	 IF (nv.EQ.0) THEN
	  WRITE(9,'((F25.15))') theta0(1:nt)
	 ELSE
	  WRITE(9,'((F25.15))') theta0(1:nt),psi0(1:np(1))
	  WRITE(11,'(<1>(I3))') Z(:)
	 ENDIF	 
	 WRITE(10,'(<nx>(F20.10))') (STATE(I,1:nx),I=1,nobs) 
	 WRITE(15,'(<ny>(F20.10))') (yk(I,1:ny),I=1,nobs)
	 CLOSE(9)  
	 CLOSE(10)  
	 CLOSE(11)  
	 CLOSE(15)  
	 GOTO 7777 
	ENDIF

C MAXIMUM LIKELIHOOD ESTIMATION      
	IF ((estimation.EQ.'ML').OR.(estimation.EQ.'ml').OR.
     &    (estimation.EQ.'Ml').OR.(estimation.EQ.'mL')) THEN
       CALL OPENFILES(estimation,seed,nv,nf,0,datasim,marglik,
	1                path,nmlname)
       ALLOCATE(HESS((nt+np(1))*(nt+np(1)+1)/2))
	 CALL ML(nobs,d,ny,nz,nx,nu,nt,nv,ns,np(1),INFOS,pdll,INDT,yk,IYK,S,
	1         thetaprior,theta0,psi0,IMSVAR,HESS,AUX)       
	 ALLOCATE(THETASE(nt),AKMSE(nobs,nx),INN(nobs,ny))
	 IF (nv.EQ.0) THEN        
        CALL OPG(nobs,d,ny,nz,nx,nu,nt,ns,pdll,yk,IYK,S,
	1           theta0,thetaprior,HESS,thetase,STATE,AKMSE,INN,IFAIL)
        WRITE(9,'(<2>(F25.15))') (theta0(I),thetase(I),I=1,nt)
	  WRITE(9,'(<2>(F25.15))') AUX,IFAIL
	  WRITE(10,'(<nx>(F20.10))') (STATE(I,1:nx),I=1,nobs)
        WRITE(10,'(<nx>(F20.10))') (AKMSE(I,1:nx),I=1,nobs)
        WRITE(12,'(<ny>(F20.10))') (INN(I,1:ny),I=1,nobs)
       ELSE           
        ALLOCATE(psise(np(1)),SSMOOTH(nobs,nstot))   
        IF(IMSVAR.EQ.1)THEN
         CALL OPGH(nobs,ny,nz,nx,nu,nt,nv,ns,nstot,np(1),pdll,yk,IYK,
     1             INFOS,theta0,psi0,thetaprior,HESS,thetase,psise,
     1             SSMOOTH,INN,IFAIL)
        ELSE            
         CALL OPGKIM(nobs,d,ny,nz,nx,nu,nt,nv,ns,nstot,np(1),pdll,
     1               yk,IYK,INFOS,theta0,psi0,thetaprior,HESS,
     1               thetase,psise,STATE,AKMSE,SSMOOTH,INN,IFAIL)
	   WRITE(10,'(<nx>(F20.10))') (STATE(I,1:nx),I=1,nobs)
         WRITE(10,'(<nx>(F20.10))') (AKMSE(I,1:nx),I=1,nobs)
        ENDIF        
        WRITE(9,'(<2>(F25.15))') (theta0(I),thetase(I),I=1,nt)
        WRITE(9,'(<2>(F25.15))') (psi0(I),psise(I),I=1,np(1))
	  WRITE(9,'(<2>(F25.15))') AUX,IFAIL
        WRITE(11,'(<nstot>(F20.10))') (SSMOOTH(I,1:nstot),I=1,nobs)
        WRITE(12,'(<ny>(F20.10))') (INN(I,1:ny),I=1,nobs)
        CLOSE(11)         
        DEALLOCATE(PSISE,SSMOOTH,HESS)
       ENDIF
	 CLOSE(9)
       CLOSE(10)
       CLOSE(12)
	 DEALLOCATE(THETASE,AKMSE,INN)
     	 GOTO 7777 
	ENDIF
 
C MCMC BURN-IN
	IF ((nv.GT.0).AND.(HBL.GT.1)) THEN
	 ALLOCATE(PTR(nobs,nstot,nstot),PM(nobs,nstot),ACCRATE(nobs))
	 PM(:,:)    = 1.D0/DFLOAT(nstot)
       PTR(:,:,:) = 1.D0/DFLOAT(nstot)
	 ACCRATE(:) = 0
	ENDIF
	CUMN(1:nt) = 0
	NEVAL(1:nt)= 0
	IND   = 100
	Z(:)  = 1
	ZW(:) = 1
      L1    = 0
	IF (nmis.GT.0) THEN  ! MISSINGS
	 DO jjj = 1,burnin
	  IF (nv.GT.0) THEN		    
	   CALL GCK(nobs,d,ny,nz,nx,nu,nv,ns,nstot,nt,np(1),
     1            yk(1:nobs,:),IYK(1:nobs,:),theta0,psi0,
     2            INFOS,pdll,Z,S)
	   IF (HBL.GT.1) THEN
		CALL RECPR(jjj,nstot,nobs,Z,ZW,PM,PTR)
	   ENDIF
	   CALL DRAWPSI(nobs,nv,np,INFOS,Z,psiprior,psi0,psi)	  
	  ENDIF
	  DO it = 1,nt
	   IF (thetaprior(it,3).LT.thetaprior(it,4)) THEN
	    CALL SLICE(it,nobs,d,ny,nz,nx,nu,ns,nt,S,
	1               yk(1:nobs,:),IYK(1:nobs,:),theta0,
	2               thetaprior(it,:),pdftheta(it),pdll,
     3               NEVAL(it),theta(it))
          theta0(it) = theta(it) 
         ENDIF  
	  END DO  
	  CUMN = CUMN + NEVAL
	  IF (jjj/IND*IND.EQ.jjj) THEN
	   IMIN    = MINLOC(CUMN(INDT(1:ntf)))
	   IMIN(1) = CUMN(INDT(IMIN(1))) !CUMN(IMIN(1))
	   IMAX    = MAXLOC(CUMN(INDT(1:ntf)))
	   IMAX(1) = CUMN(INDT(IMAX(1))) !CUMN(IMAX(1))
	   CALL system('cls')
	   WRITE(6,1113) jjj,ntf,IMIN(1)/dfloat(jjj),
     #	             IMAX(1)/dfloat(jjj)
        ENDIF
       ENDDO
	ELSE  ! NO MISSING 
	 DO jjj = 1,burnin
	  IF (nv.GT.0) THEN
	   CALL GCK2(nobs,d,ny,nz,nx,nu,nv,ns,nstot,nt,np(1),
	1             yk(1:nobs,:),theta0,psi0,INFOS,pdll,Z,S)
	   IF (HBL.GT.1) THEN
	    CALL RECPR(jjj,nstot,nobs,Z,ZW,PM,PTR)
	   ENDIF
   	   CALL DRAWPSI(nobs,nv,np,INFOS,Z,psiprior,psi0,psi)	  
	  ENDIF
        DO it = 1,nt
	   IF (thetaprior(it,3).LT.thetaprior(it,4)) THEN
	    CALL SLICE2(it,nobs,d,ny,nz,nx,nu,ns,nt,S,yk(1:nobs,:),
	1              	theta0,thetaprior(it,:),pdftheta(it),pdll,
     2                NEVAL(it),theta(it))
          theta0(it) = theta(it) 
         ENDIF
	  END DO
	  CUMN = CUMN + NEVAL
	  IF (jjj/IND*IND.EQ.jjj) THEN
	   IMIN    = MINLOC(CUMN(INDT(1:ntf)))
	   IMIN(1) = CUMN(INDT(IMIN(1))) !CUMN(IMIN(1))
	   IMAX    = MAXLOC(CUMN(INDT(1:ntf)))
	   IMAX(1) = CUMN(INDT(IMAX(1))) !CUMN(IMAX(1))
	   CALL system('cls')
	   WRITE(6,1113) jjj,ntf,IMIN(1)/dfloat(jjj),
     #	             IMAX(1)/dfloat(jjj)
        ENDIF
       ENDDO
	ENDIF
	lastl = IMIN(1)/DFLOAT(burnin)
	lasth = IMAX(1)/DFLOAT(burnin)
	CUMN(1:nt) = 0
	NEVAL(:)   = 0 
	
C OPEN OUTPUT FILES
C 9 '.PAR', 10 '.UNB', 11 '.DIS', 12 '.INN', 13 '.FST', 14 '.MIS', 15 '.ML' o '.DAT'
      CALL OPENFILES(estimation,seed,nv,nf,nmis*INDMIS,datasim,Marglik,
	1               PATH,NMLNAME)

C MCMC RECORDING phase
	ALLOCATE(INN(nobs,ny))
	IF (nf.GT.0) THEN
       ALLOCATE(FORE(nf,ny+nx+1))
	ENDIF
	IF (indmis*NMIS.GE.1) THEN 
	 ALLOCATE(ykmis(nmis))
	ENDIF
	IF ((MargLik.EQ.'Y').OR.(MargLik.EQ.'y')) THEN
	 ALLOCATE(gibtheta(GGG,nt+np(1)),gibZ(GGG,nobs),MLHM(11,2),
	1          MLMW(2,2))
	ENDIF
	IF (nmis.GT.0) THEN ! MISSINGS 
	 DO jjj = 1,GGG*thin 
	  IF (nv.GT.0) THEN
	   IF (HBL.EQ.1) THEN
	    CALL GCK(nobs,d,ny,nz,nx,nu,nv,ns,nstot,nt,np(1),
	1             yk(1:nobs,:),IYK(1:nobs,:),theta0,psi0,
     2             INFOS,pdll,Z,S)
	   ELSE
	    CALL AMH(HBL,nobs,d,ny,nz,nx,nu,nv,ns,nstot,nt,np(1),
	1             yk(1:nobs,:),IYK(1:nobs,:),theta0,psi0,
	2    		 PTR,PM,INFOS,pdll,Z,S,ACCRATE)
		CALL RECPR(jjj+burnin,nstot,nobs,Z,ZW,PM,PTR)
	   ENDIF
	   CALL DRAWPSI(nobs,nv,np,INFOS,Z,psiprior,psi0,psi)
	  ENDIF 
	  DO it = 1,nt
	   IF (thetaprior(it,3).LT.thetaprior(it,4)) THEN
	    CALL SLICE(it,nobs,d,ny,nz,nx,nu,ns,nt,S,yk(1:nobs,:),
	1               IYK(1:nobs,:),theta0,thetaprior(it,:),
     2               pdftheta(it),pdll,NEVAL(it),theta(it))
          theta0(it) = theta(it) 
         ENDIF
	  END DO
	  CUMN = CUMN+NEVAL 
	  CALL SIMSTATE(nobs,d,ny,nz,nx,nu,ns,nt,yk(1:nobs,:),
	1                IYK(1:nobs,:),theta,S,pdll,STATE)
	  CALL INNOV(nobs,d,ny,nz,nx,nu,ns,nt,S,
	1             yk(1:nobs,:),IYK(1:nobs,:),theta,pdll,INN)
	  IF (nf.GT.0) THEN 
	   CALL FORECAST(yk(nobs+1:nobs+nf,ny+1:ny+nz),nf,ny,nz,nx,nu,nv,
	1                 ns,nstot,nt,np,theta,psi,INFOS,Z(nobs),
     2                 STATE(nobs,:),pdll,FORE)
	  ENDIF
	  IF (INDMIS*nmis.GE.1) THEN
	   J = 1
	   DO I = 1,nobs
	    IF (IYK(I,ny+1).LT.ny) THEN
		  K = ny-IYK(I,ny+1)	    
	      CALL MISSING(yk(I,:),ny,nz,nx,nu,nv,ns,nt,K,theta,INFOS,
	1                  S(I,1:6),STATE(I,:),pdll,ykmis(J:J+K-1))
	     J = J+K
	    ENDIF 
         ENDDO 
	  ENDIF
	  IF (jjj/IND*IND.EQ.jjj) THEN
	   IMIN    = MINLOC(CUMN(INDT(L1+1:ntf)))
	   IMIN(1) = CUMN(INDT(L1+IMIN(1)))
	   IMAX    = MAXLOC(CUMN(INDT(L1+1:ntf)))
	   IMAX(1) = CUMN(INDT(L1+IMAX(1)))
	   CALL system('cls')
	   WRITE(6,1113) BURNIN,ntf,lastl,lasth
	   IF ((HBL.EQ.1).OR.(nv.EQ.0)) THEN
	    WRITE(6,1114) jjj,ntf,IMIN(1)/dfloat(jjj),
     #	              IMAX(1)/dfloat(jjj)
         ELSEIF ((HBL.GT.1).AND.(nv.GT.0)) THEN
	    WRITE(6,1115) jjj,ntf,IMIN(1)/dfloat(jjj),
     #           IMAX(1)/dfloat(jjj),
     #           SUM(1.D0-ACCRATE(1:nobs)/DFLOAT(jjj))/DFLOAT(nobs) 
	   ENDIF 
        ENDIF
        IF (jjj/thin*thin.EQ.jjj) THEN	   
	   WRITE(12,'(<nobs*ny>(F20.10))') (INN(1:nobs,I),I=1,ny)
	   WRITE(10,'(<nobs*nx>(F20.10))') (STATE(1:nobs,I),I=1,nx)
	   IF ((MargLik.EQ.'Y').OR.(MargLik.EQ.'y')) THEN
	    gibtheta(jjj/thin,1:nt) = theta(1:nt)	    
	   ENDIF
	   IF (nv.EQ.0) THEN
	    WRITE(9,'(<nt>(F25.15))') theta(1:nt)	    
	   ELSE
	   IF ((MargLik.EQ.'Y').OR.(MargLik.EQ.'y')) THEN
	    gibtheta(jjj/thin,nt+1:nt+np(1)) = psi(1:np(1))	    
	    gibZ(jjj/thin,1:nobs) = Z(1:nobs)
	   ENDIF
	   WRITE(9,'(<nt+np(1)>(F25.15))') theta(1:nt),psi(1:np(1))
	   WRITE(11,'(<nobs>(I3))') Z(:)
	  ENDIF
	  IF (nf.GT.0) THEN
	   J = min(nv,1)
	   WRITE(13,'(<nf*(nx+ny+J)>(F20.10))') (FORE(1:nf,I),I=1,nx+ny+J)	  
	  ENDIF
	  IF (INDMIS*nmis.GE.1) WRITE(14,'(<nmis>(F20.10))') ykmis(1:nmis)	  
	  ENDIF
	 ENDDO
	ELSE  ! NO MISSINGS
	 DO jjj = 1,GGG*thin 
        IF (nv.GT.0) THEN
	   IF (HBL.EQ.1) THEN
	    CALL GCK2(nobs,d,ny,nz,nx,nu,nv,ns,nstot,nt,np(1),
	1              yk(1:nobs,:),theta0,psi0,INFOS,pdll,Z,S)
	   ELSE
		CALL AMH2(hbl,nobs,d,ny,nz,nx,nu,nv,ns,nstot,nt,np(1),
	1              yk(1:nobs,:),theta0,psi0,
	2    		  PTR,PM,INFOS,pdll,Z,S,ACCRATE)
		CALL RECPR(jjj+burnin,nstot,nobs,Z,ZW,PM,PTR)
	   ENDIF
	   CALL DRAWPSI(nobs,nv,np,INFOS,Z,psiprior,psi0,psi)
	  ENDIF 
	  DO it = 1,nt
	   IF (thetaprior(it,3).LT.thetaprior(it,4)) THEN
	    CALL SLICE2(it,nobs,d,ny,nz,nx,nu,ns,nt,S,yk(1:nobs,:),
	1             	theta0,thetaprior(it,:),pdftheta(it),pdll,
     2                NEVAL(it),theta(it))
          theta0(it) = theta(it) 
         ENDIF
	  END DO
	  CUMN = CUMN+NEVAL 
	  CALL SIMSTATE2(nobs,d,ny,nz,nx,nu,ns,nt,yk(1:nobs,:),
	1                 theta,S,pdll,STATE)
	  CALL INNOV2(nobs,d,ny,nz,nx,nu,ns,nt,S,
	1              yk(1:nobs,:),theta,pdll,INN)
	  IF (nf.GT.0) THEN 
	   CALL FORECAST(yk(nobs+1:nobs+nf,ny+1:ny+nz),nf,ny,nz,nx,nu,nv,
	1                 ns,nstot,nt,np,theta,psi,INFOS,Z(nobs),
     2                 STATE(nobs,:),pdll,FORE)
	  ENDIF
	  IF (jjj/IND*IND.EQ.jjj) THEN
	   IMIN    = MINLOC(CUMN(INDT(L1+1:ntf)))
	   IMIN(1) = CUMN(INDT(L1+IMIN(1)))
	   IMAX    = MAXLOC(CUMN(INDT(L1+1:ntf)))
	   IMAX(1) = CUMN(INDT(L1+IMAX(1)))
	   CALL system('cls')
	   WRITE(6,1113) BURNIN,ntf,lastl,lasth
	   IF ((HBL.EQ.1).OR.(nv.EQ.0)) THEN
	    WRITE(6,1114) jjj,ntf,IMIN(1)/dfloat(jjj),
     #	              IMAX(1)/dfloat(jjj)
         ELSEIF ((HBL.GT.1).AND.(nv.GT.0)) THEN
	    WRITE(6,1115) jjj,ntf,IMIN(1)/dfloat(jjj),
     #           IMAX(1)/dfloat(jjj),
     #           SUM(1.D0-ACCRATE(1:nobs)/DFLOAT(jjj))/DFLOAT(nobs) 
	   ENDIF         
        ENDIF
        IF (jjj/thin*thin.EQ.jjj) THEN	   
	   IF ((MargLik.EQ.'Y').OR.(MargLik.EQ.'y')) THEN
	    gibtheta(jjj/thin,1:nt) = theta(1:nt)	    
	   ENDIF 
	   WRITE(12,'(<nobs*ny>(F20.10))') (INN(1:nobs,I),I=1,ny)
	   WRITE(10,'(<nobs*nx>(F20.10))') (STATE(1:nobs,I),I=1,nx)
	   IF (nv.EQ.0) THEN	    
	    WRITE(9,'(<nt>(F25.15))') theta(1:nt)
	   ELSE
	    IF ((MargLik.EQ.'Y').OR.(MargLik.EQ.'y')) THEN
	     gibtheta(jjj/thin,nt+1:nt+np(1)) = psi(1:np(1))	    
	     gibZ(jjj/thin,1:nobs) = Z(1:nobs)
	    ENDIF
	    WRITE(9,'(<nt+np(1)>(F25.15))') theta(1:nt),psi(1:np(1))
	    WRITE(11,'(<nobs>(I3))') Z(:)
	   ENDIF
	   IF (nf.GT.0) THEN
	   J = min(nv,1)
	   WRITE(13,'(<nf*(nx+ny+J)>(F20.10))') (FORE(1:nf,I),I=1,nx+ny+J)
	   ENDIF
	  ENDIF
	 ENDDO
	ENDIF
	CLOSE(9)
	CLOSE(10)
	IF (nv.GT.0) CLOSE(11)
	CLOSE(12)
	IF (nf.GT.0) CLOSE(13)
	IF (indmis*nmis.GE.1) THEN 
	 CLOSE(14)
	 DEALLOCATE(ykmis)
	ENDIF
	DEALLOCATE(INN)
	IF ((nv.GT.0).AND.(HBL.GT.1)) THEN
	 DEALLOCATE(PTR,PM,ACCRATE)
	ENDIF

C MARGINAL LIKELIHOOD
	IF ((MargLik.EQ.'Y').OR.(MargLik.EQ.'y')) THEN
	 WRITE(*,*) ' '
       WRITE(*,*) 'Computing the marginal likelihood. Please wait ...' 
	 IF (nmis.GT.0) THEN
	  CALL HARMONIC(GGG,nobs,d,ny,nz,nx,nu,nv,ns,nstot,nt,np,
	1                INFOS,yk(1:nobs,:),IYK(1:nobs,:),gibtheta,gibZ,
     2                thetaprior,psiprior,pdftheta,pdll,MLHM)
	  WRITE(*,*) 'Modified harmonic mean: done!'	   
	  CALL MENGWONG(GGG,nobs,d,ny,nz,nx,nu,nv,ns,nstot,nt,np,
	1                INFOS,yk(1:nobs,:),IYK(1:nobs,:),gibtheta,gibZ,
     2                thetaprior,psiprior,pdftheta,pdll,MLHM(5,1),MLMW)
        WRITE(*,*) 'Bridge sampling: done!'
	  WRITE(*,*) ' '
	 ELSE
	  CALL HARMONIC2(GGG,nobs,d,ny,nz,nx,nu,nv,ns,nstot,nt,np,
	1                 INFOS,yk(1:nobs,:),gibtheta,gibZ,thetaprior,
     2                 psiprior,pdftheta,pdll,MLHM)
	  WRITE(*,*) 'Modified harmonic mean: done!'	   
	  CALL MENGWONG2(GGG,nobs,d,ny,nz,nx,nu,nv,ns,nstot,nt,np,
	1                 INFOS,yk(1:nobs,:),gibtheta,gibZ,thetaprior,
     2                 psiprior,pdftheta,pdll,MLHM(5,1),MLMW)
        WRITE(*,*) 'Bridge sampling: done!'
	  WRITE(*,*) ' '
	 ENDIF
	 WRITE(15,*) 'Modified Harmonic mean (ML and Var)'
	 WRITE(15,'(<2>(F20.10))') (MLHM(I,:),I=1,11)
	 WRITE(15,*) 'Bridge Sampling'
	 WRITE(15,'(<2>(F20.10))') (MLMW(I,:),I=1,2)
	 CLOSE(15)
	 DEALLOCATE(gibtheta,gibZ,MLHM,MLMW)
      ENDIF     
      
7777	DEALLOCATE(yk,STATE,Z,S,theta0,theta,thetaprior,pdftheta,NEVAL,
	1           CUMN,IYK,INDT)
      IF (nv.GT.0) THEN
	 DEALLOCATE(psi0,psi,psiprior,ZW)
	ENDIF
	
	STATUS = freelibrary(pdll) !libero la DLL dalla memoria alla fine del programma
	IF (TRIM(PATH).EQ.'') THEN
	 STATUS = getcwd(PATH) ! get current directory
      ENDIF

      CALL DATE_AND_TIME(REAL_CLOCK(1),REAL_CLOCK(2),REAL_CLOCK(3),
     1                   DATE_ITIME)
      IT2(1:3) = DATE_ITIME(1:3)
      IT2(4:7) = DATE_ITIME(5:8)
	IT=(IT2(4)-IT1(4))*3600+(IT2(5)-IT1(5))*60+(IT2(6)-IT1(6))
      IF ((check.EQ.'Y').OR.(check.EQ.'y')) THEN
        WRITE(6,1117) TRIM(PATH)
      ELSE
          IF ((datasim.EQ.'Y').OR.(datasim.EQ.'y')) THEN
            WRITE(6,1118) TRIM(PATH)
          ELSE
            IF ((estimation.EQ.'ML').OR.(estimation.EQ.'ml').OR.
     &          (estimation.EQ.'Ml').OR.(estimation.EQ.'mL')) THEN
              WRITE(6,1119) IT,TRIM(PATH)    
            ELSE
              WRITE(6,1116) IT,TRIM(PATH)
            ENDIF              
          ENDIF
      ENDIF
      DEALLOCATE(np,ns,INFOS,IT1,IT2,DATE_ITIME,REAL_CLOCK)       

1111  FORMAT((<4>(F25.12)), '  ',A2)
1112  FORMAT(I10,(<np(3)>(F25.12)), '  ',I2)
1113  FORMAT(/,' Burn-in draws = ',I8,
     #       /,' Parameters sampled by SLICE ',I5,
     #       /,' SLICE likelihood eval. Min/Max = ',F6.2, ' / ',F6.2)
1114  FORMAT(/,' Recording draws = ',I8,
     #       /,' Parameters sampled by SLICE ',I5,
     #       /,' SLICE likelihood eval. Min/Max = ',F6.2, ' / ',F6.2)
1115  FORMAT(/,' Recording draws = ',I8,
     #       /,' Parameters sampled by SLICE ',I5,
     #       /,' SLICE likelihood eval. Min/Max = ',F6.2, ' / ',F6.2,
     #       /,' Adaptive MH accettance rate = ',F6.2)
1116  FORMAT(/,' MCMC completed',
     #       /,' CPU-time (sec)=', I10,
     #       /,' Output printed in ',A)
1117  FORMAT(/,' Check completed',
     #       /,' Output printed in ',A)
1118  FORMAT(/,' Data simulation completed',
     #       /,' Output printed in 'A)
1119  FORMAT(/,' Maximum Likelihood completed',
     #       /,' CPU-time (sec)=', I10,
     #       /,' Output printed in ',A)

	STOP
      END