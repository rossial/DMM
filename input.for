C -----------------------------------------------------------------------
C INPUT reads the input file *.nml
C Developed by A.Rossi, C.Planas and G.Fiorentini
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
C -----------------------------------------------------------------------
	SUBROUTINE input(FILEIN,NMLNAME,PATH,ny,nz,nx,nu,d,nv,ns,
	1 nstot,np,nf,INFOS,seed,thin,burnin,simulrec,sampler,
     2 datasim,dllname,check,estimation,nt,pdftheta,hyptheta,
     3 hypS,T,obs,Ssampler,hbl,MargLik)

C	INCLUDE 'iosdef.for'
	INTEGER IERR
C INPUT
	CHARACTER*200 FILEIN
C OUTPUT
	NAMELIST /ssm/ nx,nu,d,nv,dllname,check,estimation
	INTEGER nx,nu,d(2),nv
	CHARACTER*200 dllname
	CHARACTER*1 check
	CHARACTER*2 estimation

	NAMELIST /S1/ dynS1,matS1,ns1,hypS1
	NAMELIST /S2/ dynS2,matS2,ns2,hypS2
	NAMELIST /S3/ dynS3,matS3,ns3,hypS3
	NAMELIST /S4/ dynS4,matS4,ns4,hypS4
	NAMELIST /S5/ dynS5,matS5,ns5,hypS5
	NAMELIST /S6/ dynS6,matS6,ns6,hypS6
	CHARACTER*1 dynS1,dynS2,dynS3,dynS4,dynS5,dynS6,
	1 matS1(6),matS2(6),matS3(6),matS4(6),matS5(6),matS6(6)
	INTEGER ns1,ns2,ns3,ns4,ns5,ns6
	DOUBLE PRECISION hypS1(50,50),hypS2(50,50),hypS3(50,50),
	1 hypS4(50,50),hypS5(50,50),hypS6(50,50),hypS(50,50,6)

	NAMELIST /mcmc/ seed,thin,burnin,simulrec,sampler,
	1 Ssampler,hbl,marglik
	INTEGER seed,thin,burnin,simulrec,hbl
	CHARACTER*1 MargLik
	CHARACTER*2 sampler
      CHARACTER*3 Ssampler

	NAMELIST /prior/ nt,pdftheta,hyptheta
	INTEGER nt
	DOUBLE PRECISION hyptheta(4,200)
	CHARACTER*2  pdftheta(200)

	NAMELIST /dataset/ T,ny,nz,nf,datasim,obs
	INTEGER T,ny,nz,nf
	DOUBLE PRECISION obs(30000)
	CHARACTER*1 datasim

C OUTPUT not in the namelist
	INTEGER ns(6),INFOS(9,6),IMAX(1),np(3),nstot,IND
	CHARACTER*200 NMLNAME,PATH

C LOCALS
	CHARACTER*3 IC
	CHARACTER*200 STR
	INTEGER I,J,K,IFAIL
	INTEGER, ALLOCATABLE:: nstateS(:)
	CHARACTER*1, ALLOCATABLE:: matS(:,:),dynS(:)

C IDENTIFY the PATH and NAME of the .NML INPUT FILE
	I = SCAN(FILEIN,'\', BACK = .TRUE.)
	IF((I.LE.0).OR.(I.GE.200)) THEN
	 NMLNAME = FILEIN
	 PATH    = ''
	ELSE
	 NMLNAME = FILEIN(I+1:200)
	 PATH    = FILEIN(1:I)
	ENDIF
	I = SCAN(NMLNAME,'.', BACK = .TRUE.)
      NMLNAME = NMLNAME(1:I-1)

C FIND namelist ssm
	OPEN(1,File=TRIM(FILEIN), ACCESS='SEQUENTIAL',
	1     STATUS='OLD',IOSTAT=IERR, ERR=5000)
	IFAIL = -1
	DO WHILE (.NOT.EOF(1))
	 READ(1,'(A)') STR
	 IF (INDEX(STR,'&ssm').GT.0) IFAIL = 0
	ENDDO
	CLOSE(1)
	IF (IFAIL.EQ.-1) THEN
	 TYPE *, ' Namelist ssm not found'
	 TYPE *, ' Program aborting'
	 PAUSE
	 RETURN
	ENDIF

C READ namelist ssm
      ns1=0
      ns2=0
      ns3=0
      ns4=0
      ns5=0
      ns6=0
	OPEN(1,File=TRIM(FILEIN), ACCESS='SEQUENTIAL')
	nx      = -1
	nu      = -1
	d(:)    = -1
	nv      =  0
	dllname = ''
	check   = 'N'
	estimation = 'BA'
	READ(1,NML=ssm,END=5001,ERR=5001)
	IF (nx.LE.0) THEN
	 TYPE *, ' Check nx in namelist ssm'
	 TYPE *, ' Program aborting'
	 PAUSE
	 STOP
	ENDIF
	IF(nu.LE.0) THEN
	 TYPE *, ' Check nu in namelist ssm'
	 TYPE *, ' Program aborting'
	 PAUSE
	 STOP
	ENDIF
	IF((d(1).LT.0).OR.(d(2).LT.0).OR.(d(2).GT.nx)) THEN
	 TYPE *, ' Check d in namelist ssm'
	 TYPE *, ' Program aborting'
	 PAUSE
	 STOP
	ENDIF
	IF((nv.LT.0).OR.(nv.GT.6)) THEN
	 TYPE *, ' Check nv in namelist ssm'
	 TYPE *, ' Program aborting'
	 PAUSE
	 STOP
      ENDIF
      IF(dllname.EQ.'') THEN
	 TYPE *, ' Check dllname in namelist ssm'
	 TYPE *, ' Program aborting'
	 PAUSE
	 STOP
	ENDIF
	CLOSE(1)
	IF ((estimation.NE.'ML').AND.(estimation.NE.'ml').AND.
     &    (estimation.NE.'Ml').AND.(estimation.NE.'mL').AND.
     &    (estimation.NE.'BA').AND.(estimation.NE.'ba').AND.
     &    (estimation.NE.'Ba').AND.(estimation.NE.'bA')) THEN
	 estimation = 'BA'
	ENDIF

	IF (nv.GT.0) THEN
C FIND namelist S1
	 OPEN(1,File=TRIM(FILEIN), ACCESS='SEQUENTIAL')
	 IFAIL = -1
	 DO WHILE (.NOT.EOF(1))
	  READ(1,'(A)') STR
	  IF (INDEX(STR,'&S1').GT.0) IFAIL = 0
	 ENDDO
	 CLOSE(1)
	 IF (IFAIL.EQ.-1) THEN
	  TYPE *, ' Namelist S1 not found'
	  TYPE *, ' Program aborting'
	  PAUSE
	  RETURN
	 ENDIF

C READ namelist S1
	 OPEN(1,File=TRIM(FILEIN), ACCESS='SEQUENTIAL')
	 dynS1      = '-'
	 ns1        = -1
	 hypS1(:,:) = 1
	 matS1(:)   = '-'
	 READ(1,NML=S1,END=5002,ERR=5002)

	 IF ((dynS1.NE.'I').AND.(dynS1.NE.'M'))THEN
	  TYPE *, ' Check dynS1 in namelist S1'
	  TYPE *, ' Program aborting'
	  PAUSE
	  STOP
	 ENDIF

	 IF (ns1.LT.2) THEN
	  TYPE *, ' Check ns1 in namelist S1'
	  TYPE *, ' Program aborting'
	  PAUSE
	  STOP
	 ENDIF

	 IF (dynS1.EQ.'I') THEN
        DO J = 1,ns1
	   IF (hypS1(J,1).LE.0) THEN
	    TYPE *, ' Check hypS1 in namelist S1'
		TYPE *, ' Program aborting'
	    PAUSE
		STOP
	   ENDIF
        ENDDO
	 ELSEIF (dynS1.EQ.'M') THEN
	  DO J = 1,ns1
	  DO K = 1,ns1
	   IF (hypS1(J,K).LE.0) THEN
	    TYPE *, ' Check hypS1 in namelist S1'
		TYPE *, ' Program aborting'
	    PAUSE
		STOP
	   ENDIF
        ENDDO
	  ENDDO
	 ENDIF

	 IF ((matS1(1).EQ.'-').OR.((matS1(1).NE.'a')
     #	  .AND.(matS1(1).NE.'H').AND.(matS1(1).NE.'G')
     #	  .AND.(matS1(1).NE.'c').AND.(matS1(1).NE.'F')
     #	  .AND.(matS1(1).NE.'R'))) THEN
	   TYPE *, ' Check matS1 in namelist S1'
	   TYPE *, ' Program aborting'
	   PAUSE
	   STOP
       ENDIF
	 DO I = 2,6
	  IF ((matS1(I).NE.'-').AND.(matS1(I).NE.'a')
     #	  .AND.(matS1(I).NE.'H').AND.(matS1(I).NE.'G')
     #	  .AND.(matS1(I).NE.'c').AND.(matS1(I).NE.'F')
     #	  .AND.(matS1(I).NE.'R')) THEN
	   TYPE *, ' Check matS1 in namelist S1'
	   TYPE *, ' Program aborting'
	   PAUSE
	   STOP
	  ENDIF
	 ENDDO
	ENDIF
	CLOSE(1)

	IF (nv.GT.1) THEN
C FIND namelist S2
	 OPEN(1,File=TRIM(FILEIN), ACCESS='SEQUENTIAL')
	 IFAIL = -1
	 DO WHILE (.NOT.EOF(1))
	  READ(1,'(A)') STR
	  IF (INDEX(STR,'&S2').GT.0) IFAIL = 0
	 ENDDO
	 CLOSE(1)
	 IF (IFAIL.EQ.-1) THEN
	  TYPE *, ' Namelist S2 not found'
	  TYPE *, ' Program aborting'
	  PAUSE
	  RETURN
	 ENDIF
C READ namelist S2
	 OPEN(1,File=TRIM(FILEIN), ACCESS='SEQUENTIAL')
	 dynS2      = '-'
	 ns2   = -1
	 hypS2(:,:) = 1  ! Uniform
	 matS2(:)   = '-'
	 READ(1,NML=S2,END=5003,ERR=5003)

	 IF ((dynS2.NE.'I').AND.(dynS2.NE.'M'))THEN
	  TYPE *, ' Check dynS2 in namelist S2'
	  TYPE *, ' Program aborting'
	  PAUSE
	  STOP
	 ENDIF

	 IF (ns2.LT.2) THEN
	  TYPE *, ' Check ns2 in namelist S2'
	  TYPE *, ' Program aborting'
	  PAUSE
	  STOP
	 ENDIF

	 IF (dynS2.EQ.'I') THEN
        DO J = 1,ns2
	   IF (hypS2(J,1).LE.0) THEN
	    TYPE *, ' Check hypS2 in namelist S2'
	    TYPE *, ' Program aborting'
	    PAUSE
		STOP
	   ENDIF
        ENDDO
	 ELSEIF (dynS2.EQ.'M') THEN
	  DO J = 1,ns2
	  DO K = 1,ns2
	   IF (hypS2(J,K).LE.0) THEN
	    TYPE *, ' Check hypS2 in namelist S2'
	    TYPE *, ' Program aborting'
	    PAUSE
		STOP
	   ENDIF
        ENDDO
	  ENDDO
	 ENDIF
	 I = 1
	 IF ((matS2(I).EQ.'-').OR.((matS2(I).NE.'a')
     #	  .AND.(matS2(I).NE.'H').AND.(matS2(I).NE.'G')
     #	  .AND.(matS2(I).NE.'c').AND.(matS2(I).NE.'F')
     #	  .AND.(matS2(I).NE.'R'))) THEN
	   TYPE *, ' Check matS2 in namelist S2'
	   TYPE *, ' Program aborting'
	   PAUSE
	   STOP
       ENDIF
	 DO I = 2,6
	  IF ((matS2(I).NE.'-').AND.(matS2(I).NE.'a')
     #	  .AND.(matS2(I).NE.'H').AND.(matS2(I).NE.'G')
     #	  .AND.(matS2(I).NE.'c').AND.(matS2(I).NE.'F')
     #	  .AND.(matS2(I).NE.'R')) THEN
	   PAUSE
	   TYPE *, ' Check matS2 in namelist S2'
	   TYPE *, ' Program aborting'
	   STOP
	  ENDIF
	 ENDDO
	ENDIF
	CLOSE(1)

	IF (nv.GT.2) THEN
C FIND namelist S3
	 OPEN(1,File=TRIM(FILEIN), ACCESS='SEQUENTIAL')
	 IFAIL = -1
	 DO WHILE (.NOT.EOF(1))
	  READ(1,'(A)') STR
	  IF (INDEX(STR,'&S3').GT.0) IFAIL = 0
	 ENDDO
	 CLOSE(1)
	 IF (IFAIL.EQ.-1) THEN
	  TYPE *, ' Namelist S3 not found'
	  TYPE *, ' Program aborting'
	  PAUSE
	  RETURN
	 ENDIF

C READ namelist S3
	 OPEN(1,File=TRIM(FILEIN), ACCESS='SEQUENTIAL')
	 dynS3      = '-'
	 ns3   = -1
	 hypS3(:,:) = 1
	 matS3(:)   = '-'
	 READ(1,NML=S3,END=5004,ERR=5004)

	 IF ((dynS3.NE.'I').AND.(dynS3.NE.'M'))THEN
	  TYPE *, ' Check dynS3 in namelist S3'
	  TYPE *, ' Program aborting'
	  PAUSE
	  STOP
	 ENDIF

	 IF (ns3.LT.2) THEN
	  TYPE *, ' Check ns3 in namelist S3'
	  TYPE *, ' Program aborting'
	  PAUSE
	  STOP
	 ENDIF

	 IF (dynS3.EQ.'I') THEN
        DO J = 1,ns3
	   IF (hypS3(J,1).LE.0) THEN
	    TYPE *, ' Check hypS3 in namelist S3'
	    TYPE *, ' Program aborting'
		PAUSE
		STOP
	   ENDIF
        ENDDO
	 ELSEIF (dynS3.EQ.'M') THEN
	  DO J = 1,ns3
	  DO K = 1,ns3
	   IF (hypS3(J,K).LE.0) THEN
	    TYPE *, ' Check hypS3 in namelist S3'
	    TYPE *, ' Program aborting'
		PAUSE
		STOP
	   ENDIF
        ENDDO
	  ENDDO
	 ENDIF

	 IF ((matS3(1).EQ.'-').OR.((matS3(1).NE.'a')
     #	  .AND.(matS3(1).NE.'H').AND.(matS3(1).NE.'G')
     #	  .AND.(matS3(1).NE.'c').AND.(matS3(1).NE.'F')
     #	  .AND.(matS3(1).NE.'R'))) THEN
	   TYPE *, ' Check matS3 in namelist S3'
	   TYPE *, ' Program aborting'
	   PAUSE
	   STOP
       ENDIF
	 DO I = 2,6
	  IF ((matS3(I).NE.'-').AND.(matS3(I).NE.'a')
     #	  .AND.(matS3(I).NE.'H').AND.(matS3(I).NE.'G')
     #	  .AND.(matS3(I).NE.'c').AND.(matS3(I).NE.'F')
     #	  .AND.(matS3(I).NE.'R')) THEN
	   TYPE *, ' Check matS3 in namelist S3'
	   TYPE *, ' Program aborting'
	   PAUSE
	   STOP
	  ENDIF
	 ENDDO
	ENDIF
	CLOSE(1)

	IF (nv.GT.3) THEN
C FIND namelist S4
	 OPEN(1,File=TRIM(FILEIN), ACCESS='SEQUENTIAL')
	 IFAIL = -1
	 DO WHILE (.NOT.EOF(1))
	  READ(1,'(A)') STR
	  IF (INDEX(STR,'&S4').GT.0) IFAIL = 0
	 ENDDO
	 CLOSE(1)
	 IF (IFAIL.EQ.-1) THEN
	  TYPE *, ' Namlist S4 not found'
	  TYPE *, ' Program aborting'
	  PAUSE
	  RETURN
	 ENDIF

C READ namelist S4
	 OPEN(1,File=TRIM(FILEIN), ACCESS='SEQUENTIAL')
	 dynS4      = '-'
	 ns4   = -1
	 hypS4(:,:) = 1
	 matS4(:)   = '-'
	 READ(1,NML=S4,END=5005,ERR=5005)

	 IF ((dynS4.NE.'I').AND.(dynS4.NE.'M'))THEN
	  TYPE *, ' Check dynS4 in namelist S4'
	  TYPE *, ' Program aborting'
	  PAUSE
	  STOP
	 ENDIF

	 IF (ns4.LT.2) THEN
	  TYPE *, ' Check ns4 in namelist S4'
	  TYPE *, ' Program aborting'
	  PAUSE
	  STOP
	 ENDIF

	 IF (dynS4.EQ.'I') THEN
        DO J = 1,ns4
	   IF (hypS4(J,1).LE.0) THEN
	    TYPE *, ' Check hypS4 in namelist S4'
	    TYPE *, ' Program aborting'
		PAUSE
		STOP
	   ENDIF
        ENDDO
	 ELSEIF (dynS4.EQ.'M') THEN
	  DO J = 1,ns4
	  DO K = 1,ns4
	   IF (hypS4(J,K).LE.0) THEN
	    TYPE *, ' Check hypS4 in namelist S4'
	    TYPE *, ' Program aborting'
		PAUSE
		STOP
	   ENDIF
        ENDDO
	  ENDDO
	 ENDIF

	 IF ((matS4(1).EQ.'-').OR.((matS4(1).NE.'a')
     #	  .AND.(matS4(1).NE.'H').AND.(matS4(1).NE.'G')
     #	  .AND.(matS4(1).NE.'c').AND.(matS4(1).NE.'F')
     #	  .AND.(matS4(1).NE.'R'))) THEN
	   TYPE *, ' Check matS4 in namelist S4'
	   TYPE *, ' Program aborting'
	   PAUSE
	   STOP
       ENDIF
	 DO I = 2,6
	  IF ((matS4(I).NE.'-').AND.(matS4(I).NE.'a')
     #	  .AND.(matS4(I).NE.'H').AND.(matS4(I).NE.'G')
     #	  .AND.(matS4(I).NE.'c').AND.(matS4(I).NE.'F')
     #	  .AND.(matS4(I).NE.'R')) THEN
	   TYPE *, ' Check matS4 in namelist S4'
	   TYPE *, ' Program aborting'
	   PAUSE
	   STOP
	  ENDIF
	 ENDDO
	ENDIF
	CLOSE(1)

	IF (nv.GT.4) THEN
C FIND namelist S5
	 OPEN(1,File=TRIM(FILEIN), ACCESS='SEQUENTIAL')
	 IFAIL = -1
	 DO WHILE (.NOT.EOF(1))
	  READ(1,'(A)') STR
	  IF (INDEX(STR,'&S5').GT.0) IFAIL = 0
	 ENDDO
	 CLOSE(1)
	 IF (IFAIL.EQ.-1) THEN
	  TYPE *, ' Namlist S5 not found'
	  TYPE *, ' Program aborting'
	  PAUSE
	  RETURN
	 ENDIF

C READ namelist S5
	 OPEN(1,File=TRIM(FILEIN), ACCESS='SEQUENTIAL')
	 dynS5      = '-'
	 ns5   = -1
	 hypS5(:,:) = 1
	 matS5(:)   = '-'
	 READ(1,NML=S5,END=5006,ERR=5006)

	 IF ((dynS5.NE.'I').AND.(dynS5.NE.'M'))THEN
	  TYPE *, ' Check dynS5 in namelist S5'
	  TYPE *, ' Program aborting'
	  PAUSE
	  STOP
	 ENDIF

	 IF (ns5.LT.2) THEN
	  TYPE *, ' Check ns5 in namelist S5'
	  TYPE *, ' Program aborting'
	  PAUSE
	  STOP
	 ENDIF

	 IF (dynS5.EQ.'I') THEN
        DO J = 1,ns5
	   IF (hypS5(J,1).LE.0) THEN
	    TYPE *, ' Check hypS5 in namelist S5'
	    TYPE *, ' Program aborting'
		PAUSE
		STOP
	   ENDIF
        ENDDO
	 ELSEIF (dynS5.EQ.'M') THEN
	  DO J = 1,ns5
	  DO K = 1,ns5
	   IF (hypS5(J,K).LE.0) THEN
	    TYPE *, ' Check hypS5 in namelist S5'
	    TYPE *, ' Program aborting'
		PAUSE
		STOP
	   ENDIF
        ENDDO
	  ENDDO
	 ENDIF

	 IF ((matS5(1).EQ.'-').OR.((matS5(1).NE.'a')
     #	  .AND.(matS5(1).NE.'H').AND.(matS5(1).NE.'G')
     #	  .AND.(matS5(1).NE.'c').AND.(matS5(1).NE.'F')
     #	  .AND.(matS5(1).NE.'R'))) THEN
	   TYPE *, ' Check matS5 in namelist S5'
	   TYPE *, ' Program aborting'
	   PAUSE
	   STOP
       ENDIF
	 DO I = 2,6
	  IF ((matS5(I).NE.'-').AND.(matS5(I).NE.'a')
     #	  .AND.(matS5(I).NE.'H').AND.(matS5(I).NE.'G')
     #	  .AND.(matS5(I).NE.'c').AND.(matS5(I).NE.'F')
     #	  .AND.(matS5(I).NE.'R')) THEN
	   TYPE *, ' Check matS5 in namelist S5'
	   TYPE *, ' Program aborting'
	   PAUSE
	   STOP
	  ENDIF
	 ENDDO
	ENDIF
	CLOSE(1)

	IF (nv.GT.5) THEN
C FIND namelist S6
	 OPEN(1,File=TRIM(FILEIN), ACCESS='SEQUENTIAL')
	 IFAIL = -1
	 DO WHILE (.NOT.EOF(1))
	  READ(1,'(A)') STR
	  IF (INDEX(STR,'&S6').GT.0) IFAIL = 0
	 ENDDO
	 CLOSE(1)
	 IF (IFAIL.EQ.-1) THEN
	  TYPE *, ' Namlist S6 not found'
	  TYPE *, ' Program aborting'
	  PAUSE
	  RETURN
	 ENDIF

C READ namelist S6
	 OPEN(1,File=TRIM(FILEIN), ACCESS='SEQUENTIAL')
	 dynS6      = '-'
	 ns6   = -1
	 hypS6(:,:) = 1
	 matS6(:)   = '-'
	 READ(1,NML=S6,END=5007,ERR=5007)

	 IF ((dynS6.NE.'I').AND.(dynS6.NE.'M'))THEN
	  TYPE *, ' Check dynS6 in namelist S6'
	  TYPE *, ' Program aborting'
	  PAUSE
	  STOP
	 ENDIF

	 IF (ns6.LT.2) THEN
	  TYPE *, ' Check ns6 in namelist S6'
	  TYPE *, ' Program aborting'
	  PAUSE
	  STOP
	 ENDIF

	 IF (dynS6.EQ.'I') THEN
        DO J = 1,ns6
	   IF (hypS6(J,1).LE.0) THEN
	    TYPE *, ' Check hypS6 in namelist S6'
	    TYPE *, ' Program aborting'
		PAUSE
		STOP
	   ENDIF
        ENDDO
	 ELSEIF (dynS6.EQ.'M') THEN
	  DO J = 1,ns6
	  DO K = 1,ns6
	   IF (hypS6(J,K).LE.0) THEN
	    TYPE *, ' Check hypS6 in namelist S6'
	    TYPE *, ' Program aborting'
		PAUSE
		STOP
	   ENDIF
        ENDDO
	  ENDDO
	 ENDIF

	 IF ((matS6(1).EQ.'-').OR.((matS6(1).NE.'a')
     #	  .AND.(matS6(1).NE.'H').AND.(matS6(1).NE.'G')
     #	  .AND.(matS6(1).NE.'c').AND.(matS6(1).NE.'F')
     #	  .AND.(matS6(1).NE.'R'))) THEN
	   TYPE *, ' Check matS6 in namelist S6'
	   TYPE *, ' Program aborting'
	   PAUSE
	   STOP
       ENDIF
	 DO I = 2,6
	  IF ((matS6(I).NE.'-').AND.(matS6(I).NE.'a')
     #	  .AND.(matS6(I).NE.'H').AND.(matS6(I).NE.'G')
     #	  .AND.(matS6(I).NE.'c').AND.(matS6(I).NE.'F')
     #	  .AND.(matS6(I).NE.'R')) THEN
	   TYPE *, ' Check matS6 in namelist S6'
	   TYPE *, ' Program aborting'
	   PAUSE
	   STOP
	  ENDIF
	 ENDDO
	ENDIF
	CLOSE(1)

C FIND namelist prior
	OPEN(1,File=TRIM(FILEIN), ACCESS='SEQUENTIAL')
	IFAIL = -1
	DO WHILE (.NOT.EOF(1))
	 READ(1,'(A)') STR
	 IF (INDEX(STR,'&prior').GT.0) IFAIL = 0
	ENDDO
	CLOSE(1)
	IF (IFAIL.EQ.-1) THEN
	 TYPE *, ' Namelist prior not found'
	 TYPE *, ' Program aborting'
	 PAUSE
	 RETURN
	ENDIF
C READ namelist prior
	OPEN(1,File=TRIM(FILEIN), ACCESS='SEQUENTIAL')
	nt            = 0
	pdftheta(:)   = '  '
	hyptheta(:,:) = -1
	READ(1,NML = prior,END=5008,ERR=5008)
	IF (nt.LE.0) THEN
	 TYPE *, ' Check nt in namelist prior'
	 TYPE *, ' Program aborting'
	 PAUSE
	 STOP
      ENDIF
	IF (nt.GT.200) THEN
	 TYPE *, ' nt is too large '
	 TYPE *, ' Program aborting'
	 PAUSE
	 STOP
	ENDIF
	IF (estimation.EQ.'BA') THEN
       DO I = 1,nt
	  WRITE(IC,'(I3)') I
	  IF ((pdftheta(I).NE.'BE').AND.(pdftheta(I).NE.'NT').AND.
     #     (pdftheta(I).NE.'IG')) THEN
	   TYPE *, ' Check pdftheta('//IC//') in namelist prior'
	   TYPE *, ' Program aborting'
	   PAUSE
	   STOP
	  ENDIF
	  IF (hyptheta(3,I).GT.hyptheta(4,I)) THEN
	   TYPE *, ' Check hyptheta('//IC//') in namelist prior'
	   TYPE *, ' Program aborting'
	   PAUSE
	   STOP
	  ENDIF
	  IF (pdftheta(I).EQ.'BE') THEN
	   IF (hyptheta(3,I).LT.hyptheta(4,I)) THEN
	    IF ((hyptheta(1,I).LE.0.).OR.(hyptheta(2,I).LE.0.).OR.
     #       (hyptheta(3,I).GT.hyptheta(4,I))) THEN
	     TYPE *, ' Check hyptheta('//IC//') in namelist prior'
	     TYPE *, ' Program aborting'
	     PAUSE
	     STOP
	    ENDIF
	   ENDIF
	  ELSEIF (pdftheta(I).EQ.'NT') THEN
	   IF (hyptheta(3,I).LT.hyptheta(4,I)) THEN
	    IF ((hyptheta(2,I).LE.0.).OR.(hyptheta(3,I).GT.hyptheta(4,I)))
     #     THEN
	     TYPE *, ' Check hyptheta('//IC//') in namelist prior'
	     TYPE *, ' Program aborting'
	     PAUSE
	     STOP
	    ENDIF
	   ENDIF
	  ELSEIF (pdftheta(I).EQ.'IG') THEN
 	   IF (hyptheta(3,I).LT.hyptheta(4,I)) THEN
	    IF ((hyptheta(1,I).LE.0.).OR.(hyptheta(2,I).LE.0.).OR.
     #    (hyptheta(3,I).GT.hyptheta(4,I)).OR.(hyptheta(3,I).LT.0.))THEN
	     TYPE *, ' Check hyptheta('//IC//') in namelist prior'
	     TYPE *, ' Program aborting'
	     PAUSE
	     STOP
	    ENDIF
	   ENDIF
	  ENDIF
       ENDDO
      ELSE
       DO I = 1,nt  ! ML check
        WRITE(IC,'(I3)') I
        IF (hyptheta(3,I).GT.hyptheta(4,I)) THEN
	   TYPE *, ' Check hyptheta('//IC//') in namelist prior'
	   TYPE *, ' Program aborting'
	   PAUSE
	   STOP
	  ENDIF
       ENDDO
      ENDIF
	CLOSE(1)

C FIND namelist mcmc
	OPEN(1,File=TRIM(FILEIN), ACCESS='SEQUENTIAL')
	IFAIL = -1
	DO WHILE (.NOT.EOF(1))
	 READ(1,'(A)') STR
	 IF (INDEX(STR,'&mcmc').GT.0) IFAIL = 0
	ENDDO
	CLOSE(1)
	IF (IFAIL.EQ.-1) THEN
	 TYPE *, ' Namelist mcmc not found'
	 TYPE *, ' Program aborting'
	 PAUSE
	 RETURN
	ENDIF

C READ namelist mcmc
	OPEN(1,File=TRIM(FILEIN), ACCESS='SEQUENTIAL')
	seed         = 0
	thin         = 1
	burnin       = 1000
	simulrec     = 5000
	sampler      = 'SL'
	Ssampler     = 'GCK'
	hbl          = 1
	MargLik      = 'N'
	READ(1,NML=mcmc,END=5009,ERR=5009)
	IF ((seed.LT.0).OR.(seed.GT.999)) THEN
	 TYPE *, ' Check seed in namelist mcmc'
	 TYPE *, ' Program aborting'
	 PAUSE
	 STOP
	ENDIF
	IF (thin.LT.1) THEN
	 TYPE *, ' Check thin in namelist mcmc'
	 TYPE *, ' Program aborting'
	 PAUSE
	 STOP
	ENDIF
	IF (burnin.LE.0) THEN
	 TYPE *, ' Check burnin in namelist mcmc'
	 TYPE *, ' Program aborting'
	 PAUSE
	 STOP
	ENDIF
	IF (simulrec.LE.1) THEN
	 TYPE *, ' Check simulrec in namelist mcmc'
	 TYPE *, ' Program aborting'
	 PAUSE
	 STOP
	ENDIF
	IF ((sampler.NE.'SL').AND.(sampler.NE.'MH')) THEN
	 TYPE *, ' Check sampler in namelist mcmc'
	 TYPE *, ' Program aborting'
	 PAUSE
	 STOP
	ENDIF

c	IF ((Ssampler.NE.'GCK').AND.(Ssampler.NE.'AMH')
c     #	.AND.(Ssampler.NE.'MH ')) THEN
c	 TYPE *, ' Check Ssampler in namelist mcmc'
c	 TYPE *, ' Program aborting'
c	 PAUSE
c	 STOP
c	ENDIF
c	IF ((hbl.GT.1).AND.(Ssampler.EQ.'GCK')) THEN
c	 TYPE *, ' Check hbl in namelist mcmc'
c	 TYPE *, ' Program aborting'
c	 PAUSE
c	 STOP
c	ENDIF

	CLOSE(1)

C FIND namelist dataset
	OPEN(1,File=TRIM(FILEIN), ACCESS='SEQUENTIAL')
	IFAIL = -1
	DO WHILE (.NOT.EOF(1))
	 READ(1,'(A)') STR
	 IF (INDEX(STR,'&dataset').GT.0) IFAIL = 0
	ENDDO
	CLOSE(1)
	IF (IFAIL.EQ.-1) THEN
	 TYPE *, ' Namelist dataset not found'
	 TYPE *, ' Program aborting'
	 PAUSE
	 RETURN
	ENDIF
C READ namelist dataset
	OPEN(1,File=TRIM(FILEIN), ACCESS='SEQUENTIAL')
	ny      = -1
	nz      = -1
	nf      =  0
	datasim = 'N'
	READ(1,NML=dataset,END=5010,ERR=5010)
	IF ((T.LE.0).OR.(T.GT.3000)) THEN
	 TYPE *, ' Check T in namelist dataset (T<=3000)'
	 TYPE *, ' Program aborting'
	 PAUSE
	 STOP
	ENDIF
	IF (ny.LE.0) THEN
	 TYPE *, ' Check ny in namelist dataset'
	 TYPE *, ' Program aborting'
	 PAUSE
	 STOP
	ENDIF
	IF (nz.LT.0) THEN
	 TYPE *, ' Check nz in namelist dataset'
	 TYPE *, ' Program aborting'
	 PAUSE
	 STOP
	ENDIF
	IF (nf.LT.0) THEN
	 TYPE *, ' Check nf in namelist dataset'
	 TYPE *, ' Program aborting'
	 PAUSE
	 STOP
	ENDIF
	IF (T.LT.hbl) THEN
	 TYPE *, ' Check hbl in namelist mcmc (hbl > T)'
	 TYPE *, ' Program aborting'
	 PAUSE
	 STOP
	ENDIF
	IF ((datasim.NE.'N').AND.(datasim.NE.'n').AND.
     &    (datasim.NE.'y').AND.(datasim.NE.'Y')) THEN
	  datasim = 'N'
	ENDIF

	CLOSE(1)

C -----------------------------------------------------------------------
C ASSIGN discrete latent variables: INFOS (9 x nv)
C by cols: S1,S2,...,SNV; with nv <=6
C by row: the 1st contains the # of matrices affected by Si
C         the 2nd-3rd etc point to c (1),H (2),G (3),a (4),F (5),R (6)
C         the 8-th  row contains the # of states
C         the 9-th  row spec. the dynamic of Sj (0-deterministic,1=Indep,2=Markov)
C nstot:  total # of states i.e. ns1 x ns2 x ...x nsv
C np(1):  total # of psi parameters
C np(2):  total # of idependent PSI-Dirichlet vectors
C np(3):  max # of hyperparameters for psi
C -----------------------------------------------------------------------
	INFOS(:,:) = 0
	INFOS(8,:) = 1  ! number of states
	ns(:)      = 1
	np(:)      = 0
	IF (nv.GT.0) THEN
	 ALLOCATE(matS(6,6),dynS(6),nstateS(6))
	 matS(:,1) = matS1(:)
	 matS(:,2) = matS2(:)
	 matS(:,3) = matS3(:)
	 matS(:,4) = matS4(:)
	 matS(:,5) = matS5(:)
	 matS(:,6) = matS6(:)
	 hypS(:,:,1) = hypS1(:,:)
	 hypS(:,:,2) = hypS2(:,:)
	 hypS(:,:,3) = hypS3(:,:)
	 hypS(:,:,4) = hypS4(:,:)
	 hypS(:,:,5) = hypS5(:,:)
	 hypS(:,:,6) = hypS6(:,:)
	 dynS(1) = dynS1
	 dynS(2) = dynS2
	 dynS(3) = dynS3
	 dynS(4) = dynS4
	 dynS(5) = dynS5
	 dynS(6) = dynS6
	 nstateS(1) = ns1
       nstateS(2) = ns2
	 nstateS(3) = ns3
	 nstateS(4) = ns4
	 nstateS(5) = ns5
	 nstateS(6) = ns6
	 IMAX  = MAXLOC(nstateS(1:nv))
	 np(2) = 0
	 np(3) = nstateS(IMAX(1))
	 DO 50 J = 1,nv
	  K = 0
	  DO 40 I = 1,6
	   IF(matS(I,J).EQ.'c') THEN
	     K = K + 1
	     IND = 1
  	     INFOS(K+1,J) = IND          ! Matrices affected by SJ
	   ENDIF
	   IF(matS(I,J).EQ.'H') THEN
	     K = K + 1
	     IND = 2
	     INFOS(K+1,J) = IND
	   ENDIF
	   IF(matS(I,J).EQ.'G') THEN
	     K = K + 1
	     IND = 3
	     INFOS(K+1,J) = IND
	   ENDIF
	   IF(matS(I,J).EQ.'a') THEN
	     K = K + 1
	     IND = 4
	     INFOS(K+1,J) = IND
	   ENDIF
	   IF(matS(I,J).EQ.'F') THEN
	     K = K + 1
	     IND = 5
	     INFOS(K+1,J) = IND
	   ENDIF
	   IF(matS(I,J).EQ.'R') THEN
	     K = K + 1
	     IND = 6
	     INFOS(K+1,J) = IND
	   ENDIF
40      CONTINUE
	  INFOS(1,J) = K                 ! # of matrix affected by SJ
	  INFOS(8,J) = nstateS(J)        ! # of states for SJ
	  IF (dynS(J).EQ.'I') THEN
	   INFOS(9,J) = 1                ! dynamics for Sj
	   np(2) = np(2) + 1
	   np(1) = np(1) + nstateS(J)-1
	  ELSEIF (dynS(J).EQ.'M') THEN
	   INFOS(9,J) = 2
	   np(2) = np(2) + nstateS(J)
	   np(1) = np(1) + (nstateS(J)-1)*nstateS(J)
	  ENDIF
50     CONTINUE

	 DO 60 I = 1,nv
	 DO 60 J = 1,INFOS(1,I)
60	 ns(INFOS(J+1,I)) = INFOS(8,I)
	 nstot = PRODUCT(INFOS(8,1:nv))  ! total # of states

	 DEALLOCATE(matS,dynS,nstateS)
	ENDIF

	GO TO 7777

5000  TYPE *,'Input file not found'
	TYPE *,'Program aborting'
	PAUSE
	STOP

5001  TYPE *,'Input error in namelist ssm'
	TYPE *,'Program aborting'
	PAUSE
	STOP
5002  TYPE *,'Input error in namelist S1 '
	TYPE *,'Program aborting'
	PAUSE
	STOP
5003  TYPE *,'Input error in namelist S2 '
	TYPE *,'Program aborting'
	PAUSE
	STOP
5004  TYPE *,'Input error in namelist S3 '
	TYPE *,'Program aborting'
	PAUSE
	STOP
5005  TYPE *,'Input error in namelist S4 '
	TYPE *,'Program aborting'
	PAUSE
	STOP
5006  TYPE *,'Input error in namelist S5 '
	TYPE *,'Program aborting'
	PAUSE
	STOP
5007  TYPE *,'Input error in namelist S6 '
	TYPE *,'Program aborting'
	PAUSE
	STOP
5008  TYPE *,'Input error in namelist prior'
	TYPE *,'Program aborting'
	PAUSE
	STOP
5009  TYPE *,'Input error in namelist mcmc'
	TYPE *,'Program aborting'
	PAUSE
	STOP
5010  TYPE *,'Input error in namelist dataset'
	TYPE *,'Program aborting'
	PAUSE
	STOP

7777	RETURN
      END
