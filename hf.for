C ----------------------------------------------------------------------	
C HF computes the loglikelihood, the innovations, and the Hamilton (1989)
C Smoother for a Markov Switching VAR(1)
C Developed by A.Rossi, C.Planas and G.Fiorentini     
C      
C OUTPUT: LLILIKE:  log-likelihood
C         INN:      innovations      
C         SSMOOTH:  P(s(t)|y^T), t=1,2,...,T      
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
C ----------------------------------------------------------------------	
	SUBROUTINE HF(nobs,nx,nk,nz,nu,ns,nv,np,psi,ismoo,yk,IYK,INFOS,
     1              c,a,F,R,SSMOOTH,INN,LLIKE)

C INPUT
	INTEGER nobs,nx,nk,nz,nu,nv,np,ns(6),ismoo,IYK(nobs,nx+1),INFOS(9,6)
	DOUBLE PRECISION yk(nobs,nx+nz),c(nx,max(1,nz),ns(1)),
     1 a(nx,ns(4)),F(nx,nx,ns(5)),R(nx,nu,ns(6)),psi(max(1,np))  
      
C OUTPUT
	DOUBLE PRECISION SSMOOTH(nobs,nk),INN(nobs,nx),LLIKE(nobs)

C LOCALS
	INTEGER inobs,inx,I,J,K,IMAX(1),IS(6),SEQ(1)
	DOUBLE PRECISION,ALLOCATABLE:: ZZ(:),gam(:),mu(:),SIG(:,:),FILT(:,:),
     1 P(:,:),PE(:),PP1(:,:),PP2(:,:),PP3(:,:),PP4(:,:),PP5(:,:),
     1 PP6(:,:)
         
      DOUBLE PRECISION logmvnpdf,mvnpdf,norm,Lmax
      ALLOCATE(PP1(INFOS(8,1),INFOS(8,1)),PP2(INFOS(8,2),INFOS(8,2)),
     1 PP3(INFOS(8,3),INFOS(8,3)),PP4(INFOS(8,4),INFOS(8,4)),
     1 PP5(INFOS(8,5),INFOS(8,5)),PP6(INFOS(8,6),INFOS(8,6)))
	ALLOCATE(ZZ(nk),gam(nk),mu(nx),SIG(nx,nx),FILT(nobs,nk),
     1 P(nk,nk),PE(nk))
      
      FILT(:,:) = 0.D0
	LLIKE(:)  = 0.D0      
      CALL designz(nv,np,psi,INFOS,PP1,PP2,PP3,PP4,PP5,PP6)  
      CALL pprod(nv,nk,INFOS,PP1,PP2,PP3,PP4,PP5,PP6,P)          
C Initial condition (stationary MS-VAR(1))
      CALL ergodic(nk,P,PE)
      inx = IYK(1,nx+1)
      DO J = 1,inx
	 mu(J) = a(IYK(1,J),IS(4))+
     #       + SUM(F(IYK(1,J),IYK(1,1:inx),IS(5))
     #       * yk(1,IYK(1,1:inx)))
     #       + SUM(c(IYK(1,J),1:nz,IS(1))*yk(1,nx+1:nx+nz))
      ENDDO
C SIG = R*R'
      DO  K=1,inx
       SIG(K,K) = SUM(R(IYK(1,K),1:nu,IS(6))
     #          *     R(IYK(1,K),1:nu,IS(6)))
	 DO  J=1,K-1	
        SIG(K,J) = SUM(R(IYK(1,K),1:nu,IS(6))
     #           *     R(IYK(1,J),1:nu,IS(6)))
        SIG(J,K) = SIG(K,J)         
       ENDDO
      ENDDO
      DO I = 1,nk
       gam(I) = PE(I)*mvnpdf(yk(1,IYK(1,1:inx)),mu(1:inx),
     #                       SIG(1:inx,1:inx),inx)
      ENDDO    
      LLIKE(1) = dlog(sum(gam(1:nk)))
      FILT(1,1:nk-1) = gam(1:nk-1)/sum(gam(1:nk))
      FILT(1,nk) = 1.D0-sum(FILT(1,1:nk-1)) 

C ----------------------------------------------------------
C Filtering Z(t|t)=P*Z(t-1|t-1) (.*) F(j) / SUM of numerator
C ----------------------------------------------------------
	DO 100 inobs = 2,nobs
       inx = IYK(inobs,nx+1)
	 DO I =1,nk
        CALL int2seq(I,nv,INFOS,SEQ,IS)     
	  ZZ(I) = SUM(P(I,1:nk)*FILT(inobs-1,1:nk)) 	
c y(t) = c(t)z(t) + x(t)
C x(t) = a(t)     + F(t)x(t-1) + R(t)u(t)
        DO J = 1,inx
	   mu(J) = a(IYK(inobs,J),IS(4))+
     #         + SUM(F(IYK(inobs,J),IYK(inobs,1:inx),IS(5))
     #         * yk(inobs-1,IYK(inobs-1,1:inx)))
     #         + SUM(c(IYK(inobs,J),1:nz,IS(1))*yk(inobs,nx+1:nx+nz))
        ENDDO
C SIG = R*R'
        DO  K=1,inx
         SIG(K,K) = SUM(R(IYK(inobs,K),1:nu,IS(6))
     #            *     R(IYK(inobs,K),1:nu,IS(6)))
	   DO  J=1,K-1	
          SIG(K,J) = SUM(R(IYK(inobs,K),1:nu,IS(6))
     #             *     R(IYK(inobs,J),1:nu,IS(6)))
          SIG(J,K) = SIG(K,J)         
         ENDDO
        ENDDO
        gam(I) = logmvnpdf(yk(inobs,IYK(inobs,1:inx)),mu(1:inx),
     #                     SIG(1:inx,1:inx),inx) ! log
       ENDDO 
C ---------------------------------------------------
C Compute the log-likelihood and Normalise the filter
C ---------------------------------------------------
	 IMAX = MAXLOC(gam(1:nk))
	 Lmax = gam(IMAX(1))
       gam(:) = dexp(gam(:)-Lmax)
       FILT(inobs,1:nk) = ZZ(1:nk)*gam(1:nk)
	 norm = SUM(FILT(inobs,:))
       FILT(inobs,:) = FILT(inobs,:)/norm   
	 LLIKE(inobs)  = DLOG(norm) + Lmax   
100   CONTINUE      
      
C ------------------------------------------------------------------------------------
C Hamilton smoother for a MS-VAR(1)
C Z(t|T): Hamilton (94), pp 694      
C Innovations: 
C INN(t) = y(y)-SUM(S(t),S(t-1)) P(S(t)|S(t-1)*P(S(t-1)|x^(t-1))*E(x(t)|S(t),x^(t-1))       
C ------------------------------------------------------------------------------------      
      IF (ismoo.EQ.1) THEN	 
       SSMOOTH(nobs,:) = FILT(nobs,:)
	 DO inobs = nobs-1,1,-1
	  DO J=1,nk
	   gam(J) = SUM(P(J,1:nk)*FILT(inobs,1:nk))  
        ENDDO	 
	  DO J=1,nk
	   ZZ(J) = 1.D0
	   IF (gam(J).GT.1.D-13) ZZ(J) = SSMOOTH(inobs+1,J)/gam(J) 
        ENDDO
	  DO J=1,nk
	   gam(J) = SUM(P(1:nk,J)*ZZ(1:nk))  
        ENDDO
	  SSMOOTH(inobs,1:nk) = FILT(inobs,1:nk)*gam(1:nk)
       ENDDO
       INN(:,:) = 0.D0
       DO inobs = 2,nobs
        inx = IYK(inobs,nx+1)   
	  DO I=1,nk
         CALL int2seq(I,nv,INFOS,SEQ,IS)                 
         DO K = 1,inx
	    mu(K) = a(IYK(inobs,K),IS(4))+
     #          + SUM(F(IYK(inobs,K),IYK(inobs,1:inx),IS(5))
     #          * yk(inobs-1,IYK(inobs-1,1:inx)))
     #          + SUM(c(IYK(inobs,K),1:nz,IS(1))*yk(inobs,nx+1:nx+nz))
         ENDDO          
         DO J=1,nk
          INN(inobs,IYK(inobs,1:inx)) = INN(inobs,IYK(inobs,1:inx)) 
     1                                + mu(1:inx)*P(I,J)*FILT(inobs-1,J)
         ENDDO 
        ENDDO
        INN(inobs,IYK(inobs,1:inx)) = yk(inobs,IYK(inobs,1:inx)) 
     1                              - INN(inobs,IYK(inobs,1:inx))  
       ENDDO
      ENDIF    
	
      DEALLOCATE(gam,mu,SIG)
	RETURN
	END
