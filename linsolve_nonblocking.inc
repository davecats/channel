#define req(nm,iz,j) nm((2*nz+1)*(j)+iz+nz+1)
  !--------------------------------------------------------------!
  !------------------- solve the linear system  -----------------!
  SUBROUTINE linsolve(lambda)
    real(C_DOUBLE), intent(in) :: lambda
    integer(C_INT) :: ix,iz,i
    DO ix=nx0,nxN
      DO iz=-nz,nz
        DO CONCURRENT (iy=ny0:nyN)
          D2vmat(iy,-2:2,iz)=lambda*(der(iy)%d2(-2:2)-k2(iz,ix)*der(iy)%d0(-2:2))-OS(iy,-2:2)
          etamat(iy,-2:2,iz)=lambda*der(iy)%d0(-2:2)-SQ(iy,-2:2)
        END DO
        IF (first) THEN 
          IF (ix==0 .AND. iz==0) THEN
            bc0(iz,ix)%v=0; bc0(iz,ix)%vy=0; bc0(iz,ix)%eta=dcmplx(dreal(bc0(iz,ix)%u)-dimag(bc0(iz,ix)%w),dimag(bc0(iz,ix)%u)+dreal(bc0(iz,ix)%w))
          ELSE
            bc0(iz,ix)%vy=-ialfa(ix)*bc0(iz,ix)%u-ibeta(iz)*bc0(iz,ix)%w; bc0(iz,ix)%eta=ibeta(iz)*bc0(iz,ix)%u-ialfa(ix)*bc0(iz,ix)%w
          END IF
          bc0(iz,ix)%v=bc0(iz,ix)%v-v0bc(-2)*bc0(iz,ix)%vy/v0m1bc(-2)
          CALL applybc_0(D2vmat(:,:,iz),v0bc,v0m1bc)
          V(1,iz,ix,2)=V(1,iz,ix,2)-D2vmat(1,-2,iz)*bc0(iz,ix)%vy/v0m1bc(-2)-D2vmat(1,-1,iz)*bc0(iz,ix)%v/v0bc(-1)
          V(2,iz,ix,2)=V(2,iz,ix,2)-D2vmat(2,-2,iz)*bc0(iz,ix)%v/v0bc(-1)       
          CALL applybc_0(etamat(:,:,iz),eta0bc,eta0m1bc)
          V(1,iz,ix,1)=V(1,iz,ix,1)-etamat(1,-1,iz)*bc0(iz,ix)%eta/eta0bc(-1)
          V(2,iz,ix,1)=V(2,iz,ix,1)-etamat(2,-2,iz)*bc0(iz,ix)%eta/eta0bc(-1)
        END IF
        IF (last) THEN 
          IF (ix==0 .AND. iz==0) THEN
            bcn(iz,ix)%v=0; bcn(iz,ix)%vy=0; bcn(iz,ix)%eta=dcmplx(dreal(bcn(iz,ix)%u)-dimag(bcn(iz,ix)%w),dimag(bcn(iz,ix)%u)+dreal(bcn(iz,ix)%w))
          ELSE
            bcn(iz,ix)%vy=-ialfa(ix)*bcn(iz,ix)%u-ibeta(iz)*bcn(iz,ix)%w; bcn(iz,ix)%eta=ibeta(iz)*bcn(iz,ix)%u-ialfa(ix)*bcn(iz,ix)%w
          END IF
          bcn(iz,ix)%v=bcn(iz,ix)%v-vnbc(2)*bcn(iz,ix)%vy/vnp1bc(2)
          CALL applybc_n(D2vmat(:,:,iz),vnbc,vnp1bc)
          V(ny-1,iz,ix,2)=V(ny-1,iz,ix,2)-D2vmat(ny-1,2,iz)*bcn(iz,ix)%vy/vnp1bc(2)-D2vmat(ny-1,1,iz)*bcn(iz,ix)%v/vnbc(1)
          V(ny-2,iz,ix,2)=V(ny-2,iz,ix,2)-D2vmat(ny-2,2,iz)*bcn(iz,ix)%v/vnbc(1)
          CALL applybc_n(etamat(:,:,iz),etanbc,etanp1bc)
          V(ny-1,iz,ix,1)=V(ny-1,iz,ix,1)-etamat(ny-1,1,iz)*bcn(iz,ix)%eta/etanbc(1) 
          V(ny-2,iz,ix,1)=V(ny-2,iz,ix,1)-etamat(ny-2,2,iz)*bcn(iz,ix)%eta/etanbc(1)
        END IF
        CALL LU5decompStep(D2vmat(:,:,iz),req(ReqLinSolve,iz,0)); CALL LU5decompStep(etamat(:,:,iz),req(ReqLinSolve,iz,1))
        CALL LeftLU5divStep1(V(:,iz,ix,2),D2vmat(:,:,iz),V(:,iz,ix,2),req(ReqLinSolve,iz,2))
        CALL LeftLU5divStep1(V(:,iz,ix,1),etamat(:,:,iz),V(:,iz,ix,1),req(ReqLinSolve,iz,3))
      END DO    
      IF (.NOT. first) CALL MPI_Waitall((2*nz+1)*4,ReqLinSolve,MPI_STATUSES_IGNORE)
      IF (ix==0) eta00mat=etamat(:,:,0)
      DO iz=-nz,nz
        CALL LeftLU5divStep2(D2vmat(:,:,iz),V(:,iz,ix,2),req(ReqLinSolve,iz,0))
        CALL LeftLU5divStep2(etamat(:,:,iz),V(:,iz,ix,1),req(ReqLinSolve,iz,1))
        IF (first) THEN
          V(0,iz,ix,2)=(bc0(iz,ix)%v-sum(V(1:3,iz,ix,2)*v0bc(0:2)))/v0bc(-1)
          V(-1,iz,ix,2)=(bc0(iz,ix)%vy-sum(V(0:3,iz,ix,2)*v0m1bc(-1:2)))/v0m1bc(-2)
          V(0,iz,ix,1)=(bc0(iz,ix)%eta-sum(V(1:3,iz,ix,1)*eta0bc(0:2)))/eta0bc(-1)
          V(-1,iz,ix,1)=-sum(V(0:3,iz,ix,1)*eta0m1bc(-1:2))/eta0m1bc(-2)
        END IF
        IF (last) THEN 
          V(ny,iz,ix,2)=(bcn(iz,ix)%v-sum(V(ny-3:ny-1,iz,ix,2)*vnbc(-2:0)))/vnbc(1)
          V(ny+1,iz,ix,2)=(bcn(iz,ix)%vy-sum(V(ny-3:ny,iz,ix,2)*vnp1bc(-2:1)))/vnp1bc(2)
          V(ny,iz,ix,1)=(bcn(iz,ix)%eta-sum(V(ny-3:ny-1,iz,ix,1)*etanbc(-2:0)))/etanbc(1)
          V(ny+1,iz,ix,1)=-sum(V(ny-3:ny,iz,ix,1)*etanp1bc(-2:1))/etanp1bc(2)
        END IF
      END DO
      IF (.NOT. last) CALL MPI_Waitall((2*nz+1)*2,ReqLinSolve,MPI_STATUSES_IGNORE)
    END DO
  END SUBROUTINE linsolve

  !----------------------- computeflowrate  ---------------------!
  !--------------------------------------------------------------!
  SUBROUTINE computeflowrate(lambda)
    real(C_DOUBLE), intent(in) :: lambda
    complex(C_DOUBLE_COMPLEX) :: ucor(ny0-2:nyN+2)
    real(C_DOUBLE) :: frl(1:3)
    TYPE(MPI_REQUEST) :: Rs
    IF (has_average) THEN
         ucor(ny0-2:ny0-1)=0; ucor(ny0:nyN)=1; ucor(nyN+1:nyN+2)=0
         ! XXX TODO Make computation of ucor communication-free
         CALL LeftLU5divStep1(ucor,eta00mat(:,:),ucor,Rs)
         IF (.NOT. first) CALL MPI_WAIT(Rs,MPI_STATUS_IGNORE)
         CALL LeftLU5divStep2(eta00mat(:,:),ucor,Rs)
         IF (.NOT. last)  CALL MPI_WAIT(Rs,MPI_STATUS_IGNORE)
         IF (first) THEN
           ucor(0)=-sum(ucor(1:3)*eta0bc(0:2))/eta0bc(-1)
           ucor(-1)=-sum(ucor(0:3)*eta0m1bc(-1:2))/eta0m1bc(-2)
         END IF
         IF (last) THEN 
           ucor(ny)=-sum(ucor(ny-3:ny-1)*etanbc(-2:0))/etanbc(1)
           ucor(ny+1)=-sum(ucor(ny-3:ny)*etanp1bc(-2:1))/etanp1bc(2)  
         END IF
         frl(1)=yintegr(dreal(V(:,0,0,1))); frl(2)=yintegr(dreal(V(:,0,0,3))); frl(3)=yintegr(dreal(ucor))
         CALL MPI_Allreduce(frl,fr,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_Y)
         IF (abs(meanflowx)>1.0d-7 .AND. .NOT. CPI) THEN
           corrpx = (meanflowx-fr(1))/fr(3)
           V(:,0,0,1)=dcmplx(dreal(V(:,0,0,1))+corrpx*dreal(ucor),dimag(V(:,0,0,1)))
         END IF
         IF (abs(meanflowz)>1.0d-7 .AND. .NOT. CPI) THEN
           corrpz = (meanflowz-fr(2))/fr(3)
           V(:,0,0,3)=dcmplx(dreal(V(:,0,0,3))+corrpz*dreal(ucor),dimag(V(:,0,0,3)))
         END IF
         IF (CPI) THEN
           SELECT CASE (CPI_type)
              CASE (0)
                 meanpx = (1-gamma)*6*ni/fr(1)
              CASE (1)
                 meanpx = (1.5d0/gamma)*fr(1)*ni 
              CASE DEFAULT
                 WRITE(*,*) "Wrong selection of CPI_Type"
                 STOP
           END SELECT
         END IF 
    END IF
  END SUBROUTINE computeflowrate

  !------------------------- vetaTOuvw --------------------------!
  !--------------------------------------------------------------!
  SUBROUTINE vetaTOuvw()
    complex(C_DOUBLE_COMPLEX) :: temp(ny0-2:nyN+2)
    integer(C_INT) :: ix,iz,cnt_req
    IF (has_average) THEN
      V(:,0,0,3) = dcmplx(dimag(V(:,0,0,1)),0.d0); 
      V(:,0,0,1) = dcmplx(dreal(V(:,0,0,1)),0.d0);
    END IF
    cnt_req=1
    DO ix=nx0,nxN
      DO iz=-nz,nz
          IF ( .NOT. (ix==0 .AND. iz==0) )  THEN
            CALL COMPLEXderiv(V(:,iz,ix,2),V(:,iz,ix,3),REQvetaTOuvw(cnt_req)); cnt_req=cnt_req+1
          END IF
      END DO
    END DO
    IF (.NOT. first) CALL MPI_Waitall(cnt_req-1,REQvetaTOuvw,MPI_STATUSES_IGNORE); 
    cnt_req=1
    DO ix=nx0,nxN
      DO iz=-nz,nz
          IF (.NOT. ( ix==0 .AND. iz==0) ) THEN
            CALL LeftLU5divStep2(D0mat,V(:,iz,ix,3),REQvetaTOuvw(cnt_req)); cnt_req=cnt_req+1
            temp=(ialfa(ix)*V(:,iz,ix,3)-ibeta(iz)*V(:,iz,ix,1))/k2(iz,ix)
            V(:,iz,ix,3)=(ibeta(iz)*V(:,iz,ix,3)+ialfa(ix)*V(:,iz,ix,1))/k2(iz,ix)
            V(:,iz,ix,1)=temp
          END IF
      END DO
    END DO
    IF (.NOT. last) CALL MPI_Waitall(cnt_req-1,REQvetaTOuvw,MPI_STATUSES_IGNORE)
  END SUBROUTINE vetaTOuvw

