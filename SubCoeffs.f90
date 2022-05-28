subroutine coeff_calc
    !(geometry_inp,xsec_inp,coeff,ngroup,nx,ny,B_C)
                use dataset, only : ngroup,nx,ny,nz,D,sig_R,nu_sgf,sig_abs ,BK
                 use typedata
                 implicit none
                      !TYPE (geometry), intent(in) :: geometry_inp
                      !TYPE (xsection), intent(in) :: xsec_inp
                      !TYPE (coeffs),allocatable,dimension(:),intent(out) :: coeff
                      !TYPE (boundary), intent(in) :: B_C
                      !integer, intent(out) :: ngroup,nx,ny
  !                 integer, parameter :: dp = selected_real_kind(15, 307)
              integer :: i,j,l,g,grp,n,nregx,nregy,nregz,x,xo,y,yo,z,zo,mi,ii,jj,ll,it
              real(dp) :: gL,gR,gT,gB,gF,gBa
              
              integer, allocatable, dimension(:) :: nmeshx,nmeshy,nmeshz 
             ! real(dp), allocatable, dimension(:,:,:) ::  D,sig_R,nu_sgf
              real(dp), allocatable, dimension(:,:,:) :: sig_scatt,Scatt
              real(dp), allocatable, dimension(:) :: R,dX,dY,dZ 
              integer, allocatable, dimension(:,:) :: mcore
              integer, allocatable, dimension(:,:,:) :: mcore3d
                
                   
  !               start, finish,,tol,e,k,k_prev,w,diff,toll,phi phi_prev,s_prev,so phi1,phi2 
   !********************************************************************************************!
   !-------------------------1D Coefficients calculation----------------------------------------!
   !____________________________________________________________________________________________!           
            if (geometry_inp%gmtry==1) then
             ALLOCATE (coeff(geometry_inp%group))
               nregx= geometry_inp%nreg(1)
                allocate(nmeshx(nregx))
            do i=1,nregx
              nmeshx(i)=geometry_inp%reglen(1,i)/geometry_inp%msize(1,i)
            enddo
                nx=sum(nmeshx)
                
                grp=geometry_inp%group
            allocate(dX(nx),mcore(1,nx),R(grp),scatt(grp,1,nx))
            allocate(D(grp,1,nx),sig_R(grp,1,nx),nu_sgf(grp,1,nx),sig_abs(grp,1,nx))    
             do i=1,grp
                  allocate(coeff(i)%aC(1,nx),coeff(i)%aL(1,nx),coeff(i)%aR(1,nx))
                  allocate(coeff(i)%aT(1,nx),coeff(i)%aB(1,nx),coeff(i)%F(1,nx),coeff(i)%sig_down(1,nx))
                    
                        do j=1,nx
                            allocate(coeff(i)%sig_down(1,j)%scatt(grp))
                        enddo
                    
             enddo   
              mcore=0  
                  x=0
         do j=1,nregx
             xo=x+1
             x=x+nmeshx(j)
             do i=xo,x
                 dX(i)=geometry_inp%msize(1,j)
                 mcore(1,i) =geometry_inp%core(1,1,j)
             enddo
         enddo   
             
           !****************************************************************!      
         write(11,*) 'Core Configuration'
         
             write(11,'(50I3)') (geometry_inp%core(1,1,j), j=1,nregx)
         
         write(11,*) 'Mesh wise Core Map'
         
            ! print*, (mcore(1,j), j=1,nx)
         write(11,'(200I3)') (mcore(1,j), j=1,nx)
         
         write(11,'(200a)',advance='no') ('*-*',j=1,nx)
        write(11,'(a)') ' '
        
  !*****************************************************************!
              D=0.0d0
         sig_R=0.0d0
         sig_abs=0.0d0
         R=0.0d0
       do g=1,grp
         
             do j=1,nx
                 mi=mcore(1,j)
             if (mi==0) then
                 D(g,1,j)=0.0d0
                 sig_R(g,1,j)=0.0d0
                 coeff(g)%F(1,j)=0.0d0
                 coeff(g)%sig_down(1,j)%scatt(:)=0.0d0
             else
                 D(g,1,j)=xsec_inp%D1(g,mi)
                 sig_abs(g,1,j)=xsec_inp%sig_a1(g,mi)*dX(j)
                  R(g)=0.0d0
                 do l=g+1,grp
                     R(g)=R(g)+xsec_inp%sig_scatt1(mi,g,l)
                 enddo
                 do l=1,g-1
                     R(g)=R(g)+xsec_inp%sig_scatt1(mi,g,l)
                 enddo
                ! print*,R(g)
                 sig_R(g,1,j)=(xsec_inp%sig_a1(g,mi)+R(g)) *dX(j)
                 !print*,sig_R(g,1,j)
                 coeff(g)%F(1,j)=xsec_inp%nu_sgf1(g,mi) *dX(j)
                 coeff(g)%sig_down(1,j)%scatt(:)=0.0d0
                 do l=1,g-1
                 coeff(g)%sig_down(1,j)%scatt(l)=xsec_inp%sig_scatt1(mi,l,g)*dX(j)
                 enddo
                 do l=g+1,grp
                 coeff(g)%sig_down(1,j)%scatt(l)=xsec_inp%sig_scatt1(mi,l,g)*dX(j)
                 enddo
             endif
             enddo 
       enddo
       
        !  open(unit=30, file='scatdown.txt')
        !do i=1,nx
        !    write(30,'(10f10.6)')  (coeff(2)%sig_down(1,i)%scatt(l), l=1,grp)
        !   ! print*, (coeff(g)%sig_down(1,1)%scatt(l), l=1,grp)
        !enddo
        !stop
        ! open(unit=30, file='removal.txt')
        !do g=1,grp
        !    write(30,'(10f10.6)')  sig_R(g,1,1)
        !enddo
        !stop
            do  i=1,grp
       coeff(i)%aC=0.0d0
       coeff(i)%aT=0.0d0
       coeff(i)%aB=0.0d0
       coeff(i)%aL=0.0d0
       coeff(i)%aR=0.0d0  
       enddo
!    Boundary Conditions       
        if (B_C%Left==1) then
        gL=00000.0d0
        elseif (B_C%Left==0) then
        gL=1d5
        endif 
         if (B_C%Right==1) then
        gR=00000.0d0
        elseif (B_C%Right==0) then
        gR=1d5
        endif
        
        do g=1,grp
            do j=2,nx-1
           coeff(g)%aL(1,j)=(-2*D(g,1,j-1)*D(g,1,j))/(D(g,1,j-1)*dX(j)+D(g,1,j)*dX(j-1))
           coeff(g)%aR(1,j)=(-2*D(g,1,j)*D(g,1,j+1))/(D(g,1,j)*dX(j+1)+D(g,1,j+1)*dX(j))  
           coeff(g)%aC(1,j)=sig_R(g,1,j)-coeff(g)%aL(1,j)-coeff(g)%aR(1,j)
            enddo
        
        ! Left Boundary  
        coeff(g)%aR(1,1)=coeff(g)%aL(1,2)
        coeff(g)%aC(1,1)=sig_R(g,1,1)-coeff(g)%aR(1,1) +gL*D(g,1,1)/(D(g,1,1)+gL*0.5*dX(1))
        !+(D(g,1,1)*gL)/(1+gL*dX(1)/2)
       
         ! Right Boundary
        coeff(g)%aL(1,nx)=(-2*D(g,1,nx-1)*D(g,1,nx))/(D(g,1,nx-1)*dX(nx)+D(g,1,nx)*dX(nx-1))
       
        coeff(g)%aC(1,nx)=sig_R(g,1,nx)-coeff(g)%aL(1,nx) +gR*D(g,1,nx)/(D(g,1,nx)+gR*0.5*dX(nx))  
        !+(D(g,1,nx)*gR)/(1+gR*dX(nx)/2)
                           
        enddo
        ngroup=grp
        
        open(unit=30, file='coeffs21.txt')
        do i=1,nx
           ! write(30,'(3(F10.6))') D(1,1,i),sig_R(1,1,j),coeff(1)%F(1,j)
          write(30,'(6(F10.6))') coeff(1)%aL(1,i),coeff(1)%aC(1,i),coeff(1)%aR(1,i),coeff(2)%aL(1,i),coeff(2)%aC(1,i),coeff(2)%aR(1,i)
        enddo
        close(30)
!********************************************************************************************!
!-------------------------2D Coefficients calculation----------------------------------------!
!____________________________________________________________________________________________!         
        elseif (geometry_inp%gmtry==2) then
             ALLOCATE (coeff(geometry_inp%group))         
             nregx= geometry_inp%nreg(1)
             nregy= geometry_inp%nreg(2)
             ! print*, nregx,nregy
            
           allocate(nmeshx(nregx),nmeshy(nregy))
          do i=1,nregx
              nmeshx(i)=geometry_inp%reglen(1,i)/geometry_inp%msize(1,i)
          enddo
          do i=1,nregy
              nmeshy(i)=geometry_inp%reglen(2,i)/geometry_inp%msize(2,i)
          enddo
       !   print*, 'nmeshx',nmeshx
       !   print*, 'nmeshy',nmeshy
       !stop         
              nx=sum(nmeshx)
              ny=sum(nmeshy)
     !  print*,nx,ny;stop
!       !  pause
!              n=nx*ny
!       !       print*,n
!       !       pause
              grp=geometry_inp%group
             
              allocate(dY(ny),dX(nx),mcore(ny,nx),R(grp),scatt(grp,ny,nx))
              allocate(D(grp,ny,nx),sig_R(grp,ny,nx),nu_sgf(grp,ny,nx),sig_abs(grp,ny,nx))
              
                   do i=1,grp
                  allocate(coeff(i)%aC(ny,nx),coeff(i)%aL(ny,nx),coeff(i)%aR(ny,nx))
                  allocate(coeff(i)%aT(ny,nx),coeff(i)%aB(ny,nx),coeff(i)%F(ny,nx),coeff(i)%sig_down(ny,nx))
                    do l=1,ny
                        do j=1,nx
                            allocate(coeff(i)%sig_down(l,j)%scatt(grp))
                        enddo
                    enddo
                  enddo
              x=0
         do j=1,nregx
             xo=x+1
             x=x+nmeshx(j)
             do i=xo,x
                 dX(i)=geometry_inp%msize(1,j)
             enddo
         enddo
           y=0
         do j=1,nregy
             yo=y+1
             y=y+nmeshy(j)
             do i=yo,y
                 dY(i)=geometry_inp%msize(2,j)
             enddo
         enddo
   !       print*,dY;stop
   !       pause
         mcore=0
         x=0
         y=0
         
         do i=1,nregy
               yo=y+1
                 y=y+nmeshy(i) 
             do j=1,nregx
                 if (j==1) then
                     x=0
                 endif 
                    xo=x+1
               x=x+nmeshx(j)
             !   print*,xo,yo
             !   print*,
                do jj=xo,x
                    do ii=yo,y 
                        mcore(ii,jj) =geometry_inp%core(1,nregy+1-i,j)
                        !nregy+1-i
                       !  print*,mi
                       !  pause
                       !  mcore(ii,jj)=mi
                     enddo
                 enddo
             enddo
         enddo
   !****************************************************************!      
         write(11,*) 'Core Configuration'
         do i=1,nregy
             write(11,'(100I3)') (geometry_inp%core(1,i,j), j=1,nregx)
         enddo
         write(11,*) 'Mesh wise Core Map'
         do i=1,ny
 !            print*, (mcore(i,j), j=1,nx)
         write(11,'(500I3)') (mcore(ny+1-i,j), j=1,nx)
         enddo
         write(11,'(500a)',advance='no') ('*-*',j=1,nx)
        write(11,'(a)') ' '
        
  !*****************************************************************!
         
         D=0.0d0
         sig_R=0.0d0
         sig_abs=0.0d0
         R=0.0d0
       do g=1,grp
         do i=1,ny
             do j=1,nx
                 mi=mcore(i,j)
             if (mi==0) then
                 D(g,i,j)=0.0d0
                 sig_R(g,i,j)=0.0d0
                 coeff(g)%F(i,j)=0.0d0
                 coeff(g)%sig_down(i,j)%scatt(:)=0.0d0
             else
                 D(g,i,j)=xsec_inp%D1(g,mi)
                 sig_abs(g,i,j)=xsec_inp%sig_a1(g,mi)*dX(j)*dY(i)
                  R(g)=0.0d0
                 do l=g+1,grp
                     R(g)=R(g)+xsec_inp%sig_scatt1(mi,g,l)
                 enddo
                 do l=1,g-1
                     R(g)=R(g)+xsec_inp%sig_scatt1(mi,g,l)
                 enddo
                ! print*,R(g);stop
                 sig_R(g,i,j)=(xsec_inp%sig_a1(g,mi)+R(g)+D(g,i,j)*BK)*dX(j)*dY(i)
                ! 
                 coeff(g)%F(i,j)=xsec_inp%nu_sgf1(g,mi) *dX(j)*dY(i)
                 coeff(g)%sig_down(i,j)%scatt(:)=0.0d0
                 do l=g+1,grp
                 coeff(g)%sig_down(i,j)%scatt(l)=xsec_inp%sig_scatt1(mi,l,g)*dX(j)*dY(i)
                 enddo
                 do l=1,g-1
                 coeff(g)%sig_down(i,j)%scatt(l)=xsec_inp%sig_scatt1(mi,l,g)*dX(j)*dY(i)
                 enddo
!!!!       Included up scattering !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 
             endif
             enddo 
         enddo
       enddo
       !do i=1,1
       !    do j=1,nx
       !print*, coeff(2)%sig_down(i,j)%scatt(1)
       !    enddo
       !enddo
        !open(unit=30, file='scatdown.txt')
        !do g=1,grp
        !    write(30,'(10f12.8)')  (coeff(g)%sig_down(1,1)%scatt(l), l=1,grp)
        !    print*,  (coeff(g)%sig_down(1,1)%scatt(l), l=1,grp)
        !enddo
        !stop
       !open(unit=30, file='nuf.txt')
       !    write(30,*) 'nuf1 21'
       !    do i=1,ny
       !        write(30,'(200f10.6)') (coeff(1)%F(i,j),j=1,nx)
       !    enddo
       !    write(30,*) 'nuf2 21'
       !     do i=1,ny
       !        write(30,'(200f10.6)') (coeff(2)%F(i,j),j=1,nx)
       !     enddo
       !     close(30)
       !     stop
       !stop
        !open(unit=20, file='SGR.txt')
        !   write(20,*) 'sgr1 21'
        !   do i=1,ny
        !       write(20,'(100f8.2)') (sig_R(1,i,j),j=1,nx)
        !   enddo
        !   write(20,*) 'sgr2 21'
        !    do i=1,ny
        !       write(20,'(100f8.2)') (sig_R(2,i,j),j=1,nx)
        !    enddo
        !    close(20)
        !    stop
!      write(11,*) 'nusgf'
!      do i=1,ny
!      write(11,'(50f6.3)') (nu_sgf(1,i,j), j=1,nx)
!      enddo
!   !   stop
    ! write(11,*) 'Sig12'
    !do i=1,ny
    ! write(11,'(50f6.3)') (coeff(1)%sig12(i,j), j=1,nx)
    !enddo;stop
! !    pause
       do i=1,grp
       coeff(i)%aC=0.0d0
       coeff(i)%aT=0.0d0
       coeff(i)%aB=0.0d0
       coeff(i)%aL=0.0d0
       coeff(i)%aR=0.0d0  
       
       enddo
!    Boundary Conditions       
        if (B_C%Left==1) then
        gL=00000.0d0
        elseif (B_C%Left==0) then
        gL=100000.0d0
        elseif(B_C%Left==2) then
            gL=B_C%ext_len
        endif 
         if (B_C%Right==1) then
        gR=00000.0d0
        elseif (B_C%Right==0) then
        gR=100000.0d0
        elseif(B_C%Right==2) then
            gR=B_C%ext_len
        endif
         if (B_C%Top==1) then
         gT=00000.0d0
        elseif (B_C%Top==0) then
         gT=100000.0d0
        elseif(B_C%Top==2) then
            gT=B_C%ext_len
        endif
         if (B_C%Bottom==1) then
         gB=00000.0d0
        elseif (B_C%Bottom==0) then
         gB=100000.0d0
        elseif(B_C%Bottom==2) then
            gB=B_C%ext_len
        endif
        ! print*,gL,gR,gT,gB;stop
       ii=0
        
      do g=1,grp 
               coeff(g)%aL(1,1)=0.0d0
               coeff(g)%aR(1,1)=-2*dY(1)/((dX(2)/D(g,2,1))+(dX(1)/D(g,1,1)))
               coeff(g)%aB(1,1)=0.0d0
               coeff(g)%aT(1,1)=-2*dX(1)/((dY(2)/D(g,1,2))+(dY(1)/D(g,1,1)))
               coeff(g)%aC(1,1)=sig_R(g,1,1)-(coeff(g)%aL(1,1)+coeff(g)%aR(1,1)+coeff(g)%aB(1,1)+coeff(g)%aT(1,1))  +gL*dY(1)*D(g,1,1)/(D(g,1,1)+gL*0.5*dX(1))+gB*dX(1)*D(g,1,1)/(D(g,1,1)+gB*0.5*dY(1))                 
               !+D(g,1,1)*gL*dY(1)+D(g,1,1)*gB*dX(1)             +dX(1)/(1/gB+dY(1)*0.5/D(g,1,1))
               
               coeff(g)%aL(1,nx)=0.0d0
               coeff(g)%aR(1,nx)=-2*dY(nx)/((dX(2)/D(g,2,nx))+(dX(1)/D(g,1,nx)))
               coeff(g)%aB(1,nx)=-2*dX(1)/((dY(nx-1)/D(g,1,nx-1))+(dY(nx)/D(g,1,nx)))
               coeff(g)%aT(1,nx)=0.0d0
               coeff(g)%aC(1,nx)=sig_R(g,1,nx)-(coeff(g)%aL(1,nx)+coeff(g)%aR(1,nx)+coeff(g)%aB(1,nx)+coeff(g)%aT(1,nx)) +gL*dY(1)*D(g,1,nx)/(D(g,1,nx)+gL*0.5*dX(nx)) +gT*dX(nx)*D(g,1,nx)/(D(g,1,nx)+gT*0.5*dY(1))
               !+dY(1)/(1/gL+dX(nx)*0.5/D(g,1,nx))+dX(nx)/(1/gT+dY(1)*0.5/D(g,1,nx))        !+D(g,1,nx)*gL*dY(1)+D(g,1,nx)*gT*dX(ny)
               
               coeff(g)%aL(ny,1)=-2*dY(1)/((dX(ny-1)/D(g,ny-1,1))+(dX(ny)/D(g,ny,1)))
               coeff(g)%aR(ny,1)=0.0d0
               coeff(g)%aB(ny,1)=0.0d0
               coeff(g)%aT(ny,1)=-2*dX(ny)/((dY(2)/D(g,ny,2))+(dY(1)/D(g,ny,1)))
               coeff(g)%aC(ny,1)=sig_R(g,ny,1)-(coeff(g)%aL(ny,1)+coeff(g)%aR(ny,1)+coeff(g)%aB(ny,1)+coeff(g)%aT(ny,1)) +gR*dY(ny)*D(g,ny,1)/(D(g,ny,1)+gR*0.5*dX(1))+gB*dX(1)*D(g,ny,1)/(D(g,ny,1)+gB*0.5*dY(ny))
               !+dY(ny)/(1/gR+dX(1)*0.5/D(g,ny,1))+dX(1)/(1/gB+dY(ny)*0.5/D(g,ny,1))           ! +D(g,ny,1)*gR*dY(nx)+D(g,ny,1)*gB*dX(1)
                
               coeff(g)%aL(ny,nx)=-2*dY(nx)/((dX(nx-1)/D(g,ny-1,nx))+(dX(ny)/D(g,ny,nx)))
               coeff(g)%aR(ny,nx)=0.0d0
               coeff(g)%aB(ny,nx)=-2*dX(ny)/((dY(nx-1)/D(g,ny,nx-1))+(dY(nx)/D(g,ny,nx)))
               coeff(g)%aT(ny,nx)=0.0d0
               coeff(g)%aC(ny,nx)=sig_R(g,ny,nx)-(coeff(g)%aL(ny,nx)+coeff(g)%aR(ny,nx)+coeff(g)%aB(ny,nx)+coeff(g)%aT(ny,nx)) +gR*dY(ny)*D(g,ny,nx)/(D(g,ny,nx)+gR*0.5*dX(nx)) +gT*dX(nx)*D(g,ny,nx)/(D(g,ny,nx)+gT*0.5*dY(ny))
               !+dY(ny)/(1/gR+dX(nx)*0.5/D(g,ny,nx))+dX(nx)/(1/gT+dY(ny)*0.5/D(g,ny,nx))               ! +D(g,ny,nx)*gR*dY(nx)+D(g,ny,nx)*gT*dX(ny)
       
            do j=2,nx-1   
               coeff(g)%aL(1,j)=0.0d0
               coeff(g)%aR(1,j)=-2*dY(j)/((dX(2)/D(g,2,j))+(dX(1)/D(g,1,j)))
               coeff(g)%aB(1,j)=-2*dX(1)/((dY(j-1)/D(g,1,j-1))+(dY(j)/D(g,1,j)))
               coeff(g)%aT(1,j)=-2*dX(1)/((dY(j+1)/D(g,1,j+1))+(dY(j)/D(g,1,j)))
               coeff(g)%aC(1,j)=sig_R(g,1,j)-(coeff(g)%aL(1,j)+coeff(g)%aR(1,j)+coeff(g)%aB(1,j)+coeff(g)%aT(1,j))  +gL*dY(1)*D(g,1,j)/(D(g,1,j)+gL*0.5*dX(j))
               !+dY(1)/(1/gL+dX(j)*0.5/D(g,1,j))              !+D(g,1,j)*gL*dY(1)
            enddo
            
            do i=2,ny-1
               coeff(g)%aL(i,1)=-2*dY(1)/((dX(i-1)/D(g,i-1,1))+(dX(i)/D(g,i,1)))
               coeff(g)%aR(i,1)=-2*dY(1)/((dX(i+1)/D(g,i+1,1))+(dX(i)/D(g,i,1)))
               coeff(g)%aB(i,1)=0.0d0
               coeff(g)%aT(i,1)=-2*dX(i)/((dY(2)/D(g,i,2))+(dY(1)/D(g,i,1)))
               coeff(g)%aC(i,1)=sig_R(g,i,1)-(coeff(g)%aL(i,1)+coeff(g)%aR(i,1)+coeff(g)%aB(i,1)+coeff(g)%aT(i,1))  +gB*dX(1)*D(g,i,1)/(D(g,i,1)+gB*0.5*dY(i))
               !+dX(1)/(1/gB+dY(i)*0.5/D(g,i,1))             !+D(g,i,1)*gB*dX(1)
            enddo
            
            do j=2,nx-1
               coeff(g)%aL(ny,j)=-2*dY(j)/((dX(ny-1)/D(g,ny-1,j))+(dX(ny)/D(g,ny,j)))
               coeff(g)%aR(ny,j)=0.0d0
               coeff(g)%aB(ny,j)=-2*dX(ny)/((dY(j-1)/D(g,ny,j-1))+(dY(j)/D(g,ny,j)))
               coeff(g)%aT(ny,j)=-2*dX(ny)/((dY(j+1)/D(g,ny,j+1))+(dY(j)/D(g,ny,j)))
               coeff(g)%aC(ny,j)=sig_R(g,ny,j)-(coeff(g)%aL(ny,j)+coeff(g)%aR(ny,j)+coeff(g)%aB(ny,j)+coeff(g)%aT(ny,j)) +gR*dY(ny)*D(g,ny,j)/(D(g,ny,j)+gR*0.5*dX(j))
               !+dY(ny)/(1/gR+dX(j)*0.5/D(g,ny,j))               ! +D(g,ny,j)*gR*dY(nx)
            enddo
            
            do i=2,ny-1
               coeff(g)%aL(i,nx)=-2*dY(nx)/((dX(i-1)/D(g,i-1,nx))+(dX(i)/D(g,i,nx)))
               coeff(g)%aR(i,nx)=-2*dY(nx)/((dX(i+1)/D(g,i+1,nx))+(dX(i)/D(g,i,nx)))
               coeff(g)%aB(i,nx)=-2*dX(i)/((dY(nx-1)/D(g,i,nx-1))+(dY(nx)/D(g,i,nx)))
               coeff(g)%aT(i,nx)=0.0d0
               coeff(g)%aC(i,nx)=sig_R(g,i,nx)-(coeff(g)%aL(i,nx)+coeff(g)%aR(i,nx)+coeff(g)%aB(i,nx)+coeff(g)%aT(i,nx)) +gT*dX(nx)*D(g,i,nx)/(D(g,i,nx)+gT*0.5*dY(i))
               !+dX(nx)/(1/gT+dY(i)*0.5/D(g,i,nx))               !+D(g,i,nx)*gT*dX(ny)
            enddo
            
       do i=2,ny-1
           do j=2,nx-1
               coeff(g)%aL(i,j)=-2*dY(j)/((dX(i-1)/D(g,i-1,j))+(dX(i)/D(g,i,j)))
               coeff(g)%aR(i,j)=-2*dY(j)/((dX(i+1)/D(g,i+1,j))+(dX(i)/D(g,i,j)))
               coeff(g)%aB(i,j)=-2*dX(i)/((dY(j-1)/D(g,i,j-1))+(dY(j)/D(g,i,j)))
               coeff(g)%aT(i,j)=-2*dX(i)/((dY(j+1)/D(g,i,j+1))+(dY(j)/D(g,i,j)))
               coeff(g)%aC(i,j)=sig_R(g,i,j)-(coeff(g)%aL(i,j)+coeff(g)%aR(i,j)+coeff(g)%aB(i,j)+coeff(g)%aT(i,j))
           enddo
       enddo
      enddo
      
      
      
      do g=1,grp
        do i=1,ny
           do j=1,nx
               mi=mcore(i,j)
               if (mi==0) then
               coeff(g)%aL(i,j)=0.0d0
               coeff(g)%aR(i,j)=0.0d0
               coeff(g)%aB(i,j)=0.0d0
               coeff(g)%aT(i,j)=0.0d0
               coeff(g)%aC(i,j)=1.0d0
               endif
           enddo
        enddo  
      enddo
           ngroup=grp
            
            
      !open(unit=4, File='debug1.txt')
      !write(4,*) 'Matrix aC2'
      !do i=1,ny
      !write(4,'(1600f8.2)') (coeff(1)%aC(i,j),j=1,nx)
      !enddo
      !
      ! write(4,*) 'Matrix aT2'
      !do i=1,nx
      !write(4,'(1600f8.2)') (coeff(1)%aT(i,j),j=1,nx)
      !enddo
      ! write(4,*) 'Matrix aB2'
      !do i=1,nx
      !write(4,'(1600f8.2)') (coeff(1)%aB(i,j),j=1,nx)
      !enddo
      ! write(4,*) 'Matrix aL2'
      !do i=1,nx
      !write(4,'(1600f8.2)') (coeff(1)%aL(i,j),j=1,nx)
      !enddo                  
      ! write(4,*) 'Matrix aR2'
      !do i=1,nx
      !write(4,'(1600f8.2)') (coeff(1)%aR(i,j),j=1,nx)
      !enddo
      !close(4)
      !stop
           
    !-------------------------3D Coefficients calculation----------------------------------------!
    !____________________________________________________________________________________________!         
        elseif (geometry_inp%gmtry==3) then
             ALLOCATE (coeff(geometry_inp%group))         
             nregx= geometry_inp%nreg(1)
             nregy= geometry_inp%nreg(2)
             nregz= geometry_inp%nreg(3)
             ! print*, nregx,nregy
            
           allocate(nmeshx(nregx),nmeshy(nregy),nmeshz(nregz))
          do i=1,nregx
              nmeshx(i)=geometry_inp%reglen(1,i)/geometry_inp%msize(1,i)
          enddo
          do i=1,nregy
              nmeshy(i)=geometry_inp%reglen(2,i)/geometry_inp%msize(2,i)
          enddo
          do i=1,nregz
              nmeshz(i)=geometry_inp%reglen(3,i)/geometry_inp%msize(3,i)
          enddo
       !   print*, 'nmeshx',nmeshx
       !   print*, 'nmeshy',nmeshy
       !   print*, 'nmeshz',nmeshz
       !stop         
              nx=sum(nmeshx)
              ny=sum(nmeshy)
              nz=sum(nmeshz)
       !print*,nx,ny,nz;stop
              grp=geometry_inp%group
              allocate(dY(0:ny+1),dX(0:nx+1),dZ(0:nz+1),mcore3d(nz,ny,nx),R(grp),scatt(grp,ny,nx))
             ! allocate(D(grp,ny,nx),sig_R(grp,ny,nx),nu_sgf(grp,ny,nx),sig_abs(grp,ny,nx))
              
                do g=1,grp
                  allocate(coeff(g)%aC3d(nz,ny,nx),coeff(g)%aL3d(nz,ny,nx),coeff(g)%aR3d(nz,ny,nx))
                  allocate(coeff(g)%aT3d(nz,ny,nx),coeff(g)%aB3d(nz,ny,nx),coeff(g)%F3d(nz,ny,nx),coeff(g)%sig_down3d(nz,ny,nx))
                  allocate(coeff(g)%aF3d(nz,ny,nx),coeff(g)%aBa3d(nz,ny,nx),coeff(g)%sig_R3d(nz,ny,nx))
                  allocate(coeff(g)%D3d(nz,ny,nx))
                    do i=1,ny
                        do j=1,nx
                            do l=1,nz
                            allocate(coeff(g)%sig_down3d(l,i,j)%scatt(grp))
                           enddo
                        enddo
                     enddo
                enddo
               dX=0.0d0
               dY=0.0d0
               dZ=0.0d0
              x=0
         do j=1,nregx
             xo=x+1
             x=x+nmeshx(j)
             do i=xo,x
                 dX(i)=geometry_inp%msize(1,j)
             enddo
         enddo
           y=0
         do j=1,nregy
             yo=y+1
             y=y+nmeshy(j)
             do i=yo,y
                 dY(i)=geometry_inp%msize(2,j)
             enddo
         enddo
            z=0
         do j=1,nregz
             zo=z+1
             z=z+nmeshz(j)
             do i=zo,z
                 dZ(i)=geometry_inp%msize(3,j)
             enddo
         enddo    
      ! print*, dZ; stop 
         
        mcore3d=0
         x=0
         y=0
         z=0
         do l=1,nregz
             zo=z+1
             z=z+nmeshz(l)
         do i=1,nregy
             if (i==1) then
                 y=0
             endif
             
               yo=y+1
                 y=y+nmeshy(i) 
             do j=1,nregx
                 if (j==1) then
                     x=0
                 endif 
                    xo=x+1
               x=x+nmeshx(j)
             !   print*,xo,yo
             !   print*,
            do ll=zo,z
                do jj=xo,x
                    do ii=yo,y 
                       
                        mcore3d(ll,ii,jj) =geometry_inp%core(nregz+1-l,nregy+1-i,j)
                        !nregy+1-i
                       !  print*,mi
                       !  pause
                       !  mcore(ii,jj)=mi
                     enddo
                 enddo
             enddo
         enddo    
         enddo
         
         enddo
          !****************************************************************!      
         write(11,*) 'Core Configuration'
         do l=1,nregz
         do i=1,nregy
             write(11,'(50I3)') (geometry_inp%core(l,i,j), j=1,nregx)
         enddo
         write(11,'(a)') ' '
         enddo
         
         
         write(11,*) 'Mesh wise Core Map'
         do l=1,nz
         do i=1,ny
 !            print*, (mcore(i,j), j=1,nx)
         write(11,'(200I3)') (mcore3d(nz+1-l,ny+1-i,j), j=1,nx)
         enddo
         write(11,'(a)') ' '
         enddo
         write(11,'(200a)',advance='no') ('*-*',j=1,nx)
         write(11,'(a)') ' '
        
    !*****************************************************************!  
         do g=1,grp
         coeff(g)%D3d=0.0d0
         coeff(g)%sig_R3d=0.0d0
         enddo
         !sig_abs=0.0d0
         R=0.0d0
       !  pause
       do g=1,grp
           do l=1,nz
            do i=1,ny
             do j=1,nx
                 mi=mcore3d(l,i,j)
               if (mi==0) then
                 coeff(g)%D3d(l,i,j)=0.0d0
                 coeff(g)%sig_R3d(l,i,j)=0.0d0
                 coeff(g)%F3d(l,i,j)=0.0d0
                 coeff(g)%sig_down3d(l,i,j)%scatt(:)=0.0d0
               else
                 coeff(g)%D3d(l,i,j)=xsec_inp%D1(g,mi)
                 !sig_abs(g,i,j)=xsec_inp%sig_a1(g,mi)*dX(j)*dY(i)
                  R(g)=0.0d0
                 do ll=g+1,grp
                     R(g)=R(g)+xsec_inp%sig_scatt1(mi,g,ll)
                 enddo
                 do ll=1,g-1
                     R(g)=R(g)+xsec_inp%sig_scatt1(mi,g,ll)
                 enddo
                ! print*,R(g);stop
                 coeff(g)%sig_R3d(l,i,j)=(xsec_inp%sig_a1(g,mi)+R(g))*dX(j)*dY(i)*dZ(l)
                 coeff(g)%F3d(l,i,j)=xsec_inp%nu_sgf1(g,mi) *dX(j)*dY(i)*dZ(l)
                 coeff(g)%sig_down3d(l,i,j)%scatt(:)=0.0d0
                 do ll=g+1,grp
                 coeff(g)%sig_down3d(l,i,j)%scatt(ll)=xsec_inp%sig_scatt1(mi,ll,g)*dX(j)*dY(i)*dZ(l)
                 enddo
                 do ll=1,g-1
                 coeff(g)%sig_down3d(l,i,j)%scatt(ll)=xsec_inp%sig_scatt1(mi,ll,g)*dX(j)*dY(i)*dZ(l)
                 enddo
!!!!       Included up scattering !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 
               endif
             enddo 
            enddo
           enddo 
          enddo
        
   !    Boundary Conditions       
        if (B_C%Left==1) then
        gL=00000.0d0
        elseif (B_C%Left==0) then
        gL=100000.0d0
        endif 
         if (B_C%Right==1) then
        gR=00000.0d0
        elseif (B_C%Right==0) then
        gR=100000.0d0
        endif
         if (B_C%Top==1) then
         gT=00000.0d0
        elseif (B_C%Top==0) then
         gT=100000.0d0
        endif
         if (B_C%Bottom==1) then
         gB=00000.0d0
        elseif (B_C%Bottom==0) then
         gB=100000.0d0
         endif
        if (B_C%Front==1) then
         gF=00000.0d0
        elseif (B_C%Front==0) then
         gF=100000.0d0
        endif
         if (B_C%Back==1) then
         gBa=00000.0d0
        elseif (B_C%Back==0) then
         gBa=100000.0d0
        endif
      ! print*, gL,gR,gB,gT,gBa,gF;stop 
               ii=0
        
      do g=1,grp 
             
        do l=1,nz      
          do i=1,ny
            do j=1,nx
              ! print*, l,i,j 
                
               if ((i==1).AND.(j==1).AND.(l==1)) then
                  !  print*, '1'
               coeff(g)%aL3d(l,i,j)=0.0d0
               coeff(g)%aR3d(l,i,j)=(-2*coeff(g)%D3d(l,i+1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i,j)*dX(i+1)+coeff(g)%D3d(l,i+1,j)*dX(i))
               coeff(g)%aB3d(l,i,j)=0.0d0
               coeff(g)%aT3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j+1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j)*dY(j+1)+coeff(g)%D3d(l,i,j+1)*dY(j))
               coeff(g)%aBa3d(l,i,j)=0.0d0
               coeff(g)%aF3d(l,i,j)=(-2*coeff(g)%D3d(l+1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l,i,j)*dZ(l+1)+coeff(g)%D3d(l+1,i,j)*dZ(l))
               coeff(g)%aC3d(l,i,j)=coeff(g)%sig_R3d(l,i,j)-(coeff(g)%aL3d(l,i,j)+coeff(g)%aR3d(l,i,j)+coeff(g)%aB3d(l,i,j)+coeff(g)%aT3d(l,i,j)+coeff(g)%aF3d(l,i,j)+coeff(g)%aBa3d(l,i,j))  +dY(j)*dZ(l)/(1/gL+dX(i)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dZ(l)/(1/gB+dY(j)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dY(j)/(1/gBa+dZ(l)*0.5/coeff(g)%D3d(l,i,j))
               !+coeff(g)%D3d(l,i,j)*gL*dY(j)*dZ(l)+coeff(g)%D3d(l,i,j)*gB*dX(i)*dZ(l)+coeff(g)%D3d(l,i,j)*gBa*dY(j)*dX(i) 
               
               
               else if ((i==1).AND.(j==1).AND.(l==nz)) then
                 !  print*, '2'
               coeff(g)%aL3d(l,i,j)=0.0d0
               coeff(g)%aR3d(l,i,j)=(-2*coeff(g)%D3d(l,i+1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i,j)*dX(i+1)+coeff(g)%D3d(l,i+1,j)*dX(i))
               coeff(g)%aB3d(l,i,j)=0.0d0
               coeff(g)%aT3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j+1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j)*dY(j+1)+coeff(g)%D3d(l,i,j+1)*dY(j))
               coeff(g)%aBa3d(l,i,j)=(-2*coeff(g)%D3d(l-1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l-1,i,j)*dZ(l)+coeff(g)%D3d(l,i,j)*dZ(l-1))
               coeff(g)%aF3d(l,i,j)=0.0d0
               coeff(g)%aC3d(l,i,j)=coeff(g)%sig_R3d(l,i,j)-(coeff(g)%aL3d(l,i,j)+coeff(g)%aR3d(l,i,j)+coeff(g)%aB3d(l,i,j)+coeff(g)%aT3d(l,i,j)+coeff(g)%aF3d(l,i,j)+coeff(g)%aBa3d(l,i,j))  +dY(j)*dZ(l)/(1/gL+dX(i)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dZ(l)/(1/gB+dY(j)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dY(j)/(1/gF+dZ(l)*0.5/coeff(g)%D3d(l,i,j))
              ! +coeff(g)%D3d(l,i,j)*gL*dY(j)*dZ(l)+coeff(g)%D3d(l,i,j)*gB*dX(i)*dZ(l)+coeff(g)%D3d(l,i,j)*gF*dY(j)*dX(i)
               
               else if ((i==1).AND.(j==nx).AND.(l==1)) then
                  !  print*, '3'
               coeff(g)%aL3d(l,i,j)=0.0d0
               coeff(g)%aR3d(l,i,j)=(-2*coeff(g)%D3d(l,i+1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i,j)*dX(i+1)+coeff(g)%D3d(l,i+1,j)*dX(i))
               coeff(g)%aB3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j-1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j-1)*dY(j)+coeff(g)%D3d(l,i,j)*dY(j-1))
               coeff(g)%aT3d(l,i,j)=0.0d0
               coeff(g)%aBa3d(l,i,j)=0.0d0
               coeff(g)%aF3d(l,i,j)=(-2*coeff(g)%D3d(l+1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l,i,j)*dZ(l+1)+coeff(g)%D3d(l+1,i,j)*dZ(l))
               coeff(g)%aC3d(l,i,j)=coeff(g)%sig_R3d(l,i,j)-(coeff(g)%aL3d(l,i,j)+coeff(g)%aR3d(l,i,j)+coeff(g)%aB3d(l,i,j)+coeff(g)%aT3d(l,i,j)+coeff(g)%aF3d(l,i,j)+coeff(g)%aBa3d(l,i,j)) +dY(j)*dZ(l)/(1/gL+dX(i)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dZ(l)/(1/gT+dY(j)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dY(j)/(1/gBa+dZ(l)*0.5/coeff(g)%D3d(l,i,j))
               !+coeff(g)%D3d(l,i,j)*gL*dY(j)*dZ(l)+coeff(g)%D3d(l,i,j)*gT*dX(i)*dZ(l)+coeff(g)%D3d(l,i,j)*gBa*dY(j)*dX(i)
               ! print*, '4'
               else if ((i==ny).AND.(j==1).AND.(l==1)) then
               coeff(g)%aL3d(l,i,j)=(-2*coeff(g)%D3d(l,i-1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i-1,j)*dX(i)+coeff(g)%D3d(l,i,j)*dX(i-1))
               coeff(g)%aR3d(l,i,j)=0.0d0
               coeff(g)%aB3d(l,i,j)=0.0d0
               coeff(g)%aT3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j+1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j)*dY(j+1)+coeff(g)%D3d(l,i,j+1)*dY(j))
               coeff(g)%aBa3d(l,i,j)=0.0d0
               coeff(g)%aF3d(l,i,j)=(-2*coeff(g)%D3d(l+1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l,i,j)*dZ(l+1)+coeff(g)%D3d(l+1,i,j)*dZ(l))
               coeff(g)%aC3d(l,i,j)=coeff(g)%sig_R3d(l,i,j)-(coeff(g)%aL3d(l,i,j)+coeff(g)%aR3d(l,i,j)+coeff(g)%aB3d(l,i,j)+coeff(g)%aT3d(l,i,j)+coeff(g)%aF3d(l,i,j)+coeff(g)%aBa3d(l,i,j)) +dY(j)*dZ(l)/(1/gR+dX(i)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dZ(l)/(1/gB+dY(j)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dY(j)/(1/gBa+dZ(l)*0.5/coeff(g)%D3d(l,i,j))
               !+coeff(g)%D3d(l,i,j)*gR*dY(j)*dZ(l)+coeff(g)%D3d(l,i,j)*gB*dX(i)*dZ(l)+coeff(g)%D3d(l,i,j)*gBa*dY(j)*dX(i)
                
               else if((i==1).AND.(j==nx).AND.(l==nz)) then
                  ! print*, '5'
               coeff(g)%aL3d(l,i,j)=0.0d0
               coeff(g)%aR3d(l,i,j)=(-2*coeff(g)%D3d(l,i+1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i,j)*dX(i+1)+coeff(g)%D3d(l,i+1,j)*dX(i))
               coeff(g)%aB3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j-1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j-1)*dY(j)+coeff(g)%D3d(l,i,j)*dY(j-1))
               coeff(g)%aT3d(l,i,j)=0.0d0
               coeff(g)%aBa3d(l,i,j)=(-2*coeff(g)%D3d(l-1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l-1,i,j)*dZ(l)+coeff(g)%D3d(l,i,j)*dZ(l-1))
               coeff(g)%aF3d(l,i,j)=0.0d0
               coeff(g)%aC3d(l,i,j)=coeff(g)%sig_R3d(l,i,j)-(coeff(g)%aL3d(l,i,j)+coeff(g)%aR3d(l,i,j)+coeff(g)%aB3d(l,i,j)+coeff(g)%aT3d(l,i,j)+coeff(g)%aF3d(l,i,j)+coeff(g)%aBa3d(l,i,j)) +dY(j)*dZ(l)/(1/gL+dX(i)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dZ(l)/(1/gT+dY(j)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dY(j)/(1/gF+dZ(l)*0.5/coeff(g)%D3d(l,i,j))
               !+coeff(g)%D3d(l,i,j)*gL*dY(j)*dZ(l)+coeff(g)%D3d(l,i,j)*gT*dX(i)*dZ(l)+coeff(g)%D3d(l,i,j)*gF*dY(j)*dX(i)
                  ! print*, '6'
               else if((i==ny).AND.(j==1).AND.(l==nz)) then
               coeff(g)%aL3d(l,i,j)=(-2*coeff(g)%D3d(l,i-1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i-1,j)*dX(i)+coeff(g)%D3d(l,i,j)*dX(i-1))
               coeff(g)%aR3d(l,i,j)=0.0d0
               coeff(g)%aB3d(l,i,j)=0.0d0
               coeff(g)%aT3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j+1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j)*dY(j+1)+coeff(g)%D3d(l,i,j+1)*dY(j))
               coeff(g)%aBa3d(l,i,j)=(-2*coeff(g)%D3d(l-1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l-1,i,j)*dZ(l)+coeff(g)%D3d(l,i,j)*dZ(l-1))
               coeff(g)%aF3d(l,i,j)=0.0d0
               coeff(g)%aC3d(l,i,j)=coeff(g)%sig_R3d(l,i,j)-(coeff(g)%aL3d(l,i,j)+coeff(g)%aR3d(l,i,j)+coeff(g)%aB3d(l,i,j)+coeff(g)%aT3d(l,i,j)+coeff(g)%aF3d(l,i,j)+coeff(g)%aBa3d(l,i,j)) +dY(j)*dZ(l)/(1/gR+dX(i)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dZ(l)/(1/gB+dY(j)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dY(j)/(1/gF+dZ(l)*0.5/coeff(g)%D3d(l,i,j))
               !+coeff(g)%D3d(l,i,j)*gR*dY(j)*dZ(l)+coeff(g)%D3d(l,i,j)*gB*dX(i)*dZ(l)+coeff(g)%D3d(l,i,j)*gF*dY(j)*dX(i)
               
               else if((i==ny).AND.(j==nx).AND.(l==1)) then
                  ! print*, '7'
               coeff(g)%aL3d(l,i,j)=(-2*coeff(g)%D3d(l,i-1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i-1,j)*dX(i)+coeff(g)%D3d(l,i,j)*dX(i-1))
               coeff(g)%aR3d(l,i,j)=0.0d0
               coeff(g)%aB3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j-1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j-1)*dY(j)+coeff(g)%D3d(l,i,j)*dY(j-1))
               coeff(g)%aT3d(l,i,j)=0.0d0
               coeff(g)%aBa3d(l,i,j)=0.0d0
               coeff(g)%aF3d(l,i,j)=(-2*coeff(g)%D3d(l+1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l,i,j)*dZ(l+1)+coeff(g)%D3d(l+1,i,j)*dZ(l))
               coeff(g)%aC3d(l,i,j)=coeff(g)%sig_R3d(l,i,j)-(coeff(g)%aL3d(l,i,j)+coeff(g)%aR3d(l,i,j)+coeff(g)%aB3d(l,i,j)+coeff(g)%aT3d(l,i,j)+coeff(g)%aF3d(l,i,j)+coeff(g)%aBa3d(l,i,j)) +dY(j)*dZ(l)/(1/gR+dX(i)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dZ(l)/(1/gT+dY(j)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dY(j)/(1/gBa+dZ(l)*0.5/coeff(g)%D3d(l,i,j))
               !+coeff(g)%D3d(l,i,j)*gR*dY(j)*dZ(l)+coeff(g)%D3d(l,i,j)*gT*dX(i)*dZ(l)+coeff(g)%D3d(l,i,j)*gBa*dY(j)*dX(i)
               
               else if ((i==ny).AND.(j==nx).AND.(l==nz)) then
                  ! print*, '8'
               coeff(g)%aL3d(l,i,j)=(-2*coeff(g)%D3d(l,i-1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i-1,j)*dX(i)+coeff(g)%D3d(l,i,j)*dX(i-1))
               coeff(g)%aR3d(l,i,j)=0.0d0
               coeff(g)%aB3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j-1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j-1)*dY(j)+coeff(g)%D3d(l,i,j)*dY(j-1))
               coeff(g)%aT3d(l,i,j)=0.0d0
               coeff(g)%aBa3d(l,i,j)=(-2*coeff(g)%D3d(l-1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l-1,i,j)*dZ(l)+coeff(g)%D3d(l,i,j)*dZ(l-1))
               coeff(g)%aF3d(l,i,j)=0.0d0
               coeff(g)%aC3d(l,i,j)=coeff(g)%sig_R3d(l,i,j)-(coeff(g)%aL3d(l,i,j)+coeff(g)%aR3d(l,i,j)+coeff(g)%aB3d(l,i,j)+coeff(g)%aT3d(l,i,j)+coeff(g)%aF3d(l,i,j)+coeff(g)%aBa3d(l,i,j)) +dY(j)*dZ(l)/(1/gR+dX(i)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dZ(l)/(1/gT+dY(j)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dY(j)/(1/gF+dZ(l)*0.5/coeff(g)%D3d(l,i,j))
               !+coeff(g)%D3d(l,i,j)*gR*dY(j)*dZ(l)+coeff(g)%D3d(l,i,j)*gT*dX(i)*dZ(l)+coeff(g)%D3d(l,i,j)*gF*dY(j)*dX(i)   
               
               
               else if ((i==1).AND.(j==1).AND.(l>1 .AND. l<nz)) then
                 !  print*, '9'
               coeff(g)%aL3d(l,i,j)=0.0d0
               coeff(g)%aR3d(l,i,j)=(-2*coeff(g)%D3d(l,i+1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i,j)*dX(i+1)+coeff(g)%D3d(l,i+1,j)*dX(i))
               coeff(g)%aB3d(l,i,j)=0.0d0
               coeff(g)%aT3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j+1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j)*dY(j+1)+coeff(g)%D3d(l,i,j+1)*dY(j))
               coeff(g)%aBa3d(l,i,j)=(-2*coeff(g)%D3d(l-1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l-1,i,j)*dZ(l)+coeff(g)%D3d(l,i,j)*dZ(l-1))
               coeff(g)%aF3d(l,i,j)=(-2*coeff(g)%D3d(l+1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l,i,j)*dZ(l+1)+coeff(g)%D3d(l+1,i,j)*dZ(l))
               coeff(g)%aC3d(l,i,j)=coeff(g)%sig_R3d(l,i,j)-(coeff(g)%aL3d(l,i,j)+coeff(g)%aR3d(l,i,j)+coeff(g)%aB3d(l,i,j)+coeff(g)%aT3d(l,i,j)+coeff(g)%aF3d(l,i,j)+coeff(g)%aBa3d(l,i,j)) +dY(j)*dZ(l)/(1/gL+dX(i)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dZ(l)/(1/gB+dY(j)*0.5/coeff(g)%D3d(l,i,j))
               !+coeff(g)%D3d(l,i,j)*gL*dY(j)*dZ(l)+coeff(g)%D3d(l,i,j)*gB*dX(i)*dZ(l) 
              
               else if ((i==1).AND.(j==nx).AND.(l>1 .AND. l<nz)) then
                  !  print*, '10'
               coeff(g)%aL3d(l,i,j)=0.0d0
               coeff(g)%aR3d(l,i,j)=(-2*coeff(g)%D3d(l,i+1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i,j)*dX(i+1)+coeff(g)%D3d(l,i+1,j)*dX(i))
               coeff(g)%aB3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j-1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j-1)*dY(j)+coeff(g)%D3d(l,i,j)*dY(j-1))
               coeff(g)%aT3d(l,i,j)=0.0d0
               coeff(g)%aBa3d(l,i,j)=(-2*coeff(g)%D3d(l-1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l-1,i,j)*dZ(l)+coeff(g)%D3d(l,i,j)*dZ(l-1))
               coeff(g)%aF3d(l,i,j)=(-2*coeff(g)%D3d(l+1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l,i,j)*dZ(l+1)+coeff(g)%D3d(l+1,i,j)*dZ(l))
               coeff(g)%aC3d(l,i,j)=coeff(g)%sig_R3d(l,i,j)-(coeff(g)%aL3d(l,i,j)+coeff(g)%aR3d(l,i,j)+coeff(g)%aB3d(l,i,j)+coeff(g)%aT3d(l,i,j)+coeff(g)%aF3d(l,i,j)+coeff(g)%aBa3d(l,i,j)) +dY(j)*dZ(l)/(1/gL+dX(i)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dZ(l)/(1/gT+dY(j)*0.5/coeff(g)%D3d(l,i,j))
               !+coeff(g)%D3d(l,i,j)*gL*dY(j)*dZ(l)+coeff(g)%D3d(l,i,j)*gT*dX(i)*dZ(l)
               
               else if ((i==ny).AND.(j==1).AND.(l>1 .AND. l<nz)) then
                 !  print*, '11'
               coeff(g)%aL3d(l,i,j)=(-2*coeff(g)%D3d(l,i-1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i-1,j)*dX(i)+coeff(g)%D3d(l,i,j)*dX(i-1))
               coeff(g)%aR3d(l,i,j)=0.0d0
               coeff(g)%aB3d(l,i,j)=0.0d0
               coeff(g)%aT3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j+1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j)*dY(j+1)+coeff(g)%D3d(l,i,j+1)*dY(j))
               coeff(g)%aBa3d(l,i,j)=(-2*coeff(g)%D3d(l-1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l-1,i,j)*dZ(l)+coeff(g)%D3d(l,i,j)*dZ(l-1))
               coeff(g)%aF3d(l,i,j)=(-2*coeff(g)%D3d(l+1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l,i,j)*dZ(l+1)+coeff(g)%D3d(l+1,i,j)*dZ(l))
               coeff(g)%aC3d(l,i,j)=coeff(g)%sig_R3d(l,i,j)-(coeff(g)%aL3d(l,i,j)+coeff(g)%aR3d(l,i,j)+coeff(g)%aB3d(l,i,j)+coeff(g)%aT3d(l,i,j)+coeff(g)%aF3d(l,i,j)+coeff(g)%aBa3d(l,i,j)) +dY(j)*dZ(l)/(1/gR+dX(i)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dZ(l)/(1/gB+dY(j)*0.5/coeff(g)%D3d(l,i,j))
               !+coeff(g)%D3d(l,i,j)*gR*dY(j)*dZ(l)+coeff(g)%D3d(l,i,j)*gB*dX(i)*dZ(l)
               
               else if ((i==ny).AND.(j==nx).AND.(l>1 .AND. l<nz)) then
                 !  print*, '12'
               coeff(g)%aL3d(l,i,j)=(-2*coeff(g)%D3d(l,i-1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i-1,j)*dX(i)+coeff(g)%D3d(l,i,j)*dX(i-1))
               coeff(g)%aR3d(l,i,j)=0.0d0
               coeff(g)%aB3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j-1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j-1)*dY(j)+coeff(g)%D3d(l,i,j)*dY(j-1))
               coeff(g)%aT3d(l,i,j)=0.0d0
               coeff(g)%aBa3d(l,i,j)=(-2*coeff(g)%D3d(l-1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l-1,i,j)*dZ(l)+coeff(g)%D3d(l,i,j)*dZ(l-1))
               coeff(g)%aF3d(l,i,j)=(-2*coeff(g)%D3d(l+1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l,i,j)*dZ(l+1)+coeff(g)%D3d(l+1,i,j)*dZ(l))
               coeff(g)%aC3d(l,i,j)=coeff(g)%sig_R3d(l,i,j)-(coeff(g)%aL3d(l,i,j)+coeff(g)%aR3d(l,i,j)+coeff(g)%aB3d(l,i,j)+coeff(g)%aT3d(l,i,j)+coeff(g)%aF3d(l,i,j)+coeff(g)%aBa3d(l,i,j)) +dY(j)*dZ(l)/(1/gR+dX(i)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dZ(l)/(1/gT+dY(j)*0.5/coeff(g)%D3d(l,i,j))
               !+coeff(g)%D3d(l,i,j)*gR*dY(j)*dZ(l)+coeff(g)%D3d(l,i,j)*gT*dX(i)*dZ(l)   
              
              
               else if ((i==1).AND.(j>1 .AND. j<nx).AND.(l==1)) then
                ! print*, '13'   
               coeff(g)%aL3d(l,i,j)=0.0d0
               coeff(g)%aR3d(l,i,j)=(-2*coeff(g)%D3d(l,i+1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i,j)*dX(i+1)+coeff(g)%D3d(l,i+1,j)*dX(i))
               coeff(g)%aB3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j-1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j-1)*dY(j)+coeff(g)%D3d(l,i,j)*dY(j-1))
               coeff(g)%aT3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j+1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j)*dY(j+1)+coeff(g)%D3d(l,i,j+1)*dY(j))
               coeff(g)%aBa3d(l,i,j)=0.0d0
               coeff(g)%aF3d(l,i,j)=(-2*coeff(g)%D3d(l+1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l,i,j)*dZ(l+1)+coeff(g)%D3d(l+1,i,j)*dZ(l))
               coeff(g)%aC3d(l,i,j)=coeff(g)%sig_R3d(l,i,j)-(coeff(g)%aL3d(l,i,j)+coeff(g)%aR3d(l,i,j)+coeff(g)%aB3d(l,i,j)+coeff(g)%aT3d(l,i,j)+coeff(g)%aF3d(l,i,j)+coeff(g)%aBa3d(l,i,j)) +dY(j)*dZ(l)/(1/gL+dX(i)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dY(j)/(1/gBa+dZ(l)*0.5/coeff(g)%D3d(l,i,j))
               !+coeff(g)%D3d(l,i,j)*gL*dY(j)*dZ(l)+coeff(g)%D3d(l,i,j)*gBa*dY(j)*dX(i)   
               
               else if ((i==ny).AND.(j>1 .AND. j<nx).AND.(l==1)) then
                 !  print*, '14'
               coeff(g)%aL3d(l,i,j)=(-2*coeff(g)%D3d(l,i-1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i-1,j)*dX(i)+coeff(g)%D3d(l,i,j)*dX(i-1))
               coeff(g)%aR3d(l,i,j)=0.0d0
               coeff(g)%aB3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j-1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j-1)*dY(j)+coeff(g)%D3d(l,i,j)*dY(j-1))
               coeff(g)%aT3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j+1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j)*dY(j+1)+coeff(g)%D3d(l,i,j+1)*dY(j))
               coeff(g)%aBa3d(l,i,j)=0.0d0
               coeff(g)%aF3d(l,i,j)=(-2*coeff(g)%D3d(l+1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l,i,j)*dZ(l+1)+coeff(g)%D3d(l+1,i,j)*dZ(l))
               coeff(g)%aC3d(l,i,j)=coeff(g)%sig_R3d(l,i,j)-(coeff(g)%aL3d(l,i,j)+coeff(g)%aR3d(l,i,j)+coeff(g)%aB3d(l,i,j)+coeff(g)%aT3d(l,i,j)+coeff(g)%aF3d(l,i,j)+coeff(g)%aBa3d(l,i,j)) +dY(j)*dZ(l)/(1/gR+dX(i)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dY(j)/(1/gBa+dZ(l)*0.5/coeff(g)%D3d(l,i,j))
               !+coeff(g)%D3d(l,i,j)*gR*dY(j)*dZ(l)+coeff(g)%D3d(l,i,j)*gBa*dY(j)*dX(i)
               
               else if ((i==1).AND.(j>1 .AND. j<nx).AND.(l==nz)) then
                 !  print*, '15'
               coeff(g)%aL3d(l,i,j)=0.0d0
               coeff(g)%aR3d(l,i,j)=(-2*coeff(g)%D3d(l,i+1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i,j)*dX(i+1)+coeff(g)%D3d(l,i+1,j)*dX(i))
               coeff(g)%aB3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j-1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j-1)*dY(j)+coeff(g)%D3d(l,i,j)*dY(j-1))
               coeff(g)%aT3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j+1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j)*dY(j+1)+coeff(g)%D3d(l,i,j+1)*dY(j))
               coeff(g)%aBa3d(l,i,j)=(-2*coeff(g)%D3d(l-1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l-1,i,j)*dZ(l)+coeff(g)%D3d(l,i,j)*dZ(l-1))
               coeff(g)%aF3d(l,i,j)=0.0d0
               coeff(g)%aC3d(l,i,j)=coeff(g)%sig_R3d(l,i,j)-(coeff(g)%aL3d(l,i,j)+coeff(g)%aR3d(l,i,j)+coeff(g)%aB3d(l,i,j)+coeff(g)%aT3d(l,i,j)+coeff(g)%aF3d(l,i,j)+coeff(g)%aBa3d(l,i,j)) +dY(j)*dZ(l)/(1/gL+dX(i)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dY(j)/(1/gF+dZ(l)*0.5/coeff(g)%D3d(l,i,j))
               !+coeff(g)%D3d(l,i,j)*gL*dY(j)*dZ(l)+coeff(g)%D3d(l,i,j)*gF*dY(j)*dX(i)
               
               else if ((i==ny).AND.(j>1 .AND. j<nx).AND.(l==nz)) then
                 !  print*, '16'
               coeff(g)%aL3d(l,i,j)=(-2*coeff(g)%D3d(l,i-1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i-1,j)*dX(i)+coeff(g)%D3d(l,i,j)*dX(i-1))
               coeff(g)%aR3d(l,i,j)=0.0d0
               coeff(g)%aB3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j-1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j-1)*dY(j)+coeff(g)%D3d(l,i,j)*dY(j-1))
               coeff(g)%aT3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j+1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j)*dY(j+1)+coeff(g)%D3d(l,i,j+1)*dY(j))
               coeff(g)%aBa3d(l,i,j)=(-2*coeff(g)%D3d(l-1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l-1,i,j)*dZ(l)+coeff(g)%D3d(l,i,j)*dZ(l-1))
               coeff(g)%aF3d(l,i,j)=0.0d0
               coeff(g)%aC3d(l,i,j)=coeff(g)%sig_R3d(l,i,j)-(coeff(g)%aL3d(l,i,j)+coeff(g)%aR3d(l,i,j)+coeff(g)%aB3d(l,i,j)+coeff(g)%aT3d(l,i,j)+coeff(g)%aF3d(l,i,j)+coeff(g)%aBa3d(l,i,j)) +dY(j)*dZ(l)/(1/gR+dX(i)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dY(j)/(1/gF+dZ(l)*0.5/coeff(g)%D3d(l,i,j))
               !+coeff(g)%D3d(l,i,j)*gR*dY(j)*dZ(l)+coeff(g)%D3d(l,i,j)*gF*dY(j)*dX(i)   
               
                 ! print*, '17'
               else if ((i>1 .AND. i<ny).AND.(j==1).AND.(l==1)) then
               coeff(g)%aL3d(l,i,j)=(-2*coeff(g)%D3d(l,i-1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i-1,j)*dX(i)+coeff(g)%D3d(l,i,j)*dX(i-1))
               coeff(g)%aR3d(l,i,j)=(-2*coeff(g)%D3d(l,i+1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i,j)*dX(i+1)+coeff(g)%D3d(l,i+1,j)*dX(i))
               coeff(g)%aB3d(l,i,j)=0.0d0
               coeff(g)%aT3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j+1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j)*dY(j+1)+coeff(g)%D3d(l,i,j+1)*dY(j))
               coeff(g)%aBa3d(l,i,j)=0.0d0
               coeff(g)%aF3d(l,i,j)=(-2*coeff(g)%D3d(l+1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l,i,j)*dZ(l+1)+coeff(g)%D3d(l+1,i,j)*dZ(l))
               coeff(g)%aC3d(l,i,j)=coeff(g)%sig_R3d(l,i,j)-(coeff(g)%aL3d(l,i,j)+coeff(g)%aR3d(l,i,j)+coeff(g)%aB3d(l,i,j)+coeff(g)%aT3d(l,i,j)+coeff(g)%aF3d(l,i,j)+coeff(g)%aBa3d(l,i,j))  +dX(i)*dZ(l)/(1/gB+dY(j)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dY(j)/(1/gBa+dZ(l)*0.5/coeff(g)%D3d(l,i,j))
               !+coeff(g)%D3d(l,i,j)*gB*dX(i)*dZ(l)+coeff(g)%D3d(l,i,j)*gBa*dY(j)*dX(i) 
              
               else if ((i>1 .AND. i<ny).AND.(j==nx).AND.(l==1)) then
                   ! print*, '18'
               coeff(g)%aL3d(l,i,j)=(-2*coeff(g)%D3d(l,i-1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i-1,j)*dX(i)+coeff(g)%D3d(l,i,j)*dX(i-1))
               coeff(g)%aR3d(l,i,j)=(-2*coeff(g)%D3d(l,i+1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i,j)*dX(i+1)+coeff(g)%D3d(l,i+1,j)*dX(i))
               coeff(g)%aB3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j-1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j-1)*dY(j)+coeff(g)%D3d(l,i,j)*dY(j-1))
               coeff(g)%aT3d(l,i,j)=0.0d0
               coeff(g)%aBa3d(l,i,j)=0.0d0
               coeff(g)%aF3d(l,i,j)=(-2*coeff(g)%D3d(l+1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l,i,j)*dZ(l+1)+coeff(g)%D3d(l+1,i,j)*dZ(l))
               coeff(g)%aC3d(l,i,j)=coeff(g)%sig_R3d(l,i,j)-(coeff(g)%aL3d(l,i,j)+coeff(g)%aR3d(l,i,j)+coeff(g)%aB3d(l,i,j)+coeff(g)%aT3d(l,i,j)+coeff(g)%aF3d(l,i,j)+coeff(g)%aBa3d(l,i,j))  +dX(i)*dZ(l)/(1/gT+dY(j)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dY(j)/(1/gBa+dZ(l)*0.5/coeff(g)%D3d(l,i,j))
               !+coeff(g)%D3d(l,i,j)*gT*dX(i)*dZ(l)+coeff(g)%D3d(l,i,j)*gBa*dY(j)*dX(i) 
              
               else if ((i>1 .AND. i<ny).AND.(j==1).AND.(l==nz)) then 
                  !  print*, '19'
               coeff(g)%aL3d(l,i,j)=(-2*coeff(g)%D3d(l,i-1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i-1,j)*dX(i)+coeff(g)%D3d(l,i,j)*dX(i-1))
               coeff(g)%aR3d(l,i,j)=(-2*coeff(g)%D3d(l,i+1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i,j)*dX(i+1)+coeff(g)%D3d(l,i+1,j)*dX(i))
               coeff(g)%aB3d(l,i,j)=0.0d0
               coeff(g)%aT3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j+1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j)*dY(j+1)+coeff(g)%D3d(l,i,j+1)*dY(j))
               coeff(g)%aBa3d(l,i,j)=(-2*coeff(g)%D3d(l-1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l-1,i,j)*dZ(l)+coeff(g)%D3d(l,i,j)*dZ(l-1))
               coeff(g)%aF3d(l,i,j)=0.0d0
               coeff(g)%aC3d(l,i,j)=coeff(g)%sig_R3d(l,i,j)-(coeff(g)%aL3d(l,i,j)+coeff(g)%aR3d(l,i,j)+coeff(g)%aB3d(l,i,j)+coeff(g)%aT3d(l,i,j)+coeff(g)%aF3d(l,i,j)+coeff(g)%aBa3d(l,i,j))  +dX(i)*dZ(l)/(1/gB+dY(j)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dY(j)/(1/gF+dZ(l)*0.5/coeff(g)%D3d(l,i,j))
               !+coeff(g)%D3d(l,i,j)*gB*dX(i)*dZ(l)+coeff(g)%D3d(l,i,j)*gF*dY(j)*dX(i)
              
               else if ((i>1 .AND. i<ny).AND.(j==nx).AND.(l==nz)) then  
                  !  print*, '20'
               coeff(g)%aL3d(l,i,j)=(-2*coeff(g)%D3d(l,i-1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i-1,j)*dX(i)+coeff(g)%D3d(l,i,j)*dX(i-1))
               coeff(g)%aR3d(l,i,j)=(-2*coeff(g)%D3d(l,i+1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i,j)*dX(i+1)+coeff(g)%D3d(l,i+1,j)*dX(i))
               coeff(g)%aB3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j-1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j-1)*dY(j)+coeff(g)%D3d(l,i,j)*dY(j-1))
               coeff(g)%aT3d(l,i,j)=0.0d0
               coeff(g)%aBa3d(l,i,j)=(-2*coeff(g)%D3d(l-1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l-1,i,j)*dZ(l)+coeff(g)%D3d(l,i,j)*dZ(l-1))
               coeff(g)%aF3d(l,i,j)=0.0d0
               coeff(g)%aC3d(l,i,j)=coeff(g)%sig_R3d(l,i,j)-(coeff(g)%aL3d(l,i,j)+coeff(g)%aR3d(l,i,j)+coeff(g)%aB3d(l,i,j)+coeff(g)%aT3d(l,i,j)+coeff(g)%aF3d(l,i,j)+coeff(g)%aBa3d(l,i,j))  +dX(i)*dZ(l)/(1/gT+dY(j)*0.5/coeff(g)%D3d(l,i,j))+dX(i)*dY(j)/(1/gF+dZ(l)*0.5/coeff(g)%D3d(l,i,j))
               !+coeff(g)%D3d(l,i,j)*gT*dX(i)*dZ(l)+coeff(g)%D3d(l,i,j)*gF*dY(j)*dX(i)
                   
               
               
               else if (l==1) then
                 !  print*, '21'
               coeff(g)%aL3d(l,i,j)=(-2*coeff(g)%D3d(l,i-1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i-1,j)*dX(i)+coeff(g)%D3d(l,i,j)*dX(i-1))
               coeff(g)%aR3d(l,i,j)=(-2*coeff(g)%D3d(l,i+1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i,j)*dX(i+1)+coeff(g)%D3d(l,i+1,j)*dX(i))
               coeff(g)%aB3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j-1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j-1)*dY(j)+coeff(g)%D3d(l,i,j)*dY(j-1))
               coeff(g)%aT3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j+1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j)*dY(j+1)+coeff(g)%D3d(l,i,j+1)*dY(j))
               coeff(g)%aBa3d(l,i,j)=0.0d0
               coeff(g)%aF3d(l,i,j)=(-2*coeff(g)%D3d(l+1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l,i,j)*dZ(l+1)+coeff(g)%D3d(l+1,i,j)*dZ(l))
               coeff(g)%aC3d(l,i,j)=coeff(g)%sig_R3d(l,i,j)-(coeff(g)%aL3d(l,i,j)+coeff(g)%aR3d(l,i,j)+coeff(g)%aB3d(l,i,j)+coeff(g)%aT3d(l,i,j)+coeff(g)%aF3d(l,i,j)+coeff(g)%aBa3d(l,i,j))  +dX(i)*dY(j)/(1/gBa+dZ(l)*0.5/coeff(g)%D3d(l,i,j))
               !+coeff(g)%D3d(l,i,j)*gBa*dY(j)*dX(i)
               
               else if (l==nz) then
                 !  print*, '22'
               coeff(g)%aL3d(l,i,j)=(-2*coeff(g)%D3d(l,i-1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i-1,j)*dX(i)+coeff(g)%D3d(l,i,j)*dX(i-1))
               coeff(g)%aR3d(l,i,j)=(-2*coeff(g)%D3d(l,i+1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i,j)*dX(i+1)+coeff(g)%D3d(l,i+1,j)*dX(i))
               coeff(g)%aB3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j-1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j-1)*dY(j)+coeff(g)%D3d(l,i,j)*dY(j-1))
               coeff(g)%aT3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j+1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j)*dY(j+1)+coeff(g)%D3d(l,i,j+1)*dY(j))
               coeff(g)%aBa3d(l,i,j)=(-2*coeff(g)%D3d(l-1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l-1,i,j)*dZ(l)+coeff(g)%D3d(l,i,j)*dZ(l-1))
               coeff(g)%aF3d(l,i,j)=0.0d0
               coeff(g)%aC3d(l,i,j)=coeff(g)%sig_R3d(l,i,j)-(coeff(g)%aL3d(l,i,j)+coeff(g)%aR3d(l,i,j)+coeff(g)%aB3d(l,i,j)+coeff(g)%aT3d(l,i,j)+coeff(g)%aF3d(l,i,j)+coeff(g)%aBa3d(l,i,j))  +dX(i)*dY(j)/(1/gF+dZ(l)*0.5/coeff(g)%D3d(l,i,j))
               !+coeff(g)%D3d(l,i,j)*gF*dY(j)*dX(i)
               
               else if (i==1) then
                 !  print*, '23'
               coeff(g)%aL3d(l,i,j)=0.0d0
               coeff(g)%aR3d(l,i,j)=(-2*coeff(g)%D3d(l,i+1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i,j)*dX(i+1)+coeff(g)%D3d(l,i+1,j)*dX(i))
               coeff(g)%aB3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j-1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j-1)*dY(j)+coeff(g)%D3d(l,i,j)*dY(j-1))
               coeff(g)%aT3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j+1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j)*dY(j+1)+coeff(g)%D3d(l,i,j+1)*dY(j))
               coeff(g)%aBa3d(l,i,j)=(-2*coeff(g)%D3d(l-1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l-1,i,j)*dZ(l)+coeff(g)%D3d(l,i,j)*dZ(l-1))
               coeff(g)%aF3d(l,i,j)=(-2*coeff(g)%D3d(l+1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l,i,j)*dZ(l+1)+coeff(g)%D3d(l+1,i,j)*dZ(l))
               coeff(g)%aC3d(l,i,j)=coeff(g)%sig_R3d(l,i,j)-(coeff(g)%aL3d(l,i,j)+coeff(g)%aR3d(l,i,j)+coeff(g)%aB3d(l,i,j)+coeff(g)%aT3d(l,i,j)+coeff(g)%aF3d(l,i,j)+coeff(g)%aBa3d(l,i,j))  +dY(j)*dZ(l)/(1/gL+dX(i)*0.5/coeff(g)%D3d(l,i,j))
               !+coeff(g)%D3d(l,i,j)*gL*dY(j)*dZ(l)
               
               else if (i==ny) then 
                 !  print*, '24'
               coeff(g)%aL3d(l,i,j)=(-2*coeff(g)%D3d(l,i-1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i-1,j)*dX(i)+coeff(g)%D3d(l,i,j)*dX(i-1))
               coeff(g)%aR3d(l,i,j)=0.0d0
               coeff(g)%aB3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j-1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j-1)*dY(j)+coeff(g)%D3d(l,i,j)*dY(j-1))
               coeff(g)%aT3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j+1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j)*dY(j+1)+coeff(g)%D3d(l,i,j+1)*dY(j))
               coeff(g)%aBa3d(l,i,j)=(-2*coeff(g)%D3d(l-1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l-1,i,j)*dZ(l)+coeff(g)%D3d(l,i,j)*dZ(l-1))
               coeff(g)%aF3d(l,i,j)=(-2*coeff(g)%D3d(l+1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l,i,j)*dZ(l+1)+coeff(g)%D3d(l+1,i,j)*dZ(l))
               coeff(g)%aC3d(l,i,j)=coeff(g)%sig_R3d(l,i,j)-(coeff(g)%aL3d(l,i,j)+coeff(g)%aR3d(l,i,j)+coeff(g)%aB3d(l,i,j)+coeff(g)%aT3d(l,i,j)+coeff(g)%aF3d(l,i,j)+coeff(g)%aBa3d(l,i,j))  +dY(j)*dZ(l)/(1/gR+dX(i)*0.5/coeff(g)%D3d(l,i,j))
               !+coeff(g)%D3d(l,i,j)*gR*dY(j)*dZ(l)
               
               else if (j==1) then
                 !  print*, '25'
               coeff(g)%aL3d(l,i,j)=(-2*coeff(g)%D3d(l,i-1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i-1,j)*dX(i)+coeff(g)%D3d(l,i,j)*dX(i-1))
               coeff(g)%aR3d(l,i,j)=(-2*coeff(g)%D3d(l,i+1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i,j)*dX(i+1)+coeff(g)%D3d(l,i+1,j)*dX(i))
               coeff(g)%aB3d(l,i,j)=0.0d0
               coeff(g)%aT3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j+1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j)*dY(j+1)+coeff(g)%D3d(l,i,j+1)*dY(j))
               coeff(g)%aBa3d(l,i,j)=(-2*coeff(g)%D3d(l-1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l-1,i,j)*dZ(l)+coeff(g)%D3d(l,i,j)*dZ(l-1))
               coeff(g)%aF3d(l,i,j)=(-2*coeff(g)%D3d(l+1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l,i,j)*dZ(l+1)+coeff(g)%D3d(l+1,i,j)*dZ(l))
               coeff(g)%aC3d(l,i,j)=coeff(g)%sig_R3d(l,i,j)-(coeff(g)%aL3d(l,i,j)+coeff(g)%aR3d(l,i,j)+coeff(g)%aB3d(l,i,j)+coeff(g)%aT3d(l,i,j)+coeff(g)%aF3d(l,i,j)+coeff(g)%aBa3d(l,i,j))  +dX(i)*dZ(l)/(1/gB+dY(j)*0.5/coeff(g)%D3d(l,i,j))
               !+coeff(g)%D3d(l,i,j)*gB*dX(i)*dZ(l)
              
               else if (j==nx) then
                  !  print*, '26'
               coeff(g)%aL3d(l,i,j)=(-2*coeff(g)%D3d(l,i-1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i-1,j)*dX(i)+coeff(g)%D3d(l,i,j)*dX(i-1))
               coeff(g)%aR3d(l,i,j)=(-2*coeff(g)%D3d(l,i+1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i,j)*dX(i+1)+coeff(g)%D3d(l,i+1,j)*dX(i))
               coeff(g)%aB3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j-1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j-1)*dY(j)+coeff(g)%D3d(l,i,j)*dY(j-1))
               coeff(g)%aT3d(l,i,j)=0.0d0
               coeff(g)%aBa3d(l,i,j)=(-2*coeff(g)%D3d(l-1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l-1,i,j)*dZ(l)+coeff(g)%D3d(l,i,j)*dZ(l-1))
               coeff(g)%aF3d(l,i,j)=(-2*coeff(g)%D3d(l+1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l,i,j)*dZ(l+1)+coeff(g)%D3d(l+1,i,j)*dZ(l))
               coeff(g)%aC3d(l,i,j)=coeff(g)%sig_R3d(l,i,j)-(coeff(g)%aL3d(l,i,j)+coeff(g)%aR3d(l,i,j)+coeff(g)%aB3d(l,i,j)+coeff(g)%aT3d(l,i,j)+coeff(g)%aF3d(l,i,j)+coeff(g)%aBa3d(l,i,j))  +dX(i)*dZ(l)/(1/gT+dY(j)*0.5/coeff(g)%D3d(l,i,j))
               !+coeff(g)%D3d(l,i,j)*gT*dX(i)*dZ(l)   
               
               
               else
                 ! print*, '27'
               coeff(g)%aL3d(l,i,j)=(-2*coeff(g)%D3d(l,i-1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i-1,j)*dX(i)+coeff(g)%D3d(l,i,j)*dX(i-1))
               coeff(g)%aR3d(l,i,j)=(-2*coeff(g)%D3d(l,i+1,j)*coeff(g)%D3d(l,i,j)*dY(j)*dZ(l))/(coeff(g)%D3d(l,i,j)*dX(i+1)+coeff(g)%D3d(l,i+1,j)*dX(i))
               coeff(g)%aB3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j-1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j-1)*dY(j)+coeff(g)%D3d(l,i,j)*dY(j-1))
               coeff(g)%aT3d(l,i,j)=(-2*coeff(g)%D3d(l,i,j+1)*coeff(g)%D3d(l,i,j)*dX(i)*dZ(l))/(coeff(g)%D3d(l,i,j)*dY(j+1)+coeff(g)%D3d(l,i,j+1)*dY(j))
               coeff(g)%aBa3d(l,i,j)=(-2*coeff(g)%D3d(l-1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l-1,i,j)*dZ(l)+coeff(g)%D3d(l,i,j)*dZ(l-1))
               coeff(g)%aF3d(l,i,j)=(-2*coeff(g)%D3d(l+1,i,j)*coeff(g)%D3d(l,i,j)*dX(i)*dY(j))/(coeff(g)%D3d(l,i,j)*dZ(l+1)+coeff(g)%D3d(l+1,i,j)*dZ(l))
               coeff(g)%aC3d(l,i,j)=coeff(g)%sig_R3d(l,i,j)-(coeff(g)%aL3d(l,i,j)+coeff(g)%aR3d(l,i,j)+coeff(g)%aB3d(l,i,j)+coeff(g)%aT3d(l,i,j)+coeff(g)%aF3d(l,i,j)+coeff(g)%aBa3d(l,i,j))
               endif
               
            enddo
          enddo
        enddo
      enddo
      !pause
      !do g=1,grp
      !  do l=1,nz
      !    do i=1,ny
      !      do j=1,nx
      !         mi=mcore3d(l,i,j)
      !         if (mi==0) then
      !         coeff(g)%aL3d(l,i,j)=0.0d0
      !         coeff(g)%aR3d(l,i,j)=0.0d0
      !         coeff(g)%aB3d(l,i,j)=0.0d0
      !         coeff(g)%aT3d(l,i,j)=0.0d0
      !         coeff(g)%aC3d(l,i,j)=1.0d0
      !         endif
      !     enddo
      !    enddo  
      !  enddo
      !enddo
     !open(unit=4, File='debug1.txt')
     ! write(4,*) 'Matrix aC2'
     ! do i=1,ny
     ! write(4,'(1600f8.2)') (coeff(1)%aC3d(2,i,j),j=1,nx)
     ! enddo
     ! 
     !  write(4,*) 'Matrix aT2'
     ! do i=1,ny
     ! write(4,'(1600f8.2)') (coeff(1)%aT3d(2,i,j),j=1,nx)
     ! enddo
     !  write(4,*) 'Matrix aB2'
     ! do i=1,ny
     ! write(4,'(1600f8.2)') (coeff(1)%aB3d(2,i,j),j=1,nx)
     ! enddo
     !  write(4,*) 'Matrix aL2'
     ! do i=1,ny
     ! write(4,'(1600f8.2)') (coeff(1)%aL3d(2,i,j),j=1,nx)
     ! enddo                  
     !  write(4,*) 'Matrix aR2'
     ! do i=1,ny
     ! write(4,'(1600f8.2)') (coeff(1)%aR3d(2,i,j),j=1,nx)
     ! enddo
     ! write(4,*) 'Matrix aBa2'
     ! do i=1,ny
     ! write(4,'(1600f8.2)') (coeff(1)%aBa3d(2,i,j),j=1,nx)
     ! enddo                  
     !  write(4,*) 'Matrix aF2'
     ! do i=1,ny
     ! write(4,'(1600f8.2)') (coeff(1)%aF3d(2,i,j),j=1,nx)
     ! enddo
     ! close(4)
     ! stop
           ngroup=grp
         
          endif   
                         
             end subroutine coeff_calc