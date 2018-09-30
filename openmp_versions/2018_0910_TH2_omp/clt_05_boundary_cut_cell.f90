! This CLT program is upgraded form the old version clt of swang, 
! 4th order finite difference method is employed in the R and Z directions (upgraded by wzhang), 
! while in the phi(y) direction, either finite difference or pseudo-spectrum method is used,
! in the time-advance, 4th order Runge-Kutta scheme is chosen.
! Ref: Physics of Plasmas 22, 122504 (2015); doi: 10.1063/1.4936977

! recently added features (upgrader by hwzhang):
! 1. cut-cell method is used in the boundary,
! 2. can be used for any triangularity like East or the Circle case. Fixed boundary is used,
! 3. can be used for any nprx and nprz ( grids and nprx(z) should be divisible),
! 4. can be used for including SOL region from EFIT-gfile equilibriums,
! 5. subroutines for RMP-EAST coils is added.

! main parameters and notations
! mmode=2
! hall=flase
! fdiout=5.e-2
! mxt=256,myt=32,mzt=256
! cfl=1.2
! kap0=5.e-5
! x(1-8):total of rho,p,vx,vy(v_phi),vz,bx,by(b_phi),bz
! x1(1-8):pertubation of rho,p,vx,vy(v_phi),vz,bx,by(b_phi),bz
! xint(1-8):equilibrium of rho,p,vx,vy(v_phi),vz,bx,by(b_phi),bz
! cur(1-3):pertubation of current _(x,y(phi),z)
! cint(1-3):equilibrium of current _(x,y(phi),z)
! ef(1-3):pertubation of e-field _(x,y(phi),z)

! ****coordinates***
! xx(mx),yy(my),zz(mz): coordinates(r,phi,z) in each processes; 
! xxt(mxt),yyt(myt),zzt(mzt) for total; 
! xxst(n2th+5,npsi),zzst(n2th+5,npsi) in (theta,psi) grid; 
! xxs(n2th+5,mps4:mps),zzs(n2th+5,mps4:mps) in (theta,psi) bandary grid;

! thxz(mx,mz): theta coordinates in (r,z) grid; tht(mxt,mzt) for total; 
! tpxz(mx,mz): r.z.<->s(psi).p(pol). transit angle in (r,z); tpt(mxt,mzt) for total;
! tcxz(mx,mz): tc=ing(jcb/r^2)dth 

! thst(n2th+5): theta coordinates in (theta,psi) grid;
! tpst(n2th+5,npsi): r.z.<->s(psi).p(pol). transit angle in (theta,psi); tps(n2th+5,mps4:mps) for bndry;
! usually,
! th=arc/n2th
! tp=atan2(psdz,psdx)=atan2(bx,-bz) => cos(tp)=psdx/|grps|=-bz/bp;sin(tp)=psdz/|grps|=bx/bp;




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!       !!!!!    !!!!!  !!!!      !!!!!!!!!!!!!!!!!
!!!!!!!  !!!  !!!!!  !  !!!!  !!!!  !!!!  !!!!!!!!!!!!!!!
!!!!!!!  !!!  !!!!!  !!  !!!  !!!!  !!!!!  !!!!!!!!!!!!!!
!!!!!!!      !!!!!!  !!!  !!  !!!!  !!!!!!  !!!!!!!!!!!!!
!!!!!!!  !!!  !!!!!  !!!!  !  !!!!  !!!!!  !!!!!!!!!!!!!!
!!!!!!!  !!!!  !!!!  !!!!!    !!!!  !!!!  !!!!!!!!!!!!!!!
!!!!!!!       !!!!!  !!!!!!   !!!!      !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module for boundary subroutines.
! old method written by S. WANG.
! new method (cut cell) written by H.W. Zhang during 2018 Mar-Jul.
! sortted by H.W. Zhang during 2018 AUG.
! contact changhw@zju.edu.cn if there is any problem.

! including subroutines as follows                               -|
! bndry8_x_ex(ibnd)                                               |
! bndry8_ex(f8xz,ibnd)                                            |
! bndry3_ex(f3xz,ibnd)                                            |
! valb_atlastgrid(fxz,fst,mm,kk,ibnd)                             |> old method by S. WANG
! x1_atlastgrid_r0p1_v1(kk)                                       |
! ef_atlastgrid_r1p0_v1(kk)                                       |  
! cur_atlastgrid_r0p1_v1(kk)                                      |
! smthp_traceline_v1(kk)                                         -|
! smthp_traceline_v2(kk)                                          |
! tracelineat(jxl,jzl,jyl,lx,lz,ly,sl,vol)                        |
! smthp1_traceline(kk)                                            |> old method by S. WANG
! smthp_traceline(kk)                                             |
! tracelineat2(jxl,jzl,jyl,lx,lz,ly,sl,vol)                       |  
! tracelineat2_dx(jxl,jzl,jyl,lx,lz,ly,sl,vol)                    |
! smthp2_traceline(kk)                                           -|
! smthp_traceline_5p(kk)                                          |
! tracelineat_5p(jxl,jzl,jyl,lx,lz,ly,sl,vol)                     |
! smthp2_traceline_tm(kk)                                         |> old method by S. WANG
! tracelineat_spec1(jxl,jzl,jyl,lx,lz,ly,sl,yypl,area)            |
! smthp_traceline_spec(kk)                                        |  
! smthp1_traceline_spec(kk)                                       |
! smth_x1st(ms,me,jps,kk)                                        -|
! smth_ps1(bst,kk)                                                |
! smth_ps(bst,nt,kk)                                              |
! smthef_dis(kk)                                                  |> old method by S. WANG
! smthef_d2f(kk)                                                  |
! smthef(kk)                                                      |  
! smthef_dis_v2(kk)                                               |
! smthxzy_dis_v2(ws,mm,kk)                                       -|
! smthxzy_dis(fsm,kk)                                             |
! smthxzy(fsm,kk)                                                 |
! smthe4(fsm,kk)                                                  |> old method by S. WANG
! avrgt(fsm,kk)                                                   |
! getnp2 ( px, py, x, y, nr, lcell, lnext, xmin, ymin, &          |  
! givens ( a, b, c, s )                                           |
! qs2grd ( px, py, n, x, y, f, nr, lcell, lnext, xmin, &         -|
! qshep2 ( n, x, y, f, nq, nw, nr, lcell, lnext, xmin, &          |
! qs2val                                                         |
! rotate ( n, c, s, x, y )                                        |> old method by S. WANG
! setup2 ( xk, yk, fk, xi, yi, fi, s1, s2, r, row )               |
! store2 ( n, x, y, nr, lcell, lnext, xmin, ymin, dx, dy, ier )   |  
! timestamp ( )                                                   |
! mpi_test                                                       -|
! map_xz2st(fxz,fst,mm)                                           |
! smth_st_nrk(fstsm,js,mm,kk)                                     |
! valbm_atlastgrid_v1(fxz,mm,ibnd)                                |> old method by S. WANG
! valb8_atlastgrid_r0p1_v2(f8xz)                                  |
! valb8_atlastgrid(f8xz)                                          |  
! valb8_atlastgrid_r0p1_v1(f8xz)                                 -|
! valb3_atlastgrid(f3xz)                                          |
! valb3_atlastgrid_r0p1_v1(f3xz)                                  |
! valb3_atlastgrid_r1p0_v1(f3xz)                                 -|

! bndry8_cut_cell (not used)                                     -|
! bndry8_cut_cell_v2_fixed                                        |
! decide_grd_type_in_each_proc                                    |> cut cell method by H.W. Zhang
! calculate_dnfm_coeff_for_ir_bndx_point                          |
! decide_hyperbolic_ratio (not used)                              |
! decide_hyperbolic_ratio_v2                                     -|
! smth_irpt_with_difc(kk,f_type,ms,me,avrgh0) (not used)         -|
! smth_irpt_with_difc_v2(kk,f_type,ms,me,avrgh0)                  |
! smth_irpt_with_difc_v3(kk,f_type,ms,me,avrgh0) (not used)       |> cut cell method by H.W. Zhang
! distance_to_bnd                                                 |
! data_type_weight                                                |
! find_bnd_grd                                                   -|
! find_bnd_grd_in_each_proc                                       | 
! find_ir_pnt_bndx                                                |> cut cell method by H.W. Zhang 
! find_ir_pnt_bndz                                                | 
! decide_grd_type_bndx                                            | 
! decide_grd_type_bndz                                           -| 



!swang*****************************************************************************
!swang*****************************************************************************
!subroutines for old method
!swang*****************************************************************************
!swang*****************************************************************************

!ws************************************************************
	subroutine bndry8_x_ex(ibnd)
      use declare
      integer ibnd
      include 'mpif.h'

      do 2 jy=1,my
      x1(:,:,jy,:)=x(:,:,jy,:)-xint(:,:,:)
   2  continue
      
      call bndry8_ex(x1,ibnd)        
      if(smoothp1ll) call smthp1_traceline(3)  
      do 11 jy=1,my
      x(:,:,jy,:)=x1(:,:,jy,:)+xint(:,:,:)
   11 continue
  
      return
      end

!ws***************************************
      subroutine bndry8_ex(f8xz,ibnd)
      use declare
      integer ibnd
      real*8,dimension(mx,mz,my,8) :: f8xz
      include 'mpif.h'

      select case(ibnd)
      case(0)
      call valbm_atlastgrid_v1(f8xz(:,:,:,:),8,0)
      case(1)
      call valbm_atlastgrid_v1(f8xz(:,:,:,:),8,1)         
      case(3)
      call valbm_atlastgrid_v1(f8xz(:,:,:,1:5),5,1)
      call valbm_atlastgrid_v1(f8xz(:,:,:,6:8),3,0)
      case(10)
      call valb8_atlastgrid_r0p1_v1(f8xz)
      case(20)
      call valb8_atlastgrid_r0p1_v2(f8xz)
!      call x1_atlastgrid_r0p1_v1(3)
      case default
      call valb8_atlastgrid(f8xz)
      end select  

      call mpi_transfersm(f8xz(:,:,:,:),8) 
      if(smoothx1) then 
      do m=1,8    
       call smthxzy(f8xz(:,:,:,m),1)
      enddo
      endif

      return
      end       

!ws***************************************
      subroutine bndry3_ex(f3xz,ibnd)
      use declare
      integer ibnd
      real*8,dimension(mx,mz,my,3) :: f3xz
      include 'mpif.h'

      select case(ibnd)
      case(0)
      call valbm_atlastgrid_v1(f3xz(:,:,:,:),3,0)
      case(1)
      call valbm_atlastgrid_v1(f3xz(:,:,:,:),3,1)  
             
      case(10)
      call valb3_atlastgrid_r0p1_v1(f3xz)
!      call cur_atlastgrid_r0p1_v1(3) 
      case(11)
      call valb3_atlastgrid_r1p0_v1(f3xz)
      case default
      call valb3_atlastgrid(f3xz)
      end select  
        
      call mpi_transfersm(f3xz(:,:,:,:),3) 

      if(smoothc) then 
      do m=1,3    
       call smthxzy(f3xz(:,:,:,m),1)
      enddo
      endif

      return
      end 

!ws****************************************************************
      subroutine valb_atlastgrid(fxz,fst,mm,kk,ibnd)
      use declare
      implicit real*8 (b-h,o-z)
      integer mm,kk,ibnd
      dimension fxz(mx,mz,my,mm),fst(n2th+5,mps4:mps,my,mm),f1s(mbm_nrk,mps4:mps),wst(n2th+5)
      include 'mpif.h'

!       integer status(mpi_status_size)

!
      do 2 js=mpsa-2,mpsa-nda,-1
      if(mts_nrk(js,nrank).gt.0) then
      do 21 m=1,mm
      do 21 jy=1,my      
      call interp_xz2ps(fxz(ix_first:ix_last,iz_first:iz_last,jy,m), &
                        xx(ix_first:ix_last),zz(iz_first:iz_last),ix_last-ix_first+1,iz_last-iz_first+1, &
                       fst(itsmin(js,nrank):itsmax(js,nrank),js,jy,m), &
                        xxs(itsmin(js,nrank):itsmax(js,nrank),js),zzs(itsmin(js,nrank):itsmax(js,nrank),js),mts_nrk(js,nrank))
      
   21 continue    
!!mpi
         do irecv=1,mrkb
         do isend=1,nsend(irecv,js)
         if(nrank.eq.nranksend(irecv,js,isend) .and. nrank.ne.nrkb(irecv)) then
         ltmin=ittransmin(irecv,js,isend)
         ltmax=ittransmax(irecv,js,isend)
         if (itbmin(nrkb(irecv))+n2th .le. itsmax(js,nranksend(irecv,js,isend))) then
         ltmin=ltmin+n2th
         ltmax=ltmax+n2th
         endif
         if (itbmax(nrkb(irecv))-n2th .ge. itsmin(js,nranksend(irecv,js,isend))) then
         ltmin=ltmin-n2th
         ltmax=ltmax-n2th
         endif
         do lt=ltmin,ltmax
         call mpi_send(fst(lt,js,1:my,1:mm),my*mm, mpi_double_precision, nrkb(irecv), isend,  &
		               mpi_comm_world,ierror )
         enddo
         endif
         enddo
         enddo
      endif

      if(mb_nrk(nrank).gt.0) then  
         do irecv=1,mrkb
         do isend=1,nsend(irecv,js)
         if(nrank.eq.nrkb(irecv) .and. nrank.ne.nranksend(irecv,js,isend)) then
         ltmin=ittransmin(irecv,js,isend)
         ltmax=ittransmax(irecv,js,isend)
         do lt=ltmin,ltmax
         call mpi_recv(fst(lt,js,1:my,1:mm),my*mm, mpi_double_precision, nranksend(irecv,js,isend), isend,  &
		               mpi_comm_world,status,ierror )
         enddo
         endif
         enddo
         enddo
       endif
!!mpi
    2 continue 
      
      if(mb_nrk(nrank).gt.0) then     

      do 10 js=mpsa-2,mpsa-4,-1
!      call smth_ps_nrk(fst(itbmin(nrank):itbmax(nrank),js,:,:),itbmax(nrank)-itbmin(nrank)+1,3)
!ws_smps--------------------------
      do 11 k=1,kk
      do m=1,mm
      do jy=iy_first+2,iy_last-2
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,js,jy,m)/2.+(fst(jt+1,js,jy,m)+fst(jt-1,js,jy,m))*3./16+(fst(jt+2,js,jy,m)+fst(jt-2,js,jy,m))/16.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,js,jy,m)=wst(jt)
      enddo

      enddo
      enddo

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         call mpi_send(fst(ltmin:ltmax,js,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(inrkb(nrank)+1), 1,  &
		               mpi_comm_world,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,js,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(inrkb(nrank)-1), 1,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         call mpi_send(fst(ltmin:ltmax,js,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(1), 1,  &
		               mpi_comm_world,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,js,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(mrkb), 1,  &
		               mpi_comm_world,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         call mpi_send(fst(ltmin:ltmax,js,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(inrkb(nrank)-1), 2,  &
		               mpi_comm_world,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,js,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(inrkb(nrank)+1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         call mpi_send(fst(ltmin:ltmax,js,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(mrkb), 2,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,js,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   11 continue
!ws_smps--------------------
   10 continue

      if(ibnd==1) then
      is=mpsa-1
      do 31 m=1,mm    
      do 31 jy=1,my      
      do jt=itbmin(nrank),itbmax(nrank)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))  
!       fst(jt,is)=fst(jt,is-1)                   
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   31 continue 
!      call smth_ps_nrk(fst(itbmin(nrank):itbmax(nrank),is,:,:),itbmax(nrank)-itbmin(nrank)+1,3)
     do 32 k=1,kk
      do m=1,mm
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,is,jy,m)/2.+(fst(jt+1,is,jy,m)+fst(jt-1,is,jy,m))*3./16+(fst(jt+2,is,jy,m)+fst(jt-2,is,jy,m))/16.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,is,jy,m)=wst(jt)
      enddo

      enddo
      enddo

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         call mpi_send(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(inrkb(nrank)+1), 1,  &
		               mpi_comm_world,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(inrkb(nrank)-1), 1,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         call mpi_send(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(1), 1,  &
		               mpi_comm_world,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(mrkb), 1,  &
		               mpi_comm_world,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         call mpi_send(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(inrkb(nrank)-1), 2,  &
		               mpi_comm_world,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(inrkb(nrank)+1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         call mpi_send(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(mrkb), 2,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   32 continue
!ws_smps--------------------
      is=mpsa
      do 41 m=1,mm    
      do 41 jy=1,my  
      do jt=itbmin(nrank),itbmax(nrank)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))
!       fst(jt,is)=fst(jt,is-1)                   
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   41 continue
!      call smth_ps_nrk(fst(itbmin(nrank):itbmax(nrank),is,:,:),itbmax(nrank)-itbmin(nrank)+1,3)
!ws_smps--------------------
     do 42 k=1,kk
      do m=1,mm
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,is,jy,m)/2.+(fst(jt+1,is,jy,m)+fst(jt-1,is,jy,m))*3./16+(fst(jt+2,is,jy,m)+fst(jt-2,is,jy,m))/16.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,is,jy,m)=wst(jt)
      enddo

      enddo
      enddo

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         call mpi_send(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(inrkb(nrank)+1), 1,  &
		               mpi_comm_world,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(inrkb(nrank)-1), 1,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         call mpi_send(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(1), 1,  &
		               mpi_comm_world,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(mrkb), 1,  &
		               mpi_comm_world,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         call mpi_send(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(inrkb(nrank)-1), 2,  &
		               mpi_comm_world,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(inrkb(nrank)+1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         call mpi_send(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(mrkb), 2,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   42 continue
!ws_smps--------------------
      endif !!ibnd==1
      
      if(ibnd==0) then
      is=mpsa
      do 51 m=1,mm    
      do 51 jy=1,my  
      do jt=itbmin(nrank),itbmax(nrank)
       fst(jt,is,jy,m)=0.           
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   51 continue

      is=mpsa-1
      do 52 m=1,mm    
      do 52 jy=1,my  
      do jt=itbmin(nrank),itbmax(nrank)
      call interp1d3l(fst(jt,is-3,jy,m),fst(jt,is-2,jy,m),fst(jt,is-1,jy,m),fst(jt,is+1,jy,m), &
                        ps(is-3),ps(is-2),ps(is-1),ps(is+1),ps(is),fst(jt,is,jy,m))                     
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   52 continue

      do 53 k=1,kk
      do m=1,mm
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,is,jy,m)/2.+(fst(jt+1,is,jy,m)+fst(jt-1,is,jy,m))*3./16+(fst(jt+2,is,jy,m)+fst(jt-2,is,jy,m))/16.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,is,jy,m)=wst(jt)
      enddo

      enddo
      enddo

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         call mpi_send(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(inrkb(nrank)+1), 1,  &
		               mpi_comm_world,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(inrkb(nrank)-1), 1,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         call mpi_send(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(1), 1,  &
		               mpi_comm_world,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(mrkb), 1,  &
		               mpi_comm_world,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         call mpi_send(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(inrkb(nrank)-1), 2,  &
		               mpi_comm_world,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(inrkb(nrank)+1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         call mpi_send(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(mrkb), 2,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,is,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   53 continue
      endif !!ibnd==1 

      do 3 m=1,mm    
      do 3 jy=1,my           
      do ikb=1,mb_nrk(nrank)
        i=itb_nrk(ikb,nrank)
        ib=ib_nrk(ikb,nrank)
        do js=mpsa,mpsa-3,-1     
        call interp1d3l(fst(i-2,js,jy,m),fst(i-1,js,jy,m),fst(i,js,jy,m),fst(i+1,js,jy,m), &
                        thst(i-2),thst(i-1),thst(i),thst(i+1),thb(ib),f1s(ikb,js))
        enddo
        js=mpsa
        call interp1d3l(f1s(ikb,js),f1s(ikb,js-1),f1s(ikb,js-2),f1s(ikb,js-3), &
                        ps(js),ps(js-1),ps(js-2),ps(js-3),psb(ib),fxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,m))
!        fxz(jbx(ib),jbz(ib))=(asm_b(ib)*x1s(ib,js-1)+bsm_b(ib)*x1s(ib,js-2)+csm_b(ib)*x1s(ib,js-3)) &
!                 /(asm_b(ib)+bsm_b(ib)+csm_b(ib))
      enddo
    3 continue 
      endif !!mb_nrk(nrank).gt.0
           
!      if(ibnd==2) then
!      
!      do ika1=1,ma1_nrk(nrank)
!        i=ita1_nrk(ika1,nrank)
!        ia1=ia1_nrk(ika1,nrank)
!        do js=mpsa,mpsa-3,1
!        call interp1d3l(fst(i-2,js),fst(i-1,js),fst(i,js),fst(i+1,js), &
!                        thst(i-2),thst(i-1),thst(i),thst(i+1),tha1(ia1),x1s(ika1,js,jy))
!        enddo
!        js=mpsa
!        call interp1d3l(x1s(ika1,js,jy),x1s(ika1,js-1,jy),x1s(ika1,js-2,jy),x1s(ika1,js-3,jy), &
!                        ps(js),ps(js-1),ps(js-2),ps(js-3),psa1(ia1),fxz(jxa1_nrk(ika1,nrank),jza1_nrk(ika1,nrank)))
!!        fxz(jxa1(ia1),jza1(ia1))=(asm_a1(ia1)*x1s(ia1,js-1)+bsm_a1(ia1)*x1s(ia1,js-2)+csm_a1(ia1)*x1s(ia1,js-3)) &
!!                  /(asm_a1(ia1)+bsm_a1(ia1)+csm_a1(ia1)) 
!      enddo
!      endif
      return
      end

!ws****************************************************************
      subroutine x1_atlastgrid_r0p1_v1(kk)
      use declare
      implicit real*8 (b-h,o-z)
      integer kk
      dimension fst(n2th+5,mps4:mps,my,8),f1s(mbm_nrk,mps4:mps),fsxz(mbm_nrk,8),wst(n2th+5) !fxz(mx,mz,my,mm),
      include 'mpif.h'

!       integer status(mpi_status_size)

!
      do 2 js=mpsa-2,mpsa-nda,-1
      if(mts_nrk(js,nrank).gt.0) then
      do 21 m=1,8
      do 21 jy=1,my      
      call interp_xz2ps(x1(ix_first:ix_last,iz_first:iz_last,jy,m), &
                        xx(ix_first:ix_last),zz(iz_first:iz_last),ix_last-ix_first+1,iz_last-iz_first+1, &
                       fst(itsmin(js,nrank):itsmax(js,nrank),js,jy,m), &
                        xxs(itsmin(js,nrank):itsmax(js,nrank),js),zzs(itsmin(js,nrank):itsmax(js,nrank),js),mts_nrk(js,nrank))
      
   21 continue    
!!mpi
         do irecv=1,mrkb
         do isend=1,nsend(irecv,js)
         if(nrank.eq.nranksend(irecv,js,isend) .and. nrank.ne.nrkb(irecv)) then
         ltmin=ittransmin(irecv,js,isend)
         ltmax=ittransmax(irecv,js,isend)
         if (itbmin(nrkb(irecv))+n2th .le. itsmax(js,nranksend(irecv,js,isend))) then
         ltmin=ltmin+n2th
         ltmax=ltmax+n2th
         endif
         if (itbmax(nrkb(irecv))-n2th .ge. itsmin(js,nranksend(irecv,js,isend))) then
         ltmin=ltmin-n2th
         ltmax=ltmax-n2th
         endif
         do lt=ltmin,ltmax
         call mpi_send(fst(lt,js,1:my,:),my*8, mpi_double_precision, nrkb(irecv), isend,  &
		               mpi_comm_world,ierror )
         enddo
         endif
         enddo
         enddo
      endif

      if(mb_nrk(nrank).gt.0) then  
         do irecv=1,mrkb
         do isend=1,nsend(irecv,js)
         if(nrank.eq.nrkb(irecv) .and. nrank.ne.nranksend(irecv,js,isend)) then
         ltmin=ittransmin(irecv,js,isend)
         ltmax=ittransmax(irecv,js,isend)
         do lt=ltmin,ltmax
         call mpi_recv(fst(lt,js,1:my,:),my*8, mpi_double_precision, nranksend(irecv,js,isend), isend,  &
		               mpi_comm_world,status,ierror )
         enddo
         endif
         enddo
         enddo
       endif
!!mpi
    2 continue 
           
      if(mb_nrk(nrank).gt.0) then      

      do jy=1,my
      do js=mpsa-2,mpsa-4,-1
      do jt=itbmin(nrank),itbmax(nrank)
      vx1st=fst(jt,js,jy,3)
      vz1st=fst(jt,js,jy,5)
      bx1st=fst(jt,js,jy,6)
      bz1st=fst(jt,js,jy,8)
      fst(jt,js,jy,3)=vx1st*dcos(thst(jt))+vz1st*dsin(thst(jt))
      fst(jt,js,jy,5)=-vx1st*dsin(thst(jt))+vz1st*dcos(thst(jt))
      fst(jt,js,jy,6)=bx1st*dcos(thst(jt))+bz1st*dsin(thst(jt))
      fst(jt,js,jy,8)=-bx1st*dsin(thst(jt))+bz1st*dcos(thst(jt))
       
!      vrst=fst(jt,js,jy,3)*dcos(thst(jt))+fst(jt,js,jy,5)*dsin(thst(jt))
!      vpst=-fst(jt,js,jy,3)*dsin(thst(jt))+fst(jt,js,jy,5)*dcos(thst(jt))
!      brst=fst(jt,js,jy,6)*dcos(thst(jt))+fst(jt,js,jy,8)*dsin(thst(jt))
!      bpst=-fst(jt,js,jy,6)*dsin(thst(jt))+fst(jt,js,jy,8)*dcos(thst(jt))
!      fst(jt,js,jy,3)=vrst
!      fst(jt,js,jy,5)=vpst
!      fst(jt,js,jy,6)=brst
!      fst(jt,js,jy,8)=bpst
      enddo
      enddo
      enddo    

      do 10 js=mpsa-2,mpsa-4,-1
!      call smth_ps_nrk(fst(itbmin(nrank):itbmax(nrank),js,:,:),itbmax(nrank)-itbmin(nrank)+1,3)
!ws_smps--------------------------
      do 11 k=1,kk
      do m=1,8
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,js,jy,m)/2.+((fst(jt+1,js,jy,m)+fst(jt-1,js,jy,m))*3./16+(fst(jt+2,js,jy,m)+fst(jt-2,js,jy,m))/16.&
             +(fst(jt,js,jy+1,m)+fst(jt,js,jy-1,m))*3./16+(fst(jt,js,jy+2,m)+fst(jt,js,jy-2,m))/16.)/2.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,js,jy,m)=wst(jt)
      enddo

      enddo
      enddo

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         call mpi_send(fst(ltmin:ltmax,js,1:my,:), 2*my*8, mpi_double_precision, nrkb1(inrkb(nrank)+1), 1,  &
		               mpi_comm_world,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,js,1:my,:), 2*my*8, mpi_double_precision, nrkb1(inrkb(nrank)-1), 1,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         call mpi_send(fst(ltmin:ltmax,js,1:my,:), 2*my*8, mpi_double_precision, nrkb1(1), 1,  &
		               mpi_comm_world,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,js,1:my,:), 2*my*8, mpi_double_precision, nrkb1(mrkb), 1,  &
		               mpi_comm_world,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         call mpi_send(fst(ltmin:ltmax,js,1:my,:), 2*my*8, mpi_double_precision, nrkb1(inrkb(nrank)-1), 2,  &
		               mpi_comm_world,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,js,1:my,:), 2*my*8, mpi_double_precision, nrkb1(inrkb(nrank)+1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         call mpi_send(fst(ltmin:ltmax,js,1:my,:), 2*my*8, mpi_double_precision, nrkb1(mrkb), 2,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,js,1:my,:), 2*my*8, mpi_double_precision, nrkb1(1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   11 continue
!ws_smps--------------------
   10 continue

      do 61 m=1,8

      if(m.eq.3 .or. m.eq.6) then !!(m=3,6) ibnd==0
      is=mpsa 
      do 51 jy=1,my  
      do jt=itbmin(nrank),itbmax(nrank)
       fst(jt,is,jy,m)=0.           
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   51 continue

      is=mpsa-1 
      do 52 jy=1,my  
      do jt=itbmin(nrank),itbmax(nrank)
        call interp1d3l(fst(jt,is-3,jy,m),fst(jt,is-2,jy,m),fst(jt,is-1,jy,m),fst(jt,is+1,jy,m), &
                        ps(is-3),ps(is-2),ps(is-1),ps(is+1),ps(is),fst(jt,is,jy,m))                   
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   52 continue

      do 53 k=1,kk
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
        wst(jt)=fst(jt,is,jy,m)/2.+((fst(jt+1,is,jy,m)+fst(jt-1,is,jy,m))*3./16+(fst(jt+2,is,jy,m)+fst(jt-2,is,jy,m))/16.&
             +(fst(jt,is,jy+1,m)+fst(jt,is,jy-1,m))*3./16+(fst(jt,is,jy+2,m)+fst(jt,is,jy-2,m))/16.)/2.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,is,jy,m)=wst(jt)
      enddo

      enddo

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)+1), 1,  &
		               mpi_comm_world,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)-1), 1,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(1), 1,  &
		               mpi_comm_world,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(mrkb), 1,  &
		               mpi_comm_world,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)-1), 2,  &
		               mpi_comm_world,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)+1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(mrkb), 2,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   53 continue

      else  !!(m=1,2,4,5,7,8,) bnd==1  
      is=mpsa-1    
      do 31 jy=1,my      
      do jt=itbmin(nrank),itbmax(nrank)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))  
!       fst(jt,is)=fst(jt,is-1)                   
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   31 continue 
!      call smth_ps_nrk(fst(itbmin(nrank):itbmax(nrank),is,:,:),itbmax(nrank)-itbmin(nrank)+1,3)
      do 32 k=1,kk
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
       wst(jt)=fst(jt,is,jy,m)/2.+((fst(jt+1,is,jy,m)+fst(jt-1,is,jy,m))*3./16+(fst(jt+2,is,jy,m)+fst(jt-2,is,jy,m))/16.&
             +(fst(jt,is,jy+1,m)+fst(jt,is,jy-1,m))*3./16+(fst(jt,is,jy+2,m)+fst(jt,is,jy-2,m))/16.)/2.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,is,jy,m)=wst(jt)
      enddo

      enddo
      
         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)+1), 1,  &
		               mpi_comm_world,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)-1), 1,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(1), 1,  &
		               mpi_comm_world,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(mrkb), 1,  &
		               mpi_comm_world,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)-1), 2,  &
		               mpi_comm_world,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)+1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(mrkb), 2,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   32 continue
!ws_smps--------------------
      is=mpsa
 
      do 41 jy=1,my  
      do jt=itbmin(nrank),itbmax(nrank)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))
!       fst(jt,is)=fst(jt,is-1)                   
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   41 continue
!      call smth_ps_nrk(fst(itbmin(nrank):itbmax(nrank),is,:,:),itbmax(nrank)-itbmin(nrank)+1,3)
!ws_smps--------------------
      do 42 k=1,kk
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,is,jy,m)/2.+((fst(jt+1,is,jy,m)+fst(jt-1,is,jy,m))*3./16+(fst(jt+2,is,jy,m)+fst(jt-2,is,jy,m))/16.&
             +(fst(jt,is,jy+1,m)+fst(jt,is,jy-1,m))*3./16+(fst(jt,is,jy+2,m)+fst(jt,is,jy-2,m))/16.)/2.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,is,jy,m)=wst(jt)
      enddo

      enddo

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)+1), 1,  &
		               mpi_comm_world,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)-1), 1,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(1), 1,  &
		               mpi_comm_world,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(mrkb), 1,  &
		               mpi_comm_world,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)-1), 2,  &
		               mpi_comm_world,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)+1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(mrkb), 2,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   42 continue
!ws_smps--------------------
      endif !!(m=3,6)ibnd==0
   61 continue

      do 3 jy=1,my           
      do ikb=1,mb_nrk(nrank)
        i=itb_nrk(ikb,nrank)
        ib=ib_nrk(ikb,nrank) 
        do m=1,8 
        do js=mpsa,mpsa-3,-1     
        call interp1d3l(fst(i-2,js,jy,m),fst(i-1,js,jy,m),fst(i,js,jy,m),fst(i+1,js,jy,m), &
                        thst(i-2),thst(i-1),thst(i),thst(i+1),thb(ib),f1s(ikb,js))
        enddo
        js=mpsa
        call interp1d3l(f1s(ikb,js),f1s(ikb,js-1),f1s(ikb,js-2),f1s(ikb,js-3), &
                        ps(js),ps(js-1),ps(js-2),ps(js-3),psb(ib),fsxz(ikb,m))
        enddo
!        fxz(jbx(ib),jbz(ib))=(asm_b(ib)*x1s(ib,js-1)+bsm_b(ib)*x1s(ib,js-2)+csm_b(ib)*x1s(ib,js-3)) &
!                 /(asm_b(ib)+bsm_b(ib)+csm_b(ib))
        x1(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,1)=fsxz(ikb,1)
        x1(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,2)=fsxz(ikb,2)
        x1(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,3)=fsxz(ikb,3)*dcos(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank))) &
                                                      -fsxz(ikb,5)*dsin(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank)))
        x1(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,4)=fsxz(ikb,4)
        x1(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,5)=fsxz(ikb,3)*dsin(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank))) &
                                                      +fsxz(ikb,5)*dcos(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank)))
        x1(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,6)=fsxz(ikb,6)*dcos(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank))) &
                                                      -fsxz(ikb,8)*dsin(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank)))
        x1(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,7)=fsxz(ikb,7)
        x1(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,8)=fsxz(ikb,6)*dsin(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank))) &
                                                      +fsxz(ikb,8)*dcos(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank)))
      enddo
    3 continue 
      endif !!mb_nrk(nrank).gt.0
           
!      if(ibnd==2) then
!      
!      do ika1=1,ma1_nrk(nrank)
!        i=ita1_nrk(ika1,nrank)
!        ia1=ia1_nrk(ika1,nrank)
!        do js=mpsa,mpsa-3,1
!        call interp1d3l(fst(i-2,js),fst(i-1,js),fst(i,js),fst(i+1,js), &
!                        thst(i-2),thst(i-1),thst(i),thst(i+1),tha1(ia1),x1s(ika1,js,jy))
!        enddo
!        js=mpsa
!        call interp1d3l(x1s(ika1,js,jy),x1s(ika1,js-1,jy),x1s(ika1,js-2,jy),x1s(ika1,js-3,jy), &
!                        ps(js),ps(js-1),ps(js-2),ps(js-3),psa1(ia1),fxz(jxa1_nrk(ika1,nrank),jza1_nrk(ika1,nrank)))
!!        fxz(jxa1(ia1),jza1(ia1))=(asm_a1(ia1)*x1s(ia1,js-1)+bsm_a1(ia1)*x1s(ia1,js-2)+csm_a1(ia1)*x1s(ia1,js-3)) &
!!                  /(asm_a1(ia1)+bsm_a1(ia1)+csm_a1(ia1)) 
!      enddo
!      endif
      return
      end

!ws************************************************************
      subroutine ef_atlastgrid_r1p0_v1(kk)
      use declare
      implicit real*8 (b-h,o-z)
      integer kk
      dimension fst(n2th+5,mps4:mps,my,3),f1s(mbm_nrk,mps4:mps),fsxz(mbm_nrk,3),wst(n2th+5) !fxz(mx,mz,my,mm),
      include 'mpif.h'

!       integer status(mpi_status_size)

!
      do 2 js=mpsa-2,mpsa-nda,-1
      if(mts_nrk(js,nrank).gt.0) then
      do 21 m=1,3
      do 21 jy=1,my      
      call interp_xz2ps(ef(ix_first:ix_last,iz_first:iz_last,jy,m), &
                        xx(ix_first:ix_last),zz(iz_first:iz_last),ix_last-ix_first+1,iz_last-iz_first+1, &
                       fst(itsmin(js,nrank):itsmax(js,nrank),js,jy,m), &
                        xxs(itsmin(js,nrank):itsmax(js,nrank),js),zzs(itsmin(js,nrank):itsmax(js,nrank),js),mts_nrk(js,nrank))
      
   21 continue    
!!mpi
         do irecv=1,mrkb
         do isend=1,nsend(irecv,js)
         if(nrank.eq.nranksend(irecv,js,isend) .and. nrank.ne.nrkb(irecv)) then
         ltmin=ittransmin(irecv,js,isend)
         ltmax=ittransmax(irecv,js,isend)
         if (itbmin(nrkb(irecv))+n2th .le. itsmax(js,nranksend(irecv,js,isend))) then
         ltmin=ltmin+n2th
         ltmax=ltmax+n2th
         endif
         if (itbmax(nrkb(irecv))-n2th .ge. itsmin(js,nranksend(irecv,js,isend))) then
         ltmin=ltmin-n2th
         ltmax=ltmax-n2th
         endif
         do lt=ltmin,ltmax
         call mpi_send(fst(lt,js,1:my,:),my*3, mpi_double_precision, nrkb(irecv), isend,  &
		               mpi_comm_world,ierror )
         enddo
         endif
         enddo
         enddo
      endif

      if(mb_nrk(nrank).gt.0) then  
         do irecv=1,mrkb
         do isend=1,nsend(irecv,js)
         if(nrank.eq.nrkb(irecv) .and. nrank.ne.nranksend(irecv,js,isend)) then
         ltmin=ittransmin(irecv,js,isend)
         ltmax=ittransmax(irecv,js,isend)
         do lt=ltmin,ltmax
         call mpi_recv(fst(lt,js,1:my,:),my*3, mpi_double_precision, nranksend(irecv,js,isend), isend,  &
		               mpi_comm_world,status,ierror )
         enddo
         endif
         enddo
         enddo
       endif
!!mpi
    2 continue 
           
      if(mb_nrk(nrank).gt.0) then      

      do jy=1,my
      do js=mpsa-2,mpsa-4,-1
      do jt=itbmin(nrank),itbmax(nrank)
      exst=fst(jt,js,jy,1)
      ezst=fst(jt,js,jy,3)

      fst(jt,js,jy,1)=exst*dcos(thst(jt))+ezst*dsin(thst(jt))
      fst(jt,js,jy,3)=-exst*dsin(thst(jt))+ezst*dcos(thst(jt))
       
!      vrst=fst(jt,js,jy,3)*dcos(thst(jt))+fst(jt,js,jy,5)*dsin(thst(jt))
!      vpst=-fst(jt,js,jy,3)*dsin(thst(jt))+fst(jt,js,jy,5)*dcos(thst(jt))
!      brst=fst(jt,js,jy,6)*dcos(thst(jt))+fst(jt,js,jy,8)*dsin(thst(jt))
!      bpst=-fst(jt,js,jy,6)*dsin(thst(jt))+fst(jt,js,jy,8)*dcos(thst(jt))
!      fst(jt,js,jy,3)=vrst
!      fst(jt,js,jy,5)=vpst
!      fst(jt,js,jy,6)=brst
!      fst(jt,js,jy,8)=bpst
      enddo
      enddo
      enddo    

      do 10 js=mpsa-2,mpsa-4,-1
!      call smth_ps_nrk(fst(itbmin(nrank):itbmax(nrank),js,:,:),itbmax(nrank)-itbmin(nrank)+1,3)
!ws_smps--------------------------
      do 11 k=1,kk
      do m=1,3
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,js,jy,m)/2.+((fst(jt+1,js,jy,m)+fst(jt-1,js,jy,m))*3./16+(fst(jt+2,js,jy,m)+fst(jt-2,js,jy,m))/16.&
             +(fst(jt,js,jy+1,m)+fst(jt,js,jy-1,m))*3./16+(fst(jt,js,jy+2,m)+fst(jt,is,jy-2,m))/16.)/2.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,js,jy,m)=wst(jt)
      enddo

      enddo
      enddo

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         call mpi_send(fst(ltmin:ltmax,js,1:my,:), 2*my*3, mpi_double_precision, nrkb1(inrkb(nrank)+1), 1,  &
		               mpi_comm_world,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,js,1:my,:), 2*my*3, mpi_double_precision, nrkb1(inrkb(nrank)-1), 1,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         call mpi_send(fst(ltmin:ltmax,js,1:my,:), 2*my*3, mpi_double_precision, nrkb1(1), 1,  &
		               mpi_comm_world,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,js,1:my,:), 2*my*3, mpi_double_precision, nrkb1(mrkb), 1,  &
		               mpi_comm_world,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         call mpi_send(fst(ltmin:ltmax,js,1:my,:), 2*my*3, mpi_double_precision, nrkb1(inrkb(nrank)-1), 2,  &
		               mpi_comm_world,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,js,1:my,:), 2*my*3, mpi_double_precision, nrkb1(inrkb(nrank)+1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         call mpi_send(fst(ltmin:ltmax,js,1:my,:), 2*my*3, mpi_double_precision, nrkb1(mrkb), 2,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,js,1:my,:), 2*my*3, mpi_double_precision, nrkb1(1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   11 continue
!ws_smps--------------------
   10 continue

      do 61 m=1,3

      if(m.eq.2 .or. m.eq.3) then !!(m=2,3) ibnd==0
      is=mpsa 
      do 51 jy=1,my  
      do jt=itbmin(nrank),itbmax(nrank)
       fst(jt,is,jy,m)=0.           
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   51 continue

      is=mpsa-1 
      do 52 jy=1,my  
      do jt=itbmin(nrank),itbmax(nrank)
        call interp1d3l(fst(jt,is-3,jy,m),fst(jt,is-2,jy,m),fst(jt,is-1,jy,m),fst(jt,is+1,jy,m), &
                        ps(is-3),ps(is-2),ps(is-1),ps(is+1),ps(is),fst(jt,is,jy,m))                   
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   52 continue

      do 53 k=1,kk
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,is,jy,m)/2.+((fst(jt+1,is,jy,m)+fst(jt-1,is,jy,m))*3./16+(fst(jt+2,is,jy,m)+fst(jt-2,is,jy,m))/16.&
             +(fst(jt,is,jy+1,m)+fst(jt,is,jy-1,m))*3./16+(fst(jt,is,jy+2,m)+fst(jt,is,jy-2,m))/16.)/2.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,is,jy,m)=wst(jt)
      enddo

      enddo

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)+1), 1,  &
		               mpi_comm_world,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)-1), 1,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(1), 1,  &
		               mpi_comm_world,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(mrkb), 1,  &
		               mpi_comm_world,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)-1), 2,  &
		               mpi_comm_world,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)+1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(mrkb), 2,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   53 continue

      else  !!(m=1) bnd==1  
      is=mpsa-1    
      do 31 jy=1,my      
      do jt=itbmin(nrank),itbmax(nrank)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))  
!       fst(jt,is)=fst(jt,is-1)                   
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   31 continue 
!      call smth_ps_nrk(fst(itbmin(nrank):itbmax(nrank),is,:,:),itbmax(nrank)-itbmin(nrank)+1,3)
      do 32 k=1,kk
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,is,jy,m)/2.+((fst(jt+1,is,jy,m)+fst(jt-1,is,jy,m))*3./16+(fst(jt+2,is,jy,m)+fst(jt-2,is,jy,m))/16.&
             +(fst(jt,is,jy+1,m)+fst(jt,is,jy-1,m))*3./16+(fst(jt,is,jy+2,m)+fst(jt,is,jy-2,m))/16.)/2.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,is,jy,m)=wst(jt)
      enddo

      enddo
      
         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)+1), 1,  &
		               mpi_comm_world,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)-1), 1,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(1), 1,  &
		               mpi_comm_world,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(mrkb), 1,  &
		               mpi_comm_world,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)-1), 2,  &
		               mpi_comm_world,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)+1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(mrkb), 2,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   32 continue
!ws_smps--------------------
      is=mpsa
 
      do 41 jy=1,my  
      do jt=itbmin(nrank),itbmax(nrank)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))
!       fst(jt,is)=fst(jt,is-1)                   
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   41 continue
!      call smth_ps_nrk(fst(itbmin(nrank):itbmax(nrank),is,:,:),itbmax(nrank)-itbmin(nrank)+1,3)
!ws_smps--------------------
      do 42 k=1,kk
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,is,jy,m)/2.+((fst(jt+1,is,jy,m)+fst(jt-1,is,jy,m))*3./16+(fst(jt+2,is,jy,m)+fst(jt-2,is,jy,m))/16.&
             +(fst(jt,is,jy+1,m)+fst(jt,is,jy-1,m))*3./16+(fst(jt,is,jy+2,m)+fst(jt,is,jy-2,m))/16.)/2.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,is,jy,m)=wst(jt)
      enddo

      enddo

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)+1), 1,  &
		               mpi_comm_world,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)-1), 1,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(1), 1,  &
		               mpi_comm_world,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(mrkb), 1,  &
		               mpi_comm_world,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)-1), 2,  &
		               mpi_comm_world,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)+1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(mrkb), 2,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   42 continue
!ws_smps--------------------
      endif !!(m=2,3)ibnd==0
   61 continue

      do 3 jy=1,my           
      do ikb=1,mb_nrk(nrank)
        i=itb_nrk(ikb,nrank)
        ib=ib_nrk(ikb,nrank) 
        do m=1,3 
        do js=mpsa,mpsa-3,-1     
        call interp1d3l(fst(i-2,js,jy,m),fst(i-1,js,jy,m),fst(i,js,jy,m),fst(i+1,js,jy,m), &
                        thst(i-2),thst(i-1),thst(i),thst(i+1),thb(ib),f1s(ikb,js))
        enddo
        js=mpsa
        call interp1d3l(f1s(ikb,js),f1s(ikb,js-1),f1s(ikb,js-2),f1s(ikb,js-3), &
                        ps(js),ps(js-1),ps(js-2),ps(js-3),psb(ib),fsxz(ikb,m))
        enddo
!        fxz(jbx(ib),jbz(ib))=(asm_b(ib)*x1s(ib,js-1)+bsm_b(ib)*x1s(ib,js-2)+csm_b(ib)*x1s(ib,js-3)) &
!                 /(asm_b(ib)+bsm_b(ib)+csm_b(ib))

        ef(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,1)=fsxz(ikb,1)*dcos(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank))) &
                                                      -fsxz(ikb,3)*dsin(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank)))
        ef(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,2)=fsxz(ikb,2)
        ef(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,3)=fsxz(ikb,1)*dsin(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank))) &
                                                      +fsxz(ikb,3)*dcos(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank)))

      enddo
    3 continue 
      endif !!mb_nrk(nrank).gt.0
           
!      if(ibnd==2) then
!      
!      do ika1=1,ma1_nrk(nrank)
!        i=ita1_nrk(ika1,nrank)
!        ia1=ia1_nrk(ika1,nrank)
!        do js=mpsa,mpsa-3,1
!        call interp1d3l(fst(i-2,js),fst(i-1,js),fst(i,js),fst(i+1,js), &
!                        thst(i-2),thst(i-1),thst(i),thst(i+1),tha1(ia1),x1s(ika1,js,jy))
!        enddo
!        js=mpsa
!        call interp1d3l(x1s(ika1,js,jy),x1s(ika1,js-1,jy),x1s(ika1,js-2,jy),x1s(ika1,js-3,jy), &
!                        ps(js),ps(js-1),ps(js-2),ps(js-3),psa1(ia1),fxz(jxa1_nrk(ika1,nrank),jza1_nrk(ika1,nrank)))
!!        fxz(jxa1(ia1),jza1(ia1))=(asm_a1(ia1)*x1s(ia1,js-1)+bsm_a1(ia1)*x1s(ia1,js-2)+csm_a1(ia1)*x1s(ia1,js-3)) &
!!                  /(asm_a1(ia1)+bsm_a1(ia1)+csm_a1(ia1)) 
!      enddo
!      endif
      return
      end

!ws****************************************************************
      subroutine cur_atlastgrid_r0p1_v1(kk)
      use declare
      implicit real*8 (b-h,o-z)
      integer kk
      dimension fst(n2th+5,mps4:mps,my,3),f1s(mbm_nrk,mps4:mps),fsxz(mbm_nrk,3),wst(n2th+5) !fxz(mx,mz,my,mm),
      include 'mpif.h'

!       integer status(mpi_status_size)

!
      do 2 js=mpsa-2,mpsa-nda,-1
      if(mts_nrk(js,nrank).gt.0) then
      do 21 m=1,3
      do 21 jy=1,my      
      call interp_xz2ps(cur(ix_first:ix_last,iz_first:iz_last,jy,m), &
                        xx(ix_first:ix_last),zz(iz_first:iz_last),ix_last-ix_first+1,iz_last-iz_first+1, &
                       fst(itsmin(js,nrank):itsmax(js,nrank),js,jy,m), &
                        xxs(itsmin(js,nrank):itsmax(js,nrank),js),zzs(itsmin(js,nrank):itsmax(js,nrank),js),mts_nrk(js,nrank))
      
   21 continue    
!!mpi
         do irecv=1,mrkb
         do isend=1,nsend(irecv,js)
         if(nrank.eq.nranksend(irecv,js,isend) .and. nrank.ne.nrkb(irecv)) then
         ltmin=ittransmin(irecv,js,isend)
         ltmax=ittransmax(irecv,js,isend)
         if (itbmin(nrkb(irecv))+n2th .le. itsmax(js,nranksend(irecv,js,isend))) then
         ltmin=ltmin+n2th
         ltmax=ltmax+n2th
         endif
         if (itbmax(nrkb(irecv))-n2th .ge. itsmin(js,nranksend(irecv,js,isend))) then
         ltmin=ltmin-n2th
         ltmax=ltmax-n2th
         endif
         do lt=ltmin,ltmax
         call mpi_send(fst(lt,js,1:my,:),my*3, mpi_double_precision, nrkb(irecv), isend,  &
		               mpi_comm_world,ierror )
         enddo
         endif
         enddo
         enddo
      endif

      if(mb_nrk(nrank).gt.0) then  
         do irecv=1,mrkb
         do isend=1,nsend(irecv,js)
         if(nrank.eq.nrkb(irecv) .and. nrank.ne.nranksend(irecv,js,isend)) then
         ltmin=ittransmin(irecv,js,isend)
         ltmax=ittransmax(irecv,js,isend)
         do lt=ltmin,ltmax
         call mpi_recv(fst(lt,js,1:my,:),my*3, mpi_double_precision, nranksend(irecv,js,isend), isend,  &
		               mpi_comm_world,status,ierror )
         enddo
         endif
         enddo
         enddo
       endif
!!mpi
    2 continue 
           
      if(mb_nrk(nrank).gt.0) then      

      do jy=1,my
      do js=mpsa-2,mpsa-4,-1
      do jt=itbmin(nrank),itbmax(nrank)
      cx1st=fst(jt,js,jy,1)
      cz1st=fst(jt,js,jy,3)

      fst(jt,js,jy,1)=cx1st*dcos(thst(jt))+cz1st*dsin(thst(jt))
      fst(jt,js,jy,3)=-cx1st*dsin(thst(jt))+cz1st*dcos(thst(jt))
       
!      vrst=fst(jt,js,jy,3)*dcos(thst(jt))+fst(jt,js,jy,5)*dsin(thst(jt))
!      vpst=-fst(jt,js,jy,3)*dsin(thst(jt))+fst(jt,js,jy,5)*dcos(thst(jt))
!      brst=fst(jt,js,jy,6)*dcos(thst(jt))+fst(jt,js,jy,8)*dsin(thst(jt))
!      bpst=-fst(jt,js,jy,6)*dsin(thst(jt))+fst(jt,js,jy,8)*dcos(thst(jt))
!      fst(jt,js,jy,3)=vrst
!      fst(jt,js,jy,5)=vpst
!      fst(jt,js,jy,6)=brst
!      fst(jt,js,jy,8)=bpst
      enddo
      enddo
      enddo    

      do 10 js=mpsa-2,mpsa-4,-1
!      call smth_ps_nrk(fst(itbmin(nrank):itbmax(nrank),js,:,:),itbmax(nrank)-itbmin(nrank)+1,3)
!ws_smps--------------------------
      do 11 k=1,kk
      do m=1,3
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,js,jy,m)/2.+((fst(jt+1,js,jy,m)+fst(jt-1,js,jy,m))*3./16+(fst(jt+2,js,jy,m)+fst(jt-2,js,jy,m))/16.&
             +(fst(jt,js,jy+1,m)+fst(jt,js,jy-1,m))*3./16+(fst(jt,js,jy+2,m)+fst(jt,js,jy-2,m))/16.)/2.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,js,jy,m)=wst(jt)
      enddo

      enddo
      enddo

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         call mpi_send(fst(ltmin:ltmax,js,1:my,:), 2*my*3, mpi_double_precision, nrkb1(inrkb(nrank)+1), 1,  &
		               mpi_comm_world,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,js,1:my,:), 2*my*3, mpi_double_precision, nrkb1(inrkb(nrank)-1), 1,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         call mpi_send(fst(ltmin:ltmax,js,1:my,:), 2*my*3, mpi_double_precision, nrkb1(1), 1,  &
		               mpi_comm_world,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,js,1:my,:), 2*my*3, mpi_double_precision, nrkb1(mrkb), 1,  &
		               mpi_comm_world,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         call mpi_send(fst(ltmin:ltmax,js,1:my,:), 2*my*3, mpi_double_precision, nrkb1(inrkb(nrank)-1), 2,  &
		               mpi_comm_world,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,js,1:my,:), 2*my*3, mpi_double_precision, nrkb1(inrkb(nrank)+1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         call mpi_send(fst(ltmin:ltmax,js,1:my,:), 2*my*3, mpi_double_precision, nrkb1(mrkb), 2,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,js,1:my,:), 2*my*3, mpi_double_precision, nrkb1(1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   11 continue
!ws_smps--------------------
   10 continue

      do 61 m=1,3

      if(m.eq.1) then !!(m=1) ibnd==0
      is=mpsa 
      do 51 jy=1,my  
      do jt=itbmin(nrank),itbmax(nrank)
       fst(jt,is,jy,m)=0.           
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   51 continue

      is=mpsa-1 
      do 52 jy=1,my  
      do jt=itbmin(nrank),itbmax(nrank)
        call interp1d3l(fst(jt,is-3,jy,m),fst(jt,is-2,jy,m),fst(jt,is-1,jy,m),fst(jt,is+1,jy,m), &
                        ps(is-3),ps(is-2),ps(is-1),ps(is+1),ps(is),fst(jt,is,jy,m))                   
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   52 continue

      do 53 k=1,kk
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,is,jy,m)/2.+((fst(jt+1,is,jy,m)+fst(jt-1,is,jy,m))*3./16+(fst(jt+2,is,jy,m)+fst(jt-2,is,jy,m))/16.&
             +(fst(jt,is,jy+1,m)+fst(jt,is,jy-1,m))*3./16+(fst(jt,is,jy+2,m)+fst(jt,is,jy-2,m))/16.)/2.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,is,jy,m)=wst(jt)
      enddo

      enddo

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)+1), 1,  &
		               mpi_comm_world,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)-1), 1,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(1), 1,  &
		               mpi_comm_world,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(mrkb), 1,  &
		               mpi_comm_world,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)-1), 2,  &
		               mpi_comm_world,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)+1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(mrkb), 2,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   53 continue

      else  !!(m=2,3) bnd==1  
      is=mpsa-1    
      do 31 jy=1,my      
      do jt=itbmin(nrank),itbmax(nrank)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))  
!       fst(jt,is)=fst(jt,is-1)                   
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   31 continue 
!      call smth_ps_nrk(fst(itbmin(nrank):itbmax(nrank),is,:,:),itbmax(nrank)-itbmin(nrank)+1,3)
      do 32 k=1,kk
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,is,jy,m)/2.+((fst(jt+1,is,jy,m)+fst(jt-1,is,jy,m))*3./16+(fst(jt+2,is,jy,m)+fst(jt-2,is,jy,m))/16.&
             +(fst(jt,is,jy+1,m)+fst(jt,is,jy-1,m))*3./16+(fst(jt,is,jy+2,m)+fst(jt,is,jy-2,m))/16.)/2.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,is,jy,m)=wst(jt)
      enddo

      enddo
      
         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)+1), 1,  &
		               mpi_comm_world,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)-1), 1,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(1), 1,  &
		               mpi_comm_world,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(mrkb), 1,  &
		               mpi_comm_world,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)-1), 2,  &
		               mpi_comm_world,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)+1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(mrkb), 2,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   32 continue
!ws_smps--------------------
      is=mpsa
 
      do 41 jy=1,my  
      do jt=itbmin(nrank),itbmax(nrank)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))
!       fst(jt,is)=fst(jt,is-1)                   
      enddo
!      call smth_ps(fst(itbmin(nrank):itbmax(nrank),is,jy,m),itbmax(nrank)-itbmin(nrank)+1,3)
   41 continue
!      call smth_ps_nrk(fst(itbmin(nrank):itbmax(nrank),is,:,:),itbmax(nrank)-itbmin(nrank)+1,3)
!ws_smps--------------------
      do 42 k=1,kk
      do jy=1,my
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      wst(jt)=fst(jt,is,jy,m)/2.+((fst(jt+1,is,jy,m)+fst(jt-1,is,jy,m))*3./16+(fst(jt+2,is,jy,m)+fst(jt-2,is,jy,m))/16.&
             +(fst(jt,is,jy+1,m)+fst(jt,is,jy-1,m))*3./16+(fst(jt,is,jy+2,m)+fst(jt,is,jy-2,m))/16.)/2.
      enddo
      
      do jt=itbmin(nrank)+2,itbmax(nrank)-2
      fst(jt,is,jy,m)=wst(jt)
      enddo

      enddo

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrank)+1))
         ltmax=itbmin(nrkb1(inrkb(nrank)+1))+1
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)+1), 1,  &
		               mpi_comm_world,ierror )
         endif
   
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)-1), 1,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(1), 1,  &
		               mpi_comm_world,ierror )
         endif
         
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmin(nrank)
         ltmax=itbmin(nrank)+1
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(mrkb), 1,  &
		               mpi_comm_world,status,ierror )
         endif

  !!ws  
         if(inrkb(nrank) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrank)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrank)-1))
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)-1), 2,  &
		               mpi_comm_world,ierror )
         endif

         if(inrkb(nrank) .lt. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(inrkb(nrank)+1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   
         if(inrkb(nrank) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         call mpi_send(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(mrkb), 2,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrank) .eq. mrkb) then
         ltmin=itbmax(nrank)-1
         ltmax=itbmax(nrank)
         call mpi_recv(fst(ltmin:ltmax,is,1:my,m), 2*my, mpi_double_precision, nrkb1(1), 2,  &
		               mpi_comm_world,status,ierror )
         endif
   42 continue
!ws_smps--------------------
      endif !!(m=2,3)ibnd==0
   61 continue

      do 3 jy=1,my           
      do ikb=1,mb_nrk(nrank)
        i=itb_nrk(ikb,nrank)
        ib=ib_nrk(ikb,nrank) 
        do m=1,3 
        do js=mpsa,mpsa-3,-1     
        call interp1d3l(fst(i-2,js,jy,m),fst(i-1,js,jy,m),fst(i,js,jy,m),fst(i+1,js,jy,m), &
                        thst(i-2),thst(i-1),thst(i),thst(i+1),thb(ib),f1s(ikb,js))
        enddo
        js=mpsa
        call interp1d3l(f1s(ikb,js),f1s(ikb,js-1),f1s(ikb,js-2),f1s(ikb,js-3), &
                        ps(js),ps(js-1),ps(js-2),ps(js-3),psb(ib),fsxz(ikb,m))
        enddo
!        fxz(jbx(ib),jbz(ib))=(asm_b(ib)*x1s(ib,js-1)+bsm_b(ib)*x1s(ib,js-2)+csm_b(ib)*x1s(ib,js-3)) &
!                 /(asm_b(ib)+bsm_b(ib)+csm_b(ib))

        cur(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,1)=fsxz(ikb,1)*dcos(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank))) &
                                                       -fsxz(ikb,3)*dsin(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank)))
        cur(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,2)=fsxz(ikb,2)
        cur(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank),jy,3)=fsxz(ikb,1)*dsin(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank))) &
                                                       +fsxz(ikb,3)*dcos(tpxz(jbx_nrk(ikb,nrank),jbz_nrk(ikb,nrank)))

      enddo
    3 continue 
      endif !!mb_nrk(nrank).gt.0
           
!      if(ibnd==2) then
!      
!      do ika1=1,ma1_nrk(nrank)
!        i=ita1_nrk(ika1,nrank)
!        ia1=ia1_nrk(ika1,nrank)
!        do js=mpsa,mpsa-3,1
!        call interp1d3l(fst(i-2,js),fst(i-1,js),fst(i,js),fst(i+1,js), &
!                        thst(i-2),thst(i-1),thst(i),thst(i+1),tha1(ia1),x1s(ika1,js,jy))
!        enddo
!        js=mpsa
!        call interp1d3l(x1s(ika1,js,jy),x1s(ika1,js-1,jy),x1s(ika1,js-2,jy),x1s(ika1,js-3,jy), &
!                        ps(js),ps(js-1),ps(js-2),ps(js-3),psa1(ia1),fxz(jxa1_nrk(ika1,nrank),jza1_nrk(ika1,nrank)))
!!        fxz(jxa1(ia1),jza1(ia1))=(asm_a1(ia1)*x1s(ia1,js-1)+bsm_a1(ia1)*x1s(ia1,js-2)+csm_a1(ia1)*x1s(ia1,js-3)) &
!!                  /(asm_a1(ia1)+bsm_a1(ia1)+csm_a1(ia1)) 
!      enddo
!      endif
      return
      end

!ws************************************************************
	subroutine smthp_traceline_v1(kk)
      use declare
      include 'mpif.h'
      real*8,dimension(3) :: bbb
      real*8,dimension(2) :: ppp,sl
      real*8,dimension(mx,mz,my) :: wsm
      real*8 sl0,xxp,zzp,yyp,dyyp,dy1p,dx1p,dz1p,rxp,rxm,rzp,rzm,ryp,rym
      real*8 volume,volppp,volmpp,volppm,volmpm,volpmp,volmmp,volpmm,volmmm
      integer kk,k,isgn,jjx,jjz,jjy,im,j_add
      
      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)
!      volume=dyy*dxx*dzz
!
!      dyyp=yy(jy+1)-yy(jy)
!      dyym=yy(jy-1)-yy(jy)

!      dy1p=dyyp/nn
!      dy1m=dyym/nn
      sl0=0
      do 3 k=1,kk

      do 2 jy=1,my
      do 2 jz=iz_first+2,iz_last-2
      do 2 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
      sl(:)=0

!      jjy=jy
!      jjx=jx
!      jjz=jz

      do 1 isgn=1,2

      dyyp=yy(jy+(-1)**isgn)-yy(jy)

      dy1p=dyyp/nline    
      
100   continue
      xxp=xx(jx)
      zzp=zz(jz)
      yyp=yy(jy)
      bbb(:)=x(jx,jz,jy,6:8)

      do nc=1,mcycline
      dx1p=xxp*dy1p*bbb(1)/bbb(2)
      dz1p=xxp*dy1p*bbb(3)/bbb(2)    			
	!艗脝脣茫脟掳艙酶脪禄虏艙潞贸碌脛脨脗脦禄脰脙
      xxp=xxp+dx1p
	  zzp=zzp+dz1p
      yyp=yyp+dy1p 
      if(yyp .gt. 2*pi) yyp=yyp-2*pi
      if(yyp .lt. 0) yyp=yyp+2*pi

      if( (xxp.gt.xx(jx+1)) .or. (xxp.lt.xx(jx-1)) &
      .or.(zzp.gt.zz(jz+1)) .or. (zzp.lt.zz(jz-1))) then
      if(nc.eq.1) then
      dy1p=dy1p*0.5
      goto 100
      else
      goto 200
      endif
      endif

      sl(isgn)=sl(isgn)+(-1)**isgn*sqrt(dx1p**2+dz1p**2+xxp**2*dy1p**2)
      
      jjx=jx+j_add(xxp,xx,mx,jx)
	  jjz=jz+j_add(zzp,zz,mz,jz)
	  jjy=jy+j_add(yyp,yy,my,jy)
      if(jjy.lt.1) jjy=jjy+my
      if(jjy.gt.my)jjy=jjy-my

      volume=(yy(jjy+1)-yy(jjy))*(xx(jjx+1)-xx(jjx))*(zz(jjz+1)-zz(jjz))

      rxp=xx(jjx+1)-xxp
      rxm=xxp-xx(jjx)
      rzp=zz(jjz+1)-zzp
      rzm=zzp-zz(jjz)
      ryp=yy(jjy+1)-yyp
      rym=yyp-yy(jjy)

      volppp=rxp*rzp*ryp
      volppm=rxp*rzp*rym
      volpmp=rxp*rzm*ryp
      volpmm=rxp*rzm*rym

      volmpp=rxm*rzp*ryp
      volmpm=rxm*rzp*rym
      volmmp=rxm*rzm*ryp
      volmmm=rxm*rzm*rym

      do im=6,8
      bbb(im-5)=(volppp*x(jjx,jjz,jjy,im)       +volmpp*x(jjx+1,jjz,jjy,im) &
	           +volppm*x(jjx,jjz,jjy+1,im)  +volmpm*x(jjx+1,jjz,jjy+1,im) &
	           +volpmp*x(jjx,jjz+1,jjy,im)     +volmmp*x(jjx+1,jjz+1,jjy,im)&
	           +volpmm*x(jjx,jjz+1,jjy+1,im)+volmmm*x(jjx+1,jjz+1,jjy+1,im))/volume
      enddo

      enddo

200   continue
      ppp(isgn)=(volppp*x1(jjx,jjz,jjy,2)       +volmpp*x1(jjx+1,jjz,jjy,2) &
	            +volppm*x1(jjx,jjz,jjy+1,2)  +volmpm*x1(jjx+1,jjz,jjy+1,2) &
	            +volpmp*x1(jjx,jjz+1,jjy,2)     +volmmp*x1(jjx+1,jjz+1,jjy,2)&
	            +volpmm*x1(jjx,jjz+1,jjy+1,2)+volmmm*x1(jjx+1,jjz+1,jjy+1,2))/volume         
                
1     continue

      wsm(jx,jz,jy)=csmp0*difc(ppp(1),x1(jx,jz,jy,2),ppp(2),sl(1),sl0,sl(2))
      endif
2     continue  


      x1(:,:,:,2)=x1(:,:,:,2)+wsm(:,:,:)

      call mpi_transfersm(x1(:,:,:,2),1)

3     continue 
      
      return
      end

!ws*****************************************************************
	subroutine smthp_traceline_v2(kk)
      use declare
      include 'mpif.h'
      real*8,dimension(2,mx,mz,my) :: sl,ppp,volppp,volmpp,volppm,volmpm,volpmp,volmmp,volpmm,volmmm
      real*8,dimension(mx,mz,my) :: wsm
      real*8,dimension(3) :: bbb
      integer,dimension(2,mx) :: jjx
      integer,dimension(2,mz) :: jjz
      integer,dimension(2,my) :: jjy
      real*8 sl0,xxp,zzp,yyp,dyyp,dy1p,dx1p,dz1p,rxp,rxm,rzp,rzm,ryp,rym,volume
      integer kk,isgn,sgn

      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)
!      volume=dyy*dxx*dzz
!
!      dyyp=yy(jy+1)-yy(jy)
!      dyym=yy(jy-1)-yy(jy)

!      dy1p=dyyp/nn
!      dy1m=dyym/nn
      sl0=0
!      volume=dyy*dxx*dzz

      do 2 jy=1,my
      do 2 jz=iz_first+2,iz_last-2
      do 2 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
      sl(:,jx,jz,jy)=0

      do 1 isgn=1,2

      dyyp=yy(jy+(-1)**isgn)-yy(jy)
!      
      dy1p=dyyp/nline 

100   continue

      xxp=xx(jx)
      zzp=zz(jz)
      yyp=yy(jy)
      
      bbb(:)=x(jx,jz,jy,6:8)
      
      do nc=1,mcycline
      dx1p=xxp*dy1p*bbb(1)/bbb(2)
      dz1p=xxp*dy1p*bbb(3)/bbb(2)

    			
	!艗脝脣茫脟掳艙酶脪禄虏艙潞贸碌脛脨脗脦禄脰脙
      xxp=xxp+dx1p
	  zzp=zzp+dz1p
      yyp=yyp+dy1p 
      if(yyp.gt. 2*pi) yyp=yyp-2*pi
      if(yyp.lt. 0) yyp=yyp+2*pi

      if( (xxp.gt.xx(jx+1)) .or. (xxp.lt.xx(jx-1)) &
      .or.(zzp.gt.zz(jz+1)) .or. (zzp.lt.zz(jz-1))) then
      if(nc.eq.1) then
      dy1p=dy1p*0.5
      goto 100
      else
      goto 200
      endif
      endif

      sl(isgn,jx,jz,jy)=sl(isgn,jx,jz,jy)+(-1)**isgn*sqrt(dx1p**2+dz1p**2+xxp**2*dy1p**2)
      
      jjx(isgn,jx)=jx+floor((xxp-xx(jx))/abs(xx(jx)-xx(jx+sgn(xxp-xx(jx))) ))
	  jjz(isgn,jz)=jz+floor((zzp-zz(jz))/abs(zz(jz)-zz(jz+sgn(zzp-zz(jz))) ))   
	  jjy(isgn,jy)=jy+floor((yyp-yy(jy))/abs(yy(jy)-yy(jy+sgn(yyp-zz(jy))) ))
      if(jjy(isgn,jy).lt.1) jjy(isgn,jy)=jjy(isgn,jy)+my
      if(jjy(isgn,jy).gt.my)jjy(isgn,jy)=jjy(isgn,jy)-my

      rxp=xx(jjx(isgn,jx)+1)-xxp
      rxm=xxp-xx(jjx(isgn,jx))
      rzp=zz(jjz(isgn,jz)+1)-zzp
      rzm=zzp-zz(jjz(isgn,jz))
      ryp=yy(jjy(isgn,jy)+1)-yyp
      rym=yyp-yy(jjy(isgn,jy))
      volume=(xx(jjx(isgn,jx)+1)-xx(jjx(isgn,jx)))*(zz(jjz(isgn,jz)+1)-zz(jjz(isgn,jz)))*(yy(jjy(isgn,jy)+1)-yy(jjy(isgn,jy)))

      volppp(isgn,jx,jz,jy)=rxp*rzp*ryp/volume
      volppm(isgn,jx,jz,jy)=rxp*rzp*rym/volume
      volpmp(isgn,jx,jz,jy)=rxp*rzm*ryp/volume
      volpmm(isgn,jx,jz,jy)=rxp*rzm*rym/volume

      volmpp(isgn,jx,jz,jy)=rxm*rzp*ryp/volume
      volmpm(isgn,jx,jz,jy)=rxm*rzp*rym/volume
      volmmp(isgn,jx,jz,jy)=rxm*rzm*ryp/volume
      volmmm(isgn,jx,jz,jy)=rxm*rzm*rym/volume

      do m=6,8
      bbb(m-5)=(volppp(isgn,jx,jz,jy)*x(jjx(isgn,jx),jjz(isgn,jz),jjy(isgn,jy),m) &
               +volmpp(isgn,jx,jz,jy)*x(jjx(isgn,jx)+1,jjz(isgn,jz),jjy(isgn,jy),m) &
	           +volppm(isgn,jx,jz,jy)*x(jjx(isgn,jx),jjz(isgn,jz),jyp(jjy(isgn,jy)),m) &
               +volmpm(isgn,jx,jz,jy)*x(jjx(isgn,jx)+1,jjz(isgn,jz),jyp(jjy(isgn,jy)),m) &
	           +volpmp(isgn,jx,jz,jy)*x(jjx(isgn,jx),jjz(isgn,jz)+1,jjy(isgn,jy),m) &
               +volmmp(isgn,jx,jz,jy)*x(jjx(isgn,jx)+1,jjz(isgn,jz)+1,jjy(isgn,jy),m) &
	           +volpmm(isgn,jx,jz,jy)*x(jjx(isgn,jx),jjz(isgn,jz)+1,jyp(jjy(isgn,jy)),m) &
               +volmmm(isgn,jx,jz,jy)*x(jjx(isgn,jx)+1,jjz(isgn,jz)+1,jyp(jjy(isgn,jy)),m))

      enddo

      enddo

200   continue
                    
                
1     continue


      endif
2     continue  

      do 3 k=1,kk
      
      do 21 jy=1,my
      do 21 jz=iz_first+2,iz_last-2
      do 21 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
      ppp(isgn,jx,jz,jy)=(volppp(isgn,jx,jz,jy)*x(jjx(isgn,jx),jjz(isgn,jz),jjy(isgn,jy),2) &
               +volmpp(isgn,jx,jz,jy)*x(jjx(isgn,jx)+1,jjz(isgn,jz),jjy(isgn,jy),2) &
	           +volppm(isgn,jx,jz,jy)*x(jjx(isgn,jx),jjz(isgn,jz),jyp(jjy(isgn,jy)),2) &
               +volmpm(isgn,jx,jz,jy)*x(jjx(isgn,jx)+1,jjz(isgn,jz),jyp(jjy(isgn,jy)),2) &
	           +volpmp(isgn,jx,jz,jy)*x(jjx(isgn,jx),jjz(isgn,jz)+1,jjy(isgn,jy),2) &
               +volmmp(isgn,jx,jz,jy)*x(jjx(isgn,jx)+1,jjz(isgn,jz)+1,jjy(isgn,jy),2) &
	           +volpmm(isgn,jx,jz,jy)*x(jjx(isgn,jx),jjz(isgn,jz)+1,jyp(jjy(isgn,jy)),2) &
               +volmmm(isgn,jx,jz,jy)*x(jjx(isgn,jx)+1,jjz(isgn,jz)+1,jyp(jjy(isgn,jy)),2))

      wsm(jx,jz,jy)=csmp0*difc(ppp(1,jx,jz,jy),x1(jx,jz,jy,2),ppp(2,jx,jz,jy),sl(1,jx,jz,jy),sl0,sl(2,jx,jz,jy))
      endif
21    continue 

      x1(:,:,:,2)=x1(:,:,:,2)+wsm(:,:,:)

      call mpi_transfersm(x1(:,:,:,2),1)

3     continue 
      
      return
      end

!ws*****************************************************************
	subroutine tracelineat(jxl,jzl,jyl,lx,lz,ly,sl,vol)
      use declare
      include 'mpif.h'
      integer,dimension(2) :: lx,lz,ly
      real*8,dimension(2) :: sl
      real*8,dimension(8,2) :: vol
      real*8,dimension(3) :: bbb
      real*8 sl0,xxp,zzp,yyp,dyyp,dy1p,dx1p,dz1p,rxp,rxm,rzp,rzm,ryp,rym
      real*8 volume,volppp,volmpp,volppm,volmpm,volpmp,volmmp,volpmm,volmmm
      integer isgn,jjx,jjz,jjy,jxl,jzl,jyl,j_add
       
      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)
!      volume=dyy*dxx*dzz
!
!      dyyp=yy(jy+1)-yy(jy)
!      dyym=yy(jy-1)-yy(jy)

!      dy1p=dyyp/nn
!      dy1m=dyym/nn
      sl0=0
      sl(:)=0.

      do 1 isgn=1,2


      dyyp=yy(jyl+(-1)**isgn)-yy(jyl)
      dy1p=dyyp/nline    

100   continue   
      xxp=xx(jxl)
      zzp=zz(jzl)
      yyp=yy(jyl)
      
      bbb(:)=x(jxl,jzl,jyl,6:8)

      do nc=1,mcycline
      dx1p=xxp*dy1p*bbb(1)/bbb(2)
      dz1p=xxp*dy1p*bbb(3)/bbb(2)
    			
	!艗脝脣茫脟掳艙酶脪禄虏艙潞贸碌脛脨脗脦禄脰脙
      xxp=xxp+dx1p
	  zzp=zzp+dz1p
      yyp=yyp+dy1p 
      if(yyp.gt. 2*pi) yyp=yyp-2*pi
      if(yyp.lt. 0) yyp=yyp+2*pi

!      if( (xxp.gt.xx(jxl+1)) .or. (xxp.lt.xx(jxl-1)) &
!      .or.(zzp.gt.zz(jzl+1)) .or. (zzp.lt.zz(jzl-1))) then
!      if(nc.eq.1) then
!      dy1p=dy1p*0.5
!      goto 100
!      else
!      goto 200
!      endif
!      endif

      sl(isgn)=sl(isgn)+(-1)**isgn*sqrt(dx1p**2+dz1p**2+xxp**2*dy1p**2)

      jjx=floor((xxp-xx(ix_first))/dxx)+ix_first
	  jjz=floor((zzp-zz(iz_first))/dzz)+iz_first   
	  jjy=floor((yyp-yy(iy_first))/dyy)+iy_first      
!      jjx=jxl+j_add(xxp,xx,mx,jxl)
!	  jjz=jzl+j_add(zzp,zz,mz,jzl)
!	  jjy=jyl+j_add(yyp,yy,my,jyl)
      if(jjy.lt.1) jjy=jjy+my
      if(jjy.gt.my)jjy=jjy-my
!      write(*,*) nrank,'nc',nc,jjy,xxp,jjx,zzp,jjz

      rxp=xx(jjx+1)-xxp
      rxm=xxp-xx(jjx)
      rzp=zz(jjz+1)-zzp
      rzm=zzp-zz(jjz)
      ryp=yy(jjy+1)-yyp
      rym=yyp-yy(jjy)

      volume=(xx(jjx+1)-xx(jjx))*(zz(jjz+1)-zz(jjz))*(yy(jjy+1)-yy(jjy))

      volppp=rxp*rzp*ryp/volume
      volppm=rxp*rzp*rym/volume
      volpmp=rxp*rzm*ryp/volume
      volpmm=rxp*rzm*rym/volume

      volmpp=rxm*rzp*ryp/volume
      volmpm=rxm*rzp*rym/volume
      volmmp=rxm*rzm*ryp/volume
      volmmm=rxm*rzm*rym/volume

      do m=6,8
      bbb(m-5)=(volppp*x(jjx,jjz,jjy,m)       +volmpp*x(jjx+1,jjz,jjy,m) &
	           +volppm*x(jjx,jjz,jjy+1,m)  +volmpm*x(jjx+1,jjz,jjy+1,m) &
	           +volpmp*x(jjx,jjz+1,jjy,m)     +volmmp*x(jjx+1,jjz+1,jjy,m)&
	           +volpmm*x(jjx,jjz+1,jjy+1,m)+volmmm*x(jjx+1,jjz+1,jjy+1,m))/volume
      enddo

      enddo

200   continue

      vol(1,isgn)=volppp
      vol(2,isgn)=volmpp
      vol(3,isgn)=volppm
      vol(4,isgn)=volmpm
      vol(5,isgn)=volpmp
      vol(6,isgn)=volmmp
      vol(7,isgn)=volpmm
      vol(8,isgn)=volmmm
      
      lx(isgn)=jjx
      lz(isgn)=jjz
      ly(isgn)=jjy     
                
1     continue

      return
      end

!ws********************************************
      subroutine smthp1_traceline(kk)
      use declare
      include 'mpif.h'
      integer,dimension(2) :: lx,lz,ly
      real*8,dimension(2) :: sl,ppp
      real*8,dimension(8,2) :: vol
      integer,dimension(2,mx,mz,my) :: llx,llz,lly
      real*8,dimension(2,mx,mz,my) :: ssl
      real*8,dimension(8,2,mx,mz,my) :: voll
      real*8,dimension(mx,mz,my) :: wsm
      integer kk,k,isgn,jjx,jjz,jjy
      real*8 ssl0

      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

      ssl0=0.

      do 2 jy=1,my
      do 2 jz=iz_first+2,iz_last-2
      do 2 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
      call tracelineat(jx,jz,jy,lx,lz,ly,sl,vol)

      llx(:,jx,jz,jy)=lx(:)
      llz(:,jx,jz,jy)=lz(:)
      lly(:,jx,jz,jy)=ly(:)
      ssl(:,jx,jz,jy)=sl(:)
      voll(:,:,jx,jz,jy)=vol(:,:)
      endif
2     continue  

      do 3 k=1,kk
      
      do 21 jy=1,my
      do 21 jz=iz_first+2,iz_last-2
      do 21 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
      do isgn=1,2
      jjx=llx(isgn,jx,jz,jy)
      jjz=llz(isgn,jx,jz,jy)
      jjy=lly(isgn,jx,jz,jy)

      ppp(isgn)=voll(1,isgn,jx,jz,jy)*x1(jjx,jjz,jjy,2)  &
               +voll(2,isgn,jx,jz,jy)*x1(jjx+1,jjz,jjy,2) &
	           +voll(3,isgn,jx,jz,jy)*x1(jjx,jjz,jjy+1,2) &
               +voll(4,isgn,jx,jz,jy)*x1(jjx+1,jjz,jjy+1,2) &
	           +voll(5,isgn,jx,jz,jy)*x1(jjx,jjz+1,jjy,2) &
               +voll(6,isgn,jx,jz,jy)*x1(jjx+1,jjz+1,jjy,2) &
	           +voll(7,isgn,jx,jz,jy)*x1(jjx,jjz+1,jjy+1,2) &
               +voll(8,isgn,jx,jz,jy)*x1(jjx+1,jjz+1,jjy+1,2)
      enddo

      wsm(jx,jz,jy)=csmp0*difc(ppp(1),x1(jx,jz,jy,2),ppp(2),ssl(1,jx,jz,jy),ssl0,ssl(2,jx,jz,jy))
      endif
21    continue 
      x1(:,:,:,2)=x1(:,:,:,2)+wsm(:,:,:)

      call mpi_transfersm(x1(:,:,:,2),1)
!      do jy=1,my
!      x(:,:,jy,2)=x1(:,:,jy,2)+xint(:,:,2)
!      enddo
3     continue 
      
      return
      end

!ws********************************************
      subroutine smthp_traceline(kk)
      use declare
      include 'mpif.h'
      integer,dimension(2) :: lx,lz,ly
      real*8,dimension(2) :: sl,ppp
      real*8,dimension(8,2) :: vol
      integer,dimension(2,mx,mz,my) :: llx,llz,lly
      real*8,dimension(2,mx,mz,my) :: ssl
      real*8,dimension(8,2,mx,mz,my) :: voll
      real*8,dimension(mx,mz,my) :: wsm
      integer kk,k,isgn,jjx,jjz,jjy
      real*8 ssl0,csmp

      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

      ssl0=0.

      do 2 jy=iy_first+2,iy_last-2
      do 2 jz=iz_first+2,iz_last-2
      do 2 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
!      if(psi(jx,jz).lt.psia1) then
      call tracelineat2(jx,jz,jy,lx,lz,ly,sl,vol)

      llx(:,jx,jz,jy)=lx(:)
      llz(:,jx,jz,jy)=lz(:)
      lly(:,jx,jz,jy)=ly(:)
      ssl(:,jx,jz,jy)=sl(:)
      voll(:,:,jx,jz,jy)=vol(:,:)
      endif
2     continue  

      do 3 k=1,kk
      
      do 21 jy=iy_first+2,iy_last-2
      do 21 jz=iz_first+2,iz_last-2
      do 21 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
!      if(psi(jx,jz).lt.psia1) then
      do isgn=1,2
      jjx=llx(isgn,jx,jz,jy)
      jjz=llz(isgn,jx,jz,jy)
      jjy=lly(isgn,jx,jz,jy)

      ppp(isgn)=voll(1,isgn,jx,jz,jy)*x(jjx,jjz,jjy,2)  &
               +voll(2,isgn,jx,jz,jy)*x(jjx+1,jjz,jjy,2) &
	       +voll(3,isgn,jx,jz,jy)*x(jjx,jjz,jjy+1,2) &
               +voll(4,isgn,jx,jz,jy)*x(jjx+1,jjz,jjy+1,2) &
	       +voll(5,isgn,jx,jz,jy)*x(jjx,jjz+1,jjy,2) &
               +voll(6,isgn,jx,jz,jy)*x(jjx+1,jjz+1,jjy,2) &
	       +voll(7,isgn,jx,jz,jy)*x(jjx,jjz+1,jjy+1,2) &
               +voll(8,isgn,jx,jz,jy)*x(jjx+1,jjz+1,jjy+1,2)
      enddo
      csmp=csmp0all*(1-tanh((rr(jx,jz)-0.9)*40))/2.
!      csmp=csmp0all*(1-tanh((psi(jx,jz)-pstrans)/pstransw))/2.
      wsm(jx,jz,jy)=csmp*difc(ppp(1),x(jx,jz,jy,2),ppp(2),ssl(1,jx,jz,jy),ssl0,ssl(2,jx,jz,jy))
      endif
21    continue 
      x(:,:,:,2)=x(:,:,:,2)+wsm(:,:,:)

      call mpi_transfersm(x(:,:,:,2),1)

3     continue 
      
      return
      end

!ws*****************************************************************
	subroutine tracelineat2(jxl,jzl,jyl,lx,lz,ly,sl,vol)
      use declare
      include 'mpif.h'
      integer,dimension(2) :: lx,lz,ly
      real*8,dimension(2) :: sl
      real*8,dimension(8,2) :: vol
      real*8,dimension(3) :: bbb
      real*8 sl0,xxp,zzp,yyp,dyyp,dy1p,dx1p,dz1p,rxp,rxm,rzp,rzm,ryp,rym
      real*8 xxp0,zzp0,yyp0,xxp12,zzp12,yyp12
      real*8 volume,volppp,volmpp,volppm,volmpm,volpmp,volmmp,volpmm,volmmm
      integer isgn,jjx,jjz,jjy,jxl,jzl,jyl,j_add
       
      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)
!      volume=dyy*dxx*dzz
!
!      dyyp=yy(jy+1)-yy(jy)
!      dyym=yy(jy-1)-yy(jy)

!      dy1p=dyyp/nn
!      dy1m=dyym/nn
      sl0=0
      sl(:)=0.

      do 1 isgn=1,2


      dyyp=yy(jyl+(-1)**isgn)-yy(jyl)
      dy1p=dyyp/nline    

100   continue   
      xxp0=xx(jxl)
      zzp0=zz(jzl)
      yyp0=yy(jyl)

      bbb(:)=x(jxl,jzl,jyl,6:8)

      do nc=1,mcycline
!      dy1p=dyyp/nline/xxp0*xzero
      dx1p=xxp0*dy1p*bbb(1)/bbb(2)
      dz1p=xxp0*dy1p*bbb(3)/bbb(2)
    			
	!艗脝脣茫脟掳艙酶脪禄虏艙潞贸碌脛脨脗脦禄脰脙
      xxp=xxp0+dx1p
	  zzp=zzp0+dz1p
      yyp=yyp0+dy1p 
      if(yyp.gt. 2*pi) yyp=yyp-2*pi
      if(yyp.lt. 0) yyp=yyp+2*pi

!      if( (xxp.gt.xx(jxl+1)) .or. (xxp.lt.xx(jxl-1)) &
!      .or.(zzp.gt.zz(jzl+1)) .or. (zzp.lt.zz(jzl-1))) then
!      if(nc.eq.1) then
!      dy1p=dy1p*0.5
!      goto 100
!      else
!      goto 200
!      endif
!      endif

!      sl(isgn)=sl(isgn)+(-1)**isgn*sqrt(dx1p**2+dz1p**2+xxp**2*dy1p**2)
      jjx=floor((xxp-xx(ix_first))/dxx)+ix_first
	  jjz=floor((zzp-zz(iz_first))/dzz)+iz_first   
	  jjy=floor((yyp-yy(iy_first))/dyy)+iy_first      
!      jjx=jxl+j_add(xxp,xx,mx,jxl)
!	  jjz=jzl+j_add(zzp,zz,mz,jzl)
!	  jjy=jyl+j_add(yyp,yy,my,jyl)
      if(jjy.lt.1) jjy=jjy+myt
      if(jjy.gt.myt)jjy=jjy-myt
!      write(*,*) nrank,'nc',nc,jjy,xxp,jjx,zzp,jjz

      rxp=xx(jjx+1)-xxp
      rxm=xxp-xx(jjx)
      rzp=zz(jjz+1)-zzp
      rzm=zzp-zz(jjz)
      ryp=yy(jjy+1)-yyp
      rym=yyp-yy(jjy)

      volume=(xx(jjx+1)-xx(jjx))*(zz(jjz+1)-zz(jjz))*(yy(jjy+1)-yy(jjy))

      volppp=rxp*rzp*ryp/volume
      volppm=rxp*rzp*rym/volume
      volpmp=rxp*rzm*ryp/volume
      volpmm=rxp*rzm*rym/volume

      volmpp=rxm*rzp*ryp/volume
      volmpm=rxm*rzp*rym/volume
      volmmp=rxm*rzm*ryp/volume
      volmmm=rxm*rzm*rym/volume

      do m=6,8
      bbb(m-5)=(volppp*x(jjx,jjz,jjy,m)       +volmpp*x(jjx+1,jjz,jjy,m) &
	           +volppm*x(jjx,jjz,jjy+1,m)  +volmpm*x(jjx+1,jjz,jjy+1,m) &
	           +volpmp*x(jjx,jjz+1,jjy,m)     +volmmp*x(jjx+1,jjz+1,jjy,m)&
	           +volpmm*x(jjx,jjz+1,jjy+1,m)+volmmm*x(jjx+1,jjz+1,jjy+1,m))/volume
      enddo


      xxp12=xxp0+0.5*dx1p
	  zzp12=zzp0+0.5*dz1p
      yyp12=yyp0+0.5*dy1p

!      dy1p=dyyp/nline/xxp*xzero
      dx1p=xxp*dy1p*bbb(1)/bbb(2)
      dz1p=xxp*dy1p*bbb(3)/bbb(2)

      xxp=xxp12+0.5*dx1p
      zzp=zzp12+0.5*dz1p
      yyp=yyp12+0.5*dy1p
      if(yyp.gt. 2*pi) yyp=yyp-2*pi
      if(yyp.lt. 0) yyp=yyp+2*pi

      sl(isgn)=sl(isgn)+(-1)**isgn*sqrt(dx1p**2+dz1p**2+xxp**2*dy1p**2)
      
!      jjx=jxl+j_add(xxp,xx,mx,jxl)
!	  jjz=jzl+j_add(zzp,zz,mz,jzl)
!	  jjy=jyl+j_add(yyp,yy,my,jyl)

      jjx=floor((xxp-xx(ix_first))/dxx)+ix_first
	  jjz=floor((zzp-zz(iz_first))/dzz)+iz_first   
	  jjy=floor((yyp-yy(iy_first))/dyy)+iy_first
      if(jjy.lt.1) jjy=jjy+myt
      if(jjy.gt.myt)jjy=jjy-myt
!      write(*,*) nrank,'nc',nc,jjy,xxp,jjx,zzp,jjz

      rxp=xx(jjx+1)-xxp
      rxm=xxp-xx(jjx)
      rzp=zz(jjz+1)-zzp
      rzm=zzp-zz(jjz)
      ryp=yy(jjy+1)-yyp
      rym=yyp-yy(jjy)

      volume=(xx(jjx+1)-xx(jjx))*(zz(jjz+1)-zz(jjz))*(yy(jjy+1)-yy(jjy))

      volppp=rxp*rzp*ryp/volume
      volppm=rxp*rzp*rym/volume
      volpmp=rxp*rzm*ryp/volume
      volpmm=rxp*rzm*rym/volume

      volmpp=rxm*rzp*ryp/volume
      volmpm=rxm*rzp*rym/volume
      volmmp=rxm*rzm*ryp/volume
      volmmm=rxm*rzm*rym/volume

      do m=6,8
      bbb(m-5)=(volppp*x(jjx,jjz,jjy,m)       +volmpp*x(jjx+1,jjz,jjy,m) &
	           +volppm*x(jjx,jjz,jjy+1,m)  +volmpm*x(jjx+1,jjz,jjy+1,m) &
	           +volpmp*x(jjx,jjz+1,jjy,m)     +volmmp*x(jjx+1,jjz+1,jjy,m)&
	           +volpmm*x(jjx,jjz+1,jjy+1,m)+volmmm*x(jjx+1,jjz+1,jjy+1,m))/volume
      enddo

      xxp0=xxp
      zzp0=zzp
      yyp0=yyp

      enddo

200   continue

      vol(1,isgn)=volppp
      vol(2,isgn)=volmpp
      vol(3,isgn)=volppm
      vol(4,isgn)=volmpm
      vol(5,isgn)=volpmp
      vol(6,isgn)=volmmp
      vol(7,isgn)=volpmm
      vol(8,isgn)=volmmm
      
      lx(isgn)=jjx
      lz(isgn)=jjz
      ly(isgn)=jjy     
                
1     continue

      return
      end

!ws*****************************************************************
	subroutine tracelineat2_dx(jxl,jzl,jyl,lx,lz,ly,sl,vol)
!      use declare
!      include 'mpif.h'
!      integer,dimension(2) :: lx,lz,ly
!      real*8,dimension(2) :: sl
!      real*8,dimension(8,2) :: vol
!      real*8,dimension(3) :: bbb
!      real*8 sl0,xxp,zzp,yyp,dyyp,dy1p,dx1p,dz1p,rxp,rxm,rzp,rzm,ryp,rym
!      real*8 xxp0,zzp0,yyp0,xxp12,zzp12,yyp12
!      real*8 volume,volppp,volmpp,volppm,volmpm,volpmp,volmmp,volpmm,volmmm
!      integer isgn,jjx,jjz,jjy,jxl,jzl,jyl
!       
!      difc(fm1,f0,fp1,xm1,x0,xp1)= &
!       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)
!!      volume=dyy*dxx*dzz
!!
!!      dyyp=yy(jy+1)-yy(jy)
!!      dyym=yy(jy-1)-yy(jy)
!
!!      dy1p=dyyp/nn
!!      dy1m=dyym/nn
!      sl0=0
!      sl(:)=0.
!
!      do 1 isgn=1,2
!
!
!!      dyyp=yy(jyl+(-1)**isgn)-yy(jyl)
!!      dy1p=dyyp/nline
!      
!      dx1p=xx(jxl+(-1)**isgn)-xx(jxl)    
!
!100   continue   
!      xxp0=xx(jxl)
!      zzp0=zz(jzl)
!      yyp0=yy(jyl)
!
!      bbb(:)=x(jxl,jzl,jyl,6:8)
!
!      do nc=1,mcycline
!!      dy1p=dyyp/nline/xxp0*xzero
!!      dx1p=xxp0*dy1p*bbb(1)/bbb(2)
!!      dz1p=xxp0*dy1p*bbb(3)/bbb(2)
!       if(abs(bbb(1)) .lt. 1.e-7)
!       dz1p=dx1p*bbb(3)/bbb(1)
!       dy1p=dx1p*bbb(2)/bbb(1)/xxp0		
!	!艗脝脣茫脟掳艙酶脪禄虏艙潞贸碌脛脨脗脦禄脰脙
!      xxp=xxp0+dx1p
!	  zzp=zzp0+dz1p
!      yyp=yyp0+dy1p 
!      if(yyp.gt. 2*pi) yyp=yyp-2*pi
!      if(yyp.lt. 0) yyp=yyp+2*pi
!
!!      if( (xxp.gt.xx(jxl+1)) .or. (xxp.lt.xx(jxl-1)) &
!!      .or.(zzp.gt.zz(jzl+1)) .or. (zzp.lt.zz(jzl-1))) then
!!      if(nc.eq.1) then
!!      dy1p=dy1p*0.5
!!      goto 100
!!      else
!!      goto 200
!!      endif
!!      endif
!
!!      sl(isgn)=sl(isgn)+(-1)**isgn*sqrt(dx1p**2+dz1p**2+xxp**2*dy1p**2)
!      
!      jjx=floor((xxp-xx(ix_first))/dxx)+ix_first
!	  jjz=floor((zzp-zz(iz_first))/dzz)+iz_first   
!	  jjy=floor((yyp-yy(1))/dyy)+1
!      if(jjy.lt.1) jjy=jjy+my
!      if(jjy.gt.my)jjy=jjy-my
!!      write(*,*) nrank,'nc',nc,jjy,xxp,jjx,zzp,jjz
!
!      rxp=xx(jjx+1)-xxp
!      rxm=xxp-xx(jjx)
!      rzp=zz(jjz+1)-zzp
!      rzm=zzp-zz(jjz)
!      ryp=yy(jjy+1)-yyp
!      rym=yyp-yy(jjy)
!
!      volume=(xx(jjx+1)-xx(jjx))*(zz(jjz+1)-zz(jjz))*(yy(jjy+1)-yy(jjy))
!
!      volppp=rxp*rzp*ryp/volume
!      volppm=rxp*rzp*rym/volume
!      volpmp=rxp*rzm*ryp/volume
!      volpmm=rxp*rzm*rym/volume
!
!      volmpp=rxm*rzp*ryp/volume
!      volmpm=rxm*rzp*rym/volume
!      volmmp=rxm*rzm*ryp/volume
!      volmmm=rxm*rzm*rym/volume
!
!      do m=6,8
!      bbb(m-5)=(volppp*x(jjx,jjz,jjy,m)       +volmpp*x(jjx+1,jjz,jjy,m) &
!	           +volppm*x(jjx,jjz,jjy+1,m)  +volmpm*x(jjx+1,jjz,jjy+1,m) &
!	           +volpmp*x(jjx,jjz+1,jjy,m)     +volmmp*x(jjx+1,jjz+1,jjy,m)&
!	           +volpmm*x(jjx,jjz+1,jjy+1,m)+volmmm*x(jjx+1,jjz+1,jjy+1,m))/volume
!      enddo
!
!
!      xxp12=xxp0+0.5*dx1p
!	  zzp12=zzp0+0.5*dz1p
!      yyp12=yyp0+0.5*dy1p
!
!!      dy1p=dyyp/nline/xxp*xzero
!!      dx1p=xxp*dy1p*bbb(1)/bbb(2)
!!      dz1p=xxp*dy1p*bbb(3)/bbb(2)
!      dz1p=dx1p*bbb(3)/bbb(1)
!      dy1p=dx1p*bbb(2)/bbb(1)/xxp	
!      
!      xxp=xxp12+0.5*dx1p
!      zzp=zzp12+0.5*dz1p
!      yyp=yyp12+0.5*dy1p
!      if(yyp.gt. 2*pi) yyp=yyp-2*pi
!      if(yyp.lt. 0) yyp=yyp+2*pi
!
!      sl(isgn)=sl(isgn)+(-1)**isgn*sqrt(dx1p**2+dz1p**2+xxp**2*dy1p**2)
!      
!      jjx=floor((xxp-xx(ix_first))/dxx)+ix_first
!	  jjz=floor((zzp-zz(iz_first))/dzz)+iz_first   
!	  jjy=floor((yyp-yy(1))/dyy)+1
!      if(jjy.lt.1) jjy=jjy+my
!      if(jjy.gt.my)jjy=jjy-my
!!      write(*,*) nrank,'nc',nc,jjy,xxp,jjx,zzp,jjz
!
!      rxp=xx(jjx+1)-xxp
!      rxm=xxp-xx(jjx)
!      rzp=zz(jjz+1)-zzp
!      rzm=zzp-zz(jjz)
!      ryp=yy(jjy+1)-yyp
!      rym=yyp-yy(jjy)
!
!      volume=(xx(jjx+1)-xx(jjx))*(zz(jjz+1)-zz(jjz))*(yy(jjy+1)-yy(jjy))
!
!      volppp=rxp*rzp*ryp/volume
!      volppm=rxp*rzp*rym/volume
!      volpmp=rxp*rzm*ryp/volume
!      volpmm=rxp*rzm*rym/volume
!
!      volmpp=rxm*rzp*ryp/volume
!      volmpm=rxm*rzp*rym/volume
!      volmmp=rxm*rzm*ryp/volume
!      volmmm=rxm*rzm*rym/volume
!
!      do m=6,8
!      bbb(m-5)=(volppp*x(jjx,jjz,jjy,m)       +volmpp*x(jjx+1,jjz,jjy,m) &
!	           +volppm*x(jjx,jjz,jjy+1,m)  +volmpm*x(jjx+1,jjz,jjy+1,m) &
!	           +volpmp*x(jjx,jjz+1,jjy,m)     +volmmp*x(jjx+1,jjz+1,jjy,m)&
!	           +volpmm*x(jjx,jjz+1,jjy+1,m)+volmmm*x(jjx+1,jjz+1,jjy+1,m))/volume
!      enddo
!
!      xxp0=xxp
!      zzp0=zzp
!      yyp0=yyp
!
!      enddo
!
!200   continue
!
!      vol(1,isgn)=volppp
!      vol(2,isgn)=volmpp
!      vol(3,isgn)=volppm
!      vol(4,isgn)=volmpm
!      vol(5,isgn)=volpmp
!      vol(6,isgn)=volmmp
!      vol(7,isgn)=volpmm
!      vol(8,isgn)=volmmm
!      
!      lx(isgn)=jjx
!      lz(isgn)=jjz
!      ly(isgn)=jjy     
!                
!1     continue
!
!      return
      end

!ws********************************************
      subroutine smthp2_traceline(kk)
      use declare
      include 'mpif.h'
      integer,dimension(2) :: lx,lz,ly
      real*8,dimension(2) :: sl,ppp,ppb
      real*8,dimension(8,2) :: vol
      integer,dimension(2,mx,mz,my) :: llx,llz,lly
      real*8,dimension(2,mx,mz,my) :: ssl
      real*8,dimension(8,2,mx,mz,my) :: voll
      real*8,dimension(mx,mz,my) :: wsm,pb !,wsm1,wsm2
      integer kk,k,isgn,jjx,jjz,jjy
      real*8 ssl0,bb

      !  dif1= (d f / dx)*dx2 = [((x0-xm1)/(xp1-x0)*(fp1-f0)-(xp1-x0)/(x0-xm1)*(f0-fm1))/(xp1-xm1)] *(xp1-x0)*(x0-xm1) 
      dif1(fm1,f0,fp1,xm1,x0,xp1)= &
       ((x0-xm1)**2*(fp1-f0)-(xp1-x0)**2*(f0-fm1))/(xp1-xm1)


      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

      ssl0=0.

      do 2 jy=1,my
      do 2 jz=iz_first,iz_last
      do 2 jx=ix_first,ix_last
      if(psi(jx,jz).lt.ps(mpsa-2)) then
      call tracelineat2(jx,jz,jy,lx,lz,ly,sl,vol)

      llx(:,jx,jz,jy)=lx(:)
      llz(:,jx,jz,jy)=lz(:)
      lly(:,jx,jz,jy)=ly(:)
      ssl(:,jx,jz,jy)=sl(:)
      voll(:,:,jx,jz,jy)=vol(:,:)
      endif      
      bb=sqrt(x(jx,jz,jy,6)**2+x(jx,jz,jy,7)**2+x(jx,jz,jy,8)**2)
      pb(jx,jz,jy)=(x1(jx,jz,jy,6)*xint_dx(jx,jz,2)+x1(jx,jz,jy,8)*xint_dz(jx,jz,2))/bb
2     continue  

      do 3 k=1,kk
      
      do 21 jy=1,my
      do 21 jz=iz_first+2,iz_last-2
      do 21 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
      do isgn=1,2
      jjx=llx(isgn,jx,jz,jy)
      jjz=llz(isgn,jx,jz,jy)
      jjy=lly(isgn,jx,jz,jy)

      ppp(isgn)=voll(1,isgn,jx,jz,jy)*x1(jjx,jjz,jjy,2)  &
               +voll(2,isgn,jx,jz,jy)*x1(jjx+1,jjz,jjy,2) &
	           +voll(3,isgn,jx,jz,jy)*x1(jjx,jjz,jjy+1,2) &
               +voll(4,isgn,jx,jz,jy)*x1(jjx+1,jjz,jjy+1,2) &
	           +voll(5,isgn,jx,jz,jy)*x1(jjx,jjz+1,jjy,2) &
               +voll(6,isgn,jx,jz,jy)*x1(jjx+1,jjz+1,jjy,2) &
	           +voll(7,isgn,jx,jz,jy)*x1(jjx,jjz+1,jjy+1,2) &
               +voll(8,isgn,jx,jz,jy)*x1(jjx+1,jjz+1,jjy+1,2)

      ppb(isgn)=voll(1,isgn,jx,jz,jy)*pb(jjx,jjz,jjy)  &
               +voll(2,isgn,jx,jz,jy)*pb(jjx+1,jjz,jjy) &
	           +voll(3,isgn,jx,jz,jy)*pb(jjx,jjz,jjy+1) &
               +voll(4,isgn,jx,jz,jy)*pb(jjx+1,jjz,jjy+1) &
	           +voll(5,isgn,jx,jz,jy)*pb(jjx,jjz+1,jjy) &
               +voll(6,isgn,jx,jz,jy)*pb(jjx+1,jjz+1,jjy) &
	           +voll(7,isgn,jx,jz,jy)*pb(jjx,jjz+1,jjy+1) &
               +voll(8,isgn,jx,jz,jy)*pb(jjx+1,jjz+1,jjy+1)

      enddo
!       wsm1(jx,jz,jy)=difc(ppp(1),x1(jx,jz,jy,2),ppp(2),ssl(1,jx,jz,jy),ssl0,ssl(2,jx,jz,jy))
!       wsm2(jx,jz,jy)=dif1(ppb(1),pb(jx,jz,jy),ppb(2),ssl(1,jx,jz,jy),ssl0,ssl(2,jx,jz,jy))
!	wsm(jx,jz,jy)=csmp0all*(wsm1(jx,jz,jy)+wsm2(jx,jz,jy))
      wsm(jx,jz,jy)=csmp0all*(difc(ppp(1),x1(jx,jz,jy,2),ppp(2),ssl(1,jx,jz,jy),ssl0,ssl(2,jx,jz,jy)) &
                             +dif1(ppb(1),pb(jx,jz,jy),ppb(2),ssl(1,jx,jz,jy),ssl0,ssl(2,jx,jz,jy)))
      endif
21    continue 
      x(:,:,:,2)=x(:,:,:,2)+wsm(:,:,:)

      call mpi_transfersm(x(:,:,:,2),1)

3     continue 
      
      return
      end

!ws********************************************
      subroutine smthp_traceline_5p(kk)
      use declare
      include 'mpif.h'
      integer,dimension(4) :: lx,lz,ly
      real*8,dimension(4) :: sl,ppp
      real*8,dimension(8,4) :: vol
      integer,dimension(4,mx,mz,my) :: llx,llz,lly
      real*8,dimension(4,mx,mz,my) :: ssl
      real*8,dimension(8,4,mx,mz,my) :: voll
      real*8,dimension(mx,mz,my) :: wsm
      integer kk,k,isgn,jjx,jjz,jjy
      real*8 ssl0,csmp

      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

      ssl0=0.

      do 2 jy=1,my
      do 2 jz=iz_first+2,iz_last-2
      do 2 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
!      if(psi(jx,jz).lt.psia1) then
      call tracelineat_5p(jx,jz,jy,lx,lz,ly,sl,vol)

      llx(:,jx,jz,jy)=lx(:)
      llz(:,jx,jz,jy)=lz(:)
      lly(:,jx,jz,jy)=ly(:)
      ssl(:,jx,jz,jy)=sl(:)
      voll(:,:,jx,jz,jy)=vol(:,:)
      endif
2     continue  

      do 3 k=1,kk
      
      do 21 jy=1,my
      do 21 jz=iz_first+2,iz_last-2
      do 21 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
!      if(psi(jx,jz).lt.psia1) then
      do isgn=1,4
      jjx=llx(isgn,jx,jz,jy)
      jjz=llz(isgn,jx,jz,jy)
      jjy=lly(isgn,jx,jz,jy)

      ppp(isgn)=voll(1,isgn,jx,jz,jy)*x(jjx,jjz,jjy,2)  &
               +voll(2,isgn,jx,jz,jy)*x(jjx+1,jjz,jjy,2) &
	       +voll(3,isgn,jx,jz,jy)*x(jjx,jjz,jjy+1,2) &
               +voll(4,isgn,jx,jz,jy)*x(jjx+1,jjz,jjy+1,2) &
	       +voll(5,isgn,jx,jz,jy)*x(jjx,jjz+1,jjy,2) &
               +voll(6,isgn,jx,jz,jy)*x(jjx+1,jjz+1,jjy,2) &
	       +voll(7,isgn,jx,jz,jy)*x(jjx,jjz+1,jjy+1,2) &
               +voll(8,isgn,jx,jz,jy)*x(jjx+1,jjz+1,jjy+1,2)
      enddo
      csmp=csmp0all*(1-tanh((rr(jx,jz)-0.9)*40))/2.
!      csmp=csmp0all*(1-tanh((psi(jx,jz)-pstrans)/pstransw))/2.
      wsm(jx,jz,jy)=csmp*(0.75*difc(ppp(1),x(jx,jz,jy,2),ppp(2),ssl(1,jx,jz,jy),ssl0,ssl(2,jx,jz,jy)) &
                         +0.25*difc(ppp(3),x(jx,jz,jy,2),ppp(4),ssl(3,jx,jz,jy),ssl0,ssl(4,jx,jz,jy)))
      endif
21    continue 
      x(:,:,:,2)=x(:,:,:,2)+wsm(:,:,:)

      call mpi_transfersm(x(:,:,:,2),1)

3     continue 
      
      return
      end

!ws*****************************************************************
	subroutine tracelineat_5p(jxl,jzl,jyl,lx,lz,ly,sl,vol)
      use declare
      include 'mpif.h'
      integer,dimension(4) :: lx,lz,ly
      real*8,dimension(4) :: sl
      real*8,dimension(8,4) :: vol
      real*8,dimension(3) :: bbb
      real*8 sl0,xxp,zzp,yyp,dyyp,dy1p,dx1p,dz1p,rxp,rxm,rzp,rzm,ryp,rym
      real*8 xxp0,zzp0,yyp0,xxp12,zzp12,yyp12
      real*8 volume,volppp,volmpp,volppm,volmpm,volpmp,volmmp,volpmm,volmmm
      integer isgn,jjx,jjz,jjy,jxl,jzl,jyl,j_add
       
      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)
!      volume=dyy*dxx*dzz
!
!      dyyp=yy(jy+1)-yy(jy)
!      dyym=yy(jy-1)-yy(jy)

!      dy1p=dyyp/nn
!      dy1m=dyym/nn
      sl0=0
      sl(:)=0.

      do 1 isgn=1,2


      dyyp=yy(jyl+(-1)**isgn)-yy(jyl)
      dy1p=dyyp/nline    

100   continue   
      xxp0=xx(jxl)
      zzp0=zz(jzl)
      yyp0=yy(jyl)

      bbb(:)=x(jxl,jzl,jyl,6:8)
      do 2 kp=1,2
      do nc=1,mcycline
!      dy1p=dyyp/nline/xxp0*xzero
      dx1p=xxp0*dy1p*bbb(1)/bbb(2)
      dz1p=xxp0*dy1p*bbb(3)/bbb(2)
    			
	!艗脝脣茫脟掳艙酶脪禄虏艙潞贸碌脛脨脗脦禄脰脙
      xxp=xxp0+dx1p
	  zzp=zzp0+dz1p
      yyp=yyp0+dy1p 
      if(yyp.gt. 2*pi) yyp=yyp-2*pi
      if(yyp.lt. 0) yyp=yyp+2*pi

!      if( (xxp.gt.xx(jxl+1)) .or. (xxp.lt.xx(jxl-1)) &
!      .or.(zzp.gt.zz(jzl+1)) .or. (zzp.lt.zz(jzl-1))) then
!      if(nc.eq.1) then
!      dy1p=dy1p*0.5
!      goto 100
!      else
!      goto 200
!      endif
!      endif

!      sl(isgn)=sl(isgn)+(-1)**isgn*sqrt(dx1p**2+dz1p**2+xxp**2*dy1p**2)
      
      jjx=jxl+j_add(xxp,xx,mx,jxl)
	  jjz=jzl+j_add(zzp,zz,mz,jzl)
	  jjy=jyl+j_add(yyp,yy,my,jyl)
      if(jjy.lt.1) jjy=jjy+my
      if(jjy.gt.my)jjy=jjy-my
!      write(*,*) nrank,'nc',nc,jjy,xxp,jjx,zzp,jjz

      rxp=xx(jjx+1)-xxp
      rxm=xxp-xx(jjx)
      rzp=zz(jjz+1)-zzp
      rzm=zzp-zz(jjz)
      ryp=yy(jjy+1)-yyp
      rym=yyp-yy(jjy)

      volume=(xx(jjx+1)-xx(jjx))*(zz(jjz+1)-zz(jjz))*(yy(jjy+1)-yy(jjy))

      volppp=rxp*rzp*ryp/volume
      volppm=rxp*rzp*rym/volume
      volpmp=rxp*rzm*ryp/volume
      volpmm=rxp*rzm*rym/volume

      volmpp=rxm*rzp*ryp/volume
      volmpm=rxm*rzp*rym/volume
      volmmp=rxm*rzm*ryp/volume
      volmmm=rxm*rzm*rym/volume

      do m=6,8
      bbb(m-5)=(volppp*x(jjx,jjz,jjy,m)       +volmpp*x(jjx+1,jjz,jjy,m) &
	           +volppm*x(jjx,jjz,jjy+1,m)  +volmpm*x(jjx+1,jjz,jjy+1,m) &
	           +volpmp*x(jjx,jjz+1,jjy,m)     +volmmp*x(jjx+1,jjz+1,jjy,m)&
	           +volpmm*x(jjx,jjz+1,jjy+1,m)+volmmm*x(jjx+1,jjz+1,jjy+1,m))/volume
      enddo


      xxp12=xxp0+0.5*dx1p
	  zzp12=zzp0+0.5*dz1p
      yyp12=yyp0+0.5*dy1p

!      dy1p=dyyp/nline/xxp*xzero
      dx1p=xxp*dy1p*bbb(1)/bbb(2)
      dz1p=xxp*dy1p*bbb(3)/bbb(2)

      xxp=xxp12+0.5*dx1p
      zzp=zzp12+0.5*dz1p
      yyp=yyp12+0.5*dy1p
      if(yyp.gt. 2*pi) yyp=yyp-2*pi
      if(yyp.lt. 0) yyp=yyp+2*pi

      sl(isgn+(kp-1)*2)=sl(isgn+(kp-1)*2)+(-1)**isgn*sqrt(dx1p**2+dz1p**2+xxp**2*dy1p**2)
      
      jjx=jxl+j_add(xxp,xx,mx,jxl)
	  jjz=jzl+j_add(zzp,zz,mz,jzl)
	  jjy=jyl+j_add(yyp,yy,my,jyl)
      if(jjy.lt.1) jjy=jjy+my
      if(jjy.gt.my)jjy=jjy-my
!      write(*,*) nrank,'nc',nc,jjy,xxp,jjx,zzp,jjz

      rxp=xx(jjx+1)-xxp
      rxm=xxp-xx(jjx)
      rzp=zz(jjz+1)-zzp
      rzm=zzp-zz(jjz)
      ryp=yy(jjy+1)-yyp
      rym=yyp-yy(jjy)

      volume=(xx(jjx+1)-xx(jjx))*(zz(jjz+1)-zz(jjz))*(yy(jjy+1)-yy(jjy))

      volppp=rxp*rzp*ryp/volume
      volppm=rxp*rzp*rym/volume
      volpmp=rxp*rzm*ryp/volume
      volpmm=rxp*rzm*rym/volume

      volmpp=rxm*rzp*ryp/volume
      volmpm=rxm*rzp*rym/volume
      volmmp=rxm*rzm*ryp/volume
      volmmm=rxm*rzm*rym/volume

      do m=6,8
      bbb(m-5)=(volppp*x(jjx,jjz,jjy,m)       +volmpp*x(jjx+1,jjz,jjy,m) &
	           +volppm*x(jjx,jjz,jjy+1,m)  +volmpm*x(jjx+1,jjz,jjy+1,m) &
	           +volpmp*x(jjx,jjz+1,jjy,m)     +volmmp*x(jjx+1,jjz+1,jjy,m)&
	           +volpmm*x(jjx,jjz+1,jjy+1,m)+volmmm*x(jjx+1,jjz+1,jjy+1,m))/volume
      enddo

      xxp0=xxp
      zzp0=zzp
      yyp0=yyp

      enddo

200   continue

      vol(1,isgn+(kp-1)*2)=volppp
      vol(2,isgn+(kp-1)*2)=volmpp
      vol(3,isgn+(kp-1)*2)=volppm
      vol(4,isgn+(kp-1)*2)=volmpm
      vol(5,isgn+(kp-1)*2)=volpmp
      vol(6,isgn+(kp-1)*2)=volmmp
      vol(7,isgn+(kp-1)*2)=volpmm
      vol(8,isgn+(kp-1)*2)=volmmm
      
      lx(isgn+(kp-1)*2)=jjx
      lz(isgn+(kp-1)*2)=jjz
      ly(isgn+(kp-1)*2)=jjy 
      
2     continue  
                
1     continue

      return
      end

!ws********************************************
      subroutine smthp2_traceline_tm(kk)
      use declare
      include 'mpif.h'
      integer,dimension(2) :: lx,lz,ly
      real*8,dimension(2) :: sl,ppp,ppb
      real*8,dimension(8,2) :: vol
      integer,dimension(2,mx,mz,my) :: llx,llz,lly
      real*8,dimension(2,mx,mz,my) :: ssl
      real*8,dimension(8,2,mx,mz,my) :: voll
      real*8,dimension(mx,mz,my) :: wsm,pb,tm1 !,wsm1,wsm2
      integer kk,k,isgn,jjx,jjz,jjy
      real*8 ssl0,bb

      !  dif1= (d f / dx)*dx2 = [((x0-xm1)/(xp1-x0)*(fp1-f0)-(xp1-x0)/(x0-xm1)*(f0-fm1))/(xp1-xm1)] *(xp1-x0)*(x0-xm1) 
      dif1(fm1,f0,fp1,xm1,x0,xp1)= &
       ((x0-xm1)**2*(fp1-f0)-(xp1-x0)**2*(f0-fm1))/(xp1-xm1)


      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

      ssl0=0.

      do 2 jy=1,my
      do 2 jz=iz_first,iz_last
      do 2 jx=ix_first,ix_last
      if(psi(jx,jz).lt.ps(mpsa-2)) then
      call tracelineat2(jx,jz,jy,lx,lz,ly,sl,vol)

      llx(:,jx,jz,jy)=lx(:)
      llz(:,jx,jz,jy)=lz(:)
      lly(:,jx,jz,jy)=ly(:)
      ssl(:,jx,jz,jy)=sl(:)
      voll(:,:,jx,jz,jy)=vol(:,:)
      endif      
      bb=sqrt(x(jx,jz,jy,6)**2+x(jx,jz,jy,7)**2+x(jx,jz,jy,8)**2)
      pb(jx,jz,jy)=(x1(jx,jz,jy,6)*xint_dx(jx,jz,2)+x1(jx,jz,jy,8)*xint_dz(jx,jz,2))/bb/x(jx,jz,jy,1)
      tm1(jx,jz,jy)=x1(jx,jz,jy,2)/x(jx,jz,jy,1)
2     continue  

      do 3 k=1,kk
      
      do 21 jy=1,my
      do 21 jz=iz_first+2,iz_last-2
      do 21 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
      do isgn=1,2
      jjx=llx(isgn,jx,jz,jy)
      jjz=llz(isgn,jx,jz,jy)
      jjy=lly(isgn,jx,jz,jy)

      ppp(isgn)=voll(1,isgn,jx,jz,jy)*tm1(jjx,jjz,jjy)  &
               +voll(2,isgn,jx,jz,jy)*tm1(jjx+1,jjz,jjy) &
	           +voll(3,isgn,jx,jz,jy)*tm1(jjx,jjz,jjy+1) &
               +voll(4,isgn,jx,jz,jy)*tm1(jjx+1,jjz,jjy+1) &
	           +voll(5,isgn,jx,jz,jy)*tm1(jjx,jjz+1,jjy) &
               +voll(6,isgn,jx,jz,jy)*tm1(jjx+1,jjz+1,jjy) &
	           +voll(7,isgn,jx,jz,jy)*tm1(jjx,jjz+1,jjy+1) &
               +voll(8,isgn,jx,jz,jy)*tm1(jjx+1,jjz+1,jjy+1)

      ppb(isgn)=voll(1,isgn,jx,jz,jy)*pb(jjx,jjz,jjy)  &
               +voll(2,isgn,jx,jz,jy)*pb(jjx+1,jjz,jjy) &
	           +voll(3,isgn,jx,jz,jy)*pb(jjx,jjz,jjy+1) &
               +voll(4,isgn,jx,jz,jy)*pb(jjx+1,jjz,jjy+1) &
	           +voll(5,isgn,jx,jz,jy)*pb(jjx,jjz+1,jjy) &
               +voll(6,isgn,jx,jz,jy)*pb(jjx+1,jjz+1,jjy) &
	           +voll(7,isgn,jx,jz,jy)*pb(jjx,jjz+1,jjy+1) &
               +voll(8,isgn,jx,jz,jy)*pb(jjx+1,jjz+1,jjy+1)

      enddo
!       wsm1(jx,jz,jy)=difc(ppp(1),x1(jx,jz,jy,2),ppp(2),ssl(1,jx,jz,jy),ssl0,ssl(2,jx,jz,jy))
!       wsm2(jx,jz,jy)=dif1(ppb(1),pb(jx,jz,jy),ppb(2),ssl(1,jx,jz,jy),ssl0,ssl(2,jx,jz,jy))
!	wsm(jx,jz,jy)=csmp0all*(wsm1(jx,jz,jy)+wsm2(jx,jz,jy))
      wsm(jx,jz,jy)=csmp0all*(difc(ppp(1),x1(jx,jz,jy,2),ppp(2),ssl(1,jx,jz,jy),ssl0,ssl(2,jx,jz,jy)) &
                             +dif1(ppb(1),pb(jx,jz,jy),ppb(2),ssl(1,jx,jz,jy),ssl0,ssl(2,jx,jz,jy)))*x(jx,jz,jy,1)
      endif
21    continue 
      x(:,:,:,2)=x(:,:,:,2)+wsm(:,:,:)

      call mpi_transfersm(x(:,:,:,2),1)

3     continue 
      
      return
      end

!ws*****************************************************************
	subroutine tracelineat_spec1(jxl,jzl,jyl,lx,lz,ly,sl,yypl,area)
      use declare
      include 'mpif.h'
      integer,dimension(2) :: lx,lz,ly
      real*8,dimension(2) :: sl,yypl
      real*8,dimension(4,2) :: area
!      real*8,dimension(8,2) :: vol
      real*8,dimension(3) :: bbb
      real*8 bkc,bks
      real*8 sl0,xxp,zzp,yyp,dyyp,dy1p,dx1p,dz1p,rxp,rxm,rzp,rzm,ryp,rym,app,apm,amp,amm,atotal
!      real*8 volume,volppp,volmpp,volppm,volmpm,volpmp,volmmp,volpmm,volmmm
      integer isgn,jjx,jjz,jjy,jxl,jzl,jyl
       
      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)
!      volume=dyy*dxx*dzz
!
!      dyyp=yy(jy+1)-yy(jy)
!      dyym=yy(jy-1)-yy(jy)

!      dy1p=dyyp/nn
!      dy1m=dyym/nn
      sl0=0
      sl(:)=0.

      do 1 isgn=1,2


      dyyp=yy(jyl+(-1)**isgn)-yy(jyl)
      dy1p=dyyp/nline    

100   continue   
      xxp=xx(jxl)
      zzp=zz(jzl)
      yyp=yy(jyl)
      
      bbb(:)=x(jxl,jzl,jyl,6:8)

      do nc=1,mcycline
      dx1p=xxp*dy1p*bbb(1)/bbb(2)
      dz1p=xxp*dy1p*bbb(3)/bbb(2)
    			
	!艗脝脣茫脟掳艙酶脪禄虏艙潞贸碌脛脨脗脦禄脰脙
      xxp=xxp+dx1p
	  zzp=zzp+dz1p
      yyp=yyp+dy1p 
      if(yyp.gt. 2*pi) yyp=yyp-2*pi
      if(yyp.lt. 0) yyp=yyp+2*pi

!      if( (xxp.gt.xx(jxl+1)) .or. (xxp.lt.xx(jxl-1)) &
!      .or.(zzp.gt.zz(jzl+1)) .or. (zzp.lt.zz(jzl-1))) then
!      if(nc.eq.1) then
!      dy1p=dy1p*0.5
!      goto 100
!      else
!      goto 200
!      endif
!      endif

      sl(isgn)=sl(isgn)+(-1)**isgn*sqrt(dx1p**2+dz1p**2+xxp**2*dy1p**2)
      
      jjx=floor((xxp-xx(ix_first))/dxx)+ix_first
	  jjz=floor((zzp-zz(iz_first))/dzz)+iz_first   
	  jjy=floor((yyp-yy(1))/dyy)+1
!      if(jjy.lt.1) jjy=jjy+my
!      if(jjy.gt.my)jjy=jjy-my
!      write(*,*) nrank,'nc',nc,jjy,xxp,jjx,zzp,jjz

      rxp=xx(jjx+1)-xxp
      rxm=xxp-xx(jjx)
      rzp=zz(jjz+1)-zzp
      rzm=zzp-zz(jjz)
      atotal=(xx(jjx+1)-xx(jjx))*(zz(jjz+1)-zz(jjz))

      app=rxp*rzp/atotal
      apm=rxp*rzm/atotal
      amp=rxm*rzp/atotal
      amm=rxm*rzm/atotal

      
      do m=1,3
      bbb(m)=app*xint(jjx,jjz,m+5)+apm*xint(jjx,jjz+1,m+5)+amp*xint(jjx+1,jjz,m+5)+amm*xint(jjx+1,jjz+1,m+5)
      do lk=1,myt/2
      bkc=app*xkc(jjx,jjz,lk,m+5)+apm*xkc(jjx,jjz+1,lk,m+5)+amp*xkc(jjx+1,jjz,lk,m+5)+amm*xkc(jjx+1,jjz+1,lk,m+5)
      bks=app*xks(jjx,jjz,lk,m+5)+apm*xks(jjx,jjz+1,lk,m+5)+amp*xks(jjx+1,jjz,lk,m+5)+amm*xks(jjx+1,jjz+1,lk,m+5)
      bbb(m)=bbb(m)+bkc*cos((lk-1)*yyp)+bks*sin((lk-1)*yyp)
      enddo
      enddo
!      ryp=yy(jjy+1)-yyp
!      rym=yyp-yy(jjy)
!
!      volume=(xx(jjx+1)-xx(jjx))*(zz(jjz+1)-zz(jjz))*(yy(jjy+1)-yy(jjy))
!
!      volppp=rxp*rzp*ryp/volume
!      volppm=rxp*rzp*rym/volume
!      volpmp=rxp*rzm*ryp/volume
!      volpmm=rxp*rzm*rym/volume
!
!      volmpp=rxm*rzp*ryp/volume
!      volmpm=rxm*rzp*rym/volume
!      volmmp=rxm*rzm*ryp/volume
!      volmmm=rxm*rzm*rym/volume
!
!      do m=6,8
!      bbb(m-5)=(volppp*x(jjx,jjz,jjy,m)       +volmpp*x(jjx+1,jjz,jjy,m) &
!	           +volppm*x(jjx,jjz,jyp(jjy),m)  +volmpm*x(jjx+1,jjz,jyp(jjy),m) &
!	           +volpmp*x(jjx,jjz+1,jjy,m)     +volmmp*x(jjx+1,jjz+1,jjy,m)&
!	           +volpmm*x(jjx,jjz+1,jyp(jjy),m)+volmmm*x(jjx+1,jjz+1,jyp(jjy),m))/volume
!      enddo

      enddo

200   continue

!      vol(1,isgn)=volppp
!      vol(2,isgn)=volmpp
!      vol(3,isgn)=volppm
!      vol(4,isgn)=volmpm
!      vol(5,isgn)=volpmp
!      vol(6,isgn)=volmmp
!      vol(7,isgn)=volpmm
!      vol(8,isgn)=volmmm
      area(1,isgn)=app
      area(2,isgn)=apm
      area(3,isgn)=amp
      area(4,isgn)=amm
      yypl(isgn)=yyp
      lx(isgn)=jjx
      lz(isgn)=jjz
      ly(isgn)=jjy     
                
1     continue

      return
      end

!ws********************************************
      subroutine smthp_traceline_spec(kk)
      use declare
      include 'mpif.h'
      integer,dimension(2) :: lx,lz,ly
      real*8,dimension(2) :: sl,yypl,ppp
      real*8,dimension(4,2) :: area
      integer,dimension(2,mx,mz,my) :: llx,llz,lly
      real*8,dimension(2,mx,mz,my) :: ssl,yyl
!      real*8,dimension(8,2,mx,mz,my) :: voll
      real*8,dimension(4,2,mx,mz,my) :: arl
      real*8,dimension(mx,mz,my) :: wsm
      integer kk,k,isgn,jjx,jjz,jjy
      real*8 ssl0,csmp,pkc,pks

      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

      ssl0=0.

      do 2 jy=1,my
      do 2 jz=iz_first+2,iz_last-2
      do 2 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
!      if(psi(jx,jz).lt.psia1) then
      call tracelineat_spec1(jx,jz,jy,lx,lz,ly,sl,yypl,area)

      llx(:,jx,jz,jy)=lx(:)
      llz(:,jx,jz,jy)=lz(:)
      lly(:,jx,jz,jy)=ly(:)
      ssl(:,jx,jz,jy)=sl(:)
      yyl(:,jx,jz,jy)=yypl(:)
      arl(:,:,jx,jz,jy)=area(:,:)
      endif
2     continue  

      do 3 k=1,kk
      
      do 21 jy=1,my
      do 21 jz=iz_first+2,iz_last-2
      do 21 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
!      if(psi(jx,jz).lt.psia1) then
      do isgn=1,2
      jjx=llx(isgn,jx,jz,jy)
      jjz=llz(isgn,jx,jz,jy)
!      jjy=lly(isgn,jx,jz,jy)
      ppp(isgn)=arl(1,isgn,jx,jz,jy)*xint(jjx,jjz,2) &
         +arl(2,isgn,jx,jz,jy)*xint(jjx,jjz+1,2) &
         +arl(3,isgn,jx,jz,jy)*xint(jjx+1,jjz,2) &
         +arl(4,isgn,jx,jz,jy)*xint(jjx+1,jjz+1,2)
      do ky=1,myt/2
      pkc=arl(1,isgn,jx,jz,jy)*xkc(jjx,jjz,ky,2) &
         +arl(2,isgn,jx,jz,jy)*xkc(jjx,jjz+1,ky,2) &
         +arl(3,isgn,jx,jz,jy)*xkc(jjx+1,jjz,ky,2) &
         +arl(4,isgn,jx,jz,jy)*xkc(jjx+1,jjz+1,ky,2)

      pks=arl(1,isgn,jx,jz,jy)*xks(jjx,jjz,ky,2) &
         +arl(2,isgn,jx,jz,jy)*xks(jjx,jjz+1,ky,2) &
         +arl(3,isgn,jx,jz,jy)*xks(jjx+1,jjz,ky,2) &
         +arl(4,isgn,jx,jz,jy)*xks(jjx+1,jjz+1,ky,2)

      ppp(isgn)=ppp(isgn)+pkc*cos((ky-1)*yyl(isgn,jx,jz,jy))+pks*sin((ky-1)*yyl(isgn,jx,jz,jy))
      enddo
!      ppp(isgn)=voll(1,isgn,jx,jz,jy)*x(jjx,jjz,jjy,2)  &
!               +voll(2,isgn,jx,jz,jy)*x(jjx+1,jjz,jjy,2) &
!	       +voll(3,isgn,jx,jz,jy)*x(jjx,jjz,jyp(jjy),2) &
!               +voll(4,isgn,jx,jz,jy)*x(jjx+1,jjz,jyp(jjy),2) &
!	       +voll(5,isgn,jx,jz,jy)*x(jjx,jjz+1,jjy,2) &
!               +voll(6,isgn,jx,jz,jy)*x(jjx+1,jjz+1,jjy,2) &
!	       +voll(7,isgn,jx,jz,jy)*x(jjx,jjz+1,jyp(jjy),2) &
!               +voll(8,isgn,jx,jz,jy)*x(jjx+1,jjz+1,jyp(jjy),2)
      enddo
      csmp=csmp0all*(1-tanh((rr(jx,jz)-0.9)*40))/2.
!      csmp=csmp0all*(1-tanh((psi(jx,jz)-pstrans)/pstransw))/2.
      wsm(jx,jz,jy)=csmp*difc(ppp(1),x(jx,jz,jy,2),ppp(2),ssl(1,jx,jz,jy),ssl0,ssl(2,jx,jz,jy))
      endif
21    continue 
      x(:,:,:,2)=x(:,:,:,2)+wsm(:,:,:)

      call mpi_transfersm(x(:,:,:,2),1)

3     continue 
      
      return
      end

!ws********************************************
      subroutine smthp1_traceline_spec(kk)
      use declare
      include 'mpif.h'
      integer,dimension(2) :: lx,lz,ly
      real*8,dimension(2) :: sl,yypl,ppp
      real*8,dimension(4,2) :: area
      integer,dimension(2,mx,mz,my) :: llx,llz,lly
      real*8,dimension(2,mx,mz,my) :: ssl,yyl
!      real*8,dimension(8,2,mx,mz,my) :: voll
      real*8,dimension(4,2,mx,mz,my) :: arl
      real*8,dimension(mx,mz,my) :: wsm
      integer kk,k,isgn,jjx,jjz,jjy
      real*8 ssl0,csmp,pkc,pks

      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

      ssl0=0.

      do 2 jy=1,my
      do 2 jz=iz_first+2,iz_last-2
      do 2 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
!      if(psi(jx,jz).lt.psia1) then
      call tracelineat_spec1(jx,jz,jy,lx,lz,ly,sl,yypl,area)

      llx(:,jx,jz,jy)=lx(:)
      llz(:,jx,jz,jy)=lz(:)
      lly(:,jx,jz,jy)=ly(:)
      ssl(:,jx,jz,jy)=sl(:)
      yyl(:,jx,jz,jy)=yypl(:)
      arl(:,:,jx,jz,jy)=area(:,:)
      endif
2     continue  

      do 3 k=1,kk
      
      do 21 jy=1,my
      do 21 jz=iz_first+2,iz_last-2
      do 21 jx=ix_first+2,ix_last-2
      if(psi(jx,jz).lt.ps(mpsa-2)) then
!      if(psi(jx,jz).lt.psia1) then
      do isgn=1,2
      jjx=llx(isgn,jx,jz,jy)
      jjz=llz(isgn,jx,jz,jy)
!      jjy=lly(isgn,jx,jz,jy)
      ppp(isgn)=0
      do ky=1,myt/2
      pkc=arl(1,isgn,jx,jz,jy)*xkc(jjx,jjz,ky,2) &
         +arl(2,isgn,jx,jz,jy)*xkc(jjx,jjz+1,ky,2) &
         +arl(3,isgn,jx,jz,jy)*xkc(jjx+1,jjz,ky,2) &
         +arl(4,isgn,jx,jz,jy)*xkc(jjx+1,jjz+1,ky,2)

      pks=arl(1,isgn,jx,jz,jy)*xks(jjx,jjz,ky,2) &
         +arl(2,isgn,jx,jz,jy)*xks(jjx,jjz+1,ky,2) &
         +arl(3,isgn,jx,jz,jy)*xks(jjx+1,jjz,ky,2) &
         +arl(4,isgn,jx,jz,jy)*xks(jjx+1,jjz+1,ky,2)

      ppp(isgn)=ppp(isgn)+pkc*cos((ky-1)*yyl(isgn,jx,jz,jy))+pks*sin((ky-1)*yyl(isgn,jx,jz,jy))
      enddo
!      ppp(isgn)=voll(1,isgn,jx,jz,jy)*x(jjx,jjz,jjy,2)  &
!               +voll(2,isgn,jx,jz,jy)*x(jjx+1,jjz,jjy,2) &
!	       +voll(3,isgn,jx,jz,jy)*x(jjx,jjz,jyp(jjy),2) &
!               +voll(4,isgn,jx,jz,jy)*x(jjx+1,jjz,jyp(jjy),2) &
!	       +voll(5,isgn,jx,jz,jy)*x(jjx,jjz+1,jjy,2) &
!               +voll(6,isgn,jx,jz,jy)*x(jjx+1,jjz+1,jjy,2) &
!	       +voll(7,isgn,jx,jz,jy)*x(jjx,jjz+1,jyp(jjy),2) &
!               +voll(8,isgn,jx,jz,jy)*x(jjx+1,jjz+1,jyp(jjy),2)
      enddo
!      csmp=csmp0all*(1-tanh((rr(jx,jz)-0.9)*40))/2.
!      csmp=csmp0all*(1-tanh((psi(jx,jz)-pstrans)/pstransw))/2.
      wsm(jx,jz,jy)=csmp0*difc(ppp(1),x1(jx,jz,jy,2),ppp(2),ssl(1,jx,jz,jy),ssl0,ssl(2,jx,jz,jy))
      endif
21    continue 
      x1(:,:,:,2)=x1(:,:,:,2)+wsm(:,:,:)

      call mpi_transfersm(x1(:,:,:,2),1)

3     continue 
      
      return
      end

!ws**********************************************************************
      subroutine smth_x1st(ms,me,jps,kk)
      use declare
      real*8, dimension(n2th+5) :: wst
      include 'mpif.h'
!
! second-order diffusion
      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

!      difxz(fm1m1,fp1m1,fm1p1,fp1p1,xm1,x0,xp1,zm1,z0,zp1)= &
!       2.*( fp1p1-fm1p1-fp1m1-fm1m1)/((xp1-xm1)*(zp1-zm1)) &
!       *((xp1-x0)*(x0-xm1)*(zp1-z0)*(z0-zm1))**0.5
!
!      difxy(fm1,f0,fp1,xm1,x0,xp1,dyy)= &
!       (2.*((xm1-x0)/(xp1-x0)*(fp1-f0)-(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)/x0-f0/x0**2) &
!       *((xp1-x0)*(x0-xm1))**0.5*dyy
!
!      difzy(fm1,f0,fp1,zm1,z0,zp1,dyy,x0)= &
!       (2.*((zm1-z0)/(zp1-z0)*(fp1-f0)-(zp1-z0)/(zm1-z0)*(fm1-f0))/(zm1-zp1)/x0) &
!       *((zp1-z0)*(z0-zm1))**0.5*dyy


!       integer status(mpi_status_size)
      do 10 k=1,kk
!      call convt
      do 11 m=ms,me
!      if(m.eq.1.or.m.eq.2.or.m.eq.5.or.m.eq.8) then
      do 12 jy=1,my

      do 21 jt=2,n2th+4
      wst(jt)=difc(x1st(jt-1,jps,jy,m),x1st(jt,jps,jy,m),x1st(jt+1,jps,jy,m),thst(jt-1),thst(jt),thst(jt+1))
   21 continue
   

      
      do 22 jt=2,n2th+4
      x1st(jt,jps,jy,m)=x1st(jt,jps,jy,m)+1./4.*wst(jt)
   22 continue

   12 continue 
!      endif
   11 continue
   10 continue
      return
      end

!ws*******************************************
      subroutine smth_ps1(bst,kk)
      use declare
      real*8, dimension(n2th+5) :: bst,wst
      include 'mpif.h'
!
! second-order diffusion
      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

!       integer status(mpi_status_size)
      do 10 k=1,kk
      do 21 jt=2,n2th+4
      wst(jt)=difc(bst(jt-1),bst(jt),bst(jt+1),thst(jt-1),thst(jt),thst(jt+1))
   21 continue
      
      do 22 jt=2,n2th+4
      bst(jt)=bst(jt)+1./4.*wst(jt)
   22 continue
      bst(1)=bst(n2th+1)
      bst(n2th+5)=bst(5)
   10 continue
      return
      end

!ws*******************************************
      subroutine smth_ps(bst,nt,kk)
      use declare
      integer nt
      real*8, dimension(nt) :: bst,wst
      include 'mpif.h'

!       integer status(mpi_status_size)
      do 10 k=1,kk
      do 21 jt=3,nt-2
      wst(jt)=bst(jt)/2.+(bst(jt+1)+bst(jt-1))*3./16+(bst(jt+2)+bst(jt-2))/16.
   21 continue
      
      do 22 jt=3,nt-2
      bst(jt)=wst(jt)
   22 continue

   10 continue
      return
      end

!ws*******************************************
      subroutine smthef_dis(kk)
      use declare
      real*8,dimension(mx,mz,my) :: fsm,wx2,wz2,wy2
      include 'mpif.h'

      dis2(fm1,f0,fp1,xm1,x0,xp1)= &
       ((fp1-f0)/(xp1-x0)-(f0-fm1)/(x0-xm1))

      do 10 k=1,kk
      do 11 m=1,3
      do 21 jy=1,my
      do 21 jz=iz_first+1,iz_last-1
      do 21 jx=ix_first+1,ix_last-1
      if(psi(jx-1,jz).lt.psia .and. psi(jx+1,jz).lt.psia .and. psi(jx,jz-1).lt.psia .and. psi(jx,jz+1).lt.psia) then
      wx2(jx,jz,jy)=dis2(ef(jx-1,jz,jy,m),ef(jx,jz,jy,m),ef(jx+1,jz,jy,m),xx(jx-1),xx(jx),xx(jx+1))
      wz2(jx,jz,jy)=dis2(ef(jx,jz-1,jy,m),ef(jx,jz,jy,m),ef(jx,jz+1,jy,m),zz(jz-1),zz(jz),zz(jz+1))
      endif
   21 continue
   
      do 15 jz=iz_first+1,iz_last-1
      do 15 jx=ix_first+1,ix_last-1
      do 15 jy=iy_first+1,iy_last-1
      wy2(jx,jz,jy)=dis2(ef(jx,jz,jy-1,m),ef(jx,jz,jy,m),ef(jx,jz,jy+1,m),yy(jy-1),yy(jy),yy(jy+1))/xx(jx)

   15 continue

      do 13 jy=iy_first+1,iy_last-1
      do 13 jz=iz_first+1,iz_last-1
      do 13 jx=ix_first+1,ix_last-1
      ef(jx,jz,jy,m)=ef(jx,jz,jy,m)+cfsmb(jx,jz)*dxx*(wx2(jx,jz,jy)+wz2(jx,jz,jy)+wy2(jx,jz,jy))      
!           +(.5*(1.+dtanh(pi/2.-thetati(jx))))/20.*w(jx,jz,jy) &
!           +(.5*(1.+dtanh(pi/2.-thetate(jx))))/20.*w(jx,jz,jy)
   13 continue
   11 continue

      call valb3_atlastgrid_r1p0_v1(ef)
      call mpi_transfersm(ef(:,:,:,:),3)
   10 continue
      return
      end

!ws*******************************************
      subroutine smthef_d2f(kk)
      use declare
      real*8,dimension(mx,mz,my) :: fsm,wx2,wz2,wy2
      include 'mpif.h'
! second-order diffusion
      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

      dif2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)**2*(fp1-f0)-(xp1-x0)**2*(fm1-f0))/(xp1-xm1)/x0
!  d2f2= d2f / dx2  with second-order accuracy central difference
      d2f2(fm1,f0,fp1,xm1,x0,xp1)= &
       2*((fp1-f0)/(xp1-x0)-(f0-fm1)/(x0-xm1))/(xp1-xm1)

      d1f2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)/(xp1-x0)*(fp1-f0) &
         -(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)
      do 10 k=1,kk
      do 11 m=1,3
      do 21 jy=1,my
      do 21 jz=iz_first+1,iz_last-1
      do 21 jx=ix_first+1,ix_last-1
      if(psi(jx-1,jz).lt.psia .and. psi(jx+1,jz).lt.psia .and. psi(jx,jz-1).lt.psia .and. psi(jx,jz+1).lt.psia) then
      wx2(jx,jz,jy)=d2f2(ef(jx-1,jz,jy,m),ef(jx,jz,jy,m),ef(jx+1,jz,jy,m),xx(jx-1),xx(jx),xx(jx+1)) &
                   +d1f2(ef(jx-1,jz,jy,m),ef(jx,jz,jy,m),ef(jx+1,jz,jy,m),xx(jx-1),xx(jx),xx(jx+1))/xx(jx)
      wz2(jx,jz,jy)=d2f2(ef(jx,jz-1,jy,m),ef(jx,jz,jy,m),ef(jx,jz+1,jy,m),zz(jz-1),zz(jz),zz(jz+1))
      endif
   21 continue
   
      do 15 jz=iz_first+1,iz_last-1
      do 15 jx=ix_first+1,ix_last-1
      do 15 jy=iy_first+1,iy_last-1
      wy2(jx,jz,jy)=d2f2(ef(jx,jz,jy-1,m),ef(jx,jz,jy,m),ef(jx,jz,jy+1,m),yy(jy-1),yy(jy),yy(jy+1))/xx(jx)**2

   15 continue

      do 13 jy=iy_first+1,iy_last-1
      do 13 jz=iz_first+1,iz_last-1
      do 13 jx=ix_first+1,ix_last-1
      ef(jx,jz,jy,m)=ef(jx,jz,jy,m)+cfsmb(jx,jz)*1.e-4*(wx2(jx,jz,jy)+wz2(jx,jz,jy)+wy2(jx,jz,jy))      
!           +(.5*(1.+dtanh(pi/2.-thetati(jx))))/20.*w(jx,jz,jy) &
!           +(.5*(1.+dtanh(pi/2.-thetate(jx))))/20.*w(jx,jz,jy)
   13 continue
   11 continue

      call valb3_atlastgrid_r1p0_v1(ef)
      call mpi_transfersm(ef(:,:,:,:),3)
   10 continue
      return
      end

!ws*******************************************
      subroutine smthef(kk)
      use declare
      real*8,dimension(mx,mz,my) :: fsm,wx2,wz2,wy2
      include 'mpif.h'
! second-order diffusion
      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

      dif2(fm1,f0,fp1,xm1,x0,xp1)= &
        ((xm1-x0)**2*(fp1-f0)-(xp1-x0)**2*(fm1-f0))/(xp1-xm1)/x0

      do 10 k=1,kk
      do 11 m=1,3
      do 21 jy=1,my
      do 21 jz=iz_first+1,iz_last-1
      do 21 jx=ix_first+1,ix_last-1
      if(psi(jx-1,jz).lt.psia .and. psi(jx+1,jz).lt.psia .and. psi(jx,jz-1).lt.psia .and. psi(jx,jz+1).lt.psia) then
      wx2(jx,jz,jy)=difc(ef(jx-1,jz,jy,m),ef(jx,jz,jy,m),ef(jx+1,jz,jy,m),xx(jx-1),xx(jx),xx(jx+1)) &
                   +dif2(ef(jx-1,jz,jy,m),ef(jx,jz,jy,m),ef(jx+1,jz,jy,m),xx(jx-1),xx(jx),xx(jx+1))
      wz2(jx,jz,jy)=difc(ef(jx,jz-1,jy,m),ef(jx,jz,jy,m),ef(jx,jz+1,jy,m),zz(jz-1),zz(jz),zz(jz+1))
      endif
   21 continue
   
      do 15 jz=iz_first+1,iz_last-1
      do 15 jx=ix_first+1,ix_last-1
      do 15 jy=iy_first+1,iy_last-1
      wy2(jx,jz,jy)=difc(ef(jx,jz,jy-1,m),ef(jx,jz,jy,m),ef(jx,jz,jy+1,m),yy(jy-1),yy(jy),yy(jy+1))

   15 continue

      do 13 jy=iy_first+1,iy_last-1
      do 13 jz=iz_first+1,iz_last-1
      do 13 jx=ix_first+1,ix_last-1
      ef(jx,jz,jy,m)=ef(jx,jz,jy,m)+cfsmb(jx,jz)*(wx2(jx,jz,jy)+wz2(jx,jz,jy)+wy2(jx,jz,jy))      
!           +(.5*(1.+dtanh(pi/2.-thetati(jx))))/20.*w(jx,jz,jy) &
!           +(.5*(1.+dtanh(pi/2.-thetate(jx))))/20.*w(jx,jz,jy)
   13 continue
   11 continue

      call valb3_atlastgrid_r1p0_v1(ef)
      call mpi_transfersm(ef(:,:,:,:),3)
   10 continue
      return
      end

!ws*******************************************
      subroutine smthef_dis_v2(kk)
      use declare
        
        real*8 coeff_xm, coeff_xp, coeff_zm, coeff_zp, coeff_ym, coeff_yp
        real*8 coeff_total
        real*8, dimension(mx,mz,my,3) :: ef_smooth
       include 'mpif.h' 
!        coeff_smooth = 0.6!large coefficient means heavily smooth

      do 10 k=1,kk
  
      do 21 jy=iy_first+2,iy_last-2
      do 21 jz=iz_first+2,iz_last-2
      do 21 jx=ix_first+2,ix_last-2
      if(psi(jx-1,jz).lt.psia .and. psi(jx+1,jz).lt.psia .and. psi(jx,jz-1).lt.psia .and. psi(jx,jz+1).lt.psia) then
            coeff_xm = 1.0/(xx(jx)-xx(jx-1))
            coeff_xp = 1.0/(xx(jx+1)-xx(jx))
            coeff_zm = 1.0/(zz(jz)-zz(jz-1))
            coeff_zp = 1.0/(zz(jz+1)-zz(jz))
            coeff_ym = 1.0/(xx(jx)*(yy(jy)-yy(jy-1)))
            coeff_yp = 1.0/(xx(jx)*(yy(jy+1)-yy(jy)))
            coeff_total = coeff_xm+coeff_xp+coeff_zm+coeff_zp+coeff_ym+coeff_yp
        ef_smooth(jx,jz,jy,:) = (1.0 - cfsmb(jx,jz))*ef(jx,jz,jy,:) &
                              + cfsmb(jx,jz)*(coeff_xm*ef(jx-1,jz,jy,:) + coeff_xp*ef(jx+1,jz,jy,:) &
                              +               coeff_zm*ef(jx,jz-1,jy,:) + coeff_zp*ef(jx,jz+1,jy,:) &
                              +               coeff_ym*ef(jx,jz,jy-1,:) + coeff_yp*ef(jx,jz,jy+1,:))/coeff_total
      endif
   21 continue
       
        ef = ef_smooth
      call valb3_atlastgrid_r1p0_v1(ef)
      call mpi_transfersm(ef(:,:,:,:),3)   
   10 continue
      return
      end

!ws*******************************************
      subroutine smthxzy_dis_v2(ws,mm,kk)
      use declare
        integer mm
        real*8 coeff_xm, coeff_xp, coeff_zm, coeff_zp, coeff_ym, coeff_yp
        real*8 coeff_total
        real*8, dimension(mx,mz,my,mm) :: ws, ws_smooth
       include 'mpif.h' 
!        coeff_smooth = 0.6!large coefficient means heavily smooth

      do 10 k=1,kk
  
      do 21 jy=iy_first+2,iy_last-2
      do 21 jz=iz_first+2,iz_last-2
      do 21 jx=ix_first+2,ix_last-2
      if(psi(jx-1,jz).lt.psia .and. psi(jx+1,jz).lt.psia .and. psi(jx,jz-1).lt.psia .and. psi(jx,jz+1).lt.psia) then
            coeff_xm = 1.0/(xx(jx)-xx(jx-1))
            coeff_xp = 1.0/(xx(jx+1)-xx(jx))
            coeff_zm = 1.0/(zz(jz)-zz(jz-1))
            coeff_zp = 1.0/(zz(jz+1)-zz(jz))
            coeff_ym = 1.0/(xx(jx)*(yy(jy)-yy(jy-1)))
            coeff_yp = 1.0/(xx(jx)*(yy(jy+1)-yy(jy)))
            coeff_total = coeff_xm+coeff_xp+coeff_zm+coeff_zp+coeff_ym+coeff_yp
        ws_smooth(jx,jz,jy,:) = (1.0 - cfsmb(jx,jz))*ws(jx,jz,jy,:) &
                                  + cfsmb(jx,jz)*(coeff_xm*ws(jx-1,jz,jy,:) + coeff_xp*ws(jx+1,jz,jy,:) &
                                  +               coeff_zm*ws(jx,jz-1,jy,:) + coeff_zp*ws(jx,jz+1,jy,:) &
                                  +               coeff_ym*ws(jx,jz,jy-1,:) + coeff_yp*ws(jx,jz,jy+1,:))/coeff_total
      endif
   21 continue
       
        ws = ws_smooth
!      call valb3_atlastgrid_r1p0_v1(ef)
      call mpi_transfersm(ws,mm)   
   10 continue
      return
      end

!ws*******************************************
      subroutine smthxzy_dis(fsm,kk)
      use declare
      real*8,dimension(mx,mz,my) :: fsm,wx2,wz2,wy2
      include 'mpif.h'

      dis2(fm1,f0,fp1,xm1,x0,xp1)= &
       ((fp1-f0)/(xp1-x0)-(f0-fm1)/(x0-xm1))
      
      do 10 k=1,kk
      do 21 jy=1,my
      do 21 jz=iz_first+1,iz_last-1
      do 21 jx=ix_first+1,ix_last-1
      if(psi(jx-1,jz).lt.psia .and. psi(jx+1,jz).lt.psia .and. psi(jx,jz-1).lt.psia .and. psi(jx,jz+1).lt.psia) then
      wx2(jx,jz,jy)=dis2(fsm(jx-1,jz,jy),fsm(jx,jz,jy),fsm(jx+1,jz,jy),xx(jx-1),xx(jx),xx(jx+1))
      wz2(jx,jz,jy)=dis2(fsm(jx,jz-1,jy),fsm(jx,jz,jy),fsm(jx,jz+1,jy),zz(jz-1),zz(jz),zz(jz+1))
      endif
   21 continue
   
      do 15 jz=iz_first+1,iz_last-1
      do 15 jx=ix_first+1,ix_last-1
      do 15 jy=iy_first+1,iy_last-1
      wy2(jx,jz,jy)=dis2(fsm(jx,jz,jy-1),fsm(jx,jz,jy),fsm(jx,jz,jy+1),yy(jy-1),yy(jy),yy(jy+1))/xx(jx)

   15 continue

      do 13 jy=iy_first+1,iy_last-1
      do 13 jz=iz_first+1,iz_last-1
      do 13 jx=ix_first+1,ix_last-1
      fsm(jx,jz,jy)=fsm(jx,jz,jy)+cfsmb(jx,jz)*dxx*(wx2(jx,jz,jy)+wz2(jx,jz,jy)+wy2(jx,jz,jy))      
!           +(.5*(1.+dtanh(pi/2.-thetati(jx))))/20.*w(jx,jz,jy) &
!           +(.5*(1.+dtanh(pi/2.-thetate(jx))))/20.*w(jx,jz,jy)
   13 continue

   10 continue
      return
      end

!ws*******************************************
      subroutine smthxzy(fsm,kk)
      use declare
      real*8,dimension(mx,mz,my) :: fsm,wx2,wz2,wy2
      include 'mpif.h'
! second-order diffusion
      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

!      difxz(fm1m1,fp1m1,fm1p1,fp1p1,xm1,x0,xp1,zm1,z0,zp1)= &
!       2.*( fp1p1-fm1p1-fp1m1-fm1m1)/((xp1-xm1)*(zp1-zm1)) &
!       *((xp1-x0)*(x0-xm1)*(zp1-z0)*(z0-zm1))**0.5
!
!      difxy(fm1,f0,fp1,xm1,x0,xp1,dyy)= &
!       (2.*((xm1-x0)/(xp1-x0)*(fp1-f0)-(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)/x0-f0/x0**2) &
!       *((xp1-x0)*(x0-xm1))**0.5*dyy
!
!      difzy(fm1,f0,fp1,zm1,z0,zp1,dyy,x0)= &
!       (2.*((zm1-z0)/(zp1-z0)*(fp1-f0)-(zp1-z0)/(zm1-z0)*(fm1-f0))/(zm1-zp1)/x0) &
!       *((zp1-z0)*(z0-zm1))**0.5*dyy


!       integer status(mpi_status_size)
!      do 10 k=1,kk
!      do 12 jy=1,my
!      do 21 jz=jzamin,jzamax
!      do jx=jxam(jz)+1,jxap(jz)-1
!      wx2(jx,jz,jy)=difc(fsm(jx-1,jz,jy),fsm(jx,jz,jy),fsm(jx+1,jz,jy),xx(jx-1),xx(jx),xx(jx+1))
!      enddo
!   21 continue
!   
!      do 22 jx=jxamin,jxamax
!      do jz=jzam(jx)+1,jzap(jx)-1
!      wz2(jx,jz,jy)=difc(fsm(jx,jz-1,jy),fsm(jx,jz,jy),fsm(jx,jz+1,jy),zz(jz-1),zz(jz),zz(jz+1))
!      enddo
!   22 continue
!   12 continue
      do 10 k=1,kk
      do 21 jy=1,my
      do 21 jz=iz_first+1,iz_last-1
      do 21 jx=ix_first+1,ix_last-1
      if(psi(jx-1,jz).lt.psia .and. psi(jx+1,jz).lt.psia .and. psi(jx,jz-1).lt.psia .and. psi(jx,jz+1).lt.psia) then
      wx2(jx,jz,jy)=difc(fsm(jx-1,jz,jy),fsm(jx,jz,jy),fsm(jx+1,jz,jy),xx(jx-1),xx(jx),xx(jx+1))
      wz2(jx,jz,jy)=difc(fsm(jx,jz-1,jy),fsm(jx,jz,jy),fsm(jx,jz+1,jy),zz(jz-1),zz(jz),zz(jz+1))
      endif
   21 continue
   
      do 15 jz=iz_first+1,iz_last-1
      do 15 jx=ix_first+1,ix_last-1
      do 15 jy=iy_first+1,iy_last-1
      wy2(jx,jz,jy)=difc(fsm(jx,jz,jy-1),fsm(jx,jz,jy),fsm(jx,jz,jy+1),yy(jy-1),yy(jy),yy(jy+1))

   15 continue

      do 13 jy=iy_first+1,iy_last-1
      do 13 jz=iz_first+1,iz_last-1
      do 13 jx=ix_first+1,ix_last-1
      fsm(jx,jz,jy)=fsm(jx,jz,jy)+cfsmb(jx,jz)*(wx2(jx,jz,jy)+wz2(jx,jz,jy)+wy2(jx,jz,jy))      
!           +(.5*(1.+dtanh(pi/2.-thetati(jx))))/20.*w(jx,jz,jy) &
!           +(.5*(1.+dtanh(pi/2.-thetate(jx))))/20.*w(jx,jz,jy)
   13 continue

   10 continue
      return
      end

!ws*******************************************
      subroutine smthe4(fsm,kk)
      use declare
      real*8,dimension(mx,mz,my) :: fsm,wx2,wz2,wy2
      include 'mpif.h'
! second-order diffusion
      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

!      difxz(fm1m1,fp1m1,fm1p1,fp1p1,xm1,x0,xp1,zm1,z0,zp1)= &
!       2.*( fp1p1-fm1p1-fp1m1-fm1m1)/((xp1-xm1)*(zp1-zm1)) &
!       *((xp1-x0)*(x0-xm1)*(zp1-z0)*(z0-zm1))**0.5
!
!      difxy(fm1,f0,fp1,xm1,x0,xp1,dyy)= &
!       (2.*((xm1-x0)/(xp1-x0)*(fp1-f0)-(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)/x0-f0/x0**2) &
!       *((xp1-x0)*(x0-xm1))**0.5*dyy
!
!      difzy(fm1,f0,fp1,zm1,z0,zp1,dyy,x0)= &
!       (2.*((zm1-z0)/(zp1-z0)*(fp1-f0)-(zp1-z0)/(zm1-z0)*(fm1-f0))/(zm1-zp1)/x0) &
!       *((zp1-z0)*(z0-zm1))**0.5*dyy


!       integer status(mpi_status_size)
      do 10 k=1,kk
      do 12 jy=1,my
      do 21 jz=iz_first+1,iz_last-1
      do 21 jx=ix_first+1,ix_last-1
      if(psi(jx-1,jz).lt.psia .and. psi(jx+1,jz).lt.psia .and. psi(jx,jz-1).lt.psia .and. psi(jx,jz+1).lt.psia) then
      wx2(jx,jz,jy)=difc(fsm(jx-1,jz,jy),fsm(jx,jz,jy),fsm(jx+1,jz,jy),xx(jx-1),xx(jx),xx(jx+1))
      wz2(jx,jz,jy)=difc(fsm(jx,jz-1,jy),fsm(jx,jz,jy),fsm(jx,jz+1,jy),zz(jz-1),zz(jz),zz(jz+1))
      endif
   21 continue
   12 continue

      do 15 jz=iz_first+1,iz_last-1
      do 15 jx=ix_first+1,ix_last-1
      do 15 jy=iy_first+1,iy_last-1
      wy2(jx,jz,jy)=difc(fsm(jx,jz,jy-1),fsm(jx,jz,jy),fsm(jx,jz,jy+1),yy(jy-1),yy(jy),yy(jy+1))
   15 continue

      do 13 jy=iy_first+1,iy_last-1
      do 13 jz=iz_first+1,iz_last-1
      do 13 jx=ix_first+1,ix_last-1
      fsm(jx,jz,jy)=fsm(jx,jz,jy)+1./12*(wx2(jx,jz,jy)+wz2(jx,jz,jy)+wy2(jx,jz,jy))      
!           +(.5*(1.+dtanh(pi/2.-thetati(jx))))/20.*w(jx,jz,jy) &
!           +(.5*(1.+dtanh(pi/2.-thetate(jx))))/20.*w(jx,jz,jy)
   13 continue

   10 continue
      return
      end

!ws*******************************************
      subroutine avrgt(fsm,kk)
      use declare
      real*8,dimension(mx,mz,my) :: fsm,wx2,wz2,wy2
      include 'mpif.h'
! second-order diffusion
      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)

!      difxz(fm1m1,fp1m1,fm1p1,fp1p1,xm1,x0,xp1,zm1,z0,zp1)= &
!       2.*( fp1p1-fm1p1-fp1m1-fm1m1)/((xp1-xm1)*(zp1-zm1)) &
!       *((xp1-x0)*(x0-xm1)*(zp1-z0)*(z0-zm1))**0.5
!
!      difxy(fm1,f0,fp1,xm1,x0,xp1,dyy)= &
!       (2.*((xm1-x0)/(xp1-x0)*(fp1-f0)-(xp1-x0)/(xm1-x0)*(fm1-f0))/(xm1-xp1)/x0-f0/x0**2) &
!       *((xp1-x0)*(x0-xm1))**0.5*dyy
!
!      difzy(fm1,f0,fp1,zm1,z0,zp1,dyy,x0)= &
!       (2.*((zm1-z0)/(zp1-z0)*(fp1-f0)-(zp1-z0)/(zm1-z0)*(fm1-f0))/(zm1-zp1)/x0) &
!       *((zp1-z0)*(z0-zm1))**0.5*dyy


!       integer status(mpi_status_size)
      do 10 k=1,kk
      do 12 jy=1,my
      do 21 jz=iz_first+1,iz_last-1
      do 21 jx=ix_first+1,ix_last-1
      if(psi(jx-1,jz).lt.psia .and. psi(jx+1,jz).lt.psia .and. psi(jx,jz-1).lt.psia .and. psi(jx,jz+1).lt.psia) then
      wx2(jx,jz,jy)=difc(fsm(jx-1,jz,jy),fsm(jx,jz,jy),fsm(jx+1,jz,jy),xx(jx-1),xx(jx),xx(jx+1))
      wz2(jx,jz,jy)=difc(fsm(jx,jz-1,jy),fsm(jx,jz,jy),fsm(jx,jz+1,jy),zz(jz-1),zz(jz),zz(jz+1))
      endif
   21 continue
   12 continue

      do 15 jz=iz_first+1,iz_last-1
      do 15 jx=ix_first+1,ix_last-1
      do 15 jy=iy_first+1,iy_last-1
      wy2(jx,jz,jy)=difc(fsm(jx,jz,jy-1),fsm(jx,jz,jy),fsm(jx,jz,jy+1),yy(jy-1),yy(jy),yy(jy+1))
   15 continue

      do 13 jy=iy_first+1,iy_last-1
      do 13 jz=iz_first+1,iz_last-1
      do 13 jx=ix_first+1,ix_last-1
      fsm(jx,jz,jy)=fsm(jx,jz,jy)+(1-caf)/1024.0*(wx2(jx,jz,jy)+wz2(jx,jz,jy)+wy2(jx,jz,jy))      
!           +(.5*(1.+dtanh(pi/2.-thetati(jx))))/20.*w(jx,jz,jy) &
!           +(.5*(1.+dtanh(pi/2.-thetate(jx))))/20.*w(jx,jz,jy)
   13 continue

   10 continue
      return
      end

!ws********************************************************************
!ws********************************************************************
subroutine getnp2 ( px, py, x, y, nr, lcell, lnext, xmin, ymin, &
  dx, dy, np, dsq )
!
!***********************************************************************
!
!! getnp2 seeks the closest unmarked node to a point.
!
!
!  discussion:
!
!    getnp2 uses the cell method to find the closest unmarked node np
!    to a specified point p, given a set of n nodes and the data structure 
!    defined by store2.
!
!    np is then marked by negating lnext(np).  thus, the closest m nodes to
!    p may be determined by a sequence of m calls to this routine.  
!
!    if the point p is itself actually a node k, and you want to find the
!    nearest point to p that is not node k, then you must be sure to mark
!    node k before calling.
!
!    the search is begun in the cell containing or closest to p and proceeds 
!    outward in rectangular layers until all cells which contain points 
!    within distance r of p have been searched.  r is the distance from p to 
!    the first unmarked node encountered, or infinite if no unmarked nodes
!    are present.
!
!    input parameters other than lnext are not altered by this routine.  
!    with the exception of ( px, py ) and the signs of lnext elements, 
!    these parameters should be unaltered from their values on output 
!    from subroutine store2.
!
!  modified:
!
!    10 july 1999
!
!  author:
!
!    robert renka,
!    university of north texas
!
!  parameters:
!
!    input, real px, py, the (x,y) coordinates of the point p whose
!    nearest unmarked neighbor is to be found.
!
!    input, real x(n), y(n), the coordinates of the nodes at which
!    data has been supplied.
!
!    input, integer nr, the number of rows and columns in the cell grid.
!    nr must be at least 1.
!
!    input, integer lcell(nr,nr), array of nodal indices associated
!    with cells.
!
!    input/output, integer lnext(n), contains next-node indices ( or their 
!    negatives ).  on return, if the output value of np is nonzero, then
!    lnext(np) will be negative.
!
!    input, real xmin, ymin, dx, dy, the minimum nodal x, y coordinates,
!    and the x, y dimensions of a cell.  dx and dy must be positive.
!
!    output, integer np, the index into the vectors x and y of the nearest
!    unmarked node to the point p.  np will be 0 if all nodes are marked 
!    or if the values of nr, dx, dy are illegal.  lnext(np) will be less
!    than 0 if np is nonzero (this marks node np as being used now).
!
!    output, real dsq, if np is nonzero, then dsq is the squared distance
!    between p and node np.
!
!  local parameters:
!
!    first = true iff the first unmarked node has yet to be encountered,
!
!    imin,imax,jmin,jmax = cell indices defining the range of the search,
!
!    delx,dely = px-xmin and py-ymin,
!
!    i0,j0 = cell containing or closest to p,
!
!    i1,i2,j1,j2 = cell indices of the layer whose intersection with 
!    the range defined by imin,...,jmax is currently being searched.
!
  implicit none
!
  integer nr
!
  real*8 delx
  real*8 dely
  real*8 dsq
  real*8 dx
  real*8 dy
  logical first
  integer i
  integer i0
  integer i1
  integer i2
  integer imax
  integer imin
  integer j
  integer j0
  integer j1
  integer j2
  integer jmax
  integer jmin
  integer l
  integer lcell(nr,nr)
  integer lmin
  integer ln
  integer lnext(*)
  integer np
  real*8 px
  real*8 py
  real*8 r
  real*8 rsmin
  real*8 rsq
  real*8 x(*)
  real*8 xmin
  real*8 xp
  real*8 y(*)
  real*8 ymin
  real*8 yp
!
  xp = px
  yp = py
!
!  test for invalid input parameters.
!
  if ( nr < 1 .or. dx <= 0.0e+00 .or. dy <= 0.0e+00 ) then
    np = 0
    dsq = 0.0e+00
  end if
!
!  initialize parameters:
!
  first = .true.
  imin = 1
  imax = nr
  jmin = 1
  jmax = nr
  delx = xp - xmin
  dely = yp - ymin

  i0 = int ( delx / dx ) + 1
  i0 = max ( i0, 1 )
  i0 = min ( i0, nr )

  j0 = int ( dely / dy ) + 1
  j0 = max ( j0, 1 )
  j0 = min ( j0, nr )

  i1 = i0
  i2 = i0
  j1 = j0
  j2 = j0
!
!  outer loop on layers, inner loop on layer cells, excluding
!  those outside the range (imin,imax) x (jmin,jmax).
!
1 continue

  do j = j1, j2

    if ( j > jmax ) go to 7
    if ( j < jmin ) go to 6

    do i = i1, i2

      if ( i > imax ) go to 6
      if ( i < imin ) go to 5

      if ( j /= j1 .and. j /= j2 .and. i /= i1 .and. i /= i2 ) then
        go to 5
      end if
!
!  search cell (i,j) for unmarked nodes l.
!
      l = lcell(i,j)

      if ( l > 0 ) then
!
!  loop on nodes in cell (i,j).
!
2       continue

        ln = lnext(l)
!
!  node l is the first unmarked neighbor of p encountered.
!
!  initialize lmin to the current candidate for np, and
!  rsmin to the squared distance from p to lmin.  imin,
!  imax, jmin, and jmax are updated to define the smal-
!  lest rectangle containing a circle of radius r =
!  sqrt(rsmin) centered at p, and contained in (1,nr) x
!  (1,nr) (except that, if p is outside the rectangle
!  defined by the nodes, it is possible that imin .gt.
!  nr, imax < 1, jmin > nr, or jmax < 1).
!
        if ( ln >= 0 ) then

          rsq = ( x(l) - xp )**2 + ( y(l) - yp )**2

          if ( first ) then

            lmin = l
            rsmin = rsq
            r = sqrt ( rsmin )

            imin = int ( ( delx - r ) / dx ) + 1
            imin = max ( imin, 1 )

            imax = int ( ( delx + r ) / dx ) + 1
            imax = min ( imax, nr )

            jmin = int ( ( dely - r ) / dy ) + 1
            jmin = max ( jmin, 1 )

            jmax = int ( ( dely + r ) / dy ) + 1
            jmax = min ( jmax, nr )

            first = .false.

          else

            if ( rsq < rsmin ) then
              lmin = l
              rsmin = rsq
            end if
 
          end if

        end if

        if ( abs ( ln ) /= l ) then
          l = abs ( ln )
          go to 2
        end if

      end if

5     continue

    end do

6   continue

  end do
!
!  test for termination of loop on cell layers.
!
7 continue

  if ( i1 > imin .or. i2 < imax .or. j1 > jmin .or. j2 < jmax ) then

    i1 = i1 - 1
    i2 = i2 + 1
    j1 = j1 - 1
    j2 = j2 + 1
    go to 1

  end if

  if ( first ) then
    np = 0
    dsq = 0.0e+00
  else
    np = lmin
    dsq = rsmin
    lnext(lmin) = -lnext(lmin)
  end if

  return
end
subroutine givens ( a, b, c, s )
!
!***********************************************************************
!
!! givens constructs a givens plane rotation.
!
!
!  discussion:
!
!    the transformation has the form of a 2 by 2 matrix g(c,s):
!
!      (   c  s )
!      ( - s  c )
!
!    where c*c + s*s = 1, which zeroes the second entry of the
!    the column vector ( a, b ) when c and s are properly chosen.
!    a call to givens is normally followed by a call to rotate
!    which computes the product of g(c,s) with a 2 by n matrix.
!
!  modified:
!
!    10 july 1999
!
!  parameters:
!
!    input/output, real*8 a, b.
!
!    on input, a and b define the 2-vector whose second entry (b) is
!    to be annihilated by a givens rotation.
!
!    on output, a has been overwritten by a value
!      r = +/- sqrt ( a*a + b*b )
!    and b has been overwritten by a value z which allows c
!    and s to be recovered as:
!
!      if | z | <= 1, then
!        c = sqrt ( 1 - z*z ), 
!        s = z
!      else if | z | > 1 then
!        c = 1 / z, 
!        s = sqrt ( 1 - c*c ).
!
!    output, real*8 c, s, the components of the givens transformation, 
!    which may be computed by:
!
!      c = +/- a / sqrt ( a*a + b*b )
!      s = +/- b / sqrt ( a*a + b*b )
!
!  local parameters:
!
!  r =        c*a + s*b = +/-sqrt(a*a+b*b)
!  u,v =   variables used to scale a and b for computing r
!
!  abs(a) > abs(b)
!
!  note that r has the sign of a, c > 0, and s has
!  sign(a)*sign(b).
!
  implicit none
!
  real*8 a
  real*8 b
  real*8 c
  real*8 r
  real*8 s
  real*8 u
  real*8 v
!
  if ( abs ( a ) > abs ( b ) ) then

    u = 2.0e+00 * a
    v = b / u
    r = sqrt ( 0.25e+00 + v * v ) * u
    c = a / r
    s = 2.0e+00 * v * c
    b = s
    a = r
!
!  abs(a) <= abs(b)
!
!  store r in a.
!  note that r has the sign of b, s > 0, and c has sign(a)*sign(b).
!
  else if ( b /= 0.0e+00 ) then

    u = 2.0e+00 * b
    v = a / u
    a = sqrt ( 0.25e+00 + v * v ) * u
    s = b / a
    c = 2.0e+00 * v * s

    if ( c /= 0.0e+00 ) then
      b = 1.0e+00 / c
    else
      b = 1.0e+00
    end if
!
!  a = b = 0.
!
  else

    c = 1.0e+00
    s = 0.0e+00

  end if

  return
end
subroutine qs2grd ( px, py, n, x, y, f, nr, lcell, lnext, xmin, &
  ymin, dx, dy, rmax, rsq, a, q, qx, qy, ier )
!
!***********************************************************************
!
!! qs2grd evaluates the interpolant and its first spatial derivatives.
!
!
!  discussion:
!
!    qs2grd computes the value and the gradient at the point (px,py) 
!    of the interpolatory function q, defined by qshep2 for a given set
!    of scattered data.  q(x,y) is a weighted sum of quadratic
!    nodal functions.
!
!    input parameters are not altered by this subroutine.  the parameters 
!    other than px and py should be input unaltered from their values 
!    on output from qshep2.  this subroutine should not be called if a 
!    nonzero error flag was returned by qshep2.
!
!  modified:
!
!    10 july 1999
!
!  author:
!
!    robert renka,
!    university of north texas
!
!  parameters:
!
!    input, real*8 px, py, the coordinates of the point at which the
!    interpolant and its derivatives are to be evaluated.
!
!    input, integer n, the number of nodes and data values which
!    are to be interpolated.  n must be at least 6. 
!
!    input, real*8 x(n), y(n), the coordinates of the nodes at which
!    data has been supplied.
!
!    input, real*8 f(n), the data values at the nodes.
!
!    input, integer nr, the number of rows and columns in the cell 
!    grid.  refer to subroutine store2 for details.  nr must be at least 1.
!
!    input, integer lcell(nr,nr), array of nodal indices associated
!    with cells.
!
!    input, integer lnext(n), contains next-node indices.
!
!    input, real*8 xmin, ymin, dx, dy, the minimum nodal x, y coordinates,
!    and the x, y dimensions of a cell.  computed by qshep2.
!
!    input, real*8 rmax, the square root of the largest element in rsq,
!    the maximum radius of influence.  computed by qshep2.
!
!    input, real*8 rsq(n), the squared radii which enter into the weights 
!    defining the interpolant q.  computed by qshep2.
!
!    input, real*8 a(5,n), the coefficients for the nodal functions 
!    defining the interpolant q.  computed by qshep2.
!
!    output, real*8 q, qx, qy, the value of the interpolant, and its derivatives
!    with respect to x and y, at (px,py).
!
!    output, integer ier, error indicator.
!    0, if no errors were encountered.
!    1, if n, nr, dx, dy or rmax is invalid.
!    2, if no errors were encountered but (px,py) is not within the 
!       radius r(k) for any node k and thus q = qx = qy = 0.
!
  implicit none
!
  integer n
  integer nr
!
  real*8 a(5,n)
  real*8 delx
  real*8 dely
  real*8 ds
  real*8 dx
  real*8 dy
  real*8 f(n)
  integer i
  integer ier
  integer imax
  integer imin
  integer j
  integer jmax
  integer jmin
  integer k
  integer kp
  integer lcell(nr,nr)
  integer lnext(n)
  real*8 px
  real*8 py
  real*8 q
  real*8 qk
  real*8 qkx
  real*8 qky
  real*8 qx
  real*8 qy
  real*8 rd
  real*8 rds
  real*8 rmax
  real*8 rs
  real*8 rsq(n)
  real*8 sw
  real*8 swq
  real*8 swqx
  real*8 swqy
  real*8 sws
  real*8 swx
  real*8 swy
  real*8 t
  real*8 w
  real*8 wx
  real*8 wy
  real*8 x(n)
  real*8 xmin
  real*8 xp
  real*8 y(n)
  real*8 ymin
  real*8 yp
!
  xp = px
  yp = py

  if ( n < 6 ) then
    ier = 1
    return
  else if ( nr < 1 ) then
    ier = 1
    return
  else if ( dx <= 0.0e+00 ) then
    ier = 1
    return
  else if ( dy <= 0.0e+00 ) then
    ier = 1
    return
  else if ( rmax < 0.0e+00 ) then
    ier = 1
    return
  end if
!
!  set imin, imax, jmin, and jmax to cell indices defining
!  the range of the search for nodes whose radii include p.
!  the cells which must be searched are those inter-
!  sected by (or contained in) a circle of radius rmax
!  centered at p.
!
  imin = int ( ( xp - xmin - rmax ) / dx ) + 1
  imin = max ( imin, 1 )

  imax = int ( ( xp - xmin + rmax ) / dx ) + 1
  imax = min ( imax, nr )

  jmin = int ( ( yp - ymin - rmax ) / dy ) + 1
  jmin = max ( jmin, 1 )

  jmax = int ( ( yp - ymin + rmax ) / dy ) + 1
  jmax = min ( jmax, nr )
!
!  test for no cells within the circle of radius rmax.
!
  if ( imin > imax .or. jmin > jmax ) then
    q = 0.0e+00
    qx = 0.0e+00
    qy = 0.0e+00
    ier = 2
    return
  end if
!
!  q = swq/sw = sum(w(k)*q(k))/sum(w(k)) where the sum is
!  from k = 1 to n, q(k) is the quadratic nodal function,
!  and w(k) = ((r-d)+/(r*d))**2 for radius r(k) and distance d(k).  thus
!
!    qx = (swqx*sw - swq*swx)/sw**2  and
!    qy = (swqy*sw - swq*swy)/sw**2
!
!  where swqx and swx are partial derivatives with respect
!  to x of swq and sw, respectively.  swqy and swy are 
!  defined similarly.
!
  sw = 0.0e+00
  swx = 0.0e+00
  swy = 0.0e+00
  swq = 0.0e+00
  swqx = 0.0e+00
  swqy = 0.0e+00
!
!  outer loop on cells (i,j).
!
  do j = jmin, jmax

    do i = imin, imax

      k = lcell(i,j)
!
!  inner loop on nodes k.
!
      if ( k /= 0 ) then

1       continue

        delx = xp - x(k)
        dely = yp - y(k)
        ds = delx * delx + dely * dely
        rs = rsq(k)

        if ( ds == 0.0e+00 ) then
          q = f(k)
          qx = a(4,k)
          qy = a(5,k)
          ier = 0
          return
        end if

        if ( ds < rs ) then

          rds = rs * ds
          rd = sqrt ( rds )
          w = ( rs + ds - rd - rd ) / rds
          t = 2.0e+00 * ( rd - rs ) / ( ds * rds )
          wx = delx * t
          wy = dely * t
          qkx = 2.0e+00 * a(1,k) * delx + a(2,k) * dely
          qky = a(2,k) * delx + 2.0e+00 * a(3,k) * dely
          qk = ( qkx * delx + qky * dely ) / 2.0e+00
          qkx = qkx + a(4,k)
          qky = qky + a(5,k)
          qk = qk + a(4,k) * delx + a(5,k) * dely + f(k)
          sw = sw + w
          swx = swx + wx
          swy = swy + wy
          swq = swq + w * qk
          swqx = swqx + wx * qk + w * qkx
          swqy = swqy + wy * qk + w * qky

        end if

        kp = k
        k = lnext(kp)

        if ( k /= kp ) then
          go to 1
        end if

      end if

    end do

  end do
!
!  sw = 0 if and only if p is not within the radius r(k) for any node k.
!
  if ( sw /= 0.0e+00 ) then

    q = swq / sw
    sws = sw * sw
    qx = ( swqx * sw - swq * swx ) / sws
    qy = ( swqy * sw - swq * swy ) / sws
    ier = 0

  else

    q = 0.0e+00
    qx = 0.0e+00
    qy = 0.0e+00
    ier = 2

  end if

  return
end
subroutine qshep2 ( n, x, y, f, nq, nw, nr, lcell, lnext, xmin, &
  ymin, dx, dy, rmax, rsq, a, ier )
!
!***********************************************************************
!
!! qshep2 computes an interpolant to scattered data in the plane.
!
!
!  discussion:
!
!    qshep2 computes a set of parameters a and rsq defining a smooth, 
!    once continuously differentiable, bi-variate function q(x,y) which 
!    interpolates given data values f at scattered nodes (x,y).  
!
!    the interpolant function q(x,y) may be evaluated at an arbitrary point 
!    by passing the parameters a and rsq to the function qs2val.  the
!    first derivatives dqdx(x,y) and dqdy(x,y) may be evaluated by 
!    subroutine qs2grd.
!
!    the interpolation scheme is a modified quadratic shepard method:
!
!      q = ( w(1) * q(1) + w(2) * q(2) + .. + w(n) * q(n) ) 
!        / ( w(1)        + w(2)        + .. + w(n) )
!
!    for bivariate functions w(k) and q(k).  the nodal functions are given by
!
!      q(k)(x,y) = 
!          f(k)
!        + a(1,k) * ( x - x(k) )**2 
!        + a(2,k) * ( x - x(k) ) * ( y - y(k) )
!        + a(3,k) * ( y - y(k) )**2 
!        + a(4,k) * ( x - x(k) )
!        + a(5,k) * ( y - y(k) ).
!
!    thus, q(k) is a quadratic function which interpolates the
!    data value at node k.  its coefficients a(*,k) are obtained
!    by a weighted least squares fit to the closest nq data
!    points with weights similar to w(k).  note that the radius
!    of influence for the least squares fit is fixed for each
!    k, but varies with k.
!
!    the weights are taken to be
!
!      w(k)(x,y) = ( (r(k)-d(k))+ / r(k) * d(k) )**2
!
!    where (r(k)-d(k))+ = 0 if r(k) <= d(k) and d(k)(x,y) is
!    the euclidean distance between (x,y) and (x(k),y(k)).  the
!    radius of influence r(k) varies with k and is chosen so
!    that nw nodes are within the radius.  note that w(k) is
!    not defined at node (x(k),y(k)), but q(x,y) has limit f(k)
!    as (x,y) approaches (x(k),y(k)).
!
!  author:
!
!    robert renka,
!    university of north texas
!
!  parameters:
!
!    input, integer n, the number of nodes (x,y) at which data values
!    are given.  n must be at least 6.
!
!    input, real x(n), y(n), the coordinates of the nodes at which
!    data has been supplied.
!
!    input, real f(n), the data values.
!
!    input, integer nq, the number of data points to be used in the least
!    squares fit for coefficients defining the nodal functions q(k).  
!    a highly recommended value is nq = 13.  
!    nq must be at least 5, and no greater than the minimum of 40 and n-1.
!
!    input, integer nw, the number of nodes within (and defining) the radii
!    of influence r(k) which enter into the weights w(k).  for n 
!    sufficiently large, a recommended value is nw = 19.   nw must be
!    at least 1, and no greater than the minimum of 40 and n-1.
!
!    input, integer nr, the number of rows and columns in the cell grid 
!    defined in subroutine store2.  a rectangle containing the nodes 
!    is partitioned into cells in order to increase search efficiency.  
!    nr = sqrt(n/3) is recommended.  nr must be at least 1.
!
!    output, integer lcell(nr,nr), array of nodal indices associated
!    with cells.
!
!    output, integer lnext(n), contains next-node indices ( or their 
!    negatives ).
!
!    output, real xmin, ymin, dx, dy, the minimum nodal x, y coordinates,
!    and the x, y dimensions of a cell.
!
!    output, real rmax, the square root of the largest element in rsq,
!    the maximum radius of influence.
!
!    output, real rsq(n), the squared radii which enter into the weights 
!    defining the interpolant q.
!
!    output, real a(5,n), the coefficients for the nodal functions 
!    defining the interpolant q.
!
!    output, integer ier, error indicator.
!    0, if no errors were encountered.
!    1, if n, nq, nw, or nr is out of range.
!    2, if duplicate nodes were encountered.
!    3, if all nodes are collinear.
!
!  local parameters:
!
! av =        root-mean-square distance between k and the
!             nodes in the least squares system (unless
!             additional nodes are introduced for stabil-
!             ity).      the first 3 columns of the matrix
!             are scaled by 1/avsq, the last 2 by 1/av
! avsq =      av*av
! b =         transpose of the augmented regression matrix
! c =         first component of the plane rotation used to
!             zero the lower triangle of b**t -- computed
!             by subroutine givens
! ddx,ddy =   local variables for dx and dy
! dmin =      minimum of the magnitudes of the diagonal
!             elements of the regression matrix after
!             zeros are introduced below the diagonal
! dtol =      tolerance for detecting an ill-conditioned
!             system.  the system is accepted when dmin
!             >= dtol
! fk =        data value at node k -- f(k)
! i =         index for a, b, and npts
! ib =        do-loop index for back solve
! ierr =      error flag for the call to store2
! irow =      row index for b
! j =         index for a and b
! jp1 =       j+1
! k =         nodal function index and column index for a
! lmax =      maximum number of npts elements (must be con-
!             sistent with the dimension statement above)
! lnp =       current length of npts
! neq =       number of equations in the least squares fit
! nn,nnr =    local copies of n and nr
! np =        npts element
! npts =      array containing the indices of a sequence of
!             nodes to be used in the least squares fit
!             or to compute rsq.  the nodes are ordered
!             by distance from k and the last element
!             (usually indexed by lnp) is used only to
!             determine rq, or rsq(k) if nw > nq
! nqwmax =    max(nq,nw)
! rq =        radius of influence which enters into the
!             weights for q(k) (see subroutine setup2)
! rs =        squared distance between k and npts(lnp) --
!             used to compute rq and rsq(k)
! rsmx =      maximum rsq element encountered
! rsold =     squared distance between k and npts(lnp-1) --
!             used to compute a relative change in rs
!             between succeeding npts elements
! rtol =      tolerance for detecting a sufficiently large
!             relative change in rs.  if the change is
!             not greater than rtol, the nodes are
!             treated as being the same distance from k
! rws =       current value of rsq(k)
! s =         second component of the plane givens rotation
! sf =        marquardt stabilization factor used to damp
!             out the first 3 solution components (second
!             partials of the quadratic) when the system
!             is ill-conditioned.  as sf increases, the
!             fitting function approaches a linear
! sum2 =      sum of squared euclidean distances between
!             node k and the nodes used in the least
!             squares fit (unless additional nodes are
!             added for stability)
! t =         temporary variable for accumulating a scalar
!             product in the back solve
! xk,yk =     coordinates of node k -- x(k), y(k)
! xmn,ymn =   local variables for xmin and ymin
!
  implicit none
!
  integer n
  integer nr
!
  real*8 a(5,n)
  real*8 av
  real*8 avsq
  real*8 b(6,6)
  real*8 c
  real*8 ddx
  real*8 ddy
  real*8 dmin
  real*8, parameter :: dtol = 0.01e+00 !0.01e+00
  real*8 dx
  real*8 dy
  real*8 f(n)
  real*8 fk
  integer i
  integer ier
  integer ierr
  integer irow
  integer j
  integer jp1
  integer k
  integer lcell(nr,nr)
  integer lmax
  integer lnext(n)
  integer lnp
  integer neq
  integer nn
  integer nnr
  integer np
  integer npts(40)
  integer nq
  integer nqwmax
  integer nw
  real*8 rmax
  real*8 rq
  real*8 rs
  real*8 rsmx
  real*8 rsold
  real*8 rsq(n)
  real*8, parameter :: rtol = 1.0e-05
  real*8 rws
  real*8 s
  real*8, parameter :: sf = 1.0e+00
  real*8 sum2
  real*8 t
  real*8 x(n)
  real*8 xk
  real*8 xmin
  real*8 xmn
  real*8 y(n)
  real*8 yk
  real*8 ymin
  real*8 ymn
!
  nn = n
  nnr = nr
  nqwmax = max ( nq, nw )
  lmax = min ( 40, n-1 )

  if ( 5 > nq ) then
    ier = 1
    return
  else if ( 1 > nw ) then
    ier = 1
    return
  else if ( nqwmax > lmax ) then
    ier = 1
    return
  else if ( nr < 1 ) then
    ier = 1
    return
  end if
!
!  create the cell data structure, and initialize rsmx.
!
  call store2 ( nn, x, y, nnr, lcell, lnext, xmn, ymn, ddx, ddy, ierr )

  if ( ierr /= 0 ) then
    xmin = xmn
    ymin = ymn
    dx = ddx
    dy = ddy
    ier = 3
    return
  end if

  rsmx = 0.0e+00
!
!  outer loop on node k.
!
  do k = 1, nn

    xk = x(k)
    yk = y(k)
    fk = f(k)
!
!  mark node k to exclude it from the search for nearest neighbors.
!
    lnext(k) = - lnext(k)
!
!  initialize for loop on npts.
!
    rs = 0.0e+00
    sum2 = 0.0e+00
    rws = 0.0e+00
    rq = 0.0e+00
    lnp = 0
!
!  compute npts, lnp, rws, neq, rq, and avsq.
!
1   continue

    sum2 = sum2 + rs

    if ( lnp == lmax ) then
      go to 3
    end if

    lnp = lnp + 1
    rsold = rs

    call getnp2 ( xk, yk, x, y, nnr, lcell, lnext, xmn, ymn, ddx, ddy, np, rs )

    if ( rs == 0.0e+00 ) then
      ier = 2
      return
    end if

    npts(lnp) = np

    if ( ( rs - rsold ) / rs < rtol ) then
      go to 1
    end if

    if ( rws == 0.0e+00 .and. lnp > nw ) then
      rws = rs
    end if
!
!  rq = 0 (not yet computed) and lnp > nq.     
!
!  rq = sqrt(rs) is sufficiently large to (strictly) include nq nodes.  
!
!  the least squares fit will include neq = lnp - 1 equations for 
!  5 <= nq <= neq < lmax <= n-1.
!
    if ( rq == 0.0e+00 .and. lnp > nq ) then
      neq = lnp - 1
      rq = sqrt ( rs )
      avsq = sum2 / real ( neq )
    end if

    if ( lnp > nqwmax ) then
      go to 4
    else
      go to 1
    end if
!
!  all lmax nodes are included in npts.   rws and/or rq**2 is
!  (arbitrarily) taken to be 10 percent larger than the
!  distance rs to the last node included.
!
3   continue

    if ( rws == 0.0e+00 ) then
      rws = 1.1e+00 * rs
    end if

    if ( rq == 0.0e+00 ) then
      neq = lmax
      rq = sqrt ( 1.1e+00 * rs )
      avsq = sum2 / real ( neq )
    end if

4   continue
!
!  store rsq(k), update rsmx if necessary, and compute av.
!
    rsq(k) = rws
    rsmx = max ( rsmx, rws )
    av = sqrt ( avsq )
!
!  set up the augmented regression matrix (transposed) as the
!  columns of b, and zero out the lower triangle (upper
!  triangle of b) with givens rotations -- qr decomposition
!  with orthogonal matrix q not stored.
!
    i = 0

5   continue

    i = i + 1
    np = npts(i)
    irow = min ( i, 6 )

    call setup2 ( xk, yk, fk, x(np), y(np), f(np), av, avsq, rq, b(1,irow) )

    if ( i == 1 ) then
      go to 5
    end if

    do j = 1, irow-1
      jp1 = j + 1
      call givens ( b(j,j), b(j,irow), c, s )
      call rotate ( 6-j, c, s, b(jp1,j), b(jp1,irow) )
    end do

    if ( i < neq ) then
      go to 5
    end if
!
!  test the system for ill-conditioning.
!
    dmin =  min ( abs ( b(1,1) ), abs ( b(2,2) ), abs ( b(3,3) ), &
      abs ( b(4,4) ), abs ( b(5,5) ) )

    if ( dmin * rq >= dtol ) then
      go to 13
    end if

    if ( neq == lmax ) then
      go to 10
    end if
!
!  increase rq and add another equation to the system to improve conditioning.  
!  the number of npts elements is also increased if necessary.
!
7   continue

    rsold = rs
    neq = neq + 1

    if ( neq == lmax ) then
      go to 9
    end if
!
!   neq < lnp.
!
    if ( neq /= lnp ) then
      np = npts(neq+1)
      rs = ( x(np) - xk )**2 + ( y(np) - yk )**2
      if ( ( rs - rsold ) / rs < rtol ) then
        go to 7
      end if
      rq = sqrt(rs)
      go to 5
    end if
!
!  add an element to npts.
!
    lnp = lnp + 1
    call getnp2 ( xk, yk, x, y, nnr, lcell, lnext, xmn, ymn, ddx, ddy, np, rs )

    if ( np == 0 ) then
      ier = 2
      return
    end if

    npts(lnp) = np

    if ( ( rs - rsold ) / rs < rtol ) then
      go to 7
    end if

    rq = sqrt ( rs )
    go to 5

9   continue

    rq = sqrt ( 1.1e+00 * rs )
    go to 5
!
!  stabilize the system by damping second partials.  add multiples of the 
!  first three unit vectors to the first three equations.
!
10  continue

    do i = 1, 3

      b(i,6) = sf

      do j = i+1, 6
        b(j,6) = 0.0e+00
      end do

      do j = i, 5
        jp1 = j + 1
        call givens ( b(j,j), b(j,6), c, s )
        call rotate ( 6-j, c, s, b(jp1,j), b(jp1,6) )
      end do

    end do
!
!  test the stabilized system for ill-conditioning.
!
    dmin = min ( abs ( b(1,1) ), abs ( b(2,2) ), abs ( b(3,3) ), &
      abs ( b(4,4) ), abs ( b(5,5) ) )

!    if ( dmin * rq < dtol ) then
!      xmin = xmn
!      ymin = ymn
!      dx = ddx
!      dy = ddy
!      ier = 3
!      return
!    end if
!
!  solve the 5 by 5 triangular system for the coefficients.
!
13  continue

    do i = 5, 1, -1

      t = 0.0e+00

      do j = i+1, 5
        t = t + b(j,i) * a(j,k)
      end do

      a(i,k) = ( b(6,i) - t ) / b(i,i)

    end do
!
!  scale the coefficients to adjust for the column scaling.
!
    do i = 1, 3
      a(i,k) = a(i,k) / avsq
    end do

    a(4,k) = a(4,k) / av
    a(5,k) = a(5,k) / av
!
!  unmark k and the elements of npts.
!
    lnext(k) = - lnext(k)

    do i = 1, lnp
      np = npts(i)
      lnext(np) = - lnext(np)
    end do

  end do
!
!  no errors encountered.
!
  xmin = xmn
  ymin = ymn
  dx = ddx
  dy = ddy
  rmax = sqrt ( rsmx )
  ier = 0

  return
end
function qs2val ( px, py, n, x, y, f, nr, lcell, lnext, xmin, &
  ymin, dx, dy, rmax, rsq, a )
!
!***********************************************************************
!
!! qs2val evaluates the interpolant function at a point.
!
!
!  discussion:
!
!    qs2val returns the value q(px,py) where q is the weighted sum of 
!    quadratic nodal functions defined by qshep2.  if the spatial 
!    derivatives of q are also desired, call qs2grd instead.
!
!    input parameters are not altered by this function.  the
!    parameters other than px and py should be input unaltered
!    from their values on output from qshep2.  this function
!    should not be called if a nonzero error flag was returned
!    by qshep2.
!
!  modified:
!
!    10 july 1999
!
!  author:
!
!    robert renka,
!    university of north texas
!
!  parameters:
!
!    input, real px, py, the (x,y) coordinates of the point p at
!    which q is to be evaluated.
!
!    input, integer n, the number of nodes and data values to be 
!    interpolated.  n must be at least 6.
!
!    input, real x(n), y(n), the coordinates of the nodes at which
!    data has been supplied.
!
!    input, real f(n), the data values at the nodes.
!
!    input, integer nr, the number of rows and columns in the cell grid.
!    refer to subroutine store2.  nr must be at least 1.
!
!    input, integer lcell(nr,nr), the array of nodal indices associated
!    with cells.  refer to store2.
!
!    input, integer lnext(n), the next-node indices.  refer to store2.
!
!    input, real xmin, ymin, dx, dy, the minimum nodal x, y coordinates,
!    and the x, y dimensions of a cell.  computed by qshep2.
!
!    input, real rmax, the square root of the largest element in rsq,
!    the maximum radius of influence.  computed by qshep2.
!
!    input, real rsq(n), the squared radii which enter into the weights 
!    defining the interpolant q.  computed by qshep2.
!
!    input, real a(5,n), the coefficients for the nodal functions 
!    defining the interpolant q.  computed by qshep2.
!
!    output, real qs2val, the interpolated function value at (px,py).
!
  implicit none
!
  integer n
  integer nr
!
  real*8 a(5,n)
  real*8 delx
  real*8 dely
  real*8 dx
  real*8 dy
  real*8 f(n)
  integer i
  integer imax
  integer imin
  integer j
  integer jmax
  integer jmin
  real*8 ds
  integer k
  integer kp
  integer lcell(nr,nr)
  integer lnext(n)
  real*8 px
  real*8 py
  real*8 qs2val
  real*8 rd
  real*8 rds
  real*8 rmax
  real*8 rs
  real*8 rsq(n)
  real*8 sw
  real*8 swq
  real*8 w
  real*8 x(n)
  real*8 xmin
  real*8 xp
  real*8 y(n)
  real*8 ymin
  real*8 yp
!
  xp = px
  yp = py
  qs2val = 0.0e+00

  if ( n < 6  ) then
    return  
  else if ( nr < 1  ) then
    return
  else if ( dx <= 0.0e+00 ) then
    return
  else if ( dy <= 0.0e+00 ) then
    return
  else if ( rmax < 0.0e+00 ) then
    return
  end if
!
!  set imin, imax, jmin, and jmax to cell indices defining
!  the range of the search for nodes whose radii include
!  p.  the cells which must be searched are those inter-
!  sected by (or contained in) a circle of radius rmax
!  centered at p.
!
  imin = int ( ( xp - xmin - rmax ) / dx ) + 1
  imin = max ( imin, 1 )

  imax = int ( ( xp - xmin + rmax ) / dx ) + 1
  imax = min ( imax, nr )

  jmin = int ( ( yp - ymin - rmax ) / dy ) + 1
  jmin = max ( jmin, 1 )

  jmax = int ( ( yp - ymin + rmax ) / dy ) + 1
  jmax = min ( jmax, nr )
!
!  test for no cells within the circle of radius rmax.
!
  if ( imin > imax .or. jmin > jmax ) then
    qs2val = 0.0e+00
    return
  end if
!
!  accumulate weight values in sw and weighted nodal function
!  values in swq.  the weights are w(k) = ((r-d)+/(r*d))**2
!  for r**2 = rsq(k) and d = distance between p and node k.
!
  sw = 0.0e+00
  swq = 0.0e+00

  do j = jmin, jmax

    do i = imin, imax

      k = lcell(i,j)

      if ( k /= 0 ) then

1       continue

        delx = xp - x(k)
        dely = yp - y(k)
        ds = delx * delx + dely * dely
        rs = rsq(k)

        if ( ds < rs ) then

          if ( ds == 0.0e+00 ) then
            qs2val = f(k)
            return
          end if

          rds = rs * ds
          rd = sqrt ( rds )
          w = ( rs + ds - rd - rd ) / rds
          sw = sw + w

          swq = swq + w * ( f(k) + a(1,k) * delx * delx &
            + a(2,k) * delx * dely + a(3,k) * dely * dely &
            + a(4,k) * delx + a(5,k) * dely )

        end if

        kp = k
        k = lnext(kp)

        if ( k /= kp ) then
          go to 1
        end if

      end if

    end do

  end do
!
!  sw = 0 if and only if p is not within the radius r(k) for any node k.
!
  if ( sw == 0.0e+00 ) then
    qs2val = 0.0e+00
  else
    qs2val = swq / sw
  end if

  return
end
subroutine rotate ( n, c, s, x, y )
!
!***********************************************************************
!
!! rotate applies a givens rotation.
!
!
!  discussion:
!
!    the rotation has the form:
!
!      (   c  s )
!      ( - s  c )
!
!    and is essentially applied to a 2 by n matrix:
!
!      ( x(1) x(2) ... x(n) )
!      ( y(1) y(2) ... y(n) )
!
!  modified:
!
!    28 june 1999
!
!  parameters:
!
!    input, integer n, the dimension of the vectors.
!
!    input, real c, s, the cosine and sine entries of the givens
!    rotation matrix.  these may be determined by subroutine givens.
!
!    input/output, real x(n), y(n), the rotated vectors. 
!
  implicit none
!
  integer n
!
  real*8 c
  integer i
  real*8 s
  real*8 x(n)
  real*8 xi
  real*8 y(n)
  real*8 yi
!
  if ( n <= 0 ) then
    return
  else if ( c == 1.0e+00 .and. s == 0.0e+00 ) then
    return
  end if

  do i = 1, n
    xi = x(i)
    yi = y(i)
    x(i) =   c * xi + s * yi
    y(i) = - s * xi + c * yi
  end do

  return
end
subroutine setup2 ( xk, yk, fk, xi, yi, fi, s1, s2, r, row )
!
!***********************************************************************
!
!! setup2 sets up a row of the least squares regression matrix.
!
!
!  discussion:
!
!    setup2 sets up the i-th row of an augmented regression matrix for 
!    a weighted least-squares fit of a quadratic function q(x,y) to a set 
!    of data values f, where q(xk,yk) = fk.  
!
!    the first 3 columns are quadratic terms, and are scaled by 1/s2.
!    the fourth and fifth columns represent linear terms, and are scaled 
!    by 1/s1.  
!
!    if d = 0, or d >= r, the weight is
!      0,
!    else if d < r, the weight is 
!      (r-d)/(r*d), 
!    where d is the distance between nodes i and k, and r is a maximum
!    influence distance.
!
!  modified:
!
!    05 july 1999
!
!  author:
!
!    robert renka,
!    university of north texas
!
!  parameters:
!
!    input, real xk, yk, fk, the coordinates and value of the data
!    at data node k.
!
!    input, real xi, yi, fi, the coorindates and value of the data
!    at data node i.
!
!    input, real s1, s2, reciprocals of the scale factors.
!
!    input, real r, the maximum radius of influence about node k.
!
!    output, real row(6), a row of the augmented regression matrix.
!
  implicit none
!
  real*8 d
  real*8 dx
  real*8 dy
  real*8 fi
  real*8 fk
  integer i
  real*8 r
  real*8 row(6)
  real*8 s1
  real*8 s2
  real*8 w
  real*8 xi
  real*8 xk
  real*8 yi
  real*8 yk
!
  dx = xi - xk
  dy = yi - yk

  d = sqrt ( dx * dx + dy * dy )

  if ( d <= 0.0e+00 .or. d >= r ) then

    row(1:6) = 0.0e+00

  else

    w = ( r - d ) / r / d

    row(1) = dx * dx * w / s2
    row(2) = dx * dy * w / s2
    row(3) = dy * dy * w / s2
    row(4) = dx * w / s1
    row(5) = dy * w / s1
    row(6) = ( fi - fk ) * w

  end if

  return
end
subroutine store2 ( n, x, y, nr, lcell, lnext, xmin, ymin, dx, dy, ier )
!
!***********************************************************************
!
!! store2 creates a cell data structure for the scattered data.
!
!
!  discussion:
!
!    store2 is given a set of n arbitrarily distributed nodes in the 
!    plane and creates a data structure for a cell-based method of 
!    solving closest-point problems.  the smallest rectangle containing 
!    all the nodes is partitioned into an nr by nr uniform grid of cells, 
!    and nodes are associated with cells.      
!
!    in particular, the data structure stores the indices of the nodes 
!    contained in each cell.  for a uniform random distribution of nodes, 
!    the nearest node to an arbitrary point can be determined in constant
!    expected time.
!
!  modified:
!
!    05 july 1999
!
!  author:
!
!    robert renka
!    university of north texas
!
!  parameters:
!
!    input, integer n, the number of data nodes.  n must be at least 2.
!
!    input, real x(n), y(n), the coordinates of the data nodes.
!
!    input, integer nr, the number of rows and columns in the grid.  the
!    cell density, or average number of data nodes per cell, is
!      d = n / ( nr * nr ).
!    a recommended value, based on empirical evidence, is 
!      d = 3. 
!    hence, the corresponding value of nr is recommended to be about
!      nr = sqrt ( n / 3 ).  
!    nr must be at least 1.
!
!    output, integer lcell(nr,nr), an array set up so that lcell(i,j)
!    contains the index (for x and y) of the first data node (that is, the
!    data node with smallest index) in the (i,j) cell.  lcell(i,j) will be 0 if 
!    no data nodes are contained in the (i,j) cell.  the upper right corner of 
!    the (i,j) cell has coordinates 
!      ( xmin + i * dx, ymin + j * dy ).
!
!    output, integer lnext(n), an array of next-node indices.  lnext(k)
!    contains the index of the next node in the cell which contains node k, 
!    or lnext(k) = k if k is the last node in the cell.
!    the data nodes contained in a cell are ordered by their indices.
!    if, for example, cell (i,j) contains nodes 2, 3, and 5 and no others, 
!    then:
!
!      lcell(i,j) = 2, (index of the first data node)
!
!      lnext(2) = 3, 
!      lnext(3) = 5,
!      lnext(5) = 5.
!
!    output, real xmin, ymin, the x, y coordinates of the lower left
!    corner of the rectangle defined by the data nodes.  the upper right 
!    corner is ( xmax, ymax ), where
!      xmax = xmin + nr * dx,
!      ymax = ymin + nr * dy.
!
!    output, real dx, dy, the x and y dimensions of the individual cells.
!      dx = ( xmax - xmin ) / nr
!      dy = ( ymax - ymin ) / nr,
!    where xmin, xmax, ymin and ymax are the extrema of x and y.
!
!    output, integer ier, an error indicator.
!    0, if no errors were encountered.
!    1, if n < 2 or nr < 1.
!    2, if dx = 0 or dy = 0.
!
  implicit none
!
  integer n
  integer nr
!
  real*8 dx
  real*8 dy
  integer i
  integer ier
  integer j
  integer k
  integer l
  integer lcell(nr,nr)
  integer lnext(n)
  real*8 x(n)
  real*8 xmax
  real*8 xmin
  real*8 y(n)
  real*8 ymax
  real*8 ymin
!
  ier = 0

  if ( n < 2 ) then
    ier = 1
    return
  end if

  if ( nr < 1 ) then
    ier = 1
    return
  end if
!
!  compute the dimensions of the (x,y) rectangle containing all the data nodes.
!
  xmin = minval ( x(1:n) )
  xmax = maxval ( x(1:n) )
  ymin = minval ( y(1:n) )
  ymax = maxval ( y(1:n) )
!
!  compute the dimensions of a single cell.
!
  dx = ( xmax - xmin ) / real ( nr )
  dy = ( ymax - ymin ) / real ( nr )
!
!  test for zero area.
!
  if ( dx == 0.0e+00 .or. dy == 0.0e+00 ) then
    ier = 2
    return
  end if
!
!  initialize lcell.
!
  do j = 1, nr
    do i = 1, nr
      lcell(i,j) = 0
    end do
  end do
!
!  loop on nodes, storing indices in lcell and lnext.
!
  do k = n, 1, -1

    i = int ( ( x(k) - xmin ) / dx ) + 1
    i = min ( i, nr )

    j = int ( ( y(k) - ymin ) / dy ) + 1
    j = min ( j, nr )

    l = lcell(i,j)

    if ( l /= 0 ) then
      lnext(k) = l
    else
      lnext(k) = k
    end if

    lcell(i,j) = k

  end do

  return
end
subroutine timestamp ( )
!
!*******************************************************************************
!
!! timestamp prints the current ymdhms date as a time stamp.
!
!
!  example:
!
!    may 31 2001   9:45:54.872 am
!
!  modified:
!
!    31 may 2001
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    none
!
  implicit none
!
  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'january  ', 'february ', 'march    ', 'april    ', &
    'may      ', 'june     ', 'july     ', 'august   ', &
    'september', 'october  ', 'november ', 'december ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
  character ( len = 5 ) zone
!
  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'am'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'noon'
    else
      ampm = 'pm'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'pm'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'midnight'
      else
        ampm = 'am'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

!ws**********************************************
      subroutine mpi_test
      implicit real*8 (b-h,o-z)
      real*8 www
      integer mm
      dimension fxz(3,3)
      include 'mpif.h'
!       integer status(mpi_status_size)
      if(nrank.eq.0) then
      do i=1,3
      do j=1,3
      fxz(i,j)=i+j
      enddo
      enddo
      www=1.
      call mpi_send(www, 1, mpi_double_precision, 1, 1,  &
		               mpi_comm_world,ierror )
      call mpi_send(fxz(1:2,:), 2*3, mpi_double_precision, 1, 2,  &
		               mpi_comm_world,ierror )
      endif
      if(nrank.eq.1) then    
      call mpi_recv(www, 1, mpi_double_precision, 0, 1,  &
		               mpi_comm_world,status,ierror )
      call mpi_recv(fxz(1:2,:), 2*3, mpi_double_precision, 0, 2,  &
		               mpi_comm_world,status,ierror )
      endif
      write(*,*) nrank,www,fxz(1,1),fxz(2,2)

      return
      end

!ws:bndry8&3
!ws****************************************************************
      subroutine map_xz2st(fxz,fst,mm)
      use declare
      implicit real*8 (b-h,o-z)
      integer mm,im
      dimension fxz(mx,mz,my,mm),fst(n2th+5,mps4:mps,my,mm)
      include 'mpif.h'

!       integer status(mpi_status_size)

!
      do 2 js=mpsa-2,mpsa-nda,-1
      if(mts_nrk(js,nrankxz).gt.0) then
      do 21 im=1,mm
      do 21 jy=1,my      
      call interp_xz2ps(fxz(ix_first:ix_last,iz_first:iz_last,jy,im), &
                        xx(ix_first:ix_last),zz(iz_first:iz_last),ix_last-ix_first+1,iz_last-iz_first+1, &
                       fst(itsmin(js,nrankxz):itsmax(js,nrankxz),js,jy,im), &
                        xxs(itsmin(js,nrankxz):itsmax(js,nrankxz),js),zzs(itsmin(js,nrankxz):itsmax(js,nrankxz),js),mts_nrk(js,nrankxz))
      
   21 continue    
!!mpi
!	write(*,*) nrankxz,js,"map1 done" 
         do irecv=1,mrkb
         do isend=1,nsend(irecv,js)
         if(nrankxz.eq.nranksend(irecv,js,isend) .and. nrankxz.ne.nrkb(irecv)) then
         ltmin=ittransmin(irecv,js,isend)
         ltmax=ittransmax(irecv,js,isend)
         if (itbmin(nrkb(irecv))+n2th .le. itsmax(js,nranksend(irecv,js,isend))) then
         ltmin=ltmin+n2th
         ltmax=ltmax+n2th
         endif
         if (itbmax(nrkb(irecv))-n2th .ge. itsmin(js,nranksend(irecv,js,isend))) then
         ltmin=ltmin-n2th
         ltmax=ltmax-n2th
         endif
         do lt=ltmin,ltmax
         call mpi_send(fst(lt,js,1:my,1:mm),my*mm, mpi_double_precision, nrkb(irecv)+nrky(nrank)*nprxz, isend,  &
		               mpi_comm_world,ierror )
         enddo
         endif
         enddo
         enddo

!	write(*,*) nrank,js,"mapsend done" 
      endif

      if(mb_nrk(nrankxz).gt.0) then  
         do irecv=1,mrkb
         do isend=1,nsend(irecv,js)
         if(nrankxz.eq.nrkb(irecv) .and. nrankxz.ne.nranksend(irecv,js,isend)) then
         ltmin=ittransmin(irecv,js,isend)
         ltmax=ittransmax(irecv,js,isend)
         do lt=ltmin,ltmax
         call mpi_recv(fst(lt,js,1:my,1:mm),my*mm, mpi_double_precision, nranksend(irecv,js,isend)+nrky(nrank)*nprxz, isend,  &
		               mpi_comm_world,status,ierror )
         enddo
         endif
         enddo
         enddo
!	write(*,*) nrank,js,"maprecv done" 
       endif


!!mpi
    2 continue 
      return
      end

!ws****************************************************************
      subroutine smth_st_nrk(fstsm,js,mm,kk)
      use declare
      implicit real*8 (b-h,o-z)
      integer mm,js,kk,ltmin,ltmax,im
      dimension fstsm(n2th+5,mps4:mps,my,mm),wst(n2th+5)
      include 'mpif.h'
!       integer status(mpi_status_size)

      do 11 k=1,kk
      do im=1,mm
      do jy=iy_first+2,iy_last-2
      do jt=itbmin(nrankxz)+2,itbmax(nrankxz)-2
      wst(jt)=fstsm(jt,js,jy,im)/2.+((fstsm(jt+1,js,jy,im)+fstsm(jt-1,js,jy,im))*3./16+(fstsm(jt+2,js,jy,im)+fstsm(jt-2,js,jy,im))/16.&
             +(fstsm(jt,js,jy+1,im)+fstsm(jt,js,jy-1,im))*3./16+(fstsm(jt,js,jy+2,im)+fstsm(jt,js,jy-2,im))/16.)/2.
      enddo
      
      do jt=itbmin(nrankxz)+2,itbmax(nrankxz)-2
      fstsm(jt,js,jy,im)=wst(jt)
      enddo

      enddo
      enddo

         if(inrkb(nrankxz) .lt. mrkb) then
         ltmin=itbmin(nrkb1(inrkb(nrankxz)+1))
         ltmax=itbmin(nrkb1(inrkb(nrankxz)+1))+1
         call mpi_send(fstsm(ltmin:ltmax,js,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(inrkb(nrankxz)+1)+nrky(nrank)*nprxz, 1,  &
		               mpi_comm_world,ierror )
         endif
   
         if(inrkb(nrankxz) .gt. 1) then
         ltmin=itbmin(nrankxz)
         ltmax=itbmin(nrankxz)+1
         call mpi_recv(fstsm(ltmin:ltmax,js,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(inrkb(nrankxz)-1)+nrky(nrank)*nprxz, 1,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrankxz) .eq. mrkb) then
         ltmin=itbmin(nrkb1(1))+n2th
         ltmax=itbmin(nrkb1(1))+n2th+1
         call mpi_send(fstsm(ltmin:ltmax,js,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(1)+nrky(nrank)*nprxz, 1,  &
		               mpi_comm_world,ierror )
         endif
         
         if(inrkb(nrankxz) .eq. 1) then
         ltmin=itbmin(nrankxz)
         ltmax=itbmin(nrankxz)+1
         call mpi_recv(fstsm(ltmin:ltmax,js,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(mrkb)+nrky(nrank)*nprxz, 1,  &
		               mpi_comm_world,status,ierror )
         endif

  !!ws  
         if(inrkb(nrankxz) .gt. 1) then
         ltmin=itbmax(nrkb1(inrkb(nrankxz)-1))-1
         ltmax=itbmax(nrkb1(inrkb(nrankxz)-1))
         call mpi_send(fstsm(ltmin:ltmax,js,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(inrkb(nrankxz)-1)+nrky(nrank)*nprxz, 2,  &
		               mpi_comm_world,ierror )
         endif

         if(inrkb(nrankxz) .lt. mrkb) then
         ltmin=itbmax(nrankxz)-1
         ltmax=itbmax(nrankxz)
         call mpi_recv(fstsm(ltmin:ltmax,js,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(inrkb(nrankxz)+1)+nrky(nrank)*nprxz, 2,  &
		               mpi_comm_world,status,ierror )
         endif
   
         if(inrkb(nrankxz) .eq. 1) then
         ltmin=itbmax(nrkb1(mrkb))-n2th-1
         ltmax=itbmax(nrkb1(mrkb))-n2th
         call mpi_send(fstsm(ltmin:ltmax,js,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(mrkb)+nrky(nrank)*nprxz, 2,  &
		               mpi_comm_world,status,ierror )
         endif

         if(inrkb(nrankxz) .eq. mrkb) then
         ltmin=itbmax(nrankxz)-1
         ltmax=itbmax(nrankxz)
         call mpi_recv(fstsm(ltmin:ltmax,js,1:my,:), 2*my*mm, mpi_double_precision, nrkb1(1)+nrky(nrank)*nprxz, 2,  &
		               mpi_comm_world,status,ierror )
         endif
   11 continue
      return
      end

!ws****************************************************************
      subroutine valbm_atlastgrid_v1(fxz,mm,ibnd)
      use declare
      implicit real*8 (b-h,o-z)
      integer mm,ibnd
      dimension fxz(mx,mz,my,mm),fst(n2th+5,mps4:mps,my,mm),f1s(mbm_nrk,mps4:mps)
      include 'mpif.h'

!       integer status(mpi_status_size)
!
      call map_xz2st(fxz,fst,mm)
      
      if(mb_nrk(nrankxz).gt.0) then 

      do 10 js=mpsa-2,mpsa-4,-1
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),js,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------------
      call smth_st_nrk(fst,js,mm,3)
!ws_smps--------------------
   10 continue

      if(ibnd==1) then
      is=mpsa-1
      do 31 m=1,mm    
      do 31 jy=1,my      
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))  
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0

       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0
!       fst(jt,is)=fst(jt,is-1)                   
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   31 continue 
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
      call smth_st_nrk(fst,is,mm,3)
!ws_smps--------------------
      is=mpsa
      do 41 m=1,mm    
      do 41 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))     
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0

       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0
!       fst(jt,is)=fst(jt,is-1)                   
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   41 continue
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------
       call smth_st_nrk(fst,is,mm,3)
!ws_smps--------------------
      endif !!ibnd==1
      
      if(ibnd==0) then
      is=mpsa
      do 51 m=1,mm    
      do 51 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
       fst(jt,is,jy,m)=0.           
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   51 continue

      is=mpsa-1
      do 52 m=1,mm    
      do 52 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
      call interp1d3l(fst(jt,is-3,jy,m),fst(jt,is-2,jy,m),fst(jt,is-1,jy,m),fst(jt,is+1,jy,m), &
                        ps(is-3),ps(is-2),ps(is-1),ps(is+1),ps(is),fst(jt,is,jy,m))                     
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   52 continue

      call smth_st_nrk(fst,is,mm,3)
      endif !!ibnd==0 

      do 3 m=1,mm    
      do 3 jy=1,my           
      do ikb=1,mb_nrk(nrankxz)
        i=itb_nrk(ikb,nrankxz)
        ib=ib_nrk(ikb,nrankxz)
        do js=mpsa,mpsa-3,-1     
        call interp1d3l(fst(i-2,js,jy,m),fst(i-1,js,jy,m),fst(i,js,jy,m),fst(i+1,js,jy,m), &
                        thst(i-2),thst(i-1),thst(i),thst(i+1),thb(ib),f1s(ikb,js))
        enddo
        js=mpsa
        call interp1d3l(f1s(ikb,js),f1s(ikb,js-1),f1s(ikb,js-2),f1s(ikb,js-3), &
                        ps(js),ps(js-1),ps(js-2),ps(js-3),psb(ib),fxz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,m))
!        fxz(jbx(ib),jbz(ib))=(asm_b(ib)*x1s(ib,js-1)+bsm_b(ib)*x1s(ib,js-2)+csm_b(ib)*x1s(ib,js-3)) &
!                 /(asm_b(ib)+bsm_b(ib)+csm_b(ib))
      enddo
    3 continue 
      endif !!mb_nrk(nrankxz).gt.0
           
!      if(ibnd==2) then
!      
!      do ika1=1,ma1_nrk(nrankxz)
!        i=ita1_nrk(ika1,nrankxz)
!        ia1=ia1_nrk(ika1,nrankxz)
!        do js=mpsa,mpsa-3,1
!        call interp1d3l(fst(i-2,js),fst(i-1,js),fst(i,js),fst(i+1,js), &
!                        thst(i-2),thst(i-1),thst(i),thst(i+1),tha1(ia1),x1s(ika1,js,jy))
!        enddo
!        js=mpsa
!        call interp1d3l(x1s(ika1,js,jy),x1s(ika1,js-1,jy),x1s(ika1,js-2,jy),x1s(ika1,js-3,jy), &
!                        ps(js),ps(js-1),ps(js-2),ps(js-3),psa1(ia1),fxz(jxa1_nrk(ika1,nrankxz),jza1_nrk(ika1,nrankxz)))
!!        fxz(jxa1(ia1),jza1(ia1))=(asm_a1(ia1)*x1s(ia1,js-1)+bsm_a1(ia1)*x1s(ia1,js-2)+csm_a1(ia1)*x1s(ia1,js-3)) &
!!                  /(asm_a1(ia1)+bsm_a1(ia1)+csm_a1(ia1)) 
!      enddo
!      endif
      return
      end

!ws****************************************************************
      subroutine valb8_atlastgrid_r0p1_v2(f8xz)
      use declare
      implicit real*8 (b-h,o-z)
      real*8 vx1st,vz1st,bx1st,bz1st
      integer is
      dimension f8xz(mx,mz,my,8),fst(n2th+5,mps4:mps,my,8),f1s(mbm_nrk,mps4:mps),fsxz(mbm_nrk,8) !
      include 'mpif.h'

!       integer status(mpi_status_size)

      call map_xz2st(f8xz,fst,8)
         
      if(mb_nrk(nrankxz).gt.0) then      

      do jy=1,my
      do js=mpsa-2,mpsa-4,-1
      do jt=itbmin(nrankxz),itbmax(nrankxz)
      vx1st=fst(jt,js,jy,3)
      vz1st=fst(jt,js,jy,5)
      bx1st=fst(jt,js,jy,6)
      bz1st=fst(jt,js,jy,8)
      fst(jt,js,jy,3)=vx1st*wbxr(jt,js)+vz1st*wbzr(jt,js)
      fst(jt,js,jy,5)=vx1st*wbxt(jt,js)+vz1st*wbzp(jt,js)
      fst(jt,js,jy,6)=bx1st*wbxr(jt,js)+bz1st*wbzr(jt,js)
      fst(jt,js,jy,8)=bx1st*wbxt(jt,js)+bz1st*wbzp(jt,js)
       
!      vrst=fst(jt,js,jy,3)*dcos(thst(jt))+fst(jt,js,jy,5)*dsin(thst(jt))
!      vpst=-fst(jt,js,jy,3)*dsin(thst(jt))+fst(jt,js,jy,5)*dcos(thst(jt))
!      brst=fst(jt,js,jy,6)*dcos(thst(jt))+fst(jt,js,jy,8)*dsin(thst(jt))
!      bpst=-fst(jt,js,jy,6)*dsin(thst(jt))+fst(jt,js,jy,8)*dcos(thst(jt))
!      fst(jt,js,jy,3)=vrst
!      fst(jt,js,jy,5)=vpst
!      fst(jt,js,jy,6)=brst
!      fst(jt,js,jy,8)=bpst
      enddo
      enddo
      enddo    

      do 10 js=mpsa-2,mpsa-4,-1
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),js,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------------
      call smth_st_nrk(fst,js,8,3)
!ws_smps--------------------
   10 continue

      do 61 m=1,8

      if(m.eq.2 .or. m.eq.3 .or. m.eq.6) then !!(m=2,3,6) ibnd==0
      is=mpsa 
      do 51 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
       fst(jt,is,jy,m)=0.           
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   51 continue

      is=mpsa-1 
      do 52 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
        call interp1d3l(fst(jt,is-3,jy,m),fst(jt,is-2,jy,m),fst(jt,is-1,jy,m),fst(jt,is+1,jy,m), &
                        ps(is-3),ps(is-2),ps(is-1),ps(is+1),ps(is),fst(jt,is,jy,m))                   
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   52 continue

      call smth_st_nrk(fst(:,:,:,m),is,1,3)

      else  !!(m=1,4,5,7,8,) bnd==1  
      is=mpsa-1    
      do 31 jy=1,my      
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is)) 
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0
!       fst(jt,is)=fst(jt,is-1) 
       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0.                  
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   31 continue 
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
      call smth_st_nrk(fst(:,:,:,m),is,1,3)
!ws_smps--------------------
      is=mpsa
 
      do 41 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0
!       fst(jt,is)=fst(jt,is-1) 
       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0.                  
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   41 continue
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------
      call smth_st_nrk(fst(:,:,:,m),is,1,3)
!	write(*,*) nrankxz,m,"smth done"
!ws_smps--------------------
      endif !!(m=3,6)ibnd==0
   61 continue
      
      do 3 jy=1,my           
      do ikb=1,mb_nrk(nrankxz)
        i=itb_nrk(ikb,nrankxz)
        ib=ib_nrk(ikb,nrankxz) 
        do m=1,8 
        do js=mpsa,mpsa-3,-1     
        call interp1d3l(fst(i-2,js,jy,m),fst(i-1,js,jy,m),fst(i,js,jy,m),fst(i+1,js,jy,m), &
                        thst(i-2),thst(i-1),thst(i),thst(i+1),thb(ib),f1s(ikb,js))
        enddo
        js=mpsa
        call interp1d3l(f1s(ikb,js),f1s(ikb,js-1),f1s(ikb,js-2),f1s(ikb,js-3), &
                        ps(js),ps(js-1),ps(js-2),ps(js-3),psb(ib),fsxz(ikb,m))
        enddo
!        fxz(jbx(ib),jbz(ib))=(asm_b(ib)*x1s(ib,js-1)+bsm_b(ib)*x1s(ib,js-2)+csm_b(ib)*x1s(ib,js-3)) &
!                 /(asm_b(ib)+bsm_b(ib)+csm_b(ib))
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,1)=fsxz(ikb,1)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,2)=fsxz(ikb,2)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,3)=fsxz(ikb,3)*wbrx(ib)+fsxz(ikb,5)*wbpx(ib)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,4)=fsxz(ikb,4)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,5)=fsxz(ikb,3)*wbrz(ib)+fsxz(ikb,5)*wbpz(ib)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,6)=fsxz(ikb,6)*wbrx(ib)+fsxz(ikb,8)*wbpx(ib)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,7)=fsxz(ikb,7)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,8)=fsxz(ikb,6)*wbrz(ib)+fsxz(ikb,8)*wbpz(ib)
      enddo
    3 continue 
      endif !!mb_nrk(nrankxz).gt.0
           
!      if(ibnd==2) then
!      
!      do ika1=1,ma1_nrk(nrankxz)
!        i=ita1_nrk(ika1,nrankxz)
!        ia1=ia1_nrk(ika1,nrankxz)
!        do js=mpsa,mpsa-3,1
!        call interp1d3l(fst(i-2,js),fst(i-1,js),fst(i,js),fst(i+1,js), &
!                        thst(i-2),thst(i-1),thst(i),thst(i+1),tha1(ia1),x1s(ika1,js,jy))
!        enddo
!        js=mpsa
!        call interp1d3l(x1s(ika1,js,jy),x1s(ika1,js-1,jy),x1s(ika1,js-2,jy),x1s(ika1,js-3,jy), &
!                        ps(js),ps(js-1),ps(js-2),ps(js-3),psa1(ia1),fxz(jxa1_nrk(ika1,nrankxz),jza1_nrk(ika1,nrankxz)))
!!        fxz(jxa1(ia1),jza1(ia1))=(asm_a1(ia1)*x1s(ia1,js-1)+bsm_a1(ia1)*x1s(ia1,js-2)+csm_a1(ia1)*x1s(ia1,js-3)) &
!!                  /(asm_a1(ia1)+bsm_a1(ia1)+csm_a1(ia1)) 
!      enddo
!      endif
      return
      end

!ws****************************************************************
      subroutine valb8_atlastgrid(f8xz)
      use declare
      implicit real*8 (b-h,o-z)
      real*8 vx1st,vz1st,bx1st,bz1st
      integer is
      dimension f8xz(mx,mz,my,8),fst(n2th+5,mps4:mps,my,8),f1s(mbm_nrk,mps4:mps),fsxz(mbm_nrk,8) !
      include 'mpif.h'

!       integer status(mpi_status_size)

      call map_xz2st(f8xz,fst,8)
         
      if(mb_nrk(nrankxz).gt.0) then      

      do jy=1,my
      do js=mpsa-2,mpsa-4,-1
      do jt=itbmin(nrankxz),itbmax(nrankxz)
      vx1st=fst(jt,js,jy,3)
      vz1st=fst(jt,js,jy,5)
      bx1st=fst(jt,js,jy,6)
      bz1st=fst(jt,js,jy,8)
      fst(jt,js,jy,3)=vx1st*wbxr(jt,js)+vz1st*wbzr(jt,js)
      fst(jt,js,jy,5)=vx1st*wbxt(jt,js)+vz1st*wbzp(jt,js)
      fst(jt,js,jy,6)=bx1st*wbxr(jt,js)+bz1st*wbzr(jt,js)
      fst(jt,js,jy,8)=bx1st*wbxt(jt,js)+bz1st*wbzp(jt,js)
       
!      vrst=fst(jt,js,jy,3)*dcos(thst(jt))+fst(jt,js,jy,5)*dsin(thst(jt))
!      vpst=-fst(jt,js,jy,3)*dsin(thst(jt))+fst(jt,js,jy,5)*dcos(thst(jt))
!      brst=fst(jt,js,jy,6)*dcos(thst(jt))+fst(jt,js,jy,8)*dsin(thst(jt))
!      bpst=-fst(jt,js,jy,6)*dsin(thst(jt))+fst(jt,js,jy,8)*dcos(thst(jt))
!      fst(jt,js,jy,3)=vrst
!      fst(jt,js,jy,5)=vpst
!      fst(jt,js,jy,6)=brst
!      fst(jt,js,jy,8)=bpst
      enddo
      enddo
      enddo    

      do 10 js=mpsa-2,mpsa-4,-1
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),js,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------------
      call smth_st_nrk(fst,js,8,3)
!ws_smps--------------------
   10 continue

      do 61 m=1,8

      if(lbndxfix(m)) then !!(m=3,6) ibnd==0
      is=mpsa 
      do 51 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
       fst(jt,is,jy,m)=0.           
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   51 continue

      is=mpsa-1 
      do 52 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
        !call interp1d3l(fst(jt,is-3,jy,m),fst(jt,is-2,jy,m),fst(jt,is-1,jy,m),fst(jt,is+1,jy,m), &
        !                ps(is-3),ps(is-2),ps(is-1),ps(is+1),ps(is),fst(jt,is,jy,m))                   
           call interp1d2l(fst(jt,is-2,jy,m),fst(jt,is-1,jy,m),fst(jt,is+1,jy,m), &
                        ps(is-2),ps(is-1),ps(is+1),ps(is),fst(jt,is,jy,m))    
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   52 continue

      call smth_st_nrk(fst(:,:,:,m),is,1,3)

      else  !!(m=1,2,4,5,7,8,) bnd==1  
      is=mpsa-1    
      do 31 jy=1,my      
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is)) 
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0
!       fst(jt,is)=fst(jt,is-1) 
       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0.                  
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   31 continue 
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
      call smth_st_nrk(fst(:,:,:,m),is,1,3)
!ws_smps--------------------
      is=mpsa
 
      do 41 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0
!       fst(jt,is)=fst(jt,is-1) 
       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0.                  
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   41 continue
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------
      call smth_st_nrk(fst(:,:,:,m),is,1,3)
!	write(*,*) nrankxz,m,"smth done"
!ws_smps--------------------
      endif !!(m=3,6)ibnd==0
   61 continue
      
      do 3 jy=1,my           
      do ikb=1,mb_nrk(nrankxz)
        i=itb_nrk(ikb,nrankxz)
        ib=ib_nrk(ikb,nrankxz) 
        do m=1,8 
        do js=mpsa,mpsa-3,-1     
        call interp1d3l(fst(i-2,js,jy,m),fst(i-1,js,jy,m),fst(i,js,jy,m),fst(i+1,js,jy,m), &
                        thst(i-2),thst(i-1),thst(i),thst(i+1),thb(ib),f1s(ikb,js))
        enddo
        js=mpsa
        call interp1d3l(f1s(ikb,js),f1s(ikb,js-1),f1s(ikb,js-2),f1s(ikb,js-3), &
                        ps(js),ps(js-1),ps(js-2),ps(js-3),psb(ib),fsxz(ikb,m))
        enddo
!        fxz(jbx(ib),jbz(ib))=(asm_b(ib)*x1s(ib,js-1)+bsm_b(ib)*x1s(ib,js-2)+csm_b(ib)*x1s(ib,js-3)) &
!                 /(asm_b(ib)+bsm_b(ib)+csm_b(ib))
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,1)=fsxz(ikb,1)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,2)=fsxz(ikb,2)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,3)=fsxz(ikb,3)*wbrx(ib)+fsxz(ikb,5)*wbpx(ib)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,4)=fsxz(ikb,4)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,5)=fsxz(ikb,3)*wbrz(ib)+fsxz(ikb,5)*wbpz(ib)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,6)=fsxz(ikb,6)*wbrx(ib)+fsxz(ikb,8)*wbpx(ib)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,7)=fsxz(ikb,7)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,8)=fsxz(ikb,6)*wbrz(ib)+fsxz(ikb,8)*wbpz(ib)
      enddo
    3 continue 
      endif !!mb_nrk(nrankxz).gt.0
           
!      if(ibnd==2) then
!      
!      do ika1=1,ma1_nrk(nrankxz)
!        i=ita1_nrk(ika1,nrankxz)
!        ia1=ia1_nrk(ika1,nrankxz)
!        do js=mpsa,mpsa-3,1
!        call interp1d3l(fst(i-2,js),fst(i-1,js),fst(i,js),fst(i+1,js), &
!                        thst(i-2),thst(i-1),thst(i),thst(i+1),tha1(ia1),x1s(ika1,js,jy))
!        enddo
!        js=mpsa
!        call interp1d3l(x1s(ika1,js,jy),x1s(ika1,js-1,jy),x1s(ika1,js-2,jy),x1s(ika1,js-3,jy), &
!                        ps(js),ps(js-1),ps(js-2),ps(js-3),psa1(ia1),fxz(jxa1_nrk(ika1,nrankxz),jza1_nrk(ika1,nrankxz)))
!!        fxz(jxa1(ia1),jza1(ia1))=(asm_a1(ia1)*x1s(ia1,js-1)+bsm_a1(ia1)*x1s(ia1,js-2)+csm_a1(ia1)*x1s(ia1,js-3)) &
!!                  /(asm_a1(ia1)+bsm_a1(ia1)+csm_a1(ia1)) 
!      enddo
!      endif
      return
      end

!ws****************************************************************
      subroutine valb8_atlastgrid_r0p1_v1(f8xz)
      use declare
      implicit real*8 (b-h,o-z)
      real*8 vx1st,vz1st,bx1st,bz1st
      integer is
      dimension f8xz(mx,mz,my,8),fst(n2th+5,mps4:mps,my,8),f1s(mbm_nrk,mps4:mps),fsxz(mbm_nrk,8) !
      include 'mpif.h'

!       integer status(mpi_status_size)

      call map_xz2st(f8xz,fst,8)
         
      if(mb_nrk(nrankxz).gt.0) then      

      do jy=1,my
      do js=mpsa-2,mpsa-4,-1
      do jt=itbmin(nrankxz),itbmax(nrankxz)
      vx1st=fst(jt,js,jy,3)
      vz1st=fst(jt,js,jy,5)
      bx1st=fst(jt,js,jy,6)
      bz1st=fst(jt,js,jy,8)
      fst(jt,js,jy,3)=vx1st*wbxr(jt,js)+vz1st*wbzr(jt,js)
      fst(jt,js,jy,5)=vx1st*wbxt(jt,js)+vz1st*wbzp(jt,js)
      fst(jt,js,jy,6)=bx1st*wbxr(jt,js)+bz1st*wbzr(jt,js)
      fst(jt,js,jy,8)=bx1st*wbxt(jt,js)+bz1st*wbzp(jt,js)
       
!      vrst=fst(jt,js,jy,3)*dcos(thst(jt))+fst(jt,js,jy,5)*dsin(thst(jt))
!      vpst=-fst(jt,js,jy,3)*dsin(thst(jt))+fst(jt,js,jy,5)*dcos(thst(jt))
!      brst=fst(jt,js,jy,6)*dcos(thst(jt))+fst(jt,js,jy,8)*dsin(thst(jt))
!      bpst=-fst(jt,js,jy,6)*dsin(thst(jt))+fst(jt,js,jy,8)*dcos(thst(jt))
!      fst(jt,js,jy,3)=vrst
!      fst(jt,js,jy,5)=vpst
!      fst(jt,js,jy,6)=brst
!      fst(jt,js,jy,8)=bpst
      enddo
      enddo
      enddo    

      do 10 js=mpsa-2,mpsa-4,-1
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),js,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------------
      call smth_st_nrk(fst,js,8,3)
!ws_smps--------------------
   10 continue

      do 61 m=1,8

      if(m.eq.3 .or. m.eq.6) then !!(m=3,6) ibnd==0
      is=mpsa 
      do 51 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
       fst(jt,is,jy,m)=0.           
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   51 continue

      is=mpsa-1 
      do 52 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
        call interp1d3l(fst(jt,is-3,jy,m),fst(jt,is-2,jy,m),fst(jt,is-1,jy,m),fst(jt,is+1,jy,m), &
                        ps(is-3),ps(is-2),ps(is-1),ps(is+1),ps(is),fst(jt,is,jy,m))                   
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   52 continue

      call smth_st_nrk(fst(:,:,:,m),is,1,3)

      else  !!(m=1,2,4,5,7,8,) bnd==1  
      is=mpsa-1    
      do 31 jy=1,my      
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is)) 
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0
!       fst(jt,is)=fst(jt,is-1) 
       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0.                  
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   31 continue 
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
      call smth_st_nrk(fst(:,:,:,m),is,1,3)
!ws_smps--------------------
      is=mpsa
 
      do 41 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0
!       fst(jt,is)=fst(jt,is-1) 
       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0.                  
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   41 continue
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------
      call smth_st_nrk(fst(:,:,:,m),is,1,3)
!	write(*,*) nrankxz,m,"smth done"
!ws_smps--------------------
      endif !!(m=3,6)ibnd==0
   61 continue
      
      do 3 jy=1,my           
      do ikb=1,mb_nrk(nrankxz)
        i=itb_nrk(ikb,nrankxz)
        ib=ib_nrk(ikb,nrankxz) 
        do m=1,8 
        do js=mpsa,mpsa-3,-1     
        call interp1d3l(fst(i-2,js,jy,m),fst(i-1,js,jy,m),fst(i,js,jy,m),fst(i+1,js,jy,m), &
                        thst(i-2),thst(i-1),thst(i),thst(i+1),thb(ib),f1s(ikb,js))
        enddo
        js=mpsa
        call interp1d3l(f1s(ikb,js),f1s(ikb,js-1),f1s(ikb,js-2),f1s(ikb,js-3), &
                        ps(js),ps(js-1),ps(js-2),ps(js-3),psb(ib),fsxz(ikb,m))
        enddo
!        fxz(jbx(ib),jbz(ib))=(asm_b(ib)*x1s(ib,js-1)+bsm_b(ib)*x1s(ib,js-2)+csm_b(ib)*x1s(ib,js-3)) &
!                 /(asm_b(ib)+bsm_b(ib)+csm_b(ib))
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,1)=fsxz(ikb,1)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,2)=fsxz(ikb,2)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,3)=fsxz(ikb,3)*wbrx(ib)+fsxz(ikb,5)*wbpx(ib)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,4)=fsxz(ikb,4)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,5)=fsxz(ikb,3)*wbrz(ib)+fsxz(ikb,5)*wbpz(ib)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,6)=fsxz(ikb,6)*wbrx(ib)+fsxz(ikb,8)*wbpx(ib)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,7)=fsxz(ikb,7)
        f8xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,8)=fsxz(ikb,6)*wbrz(ib)+fsxz(ikb,8)*wbpz(ib)
      enddo
    3 continue 
      endif !!mb_nrk(nrankxz).gt.0
           
!      if(ibnd==2) then
!      
!      do ika1=1,ma1_nrk(nrankxz)
!        i=ita1_nrk(ika1,nrankxz)
!        ia1=ia1_nrk(ika1,nrankxz)
!        do js=mpsa,mpsa-3,1
!        call interp1d3l(fst(i-2,js),fst(i-1,js),fst(i,js),fst(i+1,js), &
!                        thst(i-2),thst(i-1),thst(i),thst(i+1),tha1(ia1),x1s(ika1,js,jy))
!        enddo
!        js=mpsa
!        call interp1d3l(x1s(ika1,js,jy),x1s(ika1,js-1,jy),x1s(ika1,js-2,jy),x1s(ika1,js-3,jy), &
!                        ps(js),ps(js-1),ps(js-2),ps(js-3),psa1(ia1),fxz(jxa1_nrk(ika1,nrankxz),jza1_nrk(ika1,nrankxz)))
!!        fxz(jxa1(ia1),jza1(ia1))=(asm_a1(ia1)*x1s(ia1,js-1)+bsm_a1(ia1)*x1s(ia1,js-2)+csm_a1(ia1)*x1s(ia1,js-3)) &
!!                  /(asm_a1(ia1)+bsm_a1(ia1)+csm_a1(ia1)) 
!      enddo
!      endif
      return
      end

!ws****************************************************************
      subroutine valb3_atlastgrid(f3xz)
      use declare
      implicit real*8 (b-h,o-z)
      real*8 cx1st,cz1st
      integer is
      dimension f3xz(mx,mz,my,3),fst(n2th+5,mps4:mps,my,3),f1s(mbm_nrk,mps4:mps),fsxz(mbm_nrk,3)
      include 'mpif.h'

!       integer status(mpi_status_size)

      call map_xz2st(f3xz,fst,3)
         
      if(mb_nrk(nrankxz).gt.0) then      

      do jy=1,my
      do js=mpsa-2,mpsa-4,-1
      do jt=itbmin(nrankxz),itbmax(nrankxz)
      cx1st=fst(jt,js,jy,1)
      cz1st=fst(jt,js,jy,3)
      fst(jt,js,jy,1)=cx1st*wbxr(jt,js)+cz1st*wbzr(jt,js)
      fst(jt,js,jy,3)=cx1st*wbxt(jt,js)+cz1st*wbzp(jt,js)
      enddo
      enddo
      enddo    

      do 10 js=mpsa-2,mpsa-4,-1
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),js,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------------
      call smth_st_nrk(fst,js,3,3)
!ws_smps--------------------
   10 continue

      do 61 m=1,3

      if(lbndcfix(m)) then !!(m=1) ibnd==0
      is=mpsa 
      do 51 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
       fst(jt,is,jy,m)=0.           
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   51 continue

      is=mpsa-1 
      do 52 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
        call interp1d3l(fst(jt,is-3,jy,m),fst(jt,is-2,jy,m),fst(jt,is-1,jy,m),fst(jt,is+1,jy,m), &
                        ps(is-3),ps(is-2),ps(is-1),ps(is+1),ps(is),fst(jt,is,jy,m))                   
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   52 continue

      call smth_st_nrk(fst(:,:,:,m),is,1,3)

      else  !!(m=2,3) ibnd==1  
      is=mpsa-1    
      do 31 jy=1,my      
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is)) 
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0
!       fst(jt,is)=fst(jt,is-1)
       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0                   
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   31 continue 
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
      call smth_st_nrk(fst(:,:,:,m),is,1,3)
!ws_smps--------------------
      is=mpsa
 
      do 41 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0
!       fst(jt,is)=fst(jt,is-1) 
       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0                  
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   41 continue
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------
      call smth_st_nrk(fst(:,:,:,m),is,1,3)
!	write(*,*) nrankxz,m,"smth done"
!ws_smps--------------------
      endif !!(m=1)ibnd==0
   61 continue
      
      do 3 jy=1,my           
      do ikb=1,mb_nrk(nrankxz)
        i=itb_nrk(ikb,nrankxz)
        ib=ib_nrk(ikb,nrankxz) 
        do m=1,3
        do js=mpsa,mpsa-3,-1     
        call interp1d3l(fst(i-2,js,jy,m),fst(i-1,js,jy,m),fst(i,js,jy,m),fst(i+1,js,jy,m), &
                        thst(i-2),thst(i-1),thst(i),thst(i+1),thb(ib),f1s(ikb,js))
        enddo
        js=mpsa
        call interp1d3l(f1s(ikb,js),f1s(ikb,js-1),f1s(ikb,js-2),f1s(ikb,js-3), &
                        ps(js),ps(js-1),ps(js-2),ps(js-3),psb(ib),fsxz(ikb,m))
        enddo
!        fxz(jbx(ib),jbz(ib))=(asm_b(ib)*x1s(ib,js-1)+bsm_b(ib)*x1s(ib,js-2)+csm_b(ib)*x1s(ib,js-3)) &
!                 /(asm_b(ib)+bsm_b(ib)+csm_b(ib))
        f3xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,1)=fsxz(ikb,1)*wbrx(ib)+fsxz(ikb,3)*wbpx(ib)
        f3xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,2)=fsxz(ikb,2)
        f3xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,3)=fsxz(ikb,1)*wbrz(ib)+fsxz(ikb,3)*wbpz(ib)

      enddo
    3 continue 
      endif !!mb_nrk(nrankxz).gt.0
           
!      if(ibnd==2) then
!      
!      do ika1=1,ma1_nrk(nrankxz)
!        i=ita1_nrk(ika1,nrankxz)
!        ia1=ia1_nrk(ika1,nrankxz)
!        do js=mpsa,mpsa-3,1
!        call interp1d3l(fst(i-2,js),fst(i-1,js),fst(i,js),fst(i+1,js), &
!                        thst(i-2),thst(i-1),thst(i),thst(i+1),tha1(ia1),x1s(ika1,js,jy))
!        enddo
!        js=mpsa
!        call interp1d3l(x1s(ika1,js,jy),x1s(ika1,js-1,jy),x1s(ika1,js-2,jy),x1s(ika1,js-3,jy), &
!                        ps(js),ps(js-1),ps(js-2),ps(js-3),psa1(ia1),fxz(jxa1_nrk(ika1,nrankxz),jza1_nrk(ika1,nrankxz)))
!!        fxz(jxa1(ia1),jza1(ia1))=(asm_a1(ia1)*x1s(ia1,js-1)+bsm_a1(ia1)*x1s(ia1,js-2)+csm_a1(ia1)*x1s(ia1,js-3)) &
!!                  /(asm_a1(ia1)+bsm_a1(ia1)+csm_a1(ia1)) 
!      enddo
!      endif
      return
      end

!ws****************************************************************
      subroutine valb3_atlastgrid_r0p1_v1(f3xz)
      use declare
      implicit real*8 (b-h,o-z)
      real*8 cx1st,cz1st
      integer is
      dimension f3xz(mx,mz,my,3),fst(n2th+5,mps4:mps,my,3),f1s(mbm_nrk,mps4:mps),fsxz(mbm_nrk,3)
      include 'mpif.h'

!       integer status(mpi_status_size)

      call map_xz2st(f3xz,fst,3)
         
      if(mb_nrk(nrankxz).gt.0) then      

      do jy=1,my
      do js=mpsa-2,mpsa-4,-1
      do jt=itbmin(nrankxz),itbmax(nrankxz)
      cx1st=fst(jt,js,jy,1)
      cz1st=fst(jt,js,jy,3)
      fst(jt,js,jy,1)=cx1st*wbxr(jt,js)+cz1st*wbzr(jt,js)
      fst(jt,js,jy,3)=cx1st*wbxt(jt,js)+cz1st*wbzp(jt,js)
      enddo
      enddo
      enddo    

      do 10 js=mpsa-2,mpsa-4,-1
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),js,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------------
      call smth_st_nrk(fst,js,3,3)
!ws_smps--------------------
   10 continue

      do 61 m=1,3

      if(m.eq.1) then !!(m=1) ibnd==0
      is=mpsa 
      do 51 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
       fst(jt,is,jy,m)=0.           
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   51 continue

      is=mpsa-1 
      do 52 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
        call interp1d3l(fst(jt,is-3,jy,m),fst(jt,is-2,jy,m),fst(jt,is-1,jy,m),fst(jt,is+1,jy,m), &
                        ps(is-3),ps(is-2),ps(is-1),ps(is+1),ps(is),fst(jt,is,jy,m))                   
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   52 continue

      call smth_st_nrk(fst(:,:,:,m),is,1,3)

      else  !!(m=2,3) ibnd==1  
      is=mpsa-1    
      do 31 jy=1,my      
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is)) 
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0
!       fst(jt,is)=fst(jt,is-1)
       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0                   
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   31 continue 
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
      call smth_st_nrk(fst(:,:,:,m),is,1,3)
!ws_smps--------------------
      is=mpsa
 
      do 41 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0
!       fst(jt,is)=fst(jt,is-1) 
       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0                  
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   41 continue
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------
      call smth_st_nrk(fst(:,:,:,m),is,1,3)
!	write(*,*) nrankxz,m,"smth done"
!ws_smps--------------------
      endif !!(m=1)ibnd==0
   61 continue
      
      do 3 jy=1,my           
      do ikb=1,mb_nrk(nrankxz)
        i=itb_nrk(ikb,nrankxz)
        ib=ib_nrk(ikb,nrankxz) 
        do m=1,3
        do js=mpsa,mpsa-3,-1     
        call interp1d3l(fst(i-2,js,jy,m),fst(i-1,js,jy,m),fst(i,js,jy,m),fst(i+1,js,jy,m), &
                        thst(i-2),thst(i-1),thst(i),thst(i+1),thb(ib),f1s(ikb,js))
        enddo
        js=mpsa
        call interp1d3l(f1s(ikb,js),f1s(ikb,js-1),f1s(ikb,js-2),f1s(ikb,js-3), &
                        ps(js),ps(js-1),ps(js-2),ps(js-3),psb(ib),fsxz(ikb,m))
        enddo
!        fxz(jbx(ib),jbz(ib))=(asm_b(ib)*x1s(ib,js-1)+bsm_b(ib)*x1s(ib,js-2)+csm_b(ib)*x1s(ib,js-3)) &
!                 /(asm_b(ib)+bsm_b(ib)+csm_b(ib))
        f3xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,1)=fsxz(ikb,1)*wbrx(ib)+fsxz(ikb,3)*wbpx(ib)
        f3xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,2)=fsxz(ikb,2)
        f3xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,3)=fsxz(ikb,1)*wbrz(ib)+fsxz(ikb,3)*wbpz(ib)

      enddo
    3 continue 
      endif !!mb_nrk(nrankxz).gt.0
           
!      if(ibnd==2) then
!      
!      do ika1=1,ma1_nrk(nrankxz)
!        i=ita1_nrk(ika1,nrankxz)
!        ia1=ia1_nrk(ika1,nrankxz)
!        do js=mpsa,mpsa-3,1
!        call interp1d3l(fst(i-2,js),fst(i-1,js),fst(i,js),fst(i+1,js), &
!                        thst(i-2),thst(i-1),thst(i),thst(i+1),tha1(ia1),x1s(ika1,js,jy))
!        enddo
!        js=mpsa
!        call interp1d3l(x1s(ika1,js,jy),x1s(ika1,js-1,jy),x1s(ika1,js-2,jy),x1s(ika1,js-3,jy), &
!                        ps(js),ps(js-1),ps(js-2),ps(js-3),psa1(ia1),fxz(jxa1_nrk(ika1,nrankxz),jza1_nrk(ika1,nrankxz)))
!!        fxz(jxa1(ia1),jza1(ia1))=(asm_a1(ia1)*x1s(ia1,js-1)+bsm_a1(ia1)*x1s(ia1,js-2)+csm_a1(ia1)*x1s(ia1,js-3)) &
!!                  /(asm_a1(ia1)+bsm_a1(ia1)+csm_a1(ia1)) 
!      enddo
!      endif
      return
      end

!ws****************************************************************
      subroutine valb3_atlastgrid_r1p0_v1(f3xz)
      use declare
      implicit real*8 (b-h,o-z)
      real*8 cx1st,cz1st
      integer is
      dimension f3xz(mx,mz,my,3),fst(n2th+5,mps4:mps,my,3),f1s(mbm_nrk,mps4:mps),fsxz(mbm_nrk,3)
      include 'mpif.h'

!       integer status(mpi_status_size)

      call map_xz2st(f3xz,fst,3)
         
      if(mb_nrk(nrankxz).gt.0) then      

      do jy=1,my
      do js=mpsa-2,mpsa-4,-1
      do jt=itbmin(nrankxz),itbmax(nrankxz)
      cx1st=fst(jt,js,jy,1)
      cz1st=fst(jt,js,jy,3)
      fst(jt,js,jy,1)=cx1st*wbxr(jt,js)+cz1st*wbzr(jt,js)
      fst(jt,js,jy,3)=cx1st*wbxt(jt,js)+cz1st*wbzp(jt,js)
      enddo
      enddo
      enddo    

      do 10 js=mpsa-2,mpsa-4,-1
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),js,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------------
      call smth_st_nrk(fst,js,3,3)
!ws_smps--------------------
   10 continue

      do 61 m=1,3

      if(m.eq.2 .or. m.eq.3) then !!(m=2,3) ibnd==0
      is=mpsa 
      do 51 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
       fst(jt,is,jy,m)=0.           
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   51 continue

      is=mpsa-1 
      do 52 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
        call interp1d3l(fst(jt,is-3,jy,m),fst(jt,is-2,jy,m),fst(jt,is-1,jy,m),fst(jt,is+1,jy,m), &
                        ps(is-3),ps(is-2),ps(is-1),ps(is+1),ps(is),fst(jt,is,jy,m))                   
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   52 continue

      call smth_st_nrk(fst(:,:,:,m),is,1,3)

      else  !!(m=1) ibnd==1  
      is=mpsa-1    
      do 31 jy=1,my      
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is)) 
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0 
!       fst(jt,is)=fst(jt,is-1) 
       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0          
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   31 continue 
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
      call smth_st_nrk(fst(:,:,:,m),is,1,3)
!ws_smps--------------------
      is=mpsa
 
      do 41 jy=1,my  
      do jt=itbmin(nrankxz),itbmax(nrankxz)
!      call extap(fst(jt,is-3),fst(jt,is-2),fst(jt,is-1),fst(jt,is))
!       fst(jt,is,jy,m)=(asm(is)*fst(jt,is-1,jy,m)+bsm(is)*fst(jt,is-2,jy,m)+csm(is)*fst(jt,is-3,jy,m))/(asm(is)+bsm(is)+csm(is))
       fst(jt,is,jy,m)=(4*fst(jt,is-1,jy,m)-fst(jt,is-2,jy,m))/3.0
!       fst(jt,is)=fst(jt,is-1)   
       if(fst(jt,is,jy,m)*fst(jt,is-1,jy,m).le.0) fst(jt,is,jy,m)=0     
      enddo
!      call smth_ps(fst(itbmin(nrankxz):itbmax(nrankxz),is,jy,m),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
   41 continue
!      call smth_ps_nrk(fst(itbmin(nrankxz):itbmax(nrankxz),is,:,:),itbmax(nrankxz)-itbmin(nrankxz)+1,3)
!ws_smps--------------------
      call smth_st_nrk(fst(:,:,:,m),is,1,3)
!	write(*,*) nrankxz,m,"smth done"
!ws_smps--------------------
      endif !!(m=2,3)ibnd==0
   61 continue
      
      do 3 jy=1,my           
      do ikb=1,mb_nrk(nrankxz)
        i=itb_nrk(ikb,nrankxz)
        ib=ib_nrk(ikb,nrankxz) 
        do m=1,3
        do js=mpsa,mpsa-3,-1     
        call interp1d3l(fst(i-2,js,jy,m),fst(i-1,js,jy,m),fst(i,js,jy,m),fst(i+1,js,jy,m), &
                        thst(i-2),thst(i-1),thst(i),thst(i+1),thb(ib),f1s(ikb,js))
        enddo
        js=mpsa
        call interp1d3l(f1s(ikb,js),f1s(ikb,js-1),f1s(ikb,js-2),f1s(ikb,js-3), &
                        ps(js),ps(js-1),ps(js-2),ps(js-3),psb(ib),fsxz(ikb,m))
        enddo
!        fxz(jbx(ib),jbz(ib))=(asm_b(ib)*x1s(ib,js-1)+bsm_b(ib)*x1s(ib,js-2)+csm_b(ib)*x1s(ib,js-3)) &
!                 /(asm_b(ib)+bsm_b(ib)+csm_b(ib))
        f3xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,1)=fsxz(ikb,1)*wbrx(ib)+fsxz(ikb,3)*wbpx(ib)
        f3xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,2)=fsxz(ikb,2)
        f3xz(jbx_nrk(ikb,nrankxz),jbz_nrk(ikb,nrankxz),jy,3)=fsxz(ikb,1)*wbrz(ib)+fsxz(ikb,3)*wbpz(ib)

      enddo
    3 continue 
      endif !!mb_nrk(nrankxz).gt.0
           
!      if(ibnd==2) then
!      
!      do ika1=1,ma1_nrk(nrankxz)
!        i=ita1_nrk(ika1,nrankxz)
!        ia1=ia1_nrk(ika1,nrankxz)
!        do js=mpsa,mpsa-3,1
!        call interp1d3l(fst(i-2,js),fst(i-1,js),fst(i,js),fst(i+1,js), &
!                        thst(i-2),thst(i-1),thst(i),thst(i+1),tha1(ia1),x1s(ika1,js,jy))
!        enddo
!        js=mpsa
!        call interp1d3l(x1s(ika1,js,jy),x1s(ika1,js-1,jy),x1s(ika1,js-2,jy),x1s(ika1,js-3,jy), &
!                        ps(js),ps(js-1),ps(js-2),ps(js-3),psa1(ia1),fxz(jxa1_nrk(ika1,nrankxz),jza1_nrk(ika1,nrankxz)))
!!        fxz(jxa1(ia1),jza1(ia1))=(asm_a1(ia1)*x1s(ia1,js-1)+bsm_a1(ia1)*x1s(ia1,js-2)+csm_a1(ia1)*x1s(ia1,js-3)) &
!!                  /(asm_a1(ia1)+bsm_a1(ia1)+csm_a1(ia1)) 
!      enddo
!      endif
      return
      end




!hwzhang*****************************************************************************
!hwzhang*****************************************************************************
!subroutines for cut cell method
!hwzhang*****************************************************************************
!hwzhang*****************************************************************************	

!hw************************************************************
!to calculate the value of x1_8bndx and x1_8bndz at the boundary points by free
!boundary conditions respectively in x and z direction
	subroutine bndry8_cut_cell
      use declare
	implicit none
      integer ibnd, itag
      include 'mpif.h'

      do 1 jy=1,my
      x1(:,:,jy,:)=x(:,:,jy,:)-xint(:,:,:)
	x1_8bndx(:,jy,:)=x_8bndx(:,jy,:)-xint_8bndx(:,:)
	x1_8bndz(:,jy,:)=x_8bndz(:,jy,:)-xint_8bndz(:,:)
   1  continue

      call mpi_transfersm(x1(:,:,:,:),8) 

	do 3 jz=iz_first,iz_last
!
    do jx=ix_first+2,ix_last
    
    itag=gdtp_ep(jx,jz,6)
    if((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz,4).eq.5)) x1(jx,jz,:,:)=0.d0

    if(gdtp_ep(jx,jz,5).eq.1) then ! jx+1 is the boundary point, use the jx, jx-1, jx-2 to calculate the bndz point
    do 11 m=1,8
    do 11 jy=1,my
!    x1_8bndz(itag,jy,m)=(axm_bndz(itag)*x1(jx,jz,jy,m)+bxm_bndz(itag)*x1(jx-1,jz,jy,m)+ &
!            cxm_bndz(itag)*x1(jx-2,jz,jy,m))/(axm_bndz(itag)+bxm_bndz(itag)+cxm_bndz(itag))
    ! fixed
    x1_8bndz(itag,jy,m)=0.d0
 11 continue

    else if((gdtp_ep(jx,jz,4).eq.5).and.(gdtp_ep(jx,jz,1).ne.5).and.(gdtp_ep(jx-1,jz,5).eq.1)) then ! jx is the inside dropped point

    do 12 m=1,8
    do 12 jy=1,my
    call interp1d2l(x1(jx-2,jz,jy,m),x1(jx-1,jz,jy,m),x1_8bndz(itag,jy,m), &
            xx(jx-2),xx(jx-1),bndz_grd(itag,1),xx(jx),x1(jx,jz,jy,m))
 12 continue

    endif
    enddo

!
    do jx=ix_first,ix_last-2

    itag=gdtp_ep(jx,jz,6)
    if((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz,4).eq.5)) x1(jx,jz,:,:)=0.d0

    if(gdtp_ep(jx,jz,5).eq.-1) then ! jx-1 is the boundary point, use the jx, jx+1, jx+2 to calculate the bndz point

    do 13 m=1,8
    do 13 jy=1,my
!    x1_8bndz(itag,jy,m)=(axp_bndz(itag)*x1(jx,jz,jy,m)+bxp_bndz(itag)*x1(jx+1,jz,jy,m)+ &
!            cxp_bndz(itag)*x1(jx+2,jz,jy,m))/(axp_bndz(itag)+bxp_bndz(itag)+cxp_bndz(itag))
    ! fixed
    x1_8bndz(itag,jy,m)=0.d0
 13 continue
    else if((gdtp_ep(jx,jz,4).eq.5).and.(gdtp_ep(jx,jz,1).ne.5).and.(gdtp_ep(jx+1,jz,5).eq.-1)) then ! jx is the inside dropped point

    do 14 m=1,8
    do 14 jy=1,my
    call interp1d2l(x1(jx+2,jz,jy,m),x1(jx+1,jz,jy,m),x1_8bndz(itag,jy,m), &
            xx(jx+2),xx(jx+1),bndz_grd(itag,1),xx(jx),x1(jx,jz,jy,m))
 14 continue
    endif
    enddo
    
    3 continue
!
    do 4 jx=ix_first,ix_last
!
    do jz=iz_first+2,iz_last

    itag=gdtp_ep(jx,jz,3)
    if((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz,4).eq.5)) x1(jx,jz,:,:)=0.d0

    if(gdtp_ep(jx,jz,2).eq.1) then ! jz+1 is the boundary point, use the jz, jz-1, jz-2 to calculate the bndx point

    do 15 m=1,8
    do 15 jy=1,my
!    x1_8bndx(itag,jy,m)=(azm_bndx(itag)*x1(jx,jz,jy,m)+bzm_bndx(itag)*x1(jx,jz-1,jy,m)+ &
!            czm_bndx(itag)*x1(jx,jz-2,jy,m))/(azm_bndx(itag)+bzm_bndx(itag)+czm_bndx(itag))
    ! fixed
    x1_8bndx(itag,jy,m)=0.d0
 15 continue
    else if((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz,4).ne.5).and.(gdtp_ep(jx,jz-1,2).eq.1)) then ! jz is the inside dropped point

    do 16 m=1,8
    do 16 jy=1,my
    call interp1d2l(x1(jx,jz-2,jy,m),x1(jx,jz-1,jy,m),x1_8bndx(itag,jy,m), &
            zz(jz-2),zz(jz-1),bndx_grd(itag,2),zz(jz),x1(jx,jz,jy,m))
 16 continue
    endif
    enddo

!
    do jz=iz_first,iz_last-2

    itag=gdtp_ep(jx,jz,3)
    if((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz,4).eq.5)) x1(jx,jz,:,:)=0.d0

    if(gdtp_ep(jx,jz,2).eq.-1) then ! jz-1 is the boundary point, use the jz, jz+1, jz+2 to calculate the bndx point

    do 17 m=1,8
    do 17 jy=1,my
!    x1_8bndx(itag,jy,m)=(azp_bndx(itag)*x1(jx,jz,jy,m)+bzp_bndx(itag)*x1(jx,jz+1,jy,m)+ &
!            czp_bndx(itag)*x1(jx,jz+2,jy,m))/(azp_bndx(itag)+bzp_bndx(itag)+czp_bndx(itag))
    ! fixed
    x1_8bndx(itag,jy,m)=0.d0
  17 continue
    else if((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz,4).ne.5).and.(gdtp_ep(jx,jz+1,2).eq.-1)) then ! jz is the inside dropped point

    do 18 m=1,5
    do 18 jy=1,my
    call interp1d2l(x1(jx,jz+2,jy,m),x1(jx,jz+1,jy,m),x1_8bndx(itag,jy,m), &
            zz(jz+2),zz(jz+1),bndx_grd(itag,2),zz(jz),x1(jx,jz,jy,m))
  18 continue
    endif
    enddo
    

    4 continue

 !     if(smoothx1) then 
 !     do m=1,8    
 !      call smthxzy(x1(:,:,:,m),1)
 !     enddo
 !     endif

!	do 10 m=1,8
!      do 10 jy=1,my
!      do 10 jz=iz_first,iz_last
!      do 10 jx=ix_first,ix_last
!      if(psi(jx,jz).lt.psia1) then
!	if((hypb_ratio(jx,jz).ge.0.d0).and.(hypb_ratio(jx,jz).le.1.d0)) then
!	x1(jx,jz,jy,m)=x1(jx,jz,jy,m)*hypb_ratio(jx,jz)
!	endif
!	endif
!   10 continue

!      call mpi_transfersm(x1(:,:,:,:),8) 


!      if(smoothp1ll) call smthp1_traceline(3)  

!	call smth_irpt_with_difc_v3(2,1,1,8,0.85d0)
!	call smth_irpt_with_difc_v2(2,1,1,8,0.85d0)

      do 19 jy=1,my
      x(:,:,jy,:)=x1(:,:,jy,:)+xint(:,:,:)
	x_8bndx(:,jy,:)=x1_8bndx(:,jy,:)+xint_8bndx(:,:)
	x_8bndz(:,jy,:)=x1_8bndz(:,jy,:)+xint_8bndz(:,:)
   19 continue
  
      return
      end

!hw************************************************************
!to calculate the value of x1_8bndx and x1_8bndz at the boundary points by free
!boundary conditions respectively in x and z direction
	subroutine bndry8_cut_cell_v2_fixed
      use declare
	implicit none
      integer ibnd, itag
	real*8 res1,res2,res3,tmp1,tmp2,tmp3,tmp
	real*8 state_func_interp1d2l
      include 'mpif.h'
! statement function for interp1d2l	
	state_func_interp1d2l(res1,res2,res3,tmp1,tmp2,tmp3,tmp) = &
	res1*(tmp-tmp2)*(tmp-tmp3)/((tmp1-tmp2)*(tmp1-tmp3)) + &
	res2*(tmp-tmp3)*(tmp-tmp1)/((tmp2-tmp3)*(tmp2-tmp1)) + &
	res3*(tmp-tmp1)*(tmp-tmp2)/((tmp3-tmp1)*(tmp3-tmp2))



      do 1 jy=1,my
	if(rmp_east) then
      x1(:,:,jy,1:5)=x(:,:,jy,1:5)-xint(:,:,1:5)
	x1(:,:,jy,6:8)=x(:,:,jy,6:8)-xint(:,:,6:8)-b_rmp_out(:,:,jy,1:3)*ft_rmp(2)
	x1_8bndx(:,jy,1:5)=x_8bndx(:,jy,1:5)-xint_8bndx(:,1:5)
	x1_8bndx(:,jy,6:8)=x_8bndx(:,jy,6:8)-xint_8bndx(:,6:8)-b_rmp_bndx(:,jy,1:3)*ft_rmp(2)
	x1_8bndz(:,jy,1:5)=x_8bndz(:,jy,1:5)-xint_8bndz(:,1:5)
	x1_8bndz(:,jy,6:8)=x_8bndz(:,jy,6:8)-xint_8bndz(:,6:8)-b_rmp_bndz(:,jy,1:3)*ft_rmp(2)
	else
      x1(:,:,jy,:)=x(:,:,jy,:)-xint(:,:,:)
	x1_8bndx(:,jy,:)=x_8bndx(:,jy,:)-xint_8bndx(:,:)
	x1_8bndz(:,jy,:)=x_8bndz(:,jy,:)-xint_8bndz(:,:)
	endif
   1  continue

      call mpi_transfersm(x1(:,:,:,:),8) 

! fixed boundary
    x1_8bndx(:,:,:)=0.d0
    x1_8bndz(:,:,:)=0.d0

! nbndx,nbndz=124 for mx=mz=64

	!$acc parallel loop collapse(2) &
	!$acc local(jx,jz,jy,itag,m) &
	!$acc copyin(xx,zz,gdtp_ep,bndz_grd,bndx_grd,iz_first,iz_last,ix_first,ix_last,iz_first_irpt,iz_last_irpt,ix_first_irpt,ix_last_irpt,iy_first,iy_last) &
	!$acc annotate(entire(xx,zz,gdtp_ep)) &
	!$acc annotate(dimension(bndz_grd(nbndz,n7),bndx_grd(nbndx,n7)))
    do 3 jz=iz_first_irpt,iz_last_irpt
    do 3 jx=ix_first_irpt,ix_last_irpt
!    do 3 jz=iz_first,iz_last
!    do 3 jx=ix_first,ix_last

    if((gdtp_ep(jx,jz,1).ge.4).and.(gdtp_ep(jx,jz,4).ge.4)) x1(jx,jz,:,:)=0.d0

    itag=gdtp_ep(jx,jz,6)
!    if((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz,4).eq.5)) x1(jx,jz,:,:)=0.d0

!    if((gdtp_ep(jx,jz,5).eq.1).or.(gdtp_ep(jx,jz,5).eq.-1)) then ! jx+-1 is the boundary point, use the jx, jx-1, jx-2 to calculate the bndz point
!    x1_8bndz(itag,:,:)=0.d0

!    else if((gdtp_ep(jx,jz,4).eq.5).and.(gdtp_ep(jx,jz,1).ne.5).and.(gdtp_ep(jx-1,jz,5).eq.1).and.(jx.gt.ix_first+1)) then ! jx is the inside dropped point

    if((gdtp_ep(jx,jz,4).eq.5).and.(gdtp_ep(jx,jz,1).ne.5).and.(gdtp_ep(jx-1,jz,5).eq.1).and.(jx.gt.ix_first+1)) then ! jx is the inside dropped point

    do m=1,8
	!$acc data copyin(x1(jx-2:jx-1,jz,*,m),x1_8bndz(itag,*,m)) &
	!$acc copyout(x1(jx,jz,*,m)) &
	!$acc annotate(dimension(x1_8bndz(nbndz,my,8))) &
	!$acc present(xx,zz,gdtp_ep,bndz_grd,bndx_grd,iz_first,iz_last,ix_first,ix_last,iz_first_irpt,iz_last_irpt,ix_first_irpt,ix_last_irpt,iy_first,iy_last,jx,jz,jy,itag,m)
    do jy=iy_first,iy_last
    x1(jx,jz,jy,m)=state_func_interp1d2l(x1(jx-2,jz,jy,m),x1(jx-1,jz,jy,m),x1_8bndz(itag,jy,m), &
            xx(jx-2),xx(jx-1),bndz_grd(itag,1),xx(jx))
!    call interp1d2l(x1(jx-2,jz,jy,m),x1(jx-1,jz,jy,m),x1_8bndz(itag,jy,m), &
!!    call interp1d2l(x1(jx-2,jz,jy,m),x1(jx-1,jz,jy,m),0.d0, &
!            xx(jx-2),xx(jx-1),bndz_grd(itag,1),xx(jx),x1(jx,jz,jy,m))
    enddo
 	!$acc end data
    enddo

    else if((gdtp_ep(jx,jz,4).eq.5).and.(gdtp_ep(jx,jz,1).ne.5).and.(gdtp_ep(jx+1,jz,5).eq.-1).and.(jx.lt.ix_last-1)) then ! jx is the inside dropped point

    do m=1,8
	!$acc data copyin(x1(jx+1:jx+2,jz,*,m),x1_8bndz(itag,*,m)) &
	!$acc copyout(x1(jx,jz,*,m)) &
	!$acc annotate(dimension(x1_8bndz(nbndz,my,8))) &
	!$acc present(xx,zz,gdtp_ep,bndz_grd,bndx_grd,iz_first,iz_last,ix_first,ix_last,iz_first_irpt,iz_last_irpt,ix_first_irpt,ix_last_irpt,iy_first,iy_last,jx,jz,jy,itag,m)
    do jy=iy_first,iy_last
    x1(jx,jz,jy,m)=state_func_interp1d2l(x1(jx+2,jz,jy,m),x1(jx+1,jz,jy,m),x1_8bndz(itag,jy,m), &
            xx(jx+2),xx(jx+1),bndz_grd(itag,1),xx(jx))
!    call interp1d2l(x1(jx+2,jz,jy,m),x1(jx+1,jz,jy,m),x1_8bndz(itag,jy,m), &
!!    call interp1d2l(x1(jx+2,jz,jy,m),x1(jx+1,jz,jy,m),0.d0, &
!            xx(jx+2),xx(jx+1),bndz_grd(itag,1),xx(jx),x1(jx,jz,jy,m))
    enddo
 	!$acc end data
    enddo

    endif


!
    itag=gdtp_ep(jx,jz,3)
!    if((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz,4).eq.5)) x1(jx,jz,:,:)=0.d0

!    if((gdtp_ep(jx,jz,2).eq.1).or.(gdtp_ep(jx,jz,2).eq.-1)) then ! jz+-1 is the boundary point, use the jz, jz-1, jz-2 to calculate the bndx point
!    x1_8bndx(itag,:,:)=0.d0
!    else if((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz,4).ne.5).and.(gdtp_ep(jx,jz-1,2).eq.1).and.(jz.gt.iz_first+1)) then ! jz is the inside dropped point

    if((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz,4).ne.5).and.(gdtp_ep(jx,jz-1,2).eq.1).and.(jz.gt.iz_first+1)) then ! jz is the inside dropped point

    do m=1,8
	!$acc data copyin(x1(jx,jz-2:jz-1,*,m),x1_8bndx(itag,*,m)) &
	!$acc copyout(x1(jx,jz,*,m)) &
	!$acc annotate(dimension(x1_8bndx(nbndx,my,8))) &
	!$acc present(xx,zz,gdtp_ep,bndz_grd,bndx_grd,iz_first,iz_last,ix_first,ix_last,iz_first_irpt,iz_last_irpt,ix_first_irpt,ix_last_irpt,iy_first,iy_last,jx,jz,jy,itag,m)
    do jy=iy_first,iy_last
    x1(jx,jz,jy,m)=state_func_interp1d2l(x1(jx,jz-2,jy,m),x1(jx,jz-1,jy,m),x1_8bndx(itag,jy,m), &
            zz(jz-2),zz(jz-1),bndx_grd(itag,2),zz(jz))
!    call interp1d2l(x1(jx,jz-2,jy,m),x1(jx,jz-1,jy,m),x1_8bndx(itag,jy,m), &
!!    call interp1d2l(x1(jx,jz-2,jy,m),x1(jx,jz-1,jy,m),0.d0, &
!            zz(jz-2),zz(jz-1),bndx_grd(itag,2),zz(jz),x1(jx,jz,jy,m))
    enddo
 	!$acc end data
    enddo

    else if((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz,4).ne.5).and.(gdtp_ep(jx,jz+1,2).eq.-1).and.(jz.lt.iz_last-1)) then ! jz is the inside dropped point

    do m=1,8
	!$acc data copyin(x1(jx,jz+1:jz+2,*,m),x1_8bndx(itag,*,m)) &
	!$acc copyout(x1(jx,jz,*,m)) &
	!$acc annotate(dimension(x1_8bndx(nbndx,my,8))) &
	!$acc present(xx,zz,gdtp_ep,bndz_grd,bndx_grd,iz_first,iz_last,ix_first,ix_last,iz_first_irpt,iz_last_irpt,ix_first_irpt,ix_last_irpt,iy_first,iy_last,jx,jz,jy,itag,m)
    do jy=iy_first,iy_last
    x1(jx,jz,jy,m)=state_func_interp1d2l(x1(jx,jz+2,jy,m),x1(jx,jz+1,jy,m),x1_8bndx(itag,jy,m), &
            zz(jz+2),zz(jz+1),bndx_grd(itag,2),zz(jz))
!    call interp1d2l(x1(jx,jz+2,jy,m),x1(jx,jz+1,jy,m),x1_8bndx(itag,jy,m), &
!!    call interp1d2l(x1(jx,jz+2,jy,m),x1(jx,jz+1,jy,m),0.d0, &
!            zz(jz+2),zz(jz+1),bndx_grd(itag,2),zz(jz),x1(jx,jz,jy,m))
    enddo
 	!$acc end data
    enddo
    endif

    3 continue
    	!$acc end parallel loop

 !     if(smoothx1) then 
 !     do m=1,8    
 !      call smthxzy(x1(:,:,:,m),1)
 !     enddo
 !     endif

      do 10 m=1,5
      do 10 jy=1,my
      do 10 jz=iz_first,iz_last
      do 10 jx=ix_first,ix_last
!      if(psi(jx,jz).lt.psia1) then
!	if((hypb_ratio(jx,jz).ge.0.d0).and.(hypb_ratio(jx,jz).le.1.d0)) then
	x1(jx,jz,jy,m)=x1(jx,jz,jy,m)*hypb_ratio(jx,jz)
!	endif
!	endif
   10 continue

!      print*, 'success 23052', nrank
!      call mpi_transfersm(x1(:,:,:,:),8) 


!      if(smoothp1ll) call smthp1_traceline(3)  

!	call smth_irpt_with_difc(2,1,1,8,0.85d0)
        if(nstep.eq.0) x1(:,:,:,:)=0.d0 
!        call smth_irpt_with_difc_v2(2,1,1,2,0.85d0)
!	call smth_irpt_with_difc_v3(2,1,1,8,0.85d0) !open if necessary

!      print*, 'success 23062', nstep
      do 19 jy=1,my
	if(rmp_east) then
      x(:,:,jy,1:5)=x1(:,:,jy,1:5)+xint(:,:,1:5)
	x(:,:,jy,6:8)=x1(:,:,jy,6:8)+xint(:,:,6:8)+b_rmp_out(:,:,jy,1:3)*ft_rmp(2)
	x_8bndx(:,jy,1:5)=x1_8bndx(:,jy,1:5)+xint_8bndx(:,1:5)
	x_8bndx(:,jy,6:8)=x1_8bndx(:,jy,6:8)+xint_8bndx(:,6:8)+b_rmp_bndx(:,jy,1:3)*ft_rmp(2)
	x_8bndz(:,jy,1:5)=x1_8bndz(:,jy,1:5)+xint_8bndz(:,1:5)
	x_8bndz(:,jy,6:8)=x1_8bndz(:,jy,6:8)+xint_8bndz(:,6:8)+b_rmp_bndz(:,jy,1:3)*ft_rmp(2)

	x1(:,:,jy,6:8)=x1(:,:,jy,6:8)+b_rmp_out(:,:,jy,1:3)*ft_rmp(2)
	x1_8bndx(:,jy,6:8)=x1_8bndx(:,jy,6:8)+b_rmp_bndx(:,jy,1:3)*ft_rmp(2)
	x1_8bndz(:,jy,6:8)=x1_8bndz(:,jy,6:8)+b_rmp_bndz(:,jy,1:3)*ft_rmp(2)
	else
      x(:,:,jy,:)=x1(:,:,jy,:)+xint(:,:,:)
	x_8bndx(:,jy,:)=x1_8bndx(:,jy,:)+xint_8bndx(:,:)
	x_8bndz(:,jy,:)=x1_8bndz(:,jy,:)+xint_8bndz(:,:)
	endif
   19 continue
  
!      print*, 'success 23082', nstep
      return
      end

!hw**************************************************************
	subroutine decide_grd_type_in_each_proc
	use declare
	implicit none
!	integer, dimension(mxt,mzt,2) :: grd_type 
!!1-4 means regualr, irregular, boundary, dropped for bndx in z direction or for bndz in x direction
!	integer, dimension(mx,mz,6) :: gdtp_ep ! the grid type for each point in every processor

      character*9 output
      character*3 cn
	include 'mpif.h'
      

! set to dropped point at first
	gdtp_ep(:,:,1)=4
	gdtp_ep(:,:,4)=4

	do jz=iz_first,iz_last
	do jx=ix_first,ix_last
	gdtp_ep(jx,jz,1)=grd_type(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2,1)
	gdtp_ep(jx,jz,2:3)=gdtp_bndx(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2,1:2)
	gdtp_ep(jx,jz,4)=grd_type(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2,2)
	gdtp_ep(jx,jz,5:6)=gdtp_bndz(nrkx(nrank)*mxm+jx-2,nrkz(nrank)*mzm+jz-2,1:2)
	enddo
	enddo

	output='gdtpx'//cn(nrank)
	open(unit=1,file=output,status='unknown',form='formatted')
	output='gdtpz'//cn(nrank)
	open(unit=2,file=output,status='unknown',form='formatted')
	write(1,3)((xx(jx),zz(jz),gdtp_ep(jx,jz,1)*1.d0,gdtp_ep(jx,jz,2)*1.d0,gdtp_ep(jx,jz,3)*1.d0,jx=1,mx),jz=1,mz)
	write(2,3)((xx(jx),zz(jz),gdtp_ep(jx,jz,4)*1.d0,gdtp_ep(jx,jz,5)*1.d0,gdtp_ep(jx,jz,6)*1.d0,jx=1,mx),jz=1,mz)
	close(1)
	close(2)
    3 format(5(1x,e12.5)) 

	return
	end

!hw******************************************************************
	subroutine calculate_dnfm_coeff_for_ir_bndx_point
	use declare
	implicit none
	integer itag
	real*8 xzp2,xzp1,xz00,xzm1,xzm2
	real*8, dimension(4,4) :: inpA
	real*8, dimension(4) :: inpB, outX 
	real*8 dxp1,dxm1,dxp2,dxm2,dzp1,dzm1,dzp2,dzm2,dxp3,dxm3,dzp3,dzm3,ca1
	include 'mpif.h'
! coeficient ax1_irx,az1_irx,....dx1_irx, dz1_irx related to the first 
! differential with a fourth order accuracy, only use for the IR2 point with
! gdtp_bndx z .eq. +2 or -2, whose stencil is as R R R IR2 IR1 (DI) B DO DO DO DO

!      real*8, dimension(mx,mz) :: az1_irx,bz1_irx,cz1_irx,dz1_irx
!      real*8, dimension(mx,mz) :: ax1_irz,bx1_irz,cx1_irz,dx1_irz
!      real*8, dimension(mx,mz) :: azbp_irx,bzbp_irx,czbp_irx,dzbp_irx
!      real*8, dimension(mx,mz) :: azbm_irx,bzbm_irx,czbm_irx,dzbm_irx
!      real*8, dimension(mx,mz) :: axbp_irz,bxbp_irz,cxbp_irz,dxbp_irz
!      real*8, dimension(mx,mz) :: axbm_irz,bxbm_irz,cxbm_irz,dxbm_irz

!      real*8, dimension(mx,mz) :: a2zbm_irx,b2zbm_irx,c2zbm_irx !d2fbm
!      real*8, dimension(mx,mz) :: a2zbp_irx,b2zbp_irx,c2zbp_irx !d2fbp
!      real*8, dimension(mx,mz) :: a2xbm_irz,b2xbm_irz,c2xbm_irz !d2fbm
!      real*8, dimension(mx,mz) :: a2xbp_irz,b2xbp_irz,c2xbp_irz !d2fbp

!      real*8, dimension(mx,mz) :: az2_irx,bz2_irx,cz2_irx,dz2_irx
!      real*8, dimension(mx,mz) :: ax2_irz,bx2_irz,cx2_irz,dx2_irz

!	real*8, allocatable :: axm_bndz(:), bxm_bndz(:), cxm_bndz(:), axp_bndz(:), bxp_bndz(:), cxp_bndz(:)
!	real*8, allocatable :: azm_bndx(:), bzm_bndx(:), czm_bndx(:), azp_bndx(:), bzp_bndx(:), czp_bndx(:)
	allocate(axm_bndz(nbndz))
	allocate(bxm_bndz(nbndz))
	allocate(cxm_bndz(nbndz))
	allocate(axp_bndz(nbndz))
	allocate(bxp_bndz(nbndz))
	allocate(cxp_bndz(nbndz))
	allocate(azm_bndx(nbndx))
	allocate(bzm_bndx(nbndx))
	allocate(czm_bndx(nbndx))
	allocate(azp_bndx(nbndx))
	allocate(bzp_bndx(nbndx))
	allocate(czp_bndx(nbndx))
! for bndz with ax1_irz, bx1_irz
	do 1 jz=iz_first,iz_last
	do 1 jx=ix_first+2,ix_last-2
	if (gdtp_ep(jx,jz,5).eq.2) then  !means upper two xgrds is the bnd point
		  itag=gdtp_ep(jx,jz,6)  !the rank of Bnd pnt in bndz_grd
		  xzp2=bndz_grd(itag,1)
		  xzp1=xx(jx+1)
		  xz00=xx(jx)
		  xzm1=xx(jx-1)
		  xzm2=xx(jx-2)

		  dxp1=xzp1-xz00
		  dxm1=xz00-xzm1
		  dxp2=xzp2-xz00
		  dxm2=xz00-xzm2
		  f1= dxp1   +dxp1**4/dxm2**3
		  f2= dxp1**2-dxp1**4/dxm2**2
		  f3= dxp1**3+dxp1**4/dxm2
		  g1= dxm1   -dxm1**4/dxm2**3
		  g2=-dxm1**2+dxm1**4/dxm2**2
		  g3= dxm1**3-dxm1**4/dxm2
		  h1= dxp2   +dxp2**4/dxm2**3
		  h2= dxp2**2-dxp2**4/dxm2**2
		  h3= dxp2**3+dxp2**4/dxm2
		  ca1= (f1*h3-f3*h1)*(g2*h3-g3*h2)-(g1*h3-g3*h1)*(f2*h3-f3*h2)
		  ax1_irz(jx,jz)=h3*(g2*h3-g3*h2)/ca1
		  bx1_irz(jx,jz)=h3*(h2*f3-h3*f2)/ca1
		  cx1_irz(jx,jz)=h3*(f2*g3-f3*g2)/ca1
		  dx1_irz(jx,jz)=(dxp1**2*ax1_irz(jx,jz)-dxm1**2*bx1_irz(jx,jz) &
			    +dxp2**2*cx1_irz(jx,jz))/dxm2**2
	endif

	if (gdtp_ep(jx,jz,5).eq.-2) then  !means upper two xgrds is the bnd point
		  itag=gdtp_ep(jx,jz,6)  !the rank of Bnd pnt in bndz_grd
		  xzp2=xx(jx+2)
		  xzp1=xx(jx+1)
		  xz00=xx(jx)
		  xzm1=xx(jx-1)
		  xzm2=bndz_grd(itag,1)

		  dxp1=xzp1-xz00
		  dxm1=xz00-xzm1
		  dxp2=xzp2-xz00
		  dxm2=xz00-xzm2
		  f1= dxp1   +dxp1**4/dxm2**3
		  f2= dxp1**2-dxp1**4/dxm2**2
		  f3= dxp1**3+dxp1**4/dxm2
		  g1= dxm1   -dxm1**4/dxm2**3
		  g2=-dxm1**2+dxm1**4/dxm2**2
		  g3= dxm1**3-dxm1**4/dxm2
		  h1= dxp2   +dxp2**4/dxm2**3
		  h2= dxp2**2-dxp2**4/dxm2**2
		  h3= dxp2**3+dxp2**4/dxm2
		  ca1= (f1*h3-f3*h1)*(g2*h3-g3*h2)-(g1*h3-g3*h1)*(f2*h3-f3*h2)
		  ax1_irz(jx,jz)=h3*(g2*h3-g3*h2)/ca1
		  bx1_irz(jx,jz)=h3*(h2*f3-h3*f2)/ca1
		  cx1_irz(jx,jz)=h3*(f2*g3-f3*g2)/ca1
		  dx1_irz(jx,jz)=(dxp1**2*ax1_irz(jx,jz)-dxm1**2*bx1_irz(jx,jz) &
			    +dxp2**2*cx1_irz(jx,jz))/dxm2**2
	endif
    1 continue


	do 2 jx=ix_first,ix_last
      do 2 jz=iz_first+2,iz_last-2
	if (gdtp_ep(jx,jz,2).eq.2) then  !means upper two zgrds is the bnd point
		  itag=gdtp_ep(jx,jz,3)  !the rank of Bnd pnt in bndz_grd
		  xzp2=bndx_grd(itag,2)
		  xzp1=zz(jz+1)
		  xz00=zz(jz)
		  xzm1=zz(jz-1)
		  xzm2=zz(jz-2)
		  dzp1=xzp1-xz00
		  dzm1=xz00-xzm1
		  dzp2=xzp2-xz00
		  dzm2=xz00-xzm2
		  f1= dzp1   +dzp1**4/dzm2**3
		  f2= dzp1**2-dzp1**4/dzm2**2
		  f3= dzp1**3+dzp1**4/dzm2
		  g1= dzm1   -dzm1**4/dzm2**3
		  g2=-dzm1**2+dzm1**4/dzm2**2
		  g3= dzm1**3-dzm1**4/dzm2
		  h1= dzp2   +dzp2**4/dzm2**3
		  h2= dzp2**2-dzp2**4/dzm2**2
		  h3= dzp2**3+dzp2**4/dzm2
		  ca1= (f1*h3-f3*h1)*(g2*h3-g3*h2)-(g1*h3-g3*h1)*(f2*h3-f3*h2)
		  az1_irx(jx,jz)=h3*(g2*h3-g3*h2)/ca1
		  bz1_irx(jx,jz)=h3*(h2*f3-h3*f2)/ca1
		  cz1_irx(jx,jz)=h3*(f2*g3-f3*g2)/ca1
		  dz1_irx(jx,jz)=(dzp1**2*az1_irx(jx,jz)-dzm1**2*bz1_irx(jx,jz) &
			    +dzp2**2*cz1_irx(jx,jz))/dzm2**2

	endif

	if (gdtp_ep(jx,jz,2).eq.-2) then  !means lower two zgrds is the bnd point
		  itag=gdtp_ep(jx,jz,3)  !the rank of Bnd pnt in bndz_grd
		  xzp2=zz(jz+2)
		  xzp1=zz(jz+1)
		  xz00=zz(jz)
		  xzm1=zz(jz-1)
		  xzm2=bndx_grd(itag,2)
		  dzp1=xzp1-xz00
		  dzm1=xz00-xzm1
		  dzp2=xzp2-xz00
		  dzm2=xz00-xzm2
		  f1= dzp1   +dzp1**4/dzm2**3
		  f2= dzp1**2-dzp1**4/dzm2**2
		  f3= dzp1**3+dzp1**4/dzm2
		  g1= dzm1   -dzm1**4/dzm2**3
		  g2=-dzm1**2+dzm1**4/dzm2**2
		  g3= dzm1**3-dzm1**4/dzm2
		  h1= dzp2   +dzp2**4/dzm2**3
		  h2= dzp2**2-dzp2**4/dzm2**2
		  h3= dzp2**3+dzp2**4/dzm2
		  ca1= (f1*h3-f3*h1)*(g2*h3-g3*h2)-(g1*h3-g3*h1)*(f2*h3-f3*h2)
		  az1_irx(jx,jz)=h3*(g2*h3-g3*h2)/ca1
		  bz1_irx(jx,jz)=h3*(h2*f3-h3*f2)/ca1
		  cz1_irx(jx,jz)=h3*(f2*g3-f3*g2)/ca1
		  dz1_irx(jx,jz)=(dzp1**2*az1_irx(jx,jz)-dzm1**2*bz1_irx(jx,jz) &
			    +dzp2**2*cz1_irx(jx,jz))/dzm2**2

	endif
    2 continue

! for bndz with axbm_irz, bxbm_irz, ...
	do 3 jz=iz_first,iz_last
	do 3 jx=ix_first+2,ix_last-1
	if(gdtp_ep(jx,jz,5).eq.1) then !means the upper one xgrid is the bnd point. R R IR2 IR1 B
		  itag=gdtp_ep(jx,jz,6)
		  xzp1=bndz_grd(itag,1)
		  xz00=xx(jx)
		  xzm1=xx(jx-1)
		  xzm2=xx(jx-2)
		  dxp1=xzp1-xz00!xx(jx+1)-xx(jx)
		  dxm1=xz00-xzm1!xx(jx)-xx(jx-1)
		  dxm2=xz00-xzm2!xx(jx)-xx(jx-2)
		  f1=-dxp1+dxp1**3/dxm2**2
		  f2=dxp1**2+dxp1**3/dxm2
		  g1=dxm1-dxm1**3/dxm2**2
		  g2=dxm1**2-dxm1**3/dxm2
		  ca1=f1*g2-f2*g1
		  axbm_irz(jx,jz)=g2/ca1
		  bxbm_irz(jx,jz)=-f2/ca1
		  cxbm_irz(jx,jz)=(1+axbm_irz(jx,jz)*dxp1-bxbm_irz(jx,jz)*dxm1)/dxm2
	endif
    3 continue

! for bndz with axbp_irz, bxbp_irz, ...
	do 4 jz=iz_first,iz_last
	do 4 jx=ix_first+1,ix_last-2
	if(gdtp_ep(jx,jz,5).eq.-1) then !means the lower one xgrid is the bnd point. B IR1 IR2 R R
		  itag=gdtp_ep(jx,jz,6)
		  xzp2=xx(jx+2)
		  xzp1=xx(jx+1)
		  xz00=xx(jx)
		  xzm1=bndz_grd(itag,1)
		  dxm1=xz00-xzm1!xx(jx)-xx(jx-1)
		  dxp1=xzp1-xz00!xx(jx+1)-xx(jx)
		  dxp2=xzp2-xz00!xx(jx+2)-xx(jx)
		  f1=-dxm1+dxm1**3/dxp2**2
		  f2=dxm1**2+dxm1**3/dxp2
		  g1=dxp1-dxp1**3/dxp2**2
		  g2=dxp1**2-dxp1**3/dxp2
		  ca1=f1*g2-f2*g1
		  axbp_irz(jx,jz)=g2/ca1
		  bxbp_irz(jx,jz)=-f2/ca1
		  cxbp_irz(jx,jz)=(1+axbp_irz(jx,jz)*dxm1-bxbp_irz(jx,jz)*dxp1)/dxp2
	endif
    4 continue

! for bndx with azbm_irx, bzbm_irx, ...
	do 5 jx=ix_first,ix_last
	do 5 jz=iz_first+2,iz_last-1
	if(gdtp_ep(jx,jz,2).eq.1) then !means the upper one zgrid is the bnd point, R R IR2 IR1 B
		  itag=gdtp_ep(jx,jz,3)
		  xzp1=bndx_grd(itag,2)
		  xz00=zz(jz)
		  xzm1=zz(jz-1)
		  xzm2=zz(jz-2)
		  dzp1=xzp1-xz00!zz(jz+1)-zz(jz)
		  dzm1=xz00-xzm1!zz(jz)-zz(jz-1)
		  dzm2=xz00-xzm2!zz(jz)-zz(jz-2)
		  f1=-dzp1+dzp1**3/dzm2**2
		  f2=dzp1**2+dzp1**3/dzm2
		  g1=dzm1-dzm1**3/dzm2**2
		  g2=dzm1**2-dzm1**3/dzm2
		  ca1=f1*g2-f2*g1
		  azbm_irx(jx,jz)=g2/ca1
		  bzbm_irx(jx,jz)=-f2/ca1
		  czbm_irx(jx,jz)=(1+azbm_irx(jx,jz)*dzp1-bzbm_irx(jx,jz)*dzm1)/dzm2
	endif
    5 continue



! for bndx with azbp_irx, bzbp_irx, ...
	do 6 jx=ix_first,ix_last
	do 6 jz=iz_first+1,iz_last-2
	if(gdtp_ep(jx,jz,2).eq.-1) then !means the lower one zgrid is the bnd point, B IR1 IR2 R R
		  itag=gdtp_ep(jx,jz,3)
		  xzp2=zz(jz+2)
		  xzp1=zz(jz+1)
		  xz00=zz(jz)
		  xzm1=bndx_grd(itag,2)
		  dzm1=xz00-xzm1!zz(jz)-zz(jz-1)
		  dzp1=xzp1-xz00!zz(jz+1)-zz(jz)
		  dzp2=xzp2-xz00!zz(jz+2)-zz(jz)
		  f1=-dzm1+dzm1**3/dzp2**2
		  f2=dzm1**2+dzm1**3/dzp2
		  g1=dzp1-dzp1**3/dzp2**2
		  g2=dzp1**2-dzp1**3/dzp2
		  ca1=f1*g2-f2*g1
		  azbp_irx(jx,jz)=g2/ca1
		  bzbp_irx(jx,jz)=-f2/ca1
		  czbp_irx(jx,jz)=(1+azbp_irx(jx,jz)*dzm1-bzbp_irx(jx,jz)*dzp1)/dzp2
	endif
    6 continue


! for bndz with a2xbm_irz, b2xbm_irz, c2xbm_irz
	do 7 jz=iz_first,iz_last
	do 7 jx=ix_first+2,ix_last-1
	if(gdtp_ep(jx,jz,5).eq.1) then !means the upper one x grid is the bnd point, R R IR2 IR1 B
		  itag=gdtp_ep(jx,jz,6)
		  xzp1=bndz_grd(itag,1)
		  xz00=xx(jx)
		  xzm1=xx(jx-1)
		  xzm2=xx(jx-2)
		  inpA(1,:)=1.d0
		  inpA(2,1)=xzm2-xz00!xx(jx-2)-xx(jx)
		  inpA(2,2)=xzm1-xz00!xx(jx-1)-xx(jx)
		  inpA(2,3)=0.d0
		  inpA(2,4)=xzp1-xz00
		  inpA(3,1)=(xzm2-xz00)**2/2.d0
		  inpA(3,2)=(xzm1-xz00)**2/2.d0
		  inpA(3,3)=0.d0
		  inpA(3,4)=(xzp1-xz00)**2/2.d0
		  inpA(4,1)=(xzm2-xz00)**3/6.d0
		  inpA(4,2)=(xzm1-xz00)**3/6.d0
		  inpA(4,3)=0.d0
		  inpA(4,4)=(xzp1-xz00)**3/6.d0
		  inpB(:)=0.d0
		  inpB(3)=1.d0
		  call gauss_solve(inpA,inpB,outX,4)
		  a2xbm_irz(jx,jz)=outX(4)/2.d0
		  b2xbm_irz(jx,jz)=outX(2)/2.d0
		  c2xbm_irz(jx,jz)=outX(1)/2.d0
     endif
   7 continue

! for bndz with a2xbp_irz, b2xbp_irz, c2xbp_irz
	do 8 jz=iz_first,iz_last
	do 8 jx=ix_first+1,ix_last-2
	if(gdtp_ep(jx,jz,5).eq.-1) then !means the lower one xgrid is the bnd point, B IR1 IR2 R R
		  itag=gdtp_ep(jx,jz,6)
		  xzp2=xx(jx+2)
		  xzp1=xx(jx+1)
		  xz00=xx(jx)
		  xzm1=bndz_grd(itag,1)
		  inpA(1,:)=1.d0
		  inpA(2,1)=xzm1-xz00
		  inpA(2,2)=0.d0
		  inpA(2,3)=xzp1-xz00
		  inpA(2,4)=xzp2-xz00
		  inpA(3,1)=(xzm1-xz00)**2/2.d0
		  inpA(3,2)=0.d0
		  inpA(3,3)=(xzp1-xz00)**2/2.d0
		  inpA(3,4)=(xzp2-xz00)**2/2.d0
		  inpA(4,1)=(xzm1-xz00)**3/6.d0
		  inpA(4,2)=0.d0
		  inpA(4,3)=(xzp1-xz00)**3/6.d0
		  inpA(4,4)=(xzp2-xz00)**3/6.d0
		  inpB(:)=0.d0
		  inpB(3)=1.d0
		  call gauss_solve(inpA,inpB,outX,4)
		  a2xbp_irz(jx,jz)=outX(1)/2.d0
		  b2xbp_irz(jx,jz)=outX(3)/2.d0
		  c2xbp_irz(jx,jz)=outX(4)/2.d0
	endif
    8 continue
   	
! for bndx with a2zbm_irx, b2zbm_irx, c2zbm_irx
	do 9 jx=ix_first,ix_last
	do 9 jz=iz_first+2,iz_last-1
	if(gdtp_ep(jx,jz,2).eq.1) then ! means the upper one zgrid is the bnd point, R R IR2 IR1 B
		  itag=gdtp_ep(jx,jz,3)
		  xzp1=bndx_grd(itag,2)
		  xz00=zz(jz)
		  xzm1=zz(jz-1)
		  xzm2=zz(jz-2)
		  inpA(1,:)=1.d0
		  inpA(2,1)=xzm2-xz00!xx(jx-2)-xx(jx)
		  inpA(2,2)=xzm1-xz00!xx(jx-1)-xx(jx)
		  inpA(2,3)=0.d0
		  inpA(2,4)=xzp1-xz00
		  inpA(3,1)=(xzm2-xz00)**2/2.d0
		  inpA(3,2)=(xzm1-xz00)**2/2.d0
		  inpA(3,3)=0.d0
		  inpA(3,4)=(xzp1-xz00)**2/2.d0
		  inpA(4,1)=(xzm2-xz00)**3/6.d0
		  inpA(4,2)=(xzm1-xz00)**3/6.d0
		  inpA(4,3)=0.d0
		  inpA(4,4)=(xzp1-xz00)**3/6.d0
		  inpB(:)=0.d0
		  inpB(3)=1.d0
		  call gauss_solve(inpA,inpB,outX,4)
		  a2zbm_irx(jx,jz)=outX(4)/2.d0
		  b2zbm_irx(jx,jz)=outX(2)/2.d0
		  c2zbm_irx(jx,jz)=outX(1)/2.d0
	endif
    9 continue

! for bndx with a2zbp_irx, b2zbp_irx, c2zbp_irx
	do 10 jx=ix_first,ix_last
	do 10 jz=iz_first+1,iz_last-2
	if(gdtp_ep(jx,jz,2).eq.-1) then !means the lower one zgrid is the bnd point, B IR1 IR2 R R
		  itag=gdtp_ep(jx,jz,3)
		  xzp2=zz(jz+2)
		  xzp1=zz(jz+1)
		  xz00=zz(jz)
		  xzm1=bndx_grd(itag,2)
		  inpA(1,:)=1.d0
		  inpA(2,1)=xzm1-xz00
		  inpA(2,2)=0.d0
		  inpA(2,3)=xzp1-xz00
		  inpA(2,4)=xzp2-xz00
		  inpA(3,1)=(xzm1-xz00)**2/2.d0
		  inpA(3,2)=0.d0
		  inpA(3,3)=(xzp1-xz00)**2/2.d0
		  inpA(3,4)=(xzp2-xz00)**2/2.d0
		  inpA(4,1)=(xzm1-xz00)**3/6.d0
		  inpA(4,2)=0.d0
		  inpA(4,3)=(xzp1-xz00)**3/6.d0
		  inpA(4,4)=(xzp2-xz00)**3/6.d0
		  inpB(:)=0.d0
		  inpB(3)=1.d0
		  call gauss_solve(inpA,inpB,outX,4)
		  a2zbp_irx(jx,jz)=outX(1)/2.d0
		  b2zbp_irx(jx,jz)=outX(3)/2.d0
		  c2zbp_irx(jx,jz)=outX(4)/2.d0
	endif
   10 continue


! for bndz with ax2_irz,bx2_irz,...
	do 11 jz=iz_first,iz_last
	do 11 jx=ix_first+2,ix_last-2
	if(gdtp_ep(jx,jz,5).eq.2) then ! means upper two xgrid is the bnd point
		  itag=gdtp_ep(jx,jz,6)
		  xzp2=bndz_grd(itag,1)
		  xzp1=xx(jx+1)
		  xz00=xx(jx)
		  xzm1=xx(jx-1)
		  xzm2=xx(jx-2)

		  dxp1=xzp1-xz00
		  dxm1=xz00-xzm1
		  dxp2=xzp2-xz00
		  dxm2=xz00-xzm2
		  f1= dxp1**2+dxp1*dxm2
		  f2= dxp1**3-dxp1*dxm2**2
		  f3= dxp1**4+dxp1*dxm2**3
		  g1=-dxm1**2+dxm1*dxm2
		  g2= dxm1**3-dxm1*dxm2**2
		  g3=-dxm1**4+dxm1*dxm2**3
		  h1= dxp2**2+dxp2*dxm2
		  h2= dxp2**3-dxp2*dxm2**2
		  h3= dxp2**4+dxp2*dxm2**3
		  ca1= (f1*h3-f3*h1)*(g2*h3-g3*h2)-(g1*h3-g3*h1)*(f2*h3-f3*h2)
		  ax2_irz(jx,jz)=2.*h3*(g2*h3-g3*h2)/ca1
		  bx2_irz(jx,jz)=2.*h3*(h2*f3-h3*f2)/ca1
		  cx2_irz(jx,jz)=2.*h3*(f2*g3-f3*g2)/ca1
		  dx2_irz(jx,jz)=-(dxp1*ax2_irz(jx,jz)+dxm1*bx2_irz(jx,jz) &
			    +dxp2*cx2_irz(jx,jz))/dxm2
	endif

	if(gdtp_ep(jx,jz,5).eq.-2) then ! means lower two xgrid is the bnd point
		  itag=gdtp_ep(jx,jz,6)
		  xzp2=xx(jx+2)
		  xzp1=xx(jx+1)
		  xz00=xx(jx)
		  xzm1=xx(jx-1)
		  xzm2=bndz_grd(itag,1)


		  dxp1=xzp1-xz00
		  dxm1=xz00-xzm1
		  dxp2=xzp2-xz00
		  dxm2=xz00-xzm2
		  f1= dxp1**2+dxp1*dxm2
		  f2= dxp1**3-dxp1*dxm2**2
		  f3= dxp1**4+dxp1*dxm2**3
		  g1=-dxm1**2+dxm1*dxm2
		  g2= dxm1**3-dxm1*dxm2**2
		  g3=-dxm1**4+dxm1*dxm2**3
		  h1= dxp2**2+dxp2*dxm2
		  h2= dxp2**3-dxp2*dxm2**2
		  h3= dxp2**4+dxp2*dxm2**3
		  ca1= (f1*h3-f3*h1)*(g2*h3-g3*h2)-(g1*h3-g3*h1)*(f2*h3-f3*h2)
		  ax2_irz(jx,jz)=2.*h3*(g2*h3-g3*h2)/ca1
		  bx2_irz(jx,jz)=2.*h3*(h2*f3-h3*f2)/ca1
		  cx2_irz(jx,jz)=2.*h3*(f2*g3-f3*g2)/ca1
		  dx2_irz(jx,jz)=-(dxp1*ax2_irz(jx,jz)+dxm1*bx2_irz(jx,jz) &
			    +dxp2*cx2_irz(jx,jz))/dxm2
	endif

   11 continue


! for bndx with az2_irx,bz2_irx,...
	do 12 jx=ix_first,ix_last
	do 12 jz=iz_first+2,iz_last-2
	if(gdtp_ep(jx,jz,2).eq.2) then ! means upper two xgrid is the bnd point
		  itag=gdtp_ep(jx,jz,3)
		  xzp2=bndx_grd(itag,2)
		  xzp1=zz(jz+1)
		  xz00=zz(jz)
		  xzm1=zz(jz-1)
		  xzm2=zz(jz-2)

		  dzp1=xzp1-xz00
		  dzm1=xz00-xzm1
		  dzp2=xzp2-xz00
		  dzm2=xz00-xzm2
		  f1= dzp1**2+dzp1*dzm2
		  f2= dzp1**3-dzp1*dzm2**2
		  f3= dzp1**4+dzp1*dzm2**3
		  g1=-dzm1**2+dzm1*dzm2
		  g2= dzm1**3-dzm1*dzm2**2
		  g3=-dzm1**4+dzm1*dzm2**3
		  h1= dzp2**2+dzp2*dzm2
		  h2= dzp2**3-dzp2*dzm2**2
		  h3= dzp2**4+dzp2*dzm2**3
		  ca1= (f1*h3-f3*h1)*(g2*h3-g3*h2)-(g1*h3-g3*h1)*(f2*h3-f3*h2)
		  az2_irx(jx,jz)=2.*h3*(g2*h3-g3*h2)/ca1
		  bz2_irx(jx,jz)=2.*h3*(h2*f3-h3*f2)/ca1
		  cz2_irx(jx,jz)=2.*h3*(f2*g3-f3*g2)/ca1
		  dz2_irx(jx,jz)=-(dzp1*az2_irx(jx,jz)+dzm1*bz2_irx(jx,jz) &
			    +dzp2*cz2_irx(jx,jz))/dzm2
	endif

	if(gdtp_ep(jx,jz,2).eq.-2) then ! means lower two xgrid is the bnd point
		  itag=gdtp_ep(jx,jz,3)
		  xzp2=zz(jz+2)
		  xzp1=zz(jz+1)
		  xz00=zz(jz)
		  xzm1=zz(jz-1)
		  xzm2=bndx_grd(itag,2)


		  dzp1=xzp1-xz00
		  dzm1=xz00-xzm1
		  dzp2=xzp2-xz00
		  dzm2=xz00-xzm2
		  f1= dzp1**2+dzp1*dzm2
		  f2= dzp1**3-dzp1*dzm2**2
		  f3= dzp1**4+dzp1*dzm2**3
		  g1=-dzm1**2+dzm1*dzm2
		  g2= dzm1**3-dzm1*dzm2**2
		  g3=-dzm1**4+dzm1*dzm2**3
		  h1= dzp2**2+dzp2*dzm2
		  h2= dzp2**3-dzp2*dzm2**2
		  h3= dzp2**4+dzp2*dzm2**3
		  ca1= (f1*h3-f3*h1)*(g2*h3-g3*h2)-(g1*h3-g3*h1)*(f2*h3-f3*h2)
		  az2_irx(jx,jz)=2.*h3*(g2*h3-g3*h2)/ca1
		  bz2_irx(jx,jz)=2.*h3*(h2*f3-h3*f2)/ca1
		  cz2_irx(jx,jz)=2.*h3*(f2*g3-f3*g2)/ca1
		  dz2_irx(jx,jz)=-(dzp1*az2_irx(jx,jz)+dzm1*bz2_irx(jx,jz) &
			    +dzp2*cz2_irx(jx,jz))/dzm2
	endif
   12 continue



! calculate the free boundary coefficient for axm_bndz, bxm_bndz, cxm_bndz in x direction
	do 13 jz=iz_first,iz_last
	do 13 jx=ix_first+2,ix_last-2
	if(gdtp_ep(jx,jz,5).eq.1) then ! jx+1 is the boundary point, use the jx, jx-1, jx-2 to calculate the bndz point
		  itag=gdtp_ep(jx,jz,6)
		  dxm1=bndz_grd(itag,1)-xx(jx)
		  dxm2=bndz_grd(itag,1)-xx(jx-1)
		  dxm3=bndz_grd(itag,1)-xx(jx-2)
		  f1=dxm1-dxm1**3/dxm3**2
		  f2=dxm1**2-dxm1**3/dxm3
		  g1=dxm2-dxm2**3/dxm3**2
		  g2=dxm2**2-dxm2**3/dxm3
		  ca1=f1*g2-f2*g1
		  axm_bndz(itag)=g2/ca1
		  bxm_bndz(itag)=-f2/ca1
		  cxm_bndz(itag)=(1-axm_bndz(itag)*dxm1-bxm_bndz(itag)*dxm2)/dxm3
	endif
	if(gdtp_ep(jx,jz,5).eq.-1) then ! jx-1 is the boundary point, use the jx,jx+1,jx+2 to calculate the bndz point
		  itag=gdtp_ep(jx,jz,6)
		  dxp1=xx(jx)-bndz_grd(itag,1)
		  dxp2=xx(jx+1)-bndz_grd(itag,1)
		  dxp3=xx(jx+2)-bndz_grd(itag,1)
		  f1=dxp1-dxp1**3/dxp3**2
		  f2=dxp1**2-dxp1**3/dxp3
		  g1=dxp2-dxp2**3/dxp3**2
		  g2=dxp2**2-dxp2**3/dxp3
		  ca1=f1*g2-f2*g1
		  axp_bndz(itag)=g2/ca1
		  bxp_bndz(itag)=-f2/ca1
		  cxp_bndz(itag)=(1-axp_bndz(itag)*dxp1-bxp_bndz(itag)*dxp2)/dxp3
	endif
   13 continue

! calculate the free boundary coefficient for azm_bndx, bzm_bndx, czm_bndx in z direction
	do 14 jx=ix_first,ix_last
	do 14 jz=iz_first+2,iz_last-2
	if(gdtp_ep(jx,jz,2).eq.1) then ! jz+1 is the boundary point, use the jz, jz-1, jz-2 to calculate the bndx point
		  itag=gdtp_ep(jx,jz,3)
		  dzm1=bndx_grd(itag,2)-zz(jz)
		  dzm2=bndx_grd(itag,2)-zz(jz-1)
		  dzm3=bndx_grd(itag,2)-zz(jz-2)
		  f1=dzm1-dzm1**3/dzm3**2
		  f2=dzm1**2-dzm1**3/dzm3
		  g1=dzm2-dzm2**3/dzm3**2
		  g2=dzm2**2-dzm2**3/dzm3
		  ca1=f1*g2-f2*g1
		  azm_bndx(itag)=g2/ca1
		  bzm_bndx(itag)=-f2/ca1
		  czm_bndx(itag)=(1-azm_bndx(itag)*dzm1-bzm_bndx(itag)*dzm2)/dzm3
	endif
	if(gdtp_ep(jx,jz,2).eq.-1) then ! jz-1 is the boundary point, use the jz, jz+1, jz+2 to calculate the bndx point
		  itag=gdtp_ep(jx,jz,3)
		  dzp1=zz(jz)-bndx_grd(itag,2)
		  dzp2=zz(jz+1)-bndx_grd(itag,2)
		  dzp3=zz(jz+2)-bndx_grd(itag,2)
		  f1=dzp1-dzp1**3/dzp3**2
		  f2=dzp1**2-dzp1**3/dzp3
		  g1=dzp2-dzp2**3/dzp3**2
		  g2=dzp2**2-dzp2**3/dzp3
		  ca1=f1*g2-f2*g1
		  azp_bndx(itag)=g2/ca1
		  bzp_bndx(itag)=-f2/ca1
		  czp_bndx(itag)=(1-azp_bndx(itag)*dzp1-bzp_bndx(itag)*dzp2)/dzp3
	endif
   14 continue

    	return
	end
	
!hw*****************************************************************
	subroutine decide_hyperbolic_ratio
	use declare
	implicit none
	real*8 dxtt,dztt,dxzmax
    real*8 hypb_value,tmp1,tmp2
    real*8 rr_tmp
	include 'mpif.h'
	
	dztt=(maxval(zzt)-minval(zzt))/(mzt-1.d0)
	dxtt=(maxval(xxt)-minval(xxt))/(mxt-1.d0)
	dxzmax=max(dztt,dxtt)
!    hypb_value=0.98d0
    tmp1=1.d0-dztt*(2.5d0+grd_type_ratio)*2.d0/(maxval(zzt)-minval(zzt))
    tmp2=1.d0-dxtt*(2.5d0+grd_type_ratio)*2.d0/(maxval(xxt)-minval(xxt))
    hypb_value=min(tmp1,tmp2)

    if(nrank.eq.0) print*,'hypb_value=',hypb_value
	

	hypb_ratio(:,:)=0.d0
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
!      if(psi(jx,jz).lt.psia1) then
      if(gdtp_ep(jx,jz,1).ne.4) then
	if(rr(jx,jz).ge.0.d0) then
!	hypb_ratio(jx,jz)=max(0.d0,dtanh((1.d0-rr(jx,jz))*20.d0))
!	hypb_ratio(jx,jz)=max(0.d0,dtanh((hypb_value-rr(jx,jz))*30.d0))
! test only for circle case, need to expand to East case, etc.
    rr_tmp=dsqrt((xx(jx)-xzero)**2+(zz(jz)-zzero)**2) 
	hypb_ratio(jx,jz)=max(0.d0,dtanh((hypb_value-rr_tmp)*30.d0))
	else
	hypb_ratio(jx,jz)=1.d0
	endif
	endif
	enddo
	enddo
		
!	hypb_ratio(:,:)=1.d0
	if(nrank.eq.0) then
	open(unit=1,file='rrt.dat',status='unknown',form='formatted')
	write(1,2)((xxt(jx),zzt(jz),rrt(jx,jz),jx=1,mxt),jz=1,mzt)
    2 format(3(1x,e12.5)) 
    	close(1)

	open(unit=3,file='hypb.dat',status='unknown',form='formatted')
	write(3,4)((xx(jx),zz(jz),rr(jx,jz),hypb_ratio(jx,jz),x1(jx,jz,jy,7),jx=1,mx),jz=1,mz)
    4 format(5(1x,e12.5)) 
    	close(3)
	endif

	return
	end

!hw*****************************************************************
	subroutine decide_hyperbolic_ratio_v2
	use declare
	implicit none
	real*8 dxtt,dztt,dxzmax
    real*8 hypb_value,tmp1,tmp2
    real*8 rr_tmp
	include 'mpif.h'
	
    call distance_to_bnd
	dztt=(maxval(zzt)-minval(zzt))/(mzt-1.d0)
	dxtt=(maxval(xxt)-minval(xxt))/(mxt-1.d0)
	dxzmax=max(dztt,dxtt)
!    hypb_value=0.98d0
!    tmp1=dztt*(2.d0+grd_type_ratio)*2.d0/(maxval(zzt)-minval(zzt))
!    tmp2=dxtt*(2.d0+grd_type_ratio)*2.d0/(maxval(xxt)-minval(xxt))
    tmp1=dztt*(2.5d0+grd_type_ratio)/zdim
    tmp2=dxtt*(2.5d0+grd_type_ratio)/xdim
    hypb_value=max(tmp1,tmp2)

    if(nrank.eq.0) print*,'hypb_value=',hypb_value
	

	hypb_ratio(:,:)=0.d0
      do jz=iz_first,iz_last
      do jx=ix_first,ix_last
!      if(psi(jx,jz).lt.psia1) then
      if(gdtp_ep(jx,jz,1).ne.4) then
	if(rr(jx,jz).ge.0.d0) then
!	hypb_ratio(jx,jz)=max(0.d0,dtanh((1.d0-rr(jx,jz))*20.d0))
!	hypb_ratio(jx,jz)=max(0.d0,dtanh((hypb_value-rr(jx,jz))*30.d0))
! test only for circle case, need to expand to East case, etc.
!    rr_tmp=dsqrt((xx(jx)-xzero)**2+(zz(jz)-zzero)**2) 
!	hypb_ratio(jx,jz)=max(0.d0,dtanh((hypb_value-rr_tmp)*30.d0))
    hypb_ratio(jx,jz)=max(0.d0,dtanh((dist_to_bnd(jx,jz)-hypb_value)*30.d0))
    
	else
	hypb_ratio(jx,jz)=1.d0
	endif
	endif
	enddo
	enddo
		
!	hypb_ratio(:,:)=1.d0
	if(nrank.eq.0) then
	open(unit=1,file='rrt.dat',status='unknown',form='formatted')
	write(1,2)((xxt(jx),zzt(jz),rrt(jx,jz),jx=1,mxt),jz=1,mzt)
    2 format(3(1x,e12.5)) 
    	close(1)

	open(unit=3,file='hypb.dat',status='unknown',form='formatted')
	write(3,4)((xx(jx),zz(jz),rr(jx,jz),hypb_ratio(jx,jz),x1(jx,jz,jy,7),jx=1,mx),jz=1,mz)
    4 format(5(1x,e12.5)) 
    	close(3)
	endif

	return
	end

!hw**********************************************************************
	subroutine smth_irpt_with_difc(kk,f_type,ms,me,avrgh0)
      use declare
	integer itag,k,kk,f_type,ms,me
	real*8 tmp1,tmp2,difc_tmp
    real*8 avrgh0
      include 'mpif.h'
!
! second-order diffusion
      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)
! the dropped point will not be used for the ir points in its dropped direction,
! but its difc will be calculated if it is dropped only in one direction.

	do k=1,kk
	do 1 m=ms,me
	do 1 jy=iy_first,iy_last
	do 2 jz=iz_first+1,iz_last-1
	do 2 jx=ix_first+1,ix_last-1
	if((gdtp_ep(jx,jz,1).lt.4).and.(gdtp_ep(jx,jz,4).lt.4)) then
! 1st calculate the difc for x direction regular point and the inside ir point, 
! that is gdtp_ep(jx,jz,4)=1, or gdtp_ep(jx,jz,4)=2 and gdtp_ep(jx,jz,5)=-+2
!	if((gdtp_ep(jx,jz,4).eq.1).or.(gdtp_ep(jx,jz,5).eq.2).or.(gdtp_ep(jx,jz,5).eq.-2)) then
	if(f_type.eq.1) then ! for x1
	wdifcx(jx,jz)=difc(x1(jx-1,jz,jy,m),x1(jx,jz,jy,m),x1(jx+1,jz,jy,m),xx(jx-1),xx(jx),xx(jx+1))
	else if(f_type.eq.2) then ! for ef
	wdifcx(jx,jz)=difc(ef(jx-1,jz,jy,m),ef(jx,jz,jy,m),ef(jx+1,jz,jy,m),xx(jx-1),xx(jx),xx(jx+1))
	endif
!	endif

! then calculate the difc for outside ir point and inside dropped point
	if((gdtp_ep(jx,jz,5).eq.1).or.((gdtp_ep(jx,jz,4).eq.5).and.(gdtp_ep(jx-1,jz,5).eq.1))) then
	itag=gdtp_ep(jx,jz,6)

	if(f_type.eq.1) then ! for x1
	wdifcx(jx,jz)=difc(x1(jx-1,jz,jy,m),x1(jx,jz,jy,m),x1_8bndz(itag,jy,m),xx(jx-1),xx(jx),bndz_grd(itag,1))
	else if(f_type.eq.2) then ! for ef
	wdifcx(jx,jz)=difc(ef(jx-1,jz,jy,m),ef(jx,jz,jy,m),ef_3bndz(itag,jy,m),xx(jx-1),xx(jx),bndz_grd(itag,1))
	endif

	else if((gdtp_ep(jx,jz,5).eq.-1).or.((gdtp_ep(jx,jz,4).eq.5).and.(gdtp_ep(jx+1,jz,5).eq.-1))) then
	itag=gdtp_ep(jx,jz,6)

	if(f_type.eq.1) then ! for x1
	wdifcx(jx,jz)=difc(x1_8bndz(itag,jy,m),x1(jx,jz,jy,m),x1(jx+1,jz,jy,m),bndz_grd(itag,1),xx(jx),xx(jx+1))
	else if(f_type.eq.2) then ! for ef
	wdifcx(jx,jz)=difc(ef_3bndz(itag,jy,m),ef(jx,jz,jy,m),ef(jx+1,jz,jy,m),bndz_grd(itag,1),xx(jx),xx(jx+1))
	endif

	endif


! 2nd calculate the difc for z direction regular point and the inside ir point,
! that is gdtp_ep(jx,jz,1)=1, or gdtp_ep(jx,jz,4)=2 and gdtp_ep(jx,jz,5)=-+2
!	if((gdtp_ep(jx,jz,1).eq.1).or.(gdtp_ep(jx,jz,2).eq.2).or.(gdtp_ep(jx,jz,5).eq.-2)) then
	if(f_type.eq.1) then ! for x1
	wdifcz(jx,jz)=difc(x1(jx,jz-1,jy,m),x1(jx,jz,jy,m),x1(jx,jz+1,jy,m),zz(jz-1),zz(jz),zz(jz+1))
	else if(f_type.eq.2) then ! for ef
	wdifcz(jx,jz)=difc(ef(jx,jz-1,jy,m),ef(jx,jz,jy,m),ef(jx,jz+1,jy,m),zz(jz-1),zz(jz),zz(jz+1))
	endif
!	endif
	
! then calculate the difc for ouside ir point and inside dropped point
	if((gdtp_ep(jx,jz,2).eq.1).or.((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz-1,2).eq.1))) then
	itag=gdtp_ep(jx,jz,3)

	if(f_type.eq.1) then ! for x1
	wdifcz(jx,jz)=difc(x1(jx,jz-1,jy,m),x1(jx,jz,jy,m),x1_8bndx(itag,jy,m),zz(jz-1),zz(jz),bndx_grd(itag,2))
	else if(f_type.eq.2) then ! for ef
	wdifcz(jx,jz)=difc(ef(jx,jz-1,jy,m),ef(jx,jz,jy,m),ef_3bndx(itag,jy,m),zz(jz-1),zz(jz),bndx_grd(itag,2))
	endif

	else if((gdtp_ep(jx,jz,2).eq.-1).or.((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz+1,2).eq.-1))) then
	itag=gdtp_ep(jx,jz,3)

	if(f_type.eq.1) then ! for x1
	wdifcz(jx,jz)=difc(x1_8bndx(itag,jy,m),x1(jx,jz,jy,m),x1(jx,jz+1,jy,m),bndx_grd(itag,2),zz(jz),zz(jz+1))
	else if(f_type.eq.2) then ! for ef
	wdifcz(jx,jz)=difc(ef_3bndx(itag,jy,m),ef(jx,jz,jy,m),ef(jx,jz+1,jy,m),bndx_grd(itag,2),zz(jz),zz(jz+1))
	endif

	endif

!	tmp1=(xx(jx)-xzero)**2/((xx(jx)-xzero)**2+(zz(jz)-zzero)**2)
!	tmp2=(zz(jz)-zzero)**2/((xx(jx)-xzero)**2+(zz(jz)-zzero)**2)

	difc_tmp=(1.d0-avrgh0)*0.25*(1-hypb_ratio(jx,jz))	
	if(f_type.eq.1) then ! for x1
	x1(jx,jz,jy,m)=x1(jx,jz,jy,m)+difc_tmp*(wdifcx(jx,jz)+wdifcz(jx,jz))
	else if(f_type.eq.2) then ! for ef
	ef(jx,jz,jy,m)=ef(jx,jz,jy,m)+difc_tmp*(wdifcx(jx,jz)+wdifcz(jx,jz))
	endif

	endif
    2 continue
    1 continue

	if(f_type.eq.1) then ! for x1
!      call mpi_transfersm(x1(:,:,:,:),8) 
      call mpi_transfersm_one_layer(x1(:,:,:,:),8)
	else if(f_type.eq.2) then ! for ef
!	call mpi_transfersm(ef(:,:,:,:),3)
      call mpi_transfersm_one_layer(ef(:,:,:,:),3)
	endif

	enddo

	return
	end

!hw**********************************************************************
	subroutine smth_irpt_with_difc_v2(kk,f_type,ms,me,avrgh0)
      use declare
	integer itag,k,kk,f_type,ms,me
    real*8 wdifcx1,wdifcx2,wdifcx3,wdifcz1,wdifcz2,wdifcz3
	real*8 tmp1,tmp2,difc_tmp
    real*8 avrgh0
      include 'mpif.h'
!
! second-order diffusion
      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)
! the dropped point will not be used for the ir points in its dropped direction,
! but its difc will be calculated if it is dropped only in one direction.

!      print*, 'success 23846', nstep
	if(f_type.eq.1) then ! for x1
	do k=1,kk
	do 1 m=ms,me
	do 1 jy=iy_first,iy_last
	do 1 jz=iz_first+1,iz_last-1
	do 1 jx=ix_first+1,ix_last-1
	itag=gdtp_ep(jx,jz,6)

	wdifcx1=difc(x1(jx-1,jz,jy,m),x1(jx,jz,jy,m),x1(jx+1,jz,jy,m),xx(jx-1),xx(jx),xx(jx+1))
	wdifcx2=difc(x1(jx-1,jz,jy,m),x1(jx,jz,jy,m),x1_8bndz(itag,jy,m),xx(jx-1),xx(jx),bndz_grd(itag,1))
	wdifcx3=difc(x1_8bndz(itag,jy,m),x1(jx,jz,jy,m),x1(jx+1,jz,jy,m),bndz_grd(itag,1),xx(jx),xx(jx+1))
	wdifcz1=difc(x1(jx,jz-1,jy,m),x1(jx,jz,jy,m),x1(jx,jz+1,jy,m),zz(jz-1),zz(jz),zz(jz+1))
	wdifcz2=difc(x1(jx,jz-1,jy,m),x1(jx,jz,jy,m),x1_8bndx(itag,jy,m),zz(jz-1),zz(jz),bndx_grd(itag,2))
	wdifcz3=difc(x1_8bndx(itag,jy,m),x1(jx,jz,jy,m),x1(jx,jz+1,jy,m),bndx_grd(itag,2),zz(jz),zz(jz+1))

    wdifcx(jx,jz)=type_weight(jx,jz,1)*wdifcx1+type_weight(jx,jz,2)*wdifcx2+type_weight(jx,jz,3)*wdifcx3
    wdifcz(jx,jz)=type_weight(jx,jz,4)*wdifcz1+type_weight(jx,jz,5)*wdifcz2+type_weight(jx,jz,6)*wdifcz3

	difc_tmp=(1.d0-avrgh0)*0.25*(1-hypb_ratio(jx,jz))	
	x1(jx,jz,jy,m)=x1(jx,jz,jy,m)+difc_tmp*(wdifcx(jx,jz)+wdifcz(jx,jz))
    1 continue
      call mpi_transfersm_one_layer(x1(:,:,:,:),8)
	enddo

	else if(f_type.eq.2) then ! for ef

	do k=1,kk
	do 2 m=ms,me
	do 2 jy=iy_first,iy_last
	do 2 jx=ix_first+1,ix_last-1
	do 2 jz=iz_first+1,iz_last-1
	itag=gdtp_ep(jx,jz,6)

	wdifcx1=difc(ef(jx-1,jz,jy,m),ef(jx,jz,jy,m),ef(jx+1,jz,jy,m),xx(jx-1),xx(jx),xx(jx+1))
	wdifcx2=difc(ef(jx-1,jz,jy,m),ef(jx,jz,jy,m),ef_3bndz(itag,jy,m),xx(jx-1),xx(jx),bndz_grd(itag,1))
	wdifcx3=difc(ef_3bndz(itag,jy,m),ef(jx,jz,jy,m),ef(jx+1,jz,jy,m),bndz_grd(itag,1),xx(jx),xx(jx+1))
	wdifcz1=difc(ef(jx,jz-1,jy,m),ef(jx,jz,jy,m),ef(jx,jz+1,jy,m),zz(jz-1),zz(jz),zz(jz+1))
	wdifcz2=difc(ef(jx,jz-1,jy,m),ef(jx,jz,jy,m),ef_3bndx(itag,jy,m),zz(jz-1),zz(jz),bndx_grd(itag,2))
	wdifcz3=difc(ef_3bndx(itag,jy,m),ef(jx,jz,jy,m),ef(jx,jz+1,jy,m),bndx_grd(itag,2),zz(jz),zz(jz+1))

    wdifcx(jx,jz)=type_weight(jx,jz,1)*wdifcx1+type_weight(jx,jz,2)*wdifcx2+type_weight(jx,jz,3)*wdifcx3
    wdifcz(jx,jz)=type_weight(jx,jz,4)*wdifcz1+type_weight(jx,jz,5)*wdifcz2+type_weight(jx,jz,6)*wdifcz3

	difc_tmp=(1.d0-avrgh0)*0.25*(1-hypb_ratio(jx,jz))	
	ef(jx,jz,jy,m)=ef(jx,jz,jy,m)+difc_tmp*(wdifcx(jx,jz)+wdifcz(jx,jz))
    2 continue
      call mpi_transfersm_one_layer(ef(:,:,:,:),3)
    enddo
	endif

!      print*, 'success 23897', nstep

	return
	end

!hw**********************************************************************
	subroutine smth_irpt_with_difc_v3(kk,f_type,ms,me,avrgh0)
      use declare
	integer itag,k,kk,f_type,ms,me
	real*8 tmp1,tmp2,difc_tmp
	real*8 avrgh0
!	real*8, dimension(my,8) :: wdix,wdiz
      include 'mpif.h'
!
! second-order diffusion
      difc(fm1,f0,fp1,xm1,x0,xp1)= &
       2.*( (fp1-f0)*(x0-xm1)-(f0-fm1)*(xp1-x0))/(xp1-xm1)
! the dropped point will not be used for the ir points in its dropped direction,
! but its difc will be calculated if it is dropped only in one direction.

	do k=1,kk
	do 1 jz=iz_first+1,iz_last-1
	do 1 jx=ix_first+1,ix_last-1
	if((gdtp_ep(jx,jz,1).lt.4).and.(gdtp_ep(jx,jz,4).lt.4)) then
! 1st calculate the difc for x direction regular point and the inside ir point, 
! that is gdtp_ep(jx,jz,4)=1, or gdtp_ep(jx,jz,4)=2 and gdtp_ep(jx,jz,5)=-+2
!	if((gdtp_ep(jx,jz,4).eq.1).or.(gdtp_ep(jx,jz,5).eq.2).or.(gdtp_ep(jx,jz,5).eq.-2)) then
	if(f_type.eq.1) then ! for x1
	wdix(:,1:8)=2.*((x1(jx+1,jz,:,1:8)-x1(jx,jz,:,1:8))*(xx(jx)-xx(jx-1)) &
		  -(x1(jx,jz,:,1:8)-x1(jx-1,jz,:,1:8))*(xx(jx+1)-xx(jx)))/(xx(jx+1)-xx(jx-1))
!	wdifcx(jx,jz)=difc(x1(jx-1,jz,jy,m),x1(jx,jz,jy,m),x1(jx+1,jz,jy,m),xx(jx-1),xx(jx),xx(jx+1))
	else if(f_type.eq.2) then ! for ef
	wdix(:,1:3)=2.*((ef(jx+1,jz,:,1:3)-ef(jx,jz,:,1:3))*(xx(jx)-xx(jx-1)) &
		  -(ef(jx,jz,:,1:3)-ef(jx-1,jz,:,1:3))*(xx(jx+1)-xx(jx)))/(xx(jx+1)-xx(jx-1))
!	wdifcx(jx,jz)=difc(ef(jx-1,jz,jy,m),ef(jx,jz,jy,m),ef(jx+1,jz,jy,m),xx(jx-1),xx(jx),xx(jx+1))
	endif
!	endif

! then calculate the difc for outside ir point and inside dropped point
	if((gdtp_ep(jx,jz,5).eq.1).or.((gdtp_ep(jx,jz,4).eq.5).and.(gdtp_ep(jx-1,jz,5).eq.1))) then
	itag=gdtp_ep(jx,jz,6)

	if(f_type.eq.1) then ! for x1
	wdix(:,1:8)=2.*((x1_8bndz(itag,:,1:8)-x1(jx,jz,:,1:8))*(xx(jx)-xx(jx-1)) &
		  -(x1(jx,jz,:,1:8)-x1(jx-1,jz,:,1:8))*(bndz_grd(itag,1)-xx(jx)))/(bndz_grd(itag,1)-xx(jx-1))
!	wdifcx(jx,jz)=difc(x1(jx-1,jz,jy,m),x1(jx,jz,jy,m),x1_8bndz(itag,jy,m),xx(jx-1),xx(jx),bndz_grd(itag,1))
	else if(f_type.eq.2) then ! for ef
	wdix(:,1:3)=2.*((ef_3bndz(itag,:,1:3)-ef(jx,jz,:,1:3))*(xx(jx)-xx(jx-1)) &
		  -(ef(jx,jz,:,1:3)-ef(jx-1,jz,:,1:3))*(bndz_grd(itag,1)-xx(jx)))/(bndz_grd(itag,1)-xx(jx-1))
!	wdifcx(jx,jz)=difc(ef(jx-1,jz,jy,m),ef(jx,jz,jy,m),ef_3bndz(itag,jy,m),xx(jx-1),xx(jx),bndz_grd(itag,1))
	endif

	else if((gdtp_ep(jx,jz,5).eq.-1).or.((gdtp_ep(jx,jz,4).eq.5).and.(gdtp_ep(jx+1,jz,5).eq.-1))) then
	itag=gdtp_ep(jx,jz,6)

	if(f_type.eq.1) then ! for x1
	wdix(:,1:8)=2.*((x1(jx+1,jz,:,1:8)-x1(jx,jz,:,1:8))*(xx(jx)-bndz_grd(itag,1)) &
		  -(x1(jx,jz,:,1:8)-x1_8bndz(itag,:,1:8))*(xx(jx+1)-xx(jx)))/(xx(jx+1)-bndz_grd(itag,1))
!	wdifcx(jx,jz)=difc(x1_8bndz(itag,jy,m),x1(jx,jz,jy,m),x1(jx+1,jz,jy,m),bndz_grd(itag,1),xx(jx),xx(jx+1))
	else if(f_type.eq.2) then ! for ef
	wdix(:,1:3)=2.*((ef(jx+1,jz,:,1:3)-ef(jx,jz,:,1:3))*(xx(jx)-bndz_grd(itag,1)) &
		  -(ef(jx,jz,:,1:3)-ef_3bndz(itag,:,1:3))*(xx(jx+1)-xx(jx)))/(xx(jx+1)-bndz_grd(itag,1))
!	wdifcx(jx,jz)=difc(ef_3bndz(itag,jy,m),ef(jx,jz,jy,m),ef(jx+1,jz,jy,m),bndz_grd(itag,1),xx(jx),xx(jx+1))
	endif

	endif


! 2nd calculate the difc for z direction regular point and the inside ir point,
! that is gdtp_ep(jx,jz,1)=1, or gdtp_ep(jx,jz,4)=2 and gdtp_ep(jx,jz,5)=-+2
!	if((gdtp_ep(jx,jz,1).eq.1).or.(gdtp_ep(jx,jz,2).eq.2).or.(gdtp_ep(jx,jz,5).eq.-2)) then
	if(f_type.eq.1) then ! for x1
	wdiz(:,1:8)=2.*((x1(jx,jz+1,:,1:8)-x1(jx,jz,:,1:8))*(zz(jz)-zz(jz-1)) &
		  -(x1(jx,jz,:,1:8)-x1(jx,jz-1,:,1:8))*(zz(jz+1)-zz(jz)))/(zz(jz+1)-zz(jz-1))
!	wdifcz(jx,jz)=difc(x1(jx,jz-1,jy,m),x1(jx,jz,jy,m),x1(jx,jz+1,jy,m),zz(jz-1),zz(jz),zz(jz+1))
	else if(f_type.eq.2) then ! for ef
	wdiz(:,1:3)=2.*((ef(jx,jz+1,:,1:3)-ef(jx,jz,:,1:3))*(zz(jz)-zz(jz-1)) &
		  -(ef(jx,jz,:,1:3)-ef(jx,jz-1,:,1:3))*(zz(jz+1)-zz(jz)))/(zz(jz+1)-zz(jz-1))
!	wdifcz(jx,jz)=difc(ef(jx,jz-1,jy,m),ef(jx,jz,jy,m),ef(jx,jz+1,jy,m),zz(jz-1),zz(jz),zz(jz+1))
	endif
!	endif
	
! then calculate the difc for ouside ir point and inside dropped point
	if((gdtp_ep(jx,jz,2).eq.1).or.((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz-1,2).eq.1))) then
	itag=gdtp_ep(jx,jz,3)

	if(f_type.eq.1) then ! for x1
	wdiz(:,1:8)=2.*((x1_8bndx(itag,:,1:8)-x1(jx,jz,:,1:8))*(zz(jz)-zz(jz-1)) &
		  -(x1(jx,jz,:,1:8)-x1(jx,jz-1,:,1:8))*(bndx_grd(itag,2)-zz(jz)))/(bndx_grd(itag,2)-zz(jz-1))
!	wdifcz(jx,jz)=difc(x1(jx,jz-1,jy,m),x1(jx,jz,jy,m),x1_8bndx(itag,jy,m),zz(jz-1),zz(jz),bndx_grd(itag,2))
	else if(f_type.eq.2) then ! for ef
	wdiz(:,1:3)=2.*((ef_3bndx(itag,:,1:3)-ef(jx,jz,:,1:3))*(zz(jz)-zz(jz-1)) &
		  -(ef(jx,jz,:,1:3)-ef(jx,jz-1,:,1:3))*(bndx_grd(itag,2)-zz(jz)))/(bndx_grd(itag,2)-zz(jz-1))
!	wdifcz(jx,jz)=difc(ef(jx,jz-1,jy,m),ef(jx,jz,jy,m),ef_3bndx(itag,jy,m),zz(jz-1),zz(jz),bndx_grd(itag,2))
	endif

	else if((gdtp_ep(jx,jz,2).eq.-1).or.((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz+1,2).eq.-1))) then
	itag=gdtp_ep(jx,jz,3)

	if(f_type.eq.1) then ! for x1
	wdiz(:,1:8)=2.*((x1(jx,jz+1,:,1:8)-x1(jx,jz,:,1:8))*(zz(jz)-bndx_grd(itag,2)) &
		  -(x1(jx,jz,:,1:8)-x1_8bndx(itag,:,1:8))*(zz(jz+1)-zz(jz)))/(zz(jz+1)-bndx_grd(itag,2))
!	wdifcz(jx,jz)=difc(x1_8bndx(itag,jy,m),x1(jx,jz,jy,m),x1(jx,jz+1,jy,m),bndx_grd(itag,2),zz(jz),zz(jz+1))
	else if(f_type.eq.2) then ! for ef
	wdiz(:,1:3)=2.*((ef(jx,jz+1,:,1:3)-ef(jx,jz,:,1:3))*(zz(jz)-bndx_grd(itag,2)) &
		  -(ef(jx,jz,:,1:3)-ef_3bndx(itag,:,1:3))*(zz(jz+1)-zz(jz)))/(zz(jz+1)-bndx_grd(itag,2))
!	wdifcz(jx,jz)=difc(ef_3bndx(itag,jy,m),ef(jx,jz,jy,m),ef(jx,jz+1,jy,m),bndx_grd(itag,2),zz(jz),zz(jz+1))
	endif

	endif

!	tmp1=(xx(jx)-xzero)**2/((xx(jx)-xzero)**2+(zz(jz)-zzero)**2)
!	tmp2=(zz(jz)-zzero)**2/((xx(jx)-xzero)**2+(zz(jz)-zzero)**2)

	difc_tmp=(1.d0-avrgh0)*0.25*(1-hypb_ratio(jx,jz))	
	if(f_type.eq.1) then ! for x1
	x1(jx,jz,:,1:8)=x1(jx,jz,:,1:8)+difc_tmp*(wdix(:,1:8)+wdiz(:,1:8))
	else if(f_type.eq.2) then ! for ef
	ef(jx,jz,:,1:3)=ef(jx,jz,:,1:3)+difc_tmp*(wdix(:,1:3)+wdiz(:,1:3))
	endif

	endif
    1 continue

	if(f_type.eq.1) then ! for x1
!      call mpi_transfersm(x1(:,:,:,:),8) 
      call mpi_transfersm_one_layer(x1(:,:,:,:),8)
	else if(f_type.eq.2) then ! for ef
!	call mpi_transfersm(ef(:,:,:,:),3)
      call mpi_transfersm_one_layer(ef(:,:,:,:),3)
	endif

	enddo

	return
	end

!hw*****************************************************************************
      subroutine distance_to_bnd
      use declare
      implicit none
!      real*8, dimension(mx,mz) :: dist_to_bnd
!      real*8, dimension(mx,mz) :: theta_ep
!      real*8 theta_tmp
      integer id1,id2,id3,id4,i,j,js
      real*8 thet1,thet2,thet3,thet4
      real*8 rr11,rr22,rr33,rr44
      real*8, dimension(mx,mz) :: rrmax_ep
      character*9 output
      character*3 cn
      include 'mpif.h'

      do 1 jx=ix_first,ix_last
      do 1 jz=iz_first,iz_last
!      if((psi(jx,jz).lt.psia1)) then
      if(gdtp_ep(jx,jz,1).ne.4) then

      theta_ep(jx,jz)=datan2(zz(jz)-zzero,xx(jx)-xzero)
      if(theta_ep(jx,jz).lt.0.d0) theta_ep(jx,jz)=theta_ep(jx,jz)+2.d0*pi

      do i=1,ngrdb
      id1=i-1
      id2=i
      id3=i+1
      id4=i+2

      if(i.eq.1) then
      id1=ngrdb
      endif

      if(i.eq.ngrdb-1) then
      id4=1
      endif

      if(i.eq.ngrdb) then
           id3=1
           id4=2
      endif

      thet1=datan2(zb(id1)-zzero,xb(id1)-xzero)
      if(thet1.lt.0.d0) thet1=thet1+2.d0*pi
      thet2=datan2(zb(id2)-zzero,xb(id2)-xzero)
      if(thet2.lt.0.d0) thet2=thet2+2.d0*pi
      thet3=datan2(zb(id3)-zzero,xb(id3)-xzero)
      if(thet3.lt.0.d0) thet3=thet3+2.d0*pi
      thet4=datan2(zb(id4)-zzero,xb(id4)-xzero)
      if(thet4.lt.0.d0) thet4=thet4+2.d0*pi

!hw: make sure the thet1, thet2, thet3, thet4 is from small to large
      if(i.eq.1) then
!      id1=ngrdb
      if(thet1.gt.thet2) then
      thet1=thet1-2.d0*pi
      endif
      endif


      if(i.eq.ngrdb-2) then
      if(thet4.lt.thet3) then
      thet4=thet4+2.d0*pi
      endif
      endif


      if(i.eq.ngrdb-1) then
!      id4=1
      if(thet2.gt.thet3) then
      thet3=thet3+2.d0*pi
      endif
      if(thet3.gt.thet4) then
      thet4=thet4+2.d0*pi
      endif
      endif

      if(i.eq.ngrdb) then
!           id3=1
!           id4=2
      if(thet3.lt.thet2) then
      thet2=thet2-2.d0*pi
      endif
      if(thet2.lt.thet1) then
      thet1=thet1-2.d0*pi
      endif
      endif

      rr11=dsqrt(((xb(id1)-xzero)/xdim)**2+((zb(id1)-zzero)/zdim)**2)
      rr22=dsqrt(((xb(id2)-xzero)/xdim)**2+((zb(id2)-zzero)/zdim)**2)
      rr33=dsqrt(((xb(id3)-xzero)/xdim)**2+((zb(id3)-zzero)/zdim)**2)
      rr44=dsqrt(((xb(id4)-xzero)/xdim)**2+((zb(id4)-zzero)/zdim)**2)

      if(((thet2.le.theta_ep(jx,jz)).and.(thet3.ge.theta_ep(jx,jz))).or. &
              ((thet2.ge.theta_ep(jx,jz)).and.(thet3.le.theta_ep(jx,jz)))) then
      call lag_intp1d4p(thet1,rr11,thet2,rr22,thet3,rr33,thet4,rr44,theta_ep(jx,jz),rrmax_ep(jx,jz))

      dist_to_bnd(jx,jz)=dabs(dabs(rrmax_ep(jx,jz))-dsqrt(((xx(jx)-xzero)/xdim)**2+((zz(jz)-zzero)/zdim)**2))
      dist_to_bnd(jx,jz)=min(dist_to_bnd(jx,jz),0.8d0)
      exit
      endif


      enddo
    
      else
      dist_to_bnd(jx,jz)=1.d0
      endif

   1 continue


      print*,'ngrdb=',ngrdb
      output='dstan'//cn(nrank)
      open(unit=3,file=output,status='unknown',form='formatted')
      write(3,4)((xx(jx),zz(jz),dist_to_bnd(jx,jz),rrmax_ep(jx,jz),jx=1,mx),jz=1,mz)
    4 format(4(1x,e12.5)) 
      close(3)
      
      return
      end

!hw**************************************************************
	subroutine data_type_weight
      use declare
!    real*8, dimension(mx,mz,5) :: type_weight
! 1 2 3 x direction type, 4 5 6 z direction type
	integer itag,k,kk,f_type,ms,me
	real*8 tmp1,tmp2,difc_tmp
    real*8 avrgh0
      include 'mpif.h'
!
    type_weight(:,:,:)=0.d0

	do 2 jx=ix_first,ix_last
	do 2 jz=iz_first,iz_last

	if((gdtp_ep(jx,jz,1).lt.4).and.(gdtp_ep(jx,jz,4).lt.4)) then
! 1st calculate the difc for x direction regular point and the inside ir point, 
! that is gdtp_ep(jx,jz,4)=1, or gdtp_ep(jx,jz,4)=2 and gdtp_ep(jx,jz,5)=-+2
!	if((gdtp_ep(jx,jz,4).eq.1).or.(gdtp_ep(jx,jz,5).eq.2).or.(gdtp_ep(jx,jz,5).eq.-2)) then
    type_weight(jx,jz,1)=1.d0


! then calculate the difc for outside ir point and inside dropped point
	if((gdtp_ep(jx,jz,5).eq.1).or.((gdtp_ep(jx,jz,4).eq.5).and.(gdtp_ep(jx-1,jz,5).eq.1))) then
    type_weight(jx,jz,2)=1.d0

	else if((gdtp_ep(jx,jz,5).eq.-1).or.((gdtp_ep(jx,jz,4).eq.5).and.(gdtp_ep(jx+1,jz,5).eq.-1))) then

    type_weight(jx,jz,3)=1.d0

	endif


! 2nd calculate the difc for z direction regular point and the inside ir point,
! that is gdtp_ep(jx,jz,1)=1, or gdtp_ep(jx,jz,4)=2 and gdtp_ep(jx,jz,5)=-+2
!	if((gdtp_ep(jx,jz,1).eq.1).or.(gdtp_ep(jx,jz,2).eq.2).or.(gdtp_ep(jx,jz,5).eq.-2)) then
    type_weight(jx,jz,4)=1.d0
	
! then calculate the difc for ouside ir point and inside dropped point
	if((gdtp_ep(jx,jz,2).eq.1).or.((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz-1,2).eq.1))) then
    type_weight(jx,jz,5)=1.d0

	else if((gdtp_ep(jx,jz,2).eq.-1).or.((gdtp_ep(jx,jz,1).eq.5).and.(gdtp_ep(jx,jz+1,2).eq.-1))) then

    type_weight(jx,jz,6)=1.d0

	endif

!	tmp1=(xx(jx)-xzero)**2/((xx(jx)-xzero)**2+(zz(jz)-zzero)**2)
!	tmp2=(zz(jz)-zzero)**2/((xx(jx)-xzero)**2+(zz(jz)-zzero)**2)

	endif
    2 continue

	return
	end

! 2017.12.13
!hw**********************************************************************
!hw**********************************************************************
!hw**********************************************************************
!hw**********************************************************************
!hw**********************************************************************
	subroutine find_bnd_grd
	use declare
	implicit none
	integer id1,id2,id3,id4,i,j,js
      include 'mpif.h'
!	real*8, dimension(:,:), allocatable :: bnd_grd

	if(.not.initia_from_pgfile) then
	xxs5(1:2*(nthe-1)+5,:)=xxs(:,mpsa-4:mpsa)
	zzs5(1:2*(nthe-1)+5,:)=zzs(:,mpsa-4:mpsa)
	endif

!      open(unit=161,file='xxs',status='unknown',form='formatted')
!      read(161,500)((xxs5(i,js),js=mpsa5,1,-1),i=1,n2th+5)
!      close(161)
!      open(unit=162,file='zzs',status='unknown',form='formatted')
!      read(162,500)((zzs5(i,js),js=mpsa5,1,-1),i=1,n2th+5)
!      close(162)
!  500 format(5(1x,e12.5)) 



!      open(unit=201,file='gridxx.dat',status='unknown',form='formatted')
!      read(201,200)(xxt(jx),jx=1,mxt)
!      close(201)
!      open(unit=202,file='gridzz.dat',status='unknown',form='formatted')
!      read(202,200)(zzt(jz),jz=1,mzt)
!      close(202)
!  200 format(1x,e12.5)
 	
	if(nrank.eq.0) then
      open(unit=203,file='xxt.dat',status='unknown',form='formatted')
      write(203,201)(xxt(jx),jx=1,mxt)
      close(203)
      open(unit=204,file='zzt.dat',status='unknown',form='formatted')
      write(204,201)(zzt(jz),jz=1,mzt)
      close(204)
  201 format(1x,e12.5)
  	endif

 	do i=1,nxzs
	theta(i)=dasin((zzs5(i,5)-zzero)/dsqrt((zzs5(i,5)-zzero)**2+(xxs5(i,5)-xzero)**2))
	if(xxs5(i,5).le.xzero) theta(i)=pi-theta(i)
	if((xxs5(i,5).gt.xzero).and.(zzs5(i,5).le.zzero)) theta(i)=2.d0*pi+theta(i)
!	print*, theta(i)*180/pi,xxs(i,5),zzs(i,5)
	enddo
	
	theta_tmp=theta
	call shell_sort(nxzs,theta_tmp,mpsa5,xxs5)
	theta_tmp=theta
	call shell_sort(nxzs,theta_tmp,mpsa5,zzs5)
	theta=theta_tmp

! 	do i=1,205
!	print*, theta_tmp(i)*180/pi,xxs(i,5),zzs(i,5)
!	enddo

	if(nrank.eq.0) then
      open(unit=163,file='xxs.dat',status='unknown',form='formatted')
      write(163,501)((xxs(i,js),js=mpsa,mpsa-4,-1),i=1,n2th+5)
      close(163)
      open(unit=164,file='zzs.dat',status='unknown',form='formatted')
      write(164,501)((zzs(i,js),js=mpsa,mpsa-4,-1),i=1,n2th+5)
      close(164)
  501 format(5(1x,e12.5)) 
  	endif


! delete the repeated points in theta, xxs, zzs
! the left ngrdb grids are unique, xb(1:ngrdb)
	xb(:)=xxs5(:,mpsa5)
	zb(:)=zzs5(:,mpsa5)
	ngrdb=0
	do i=1,nxzs-1
!	print*,i,theta(i),theta(i+1)
	do while(theta(i).eq.theta(i+1))
	theta(i:nxzs-1)=theta(i+1:nxzs)
	xb(i:nxzs-1)=xb(i+1:nxzs)
	zb(i:nxzs-1)=zb(i+1:nxzs)
	theta(nxzs)=(ngrdb+1)*10.d0
	ngrdb=ngrdb+1
	enddo
	enddo
	ngrdb=nxzs-ngrdb
!	print*,xb,'haha',zb

	if(nrank.eq.0) then
	open(unit=191,file='xb.dat',status='unknown',form='formatted')
	write(191,502)(xb(i),i=1,ngrdb)
	close(191)
	open(unit=192,file='zb.dat',status='unknown',form='formatted')
	write(192,502)(zb(i),i=1,ngrdb)
	close(192)
  502 format(1(1x,e12.5)) 
  	endif

! the sortted bounary for theta from small to large, that is xb and zb
!	xb(:)=xxs(:,5)
!	zb(:)=zzs(:,5)
!	ngrdb=size(xb)
!	ngrdx=size(xxt)
!	ngrdz=size(zzt)
	

	nbndx=0
	nbndz=0
	do i=1,ngrdb
	id1=i-1
	id2=i
	id3=i+1
	id4=i+2
	if(i.eq.1) then
		  id1=ngrdb
!		  print*,id1,id2,id3,id4
	endif
	if(i.eq.ngrdb-1) then
		  id4=1
!		  print*,id1,id2,id3,id4
	endif
	if(i.eq.ngrdb) then
		  id3=1
		  id4=2
!		  print*,id1,id2,id3,id4
	endif
! find z boundary point, that is the intersection of xb-zb and z grid lines, zzt
	do j=1,mzt
	if(((zb(id2).le.zzt(j)).and.(zb(id3).gt.zzt(j))).or.((zb(id2).ge.zzt(j)).and.(zb(id3).lt.zzt(j)))) then
	nbndz=nbndz+1
	call lag_intp1d4p(zb(id1),xb(id1),zb(id2),xb(id2),zb(id3),xb(id3),zb(id4),xb(id4),zzt(j),bnd_z(nbndz,1))
	bnd_z(nbndz,2)=zzt(j)
	bnd_z(nbndz,3:6)=(/id1,id2,id3,id4/)
	bnd_z(nbndz,7)=j
!	exit
	endif
	enddo


! find x boundary point, that is the intersection of xb-zb and x grid lines, xxt
	do j=1,mxt
	if(((xb(id2).le.xxt(j)).and.(xb(id3).gt.xxt(j))).or.((xb(id2).ge.xxt(j)).and.(xb(id3).lt.xxt(j)))) then
	nbndx=nbndx+1
	call lag_intp1d4p(xb(id1),zb(id1),xb(id2),zb(id2),xb(id3),zb(id3),xb(id4),zb(id4),xxt(j),bnd_x(nbndx,2))
	bnd_x(nbndx,1)=xxt(j)
	bnd_x(nbndx,3:6)=(/id1,id2,id3,id4/)
	bnd_x(nbndx,7)=j
!	exit
	endif
	enddo

	enddo


	bnd_x(nbndx+1:m2xt,:)=1.d7
	bnd_tmpx(:)=bnd_x(:,7)
	call shell_sort(m2xt,bnd_tmpx,n7,bnd_x)

	bnd_z(nbndz+1:m2zt,:)=1.d7
	bnd_tmpz(:)=bnd_z(:,7)
	call shell_sort(m2zt,bnd_tmpz,n7,bnd_z)


	if(nrank.eq.0) then
      open(unit=169,file='bndx.dat',status='unknown',form='formatted')
      open(unit=170,file='bndz.dat',status='unknown',form='formatted')
      write(169,509)((bnd_x(i,js),js=1,7),i=1,nbndx)
      write(170,509)((bnd_z(i,js),js=1,7),i=1,nbndz)
      close(169)
      close(170)
  509 format(7(1x,e12.5)) 
  	endif

	call decide_grd_type_bndx	
	call decide_grd_type_bndz	
	if(nrank.eq.0) then
	open(unit=193,file='grd_typebndx.dat',status='unknown',form='formatted')
	open(unit=194,file='grd_typebndz.dat',status='unknown',form='formatted')
	write(193,511)((xxt(jx),zzt(jz),grd_type(jx,jz,1)*1.d0,jx=1,mxt),jz=1,mzt)
	write(194,511)((xxt(jx),zzt(jz),grd_type(jx,jz,2)*1.d0,jx=1,mxt),jz=1,mzt)
	close(193)
	close(194)
  511 format(3(1x,e12.5)) 
	endif
! finally we got bnd_x z and grd_type	

!	call decide_grd_type_in_each_proc


  	return
 	end


!hw**********************************************************************
	subroutine find_bnd_grd_in_each_proc
	use declare
	implicit none
	integer id1,id2,id3,id4,i,j,js
      include 'mpif.h'
!	real*8, dimension(:,:), allocatable :: bnd_grd

	if(.not.initia_from_pgfile) then
	xxs5(1:2*(nthe-1)+5,:)=xxs(:,mpsa-4:mpsa)
	zzs5(1:2*(nthe-1)+5,:)=zzs(:,mpsa-4:mpsa)
	endif

!      open(unit=161,file='xxs',status='unknown',form='formatted')
!      read(161,500)((xxs5(i,js),js=mpsa5,1,-1),i=1,n2th+5)
!      close(161)
!      open(unit=162,file='zzs',status='unknown',form='formatted')
!      read(162,500)((zzs5(i,js),js=mpsa5,1,-1),i=1,n2th+5)
!      close(162)
!  500 format(5(1x,e12.5)) 



!      open(unit=201,file='gridxx.dat',status='unknown',form='formatted')
!      read(201,200)(xxt(jx),jx=1,mxt)
!      close(201)
!      open(unit=202,file='gridzz.dat',status='unknown',form='formatted')
!      read(202,200)(zzt(jz),jz=1,mzt)
!      close(202)
!  200 format(1x,e12.5)
 	
	if(nrank.eq.0) then
      open(unit=203,file='xxn0.dat',status='unknown',form='formatted')
      write(203,201)(xx(jx),jx=1,mx)
      close(203)
      open(unit=204,file='zzn0.dat',status='unknown',form='formatted')
      write(204,201)(zz(jz),jz=1,mz)
      close(204)
  201 format(1x,e12.5)
  	endif

 	do i=1,nxzs
	theta(i)=dasin((zzs5(i,5)-zzero)/dsqrt((zzs5(i,5)-zzero)**2+(xxs5(i,5)-xzero)**2))
	if(xxs5(i,5).le.xzero) theta(i)=pi-theta(i)
	if((xxs5(i,5).gt.xzero).and.(zzs5(i,5).le.zzero)) theta(i)=2.d0*pi+theta(i)
!	print*, theta(i)*180/pi,xxs(i,5),zzs(i,5)
	enddo
	
	theta_tmp=theta
	call shell_sort(nxzs,theta_tmp,mpsa5,xxs5)
	theta_tmp=theta
	call shell_sort(nxzs,theta_tmp,mpsa5,zzs5)
	theta=theta_tmp

! 	do i=1,205
!	print*, theta_tmp(i)*180/pi,xxs(i,5),zzs(i,5)
!	enddo

	if(nrank.eq.0) then
      open(unit=163,file='xxs.dat',status='unknown',form='formatted')
      write(163,501)((xxs(i,js),js=mpsa,mpsa-4,-1),i=1,n2th+5)
      close(163)
      open(unit=164,file='zzs.dat',status='unknown',form='formatted')
      write(164,501)((zzs(i,js),js=mpsa,mpsa-4,-1),i=1,n2th+5)
      close(164)
  501 format(5(1x,e12.5)) 
  	endif


! delete the repeated points in theta, xxs, zzs
! the left ngrdb grids are unique, xb(1:ngrdb)
	xb(:)=xxs5(:,mpsa5)
	zb(:)=zzs5(:,mpsa5)
	ngrdb=0
	do i=1,nxzs-1
!	print*,i,theta(i),theta(i+1)
	do while(theta(i).eq.theta(i+1))
	theta(i:nxzs-1)=theta(i+1:nxzs)
	xb(i:nxzs-1)=xb(i+1:nxzs)
	zb(i:nxzs-1)=zb(i+1:nxzs)
	theta(nxzs)=(ngrdb+1)*10.d0
	ngrdb=ngrdb+1
	enddo
	enddo
	ngrdb=nxzs-ngrdb
!	print*,xb,'haha',zb

	if(nrank.eq.0) then
	open(unit=191,file='xb.dat',status='unknown',form='formatted')
	write(191,502)(xb(i),i=1,ngrdb)
	close(191)
	open(unit=192,file='zb.dat',status='unknown',form='formatted')
	write(192,502)(zb(i),i=1,ngrdb)
	close(192)
  502 format(1(1x,e12.5)) 
  	endif

! the sortted bounary for theta from small to large, that is xb and zb
!	xb(:)=xxs(:,5)
!	zb(:)=zzs(:,5)
!	ngrdb=size(xb)
!	ngrdx=size(xxt)
!	ngrdz=size(zzt)
	

	nbndx_ep=0
	nbndz_ep=0
	do i=1,ngrdb
	id1=i-1
	id2=i
	id3=i+1
	id4=i+2
	if(i.eq.1) then
		  id1=ngrdb
!		  print*,id1,id2,id3,id4
	endif
	if(i.eq.ngrdb-1) then
		  id4=1
!		  print*,id1,id2,id3,id4
	endif
	if(i.eq.ngrdb) then
		  id3=1
		  id4=2
!		  print*,id1,id2,id3,id4
	endif
! find z boundary point, that is the intersection of xb-zb and z grid lines, zzt
	do j=1,mz
	if(((zb(id2).le.zz(j)).and.(zb(id3).gt.zz(j))).or.((zb(id2).ge.zz(j)).and.(zb(id3).lt.zz(j)))) then
	nbndz_ep=nbndz_ep+1
	call lag_intp1d4p(zb(id1),xb(id1),zb(id2),xb(id2),zb(id3),xb(id3),zb(id4),xb(id4),zz(j),bnd_z_ep(nbndz_ep,1))
	bnd_z_ep(nbndz_ep,2)=zz(j)
	bnd_z_ep(nbndz_ep,3:6)=(/id1,id2,id3,id4/)
	bnd_z_ep(nbndz_ep,7)=j
!	exit
	endif
	enddo


! find x boundary point, that is the intersection of xb-zb and x grid lines, xxt
	do j=1,mx
	if(((xb(id2).le.xx(j)).and.(xb(id3).gt.xx(j))).or.((xb(id2).ge.xx(j)).and.(xb(id3).lt.xx(j)))) then
	nbndx_ep=nbndx_ep+1
	call lag_intp1d4p(xb(id1),zb(id1),xb(id2),zb(id2),xb(id3),zb(id3),xb(id4),zb(id4),xx(j),bnd_x_ep(nbndx_ep,2))
	bnd_x_ep(nbndx_ep,1)=xx(j)
	bnd_x_ep(nbndx_ep,3:6)=(/id1,id2,id3,id4/)
	bnd_x_ep(nbndx_ep,7)=j
!	exit
	endif
	enddo

	enddo

	bnd_x_ep(nbndx_ep+1:m2x,:)=1.d7
	bnd_tmpx_ep(:)=bnd_x_ep(:,7)
	call shell_sort(m2x,bnd_tmpx_ep,n7,bnd_x_ep)
	bnd_z_ep(nbndz_ep+1:m2z,:)=1.d7
	bnd_tmpz_ep(:)=bnd_z_ep(:,7)
	call shell_sort(m2z,bnd_tmpz_ep,n7,bnd_z_ep)

	if(nrank.eq.0) then
      open(unit=169,file='bndx_ep.dat',status='unknown',form='formatted')
      open(unit=170,file='bndz_ep.dat',status='unknown',form='formatted')
      write(169,509)((bnd_x_ep(i,js),js=1,7),i=1,nbndx_ep)
      write(170,509)((bnd_z_ep(i,js),js=1,7),i=1,nbndz_ep)
      close(169)
      close(170)
  509 format(7(1x,e12.5)) 
  	endif

!	call decide_grd_type_bndx	
!	call decide_grd_type_bndz	
!	if(nrank.eq.0) then
!	open(unit=193,file='grd_typebndx.dat',status='unknown',form='formatted')
!	open(unit=194,file='grd_typebndz.dat',status='unknown',form='formatted')
!	write(193,511)(((xxt(jx),zzt(jz),grd_type(jx,jz,1)*1.d0),jx=1,mxt),jz=1,mzt)
!	write(194,511)(((xxt(jx),zzt(jz),grd_type(jx,jz,2)*1.d0),jx=1,mxt),jz=1,mzt)
!	close(193)
!	close(194)
!  511 format(3(1x,e12.5)) 
!	endif
!! finally we got bnd_x z and grd_type	

!	call decide_grd_type_in_each_proc


  	return
 	end

!hw**********************************************************************
! find the bnd_x point for the IR point in mxt mzt in z direction
	subroutine find_ir_pnt_bndx
	use declare
	implicit none
	real*8 deltaxz,dxtt,dztt
	include 'mpif.h'
	dztt=(maxval(zzt)-minval(zzt))/(mzt-1.d0)
	deltaxz=1.0001d0*(1.d0+grd_type_ratio)*dztt
	gdtp_bndx(:,:,1)=0 ! -> -2 -1 0 +1 +2 the nearest bnd pnt for IR pnt
	gdtp_bndx(:,:,2)=0 ! -> 1 ~ nbndx corresponding to point in bndx_grd
	do jx=1,mxt
	do jz=1,mzt
	if(grd_type(jx,jz,1).eq.2) then 
		  ! its a IR point
		  if((jz.lt.mzt).and.(grd_type(jx,jz+1,1).eq.2).and.(jz.gt.1).and.(grd_type(jx,jz-1,1).gt.2)) then
		  ! its upper point is a IR point
			    do jr=1,nbndx
			    if((nint(bndx_grd(jr,7)).eq.jx).and.(dabs(zzt(jz)-bndx_grd(jr,2)).lt.deltaxz)) then
!			    if(nint(bndx_grd(jr,7)).eq.jx) then
		                  gdtp_bndx(jx,jz,1)=-1
           	                  gdtp_bndx(jx,jz+1,1)=-2
					! jr is the boundary point
					gdtp_bndx(jx,jz,2)=jr
					gdtp_bndx(jx,jz+1,2)=jr
			    endif
			    enddo
		  endif

		  if((jz.gt.1).and.(grd_type(jx,jz-1,1).eq.2).and.(jz.lt.mzt).and.(grd_type(jx,jz+1,1).gt.2)) then
		  ! its lower point is a IR point
			    do jr=1,nbndx
			    if((nint(bndx_grd(jr,7)).eq.jx).and.(dabs(bndx_grd(jr,2)-zzt(jz)).lt.deltaxz)) then
!			    if(nint(bndx_grd(jr,7)).eq.jx) then
					gdtp_bndx(jx,jz,1)=1
					gdtp_bndx(jx,jz-1,1)=2
					! jr is the boundary point
					gdtp_bndx(jx,jz,2)=jr
					gdtp_bndx(jx,jz-1,2)=jr
			    endif
			    enddo
		  endif
	endif
	enddo
	enddo

    	return
    	end

!hw**********************************************************************
! find the bnd_z point for the IR point in mxt mzt in x direction
	subroutine find_ir_pnt_bndz
	use declare
	implicit none
	real*8 deltaxz,dxtt,dztt
	include 'mpif.h'
	dxtt=(maxval(xxt)-minval(xxt))/(mxt-1.d0)
	deltaxz=1.0001d0*(1.d0+grd_type_ratio)*dxtt
	gdtp_bndz(:,:,1)=0 ! -> -2 -1 0 +1 +2 the nearest bnd pnt for IR pnt
	gdtp_bndz(:,:,2)=0 ! -> 1 ~ nbndz corresponding to point in bndz_grd
	do jz=1,mzt
	do jx=1,mxt
	if(grd_type(jx,jz,2).eq.2) then 
		  ! its a IR point
		  if((jx.lt.mxt).and.(grd_type(jx+1,jz,2).eq.2).and.(jx.gt.1).and.(grd_type(jx-1,jz,2).gt.2)) then
			    do jr=1,nbndz
			    if((nint(bndz_grd(jr,7)).eq.jz).and.(dabs(xxt(jx)-bndz_grd(jr,1)).lt.deltaxz)) then
			    		! its upper point is a IR point
					gdtp_bndz(jx,jz,1)=-1
					gdtp_bndz(jx+1,jz,1)=-2
					! jr is the boundary point
					gdtp_bndz(jx,jz,2)=jr
					gdtp_bndz(jx+1,jz,2)=jr
!					exit
			    endif
			    enddo
		  endif

		  if((jx.gt.1).and.(grd_type(jx-1,jz,2).eq.2).and.(jx.lt.mxt).and.(grd_type(jx+1,jz,2).gt.2)) then
			    do jr=1,nbndz
			    if((nint(bndz_grd(jr,7)).eq.jz).and.(dabs(bndz_grd(jr,1)-xxt(jx)).lt.deltaxz)) then
			    		! its lower point is a IR point
					gdtp_bndz(jx,jz,1)=1
					gdtp_bndz(jx-1,jz,1)=2
					! jr is the boundary point
					gdtp_bndz(jx,jz,2)=jr
					gdtp_bndz(jx-1,jz,2)=jr
!					exit
			    endif
			    enddo
		  endif
	endif
	enddo
	enddo
    	return
    	end

!hw**********************************************************************
! decide the type for bnd_x in z direction	
	subroutine decide_grd_type_bndx
	use declare
	implicit none
	integer i
	real*8 bnd,bnd1,bnd2,dxtt,dztt
      include 'mpif.h'
	dztt=(maxval(zzt)-minval(zzt))/(mzt-1.d0)
! set to dropped point at first
	grd_type(:,:,1)=4
! decide the type for bnd_x in z direction	
	do 1 jx=1,mxt
	do 2 i=1,nbndx-1

	if(((i.eq.1).and.(nint(bnd_x(i,7)).ne.nint(bnd_x(i+1,7)))).or. &
	((i.eq.nbndx-1).and.(nint(bnd_x(i,7)).ne.nint(bnd_x(i+1,7))))) then
! the first and last bndx point which is unique and tangency with the circle
	do 3 jz=1,mzt
	if((zzt(jz).eq.bnd_x(1,2)).or.(zzt(jz).eq.bnd_x(nbndx,2))) then
! boundary point, but is unique
!! boundary point, too close to the boundary, set as dropped point directly
!	grd_type(jx,jz,1)=3
!	grd_type(jx,jz,1)=4
	grd_type(jx,jz,1)=5 !inside but dropped
	endif
    3 continue


	else

	if((nint(bnd_x(i,7)).eq.jx).and. &
	(nint(bnd_x(i,7)).eq.nint(bnd_x(i+1,7)))) then
		if(bnd_x(i,2).ge.bnd_x(i+1,2)) then
			  bnd1=bnd_x(i+1,2)
			  bnd2=bnd_x(i,2)
		else
			  bnd1=bnd_x(i,2)
			  bnd2=bnd_x(i+1,2)
		endif
	do 4 jz=1,mzt
! inside the boundary, set as regular point first
	if(((zzt(jz).gt.bnd1).and.(zzt(jz).lt.bnd2)).and. &
	(min(zzt(jz)-bnd1,bnd2-zzt(jz)).ge.grd_type_ratio*dztt)) then
		  grd_type(jx,jz,1)=1
	endif
!! boundary point, too close to the boundary, set as dropped point directly
	if(((zzt(jz).gt.bnd1).and.(zzt(jz).lt.bnd2)).and. &
	(min(zzt(jz)-bnd1,bnd2-zzt(jz)).lt.grd_type_ratio*dztt)) then
		  grd_type(jx,jz,1)=5 !inside but dropped
	endif
!! boundary point, too close to the boundary, set as dropped point directly
	if((zzt(jz).eq.bnd1).or.(zzt(jz).eq.bnd2)) then
!		  grd_type(jx,jz,1)=3
!		  grd_type(jx,jz,1)=4
		  grd_type(jx,jz,1)=5 !inside but dropped
	endif

    4 continue

    	do 5 jz=2,mzt-1 ! find the irregular points,  D D B (D) I I R R R R ... R R R I I (D) B D D 
	if((grd_type(jx,jz,1).eq.1).and.(grd_type(jx,jz+1,1).eq.1).and.(grd_type(jx,jz-1,1).gt.2.99d0)) then
		  grd_type(jx,jz,1)=2
		  grd_type(jx,jz+1,1)=2
	endif
	if((grd_type(jx,jz,1).eq.1).and.(grd_type(jx,jz-1,1).eq.1).and.(grd_type(jx,jz+1,1).gt.2.99d0)) then
		  grd_type(jx,jz,1)=2
		  grd_type(jx,jz-1,1)=2
	endif
    5 continue


	endif
	endif


    2 continue
    1 continue


	
	return
	end

!hw**********************************************************************
! decide the type for bnd_z in x direction	
	subroutine decide_grd_type_bndz
	use declare
	implicit none
	integer i
	real*8 bnd,bnd1,bnd2,dxtt,dztt
      include 'mpif.h'
	dxtt=(maxval(xxt)-minval(xxt))/(mxt-1.d0)
! set to dropped point at first
	grd_type(:,:,2)=4
! decide the type for bnd_z in x direction	
	do 1 jz=1,mzt
	do 2 i=1,nbndz-1

	if(((i.eq.1).and.(nint(bnd_z(i,7)).ne.nint(bnd_z(i+1,7)))).or. &
	((i.eq.nbndz-1).and.(nint(bnd_z(i,7)).ne.nint(bnd_z(i+1,7))))) then
! the first and last bndx point which is unique and tangency with the circle
	do 3 jx=1,mxt
	if((xxt(jx).eq.bnd_z(1,1)).or.(xxt(jx).eq.bnd_z(nbndz,1))) then
! boundary point, but is unique
!! boundary point, too close to the boundary, set as dropped point directly
!	grd_type(jx,jz,2)=3
!	grd_type(jx,jz,2)=4
	grd_type(jx,jz,2)=5 !inside but dropped
	endif
    3 continue


	else

	if((nint(bnd_z(i,7)).eq.jz).and. &
	(nint(bnd_z(i,7)).eq.nint(bnd_z(i+1,7)))) then
		if(bnd_z(i,1).ge.bnd_z(i+1,1)) then
			  bnd1=bnd_z(i+1,1)
			  bnd2=bnd_z(i,1)
		else
			  bnd1=bnd_z(i,1)
			  bnd2=bnd_z(i+1,1)
		endif
	do 4 jx=1,mxt
! inside the boundary, set as regular point first
	if(((xxt(jx).gt.bnd1).and.(xxt(jx).lt.bnd2)).and. &
	(min(xxt(jx)-bnd1,bnd2-xxt(jx)).ge.grd_type_ratio*dxtt)) then
		  grd_type(jx,jz,2)=1
	endif
!! boundary point, too close to the boundary, set as dropped point directly
	if(((xxt(jx).gt.bnd1).and.(xxt(jx).lt.bnd2)).and. &
	(min(xxt(jx)-bnd1,bnd2-xxt(jx)).lt.grd_type_ratio*dxtt)) then
		  grd_type(jx,jz,2)=5 !inside but dropped
	endif
!! boundary point, too close to the boundary, set as dropped point directly
	if((xxt(jx).eq.bnd1).or.(xxt(jx).eq.bnd2)) then
!		  grd_type(jx,jz,2)=3
!		  grd_type(jx,jz,2)=4
		  grd_type(jx,jz,2)=5 !inside but dropped
	endif

    4 continue

    	do 5 jx=2,mxt-1 ! find the irregular points,  D D B (D) I I R R R R ... R R R I I (D) B D D 
	if((grd_type(jx,jz,2).eq.1).and.(grd_type(jx+1,jz,2).eq.1).and.(grd_type(jx-1,jz,2).gt.2.99d0)) then
		  grd_type(jx,jz,2)=2
		  grd_type(jx+1,jz,2)=2
	endif
	if((grd_type(jx,jz,2).eq.1).and.(grd_type(jx-1,jz,2).eq.1).and.(grd_type(jx+1,jz,2).gt.2.99d0)) then
		  grd_type(jx,jz,2)=2
		  grd_type(jx-1,jz,2)=2
	endif
    5 continue


	endif
	endif


    2 continue
    1 continue


	
	return
	end
	
