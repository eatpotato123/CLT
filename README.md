# CLT

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

! mxt=256,myt=64,mzt=256

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


