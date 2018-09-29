# CLT

CLT (tokRZ_mpi.f90) with a full name of Ci-Liu-Ti (磁流体 in Chinese, and MHD in English), is a full MHD program developed for the tokamak MHD instablity simulations.

## Table of Contents
* [Development history](#Development-istory)
* [Main parameters and coordinates](#Main-parameters-and-coordinates)
* [Terms of use](#Terms-of-use)
* [License](#License)


## Development history:

1. The 1st version is developed by Wang, Sheng @Zhejiang University. 
    * **This version code is written in 2th order finite difference method in space,**
    * **while in the time-advance, 4th order Runge-Kutta scheme is chosen,**
    * **in the phi(y) direction, either finite difference or pseudo-spectrum method is used.**
    * **the equilibrium is loaded from the psi_xz.dat, q_p_g.dat, and wch.dat(optional)**
2. The 2nd version is upgraded by Zhang, Wei @Zhejiang University mainly with 
    * **4th order finite difference method is employed in the R and Z directions.**
3. The 3rd version is improved by Zhang, Haowei @Zhejiang University, following features are added:
    * **Cut-cell method is used in the boundary,**
    * **Can be used for any triangularity like East or the Circle case. Fixed boundary is used,**
    * **Can be used for any nprx and nprz ( grids and nprx(z) should be divisible),**
    * **Can be used for including SOL region from EFIT-gfile equilibriums,**
    * **Subroutines for RMP-EAST coils is added.**
    * **Experiment equilibriums can be used (eq_pgfile_*.dat are read, transform by Matlab from p&g files).**

! Ref: Physics of Plasmas 22, 122504 (2015); doi: 10.1063/1.4936977

## Main parameters and coordinates

1. mmode=2, nmode=1, qmode=mmode/nmode=2: set the q value of the mode you want to calculate.
2. hall=true or flase: open the Hall term or not.
3. mxt=256,myt=64,mzt=256: the grids numbers in each direction.
4. cfl=1.2: decide the CFL condition number for the time step.
5. x(1-8):total of rho,p,vx,vy(v_phi),vz,bx,by(b_phi),bz.
6. x1(1-8):pertubation of rho,p,vx,vy(v_phi),vz,bx,by(b_phi),bz.
7. xint(1-8):equilibrium of rho,p,vx,vy(v_phi),vz,bx,by(b_phi),bz.
8. cur(1-3):pertubation of current (x,y(phi),z).
9. cint(1-3):equilibrium of current _(x,y(phi),z).
10. ef(1-3):pertubation of e-field _(x,y(phi),z).
11. coordinate with large cylindrical coordinate system is used (R-\phi(Y)-Z)
    * **xx(mx),yy(my),zz(mz): coordinates(r,phi,z) in each processes;**
    * **xxt(mxt),yyt(myt),zzt(mzt) for total;** 
    * **xxst(n2th+5,npsi),zzst(n2th+5,npsi) in (theta,psi) grid;**
    * **xxs(n2th+5,mps4:mps),zzs(n2th+5,mps4:mps) in (theta,psi) bandary grid;**
    * **thxz(mx,mz): theta coordinates in (r,z) grid; tht(mxt,mzt) for total;**
    * **tpxz(mx,mz): r.z.<->s(psi).p(pol). transit angle in (r,z); tpt(mxt,mzt) for total;**
    * **tcxz(mx,mz): tc=ing(jcb/r^2)dth;**
    * **thst(n2th+5): theta coordinates in (theta,psi) grid;**
    * **tpst(n2th+5,npsi): r.z.<->s(psi).p(pol). transit angle in (theta,psi); tps(n2th+5,mps4:mps) for bndry;**

## Term of use

## License

