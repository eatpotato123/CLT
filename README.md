# CLT
                                                            
CLT (tokRZ_mpi.f90) with a full name of Ci-Liu-Ti (磁流体 in Chinese, and MHD in English), is a full MHD program developed for the 
tokamak MHD instablity simulations.

## Table of Contents
* [Results presentation](#results-presentation)
* [Development history](#development-history)
* [Main parameters and coordinates](#main-parameters-and-coordinates)
* [Terms of use](#terms-of-use)
* [License](#license)

## Results presentation：
1. m/n=3/1 islands calculated with CLT:
![Image text](https://github.com/changhw/CLT/blob/master/img-folder/export_fig_out.png)

2. Parameters for East shot#52340:
![Image text](https://github.com/changhw/CLT/blob/master/img-folder/shot%23052340.03150ke.png)

## Development history:
1. The 1st version is developed by Prof. Ma, Zhiwei & Dr. Wang, Sheng @Zhejiang University. 
    * **This version code is written in 2th order finite difference method in space,**
    * **while in the time-advance, 4th order Runge-Kutta scheme is chosen,**
    * **in the phi(y) direction, either finite difference or pseudo-spectrum method is used,**
    * **the equilibrium is loaded from the psi_xz.dat, q_p_g.dat, and wch.dat(optional).**
2. The 2nd version is upgraded by Dr. Zhang, Wei @Zhejiang University mainly with 
    * **4th order finite difference method is employed in the R and Z directions.**
3. The 3rd version is improved by Zhang, Haowei @Zhejiang University, following features are added:
    * **Cut-cell method is used in the boundary,**
    * **Can be used for any triangularity like East or the Circle case. Fixed boundary is used,**
    * **Can be used for any nprx and nprz ( grids and nprx(z) should be divisible),**
    * **Can be used for including SOL region from EFIT-gfile equilibriums,**
    * **Subroutines for RMP-EAST coils is added,**
    * **Experiment equilibriums can be used (eq_pgfile_*.dat are read, transform by Matlab from p&g files).**
4. Ref: Physics of Plasmas 22, 122504 (2015); doi: 10.1063/1.4936977

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
11. mxt: total X grids number.
12. myt: total Y grids number.
13. mzt: total Z grids number.
14. npsi_pgfile: number of data rows in the eq_pgfile_1d.dat.
15. nlim: number of data rows in the eq_pgfile_rzlim.dat.
16. nbs: number of data rows in the eq_pgfile_rzbs.dat.
17. nprx: the processes number in X direction.
18. nprz: the processes number in Z direction.
19. npry: the processes number in Y direction.
20. initia_from_pgfile: true->read the initia equilibrium from Efit files; false->read the Nova files.
21. use_stepon_cut_cell: true->use the cut-cell boundary by Zhang, Haowei; false->use the old boundary by Wang, Sheng.
22. rmp_east: true->open the East-RMP coils,
23. coordinate with large cylindrical coordinate system is used (R-\phi(Y)-Z)
    * **xx(mx),yy(my),zz(mz): coordinates(r,phi,z) in each processes,**
    * **xxt(mxt),yyt(myt),zzt(mzt) for total,** 
    * **xxst(n2th+5,npsi),zzst(n2th+5,npsi) in (theta,psi) grid,**
    * **xxs(n2th+5,mps4:mps),zzs(n2th+5,mps4:mps) in (theta,psi) bandary grid,**
    * **thxz(mx,mz): theta coordinates in (r,z) grid, tht(mxt,mzt) for total,**
    * **tpxz(mx,mz): r.z.<->s(psi).p(pol). transit angle in (r,z), tpt(mxt,mzt) for total,**
    * **tcxz(mx,mz): tc=ing(jcb/r^2)dth,**
    * **thst(n2th+5): theta coordinates in (theta,psi) grid,**
    * **tpst(n2th+5,npsi): r.z.<->s(psi).p(pol). transit angle in (theta,psi), tps(n2th+5,mps4:mps) for bndry.**

## Terms of use
1. CLT is a young but powerful scientific code, and we welcome peers around the world to use and improve this code.
2. Before you start to use this CLT code, please inform us by emailing Zhang, Haowei at changhw@zju.edu.cn and Prof. Ma, Zhiwei at zwma@zju.edu.cn so that we can show you how to use this code in detail.
3. If you want to modificate the CLT source code, or publish paper with any results calculated by CLT, please inform us by emailing Zhang, Haowei and Prof. Ma, Zhiwei, and we can help to check the modification or the result before we go to the next step.
4. Thank you for your cooperation!

## License
©2015-2035 Zhejiang University. All rights reserved.
