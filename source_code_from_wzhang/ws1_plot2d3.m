% clear;
xzero=5;
mx=256;
mz=256;
my=32;
ms=5;
mr=4;
mt=326;
aa=1;
xx=load('gridxx.dat');
zz=load('gridzz.dat');
yy=load('gridyy.dat');
time=load('nstt.dat');
xx=xx*aa;
zz=zz*aa;
xmin=min(xx);
xmax=max(xx);
zmin=min(zz);
zmax=max(zz);
xs=load('xxs');
zs=load('zzs');
xb=xs(:,1);
zb=zs(:,1);
xb=xb*aa;
zb=zb*aa;
[xz_xx,xz_zz]=meshgrid(xx,zz);
%xint=load('xint.dat');

% xzint_rho=reshape(xint(:,1),mx,mz)';
% xzint_p=reshape(xint(:,2),mx,mz)';
% xzint_vx=reshape(xint(:,3),mx,mz)';
% xzint_vy=reshape(xint(:,4),mx,mz)';
% xzint_vz=reshape(xint(:,5),mx,mz)';
% xzint_bx=reshape(xint(:,6),mx,mz)';
% xzint_by=reshape(xint(:,7),mx,mz)';
% xzint_bz=reshape(xint(:,8),mx,mz)';
% xzint_br=reshape(xint(:,9),mx,mz)';
% xzint_bp=reshape(xint(:,10),mx,mz)';
% cint=load('cint.dat');
% xzint_cx=reshape(cint(:,1),mx,mz)';
% xzint_cy=reshape(cint(:,2),mx,mz)';
% xzint_cz=reshape(cint(:,3),mx,mz)';

% fint=load('fint.dat');
% fx=reshape(fint(:,1),mx,mz)';
% fy=reshape(fint(:,2),mx,mz)';
% fz=reshape(fint(:,3),mx,mz)';
% eint=load('eta.dat');
% eta=reshape(eint(:,1),mx,mz)';
% etax=reshape(eint(:,2),mx,mz)';
% etaz=reshape(eint(:,3),mx,mz)';

% xz_rr=reshape(rr,mx,mz);
% xz_th=reshape(th,mx,mz);
[xz_xx,xz_zz]=meshgrid(xx,zz);
% % xz_xx=xz_xx';
% % xz_zz=xz_zz';
% [rt_at,rt_ar]=meshgrid(at,ar);
% [rt_xx,rt_zz]=pol2cart(rt_at,rt_ar);
% xs=reshape(xz_xx,mx*mz,1);
% zs=reshape(xz_zz,mx*mz,1);
% scatter(xs,zs,1);hold on;
% scatter(xxb,zzb,5,'r');hold on;
% alpha=0:pi/160:2*pi;%½Ç¶È[0,2*pi] 
% R=min(rrb);%°ë¾¶ 
% x=xzero+R*cos(alpha); 
% y=R*sin(alpha); 
% plot(x,y,'-');hold on;
% hold off;
mkdir png1;
ndst=1;
nstartime=0;
nstart=0;
colormap(jet);
for nst=nstart:ndst:100
str_nst=num2str(nst,'%.3d');
str_nst1=num2str(nst+1,'%.3d');
str_time=num2str(time((nst-nstartime)/ndst+1,2));
w=load(['x12d' str_nst ]);
xz_rho=reshape(w(:,1),mx,mz)';
xz_p=reshape(w(:,2),mx,mz)';
xz_vx=reshape(w(:,3),mx,mz)';
xz_vy=reshape(w(:,4),mx,mz)';
xz_vz=reshape(w(:,5),mx,mz)';
xz_bx=reshape(w(:,6),mx,mz)';
xz_by=reshape(w(:,7),mx,mz)';
xz_bz=reshape(w(:,8),mx,mz)';
xz_br=reshape(w(:,9),mx,mz)';
xz_bp=reshape(w(:,10),mx,mz)';
xz_vr=reshape(w(:,11),mx,mz)';
xz_vp=reshape(w(:,12),mx,mz)';
xz_vs=reshape(w(:,13),mx,mz)';
%if(mod(nst,10)==0)
% w=load(['db' str_nst1 ]);
% db=reshape(w(:,1),mx,mz)';
% d2db=reshape(w(:,2),mx,mz)';
%end
% w=load(['cur000' str_nst ]);
% xz_cx=reshape(w(:,1),mx,mz)';
% xz_cy=reshape(w(:,2),mx,mz)';
% xz_cz=reshape(w(:,3),mx,mz)';
% 
% w=load(['Ef000' str_nst ]);
% xz_ex=reshape(w(:,1),mx,mz)';
% xz_ey=reshape(w(:,2),mx,mz)';
% xz_ez=reshape(w(:,3),mx,mz)';
% 
% wy=load(['xy_x000' str_nst ]);
% xy_rho=reshape(wy(:,1),mx,my)';
% xy_p=reshape(wy(:,2),mx,my)';
% xy_vx=reshape(wy(:,3),mx,my)';
% xy_vy=reshape(wy(:,4),mx,my)';
% xy_vz=reshape(wy(:,5),mx,my)';
% xy_bx=reshape(wy(:,6),mx,my)';
% xy_by=reshape(wy(:,7),mx,my)';
% xy_bz=reshape(wy(:,8),mx,my)';
% xy_br=reshape(wy(:,9),mx,my)';
% xy_bp=reshape(wy(:,10),mx,my)';
% xy_vr=reshape(wy(:,11),mx,my)';
% xy_vp=reshape(wy(:,12),mx,my)';
% xy_vs=reshape(wy(:,13),mx,my)';
% 
% wy=load(['xy_c000' str_nst ]);
% xy_cx=reshape(wy(:,1),mx,my)';
% xy_cy=reshape(wy(:,2),mx,my)';
% xy_cz=reshape(wy(:,3),mx,my)';

% 
% wrt=load(['xrt000' str_nst ]);
% rt_rho=reshape(wrt(:,1),mr,mt);
% rt_p=reshape(wrt(:,2),mr,mt);
% rt_vx=reshape(wrt(:,3),mr,mt);
% rt_vy=reshape(wrt(:,4),mr,mt);
% rt_vz=reshape(wrt(:,5),mr,mt);
% rt_bx=reshape(wrt(:,6),mr,mt);
% rt_by=reshape(wrt(:,7),mr,mt);
% rt_bz=reshape(wrt(:,8),mr,mt);

% wdif=load(['xdif000' str_nst ]);
% xzdif_rho=reshape(wdif(:,1),mx,mz);
% xzdif_p=reshape(wdif(:,2),mx,mz);
% xzdif_vx=reshape(wdif(:,3),mx,mz);
% xzdif_vy=reshape(wdif(:,4),mx,mz);
% xzdif_vz=reshape(wdif(:,5),mx,mz);
% xzdif_bx=reshape(wdif(:,6),mx,mz);
% xzdif_by=reshape(wdif(:,7),mx,mz);
% xzdif_bz=reshape(wdif(:,8),mx,mz);
 
% wrtdif=load(['xrtdif000' str_nst ]);
% rtdif_rho=reshape(wrtdif(:,1),mr,mt);
% rtdif_p=reshape(wrtdif(:,2),mr,mt);
% rtdif_vx=reshape(wrtdif(:,3),mr,mt);
% rtdif_vy=reshape(wrtdif(:,4),mr,mt);
% rtdif_vz=reshape(wrtdif(:,5),mr,mt);
% rtdif_bx=reshape(wrtdif(:,6),mr,mt);
% rtdif_by=reshape(wrtdif(:,7),mr,mt);
% rtdif_bz=reshape(wrtdif(:,8),mr,mt);

% w=load(['x_dx000' str_nst ]);
% xz_rhodx=reshape(w(:,1),mx,mz);
% xz_pdx=reshape(w(:,2),mx,mz);
% xz_vxdx=reshape(w(:,3),mx,mz);
% xz_vydx=reshape(w(:,4),mx,mz);
% xz_vzdx=reshape(w(:,5),mx,mz);
% xz_bxdx=reshape(w(:,6),mx,mz);
% xz_bydx=reshape(w(:,7),mx,mz);
% xz_bzdx=reshape(w(:,8),mx,mz);
% 
% w=load(['x_dz000' str_nst ]);
% xz_rhodz=reshape(w(:,1),mx,mz);
% xz_pdz=reshape(w(:,2),mx,mz);
% xz_vxdz=reshape(w(:,3),mx,mz);
% xz_vydz=reshape(w(:,4),mx,mz);
% xz_vzdz=reshape(w(:,5),mx,mz);
% xz_bxdz=reshape(w(:,6),mx,mz);
% xz_bydz=reshape(w(:,7),mx,mz);
% xz_bzdz=reshape(w(:,8),mx,mz);
% 
% w=load(['x_dy000' str_nst ]);
% xz_rhody=reshape(w(:,1),mx,mz);
% xz_pdy=reshape(w(:,2),mx,mz);
% xz_vxdy=reshape(w(:,3),mx,mz);
% xz_vydy=reshape(w(:,4),mx,mz);
% xz_vzdy=reshape(w(:,5),mx,mz);
% xz_bxdy=reshape(w(:,6),mx,mz);
% xz_bydy=reshape(w(:,7),mx,mz);
% xz_bzdy=reshape(w(:,8),mx,mz);
% 
% wrt=load(['xrt_dr000' str_nst ]);
% rt_rhodr=reshape(wrt(:,1),mr,mt);
% rt_pdr=reshape(wrt(:,2),mr,mt);
% rt_vxdr=reshape(wrt(:,3),mr,mt);
% rt_vydr=reshape(wrt(:,4),mr,mt);
% rt_vzdr=reshape(wrt(:,5),mr,mt);
% rt_bxdr=reshape(wrt(:,6),mr,mt);
% rt_bydr=reshape(wrt(:,7),mr,mt);
% rt_bzdr=reshape(wrt(:,8),mr,mt);
% 
% wrt=load(['xrt_dt000' str_nst ]);
% rt_rhodt=reshape(wrt(:,1),mr,mt);
% rt_pdt=reshape(wrt(:,2),mr,mt);
% rt_vxdt=reshape(wrt(:,3),mr,mt);
% rt_vydt=reshape(wrt(:,4),mr,mt);
% rt_vzdt=reshape(wrt(:,5),mr,mt);
% rt_bxdt=reshape(wrt(:,6),mr,mt);
% rt_bydt=reshape(wrt(:,7),mr,mt);
% rt_bzdt=reshape(wrt(:,8),mr,mt);
% 
% wrt=load(['xrt_dy000' str_nst ]);
% rt_rhody=reshape(wrt(:,1),mr,mt);
% rt_pdy=reshape(wrt(:,2),mr,mt);
% rt_vxdy=reshape(wrt(:,3),mr,mt);
% rt_vydy=reshape(wrt(:,4),mr,mt);
% rt_vzdy=reshape(wrt(:,5),mr,mt);
% rt_bxdy=reshape(wrt(:,6),mr,mt);
% rt_bydy=reshape(wrt(:,7),mr,mt);
% rt_bzdy=reshape(wrt(:,8),mr,mt);
%if(nst>0) 
% w=load(['xdif' str_nst '001' ]);
% xzdif_rho_1=reshape(w(:,1),mx,mz)';
% xzdif_p_1=reshape(w(:,2),mx,mz)';
% xzdif_vx_1=reshape(w(:,3),mx,mz)';
% xzdif_vy_1=reshape(w(:,4),mx,mz)';
% xzdif_vz_1=reshape(w(:,5),mx,mz)';
% xzdif_bx_1=reshape(w(:,6),mx,mz)';
% xzdif_by_1=reshape(w(:,7),mx,mz)';
% xzdif_bz_1=reshape(w(:,8),mx,mz)';
% xzdif_p0_1=reshape(w(:,9),mx,mz)';
% w=load(['xdif' str_nst '002' ]);
% xzdif_rho_2=reshape(w(:,1),mx,mz)';
% xzdif_p_2=reshape(w(:,2),mx,mz)';
% xzdif_vx_2=reshape(w(:,3),mx,mz)';
% xzdif_vy_2=reshape(w(:,4),mx,mz)';
% xzdif_vz_2=reshape(w(:,5),mx,mz)';
% xzdif_bx_2=reshape(w(:,6),mx,mz)';
% xzdif_by_2=reshape(w(:,7),mx,mz)';
% xzdif_bz_2=reshape(w(:,8),mx,mz)';
% xzdif_p0_2=reshape(w(:,9),mx,mz)';
% w=load(['xdif' str_nst '003' ]);
% xzdif_rho_3=reshape(w(:,1),mx,mz)';
% xzdif_p_3=reshape(w(:,2),mx,mz)';
% xzdif_vx_3=reshape(w(:,3),mx,mz)';
% xzdif_vy_3=reshape(w(:,4),mx,mz)';
% xzdif_vz_3=reshape(w(:,5),mx,mz)';
% xzdif_bx_3=reshape(w(:,6),mx,mz)';
% xzdif_by_3=reshape(w(:,7),mx,mz)';
% xzdif_bz_3=reshape(w(:,8),mx,mz)';
% xzdif_p0_3=reshape(w(:,9),mx,mz)';
% w=load(['xdif' str_nst '004' ]);
% xzdif_rho_4=reshape(w(:,1),mx,mz)';
% xzdif_p_4=reshape(w(:,2),mx,mz)';
% xzdif_vx_4=reshape(w(:,3),mx,mz)';
% xzdif_vy_4=reshape(w(:,4),mx,mz)';
% xzdif_vz_4=reshape(w(:,5),mx,mz)';
% xzdif_bx_4=reshape(w(:,6),mx,mz)';
% xzdif_by_4=reshape(w(:,7),mx,mz)';
% xzdif_bz_4=reshape(w(:,8),mx,mz)';
% xzdif_p0_4=reshape(w(:,9),mx,mz)';
% 
% %Ef
% w=load(['Efd' str_nst '001' ]);
% xzex_1=reshape(w(:,1),mx,mz)';
% xzexdx_1=reshape(w(:,2),mx,mz)';
% xzexdz_1=reshape(w(:,3),mx,mz)';
% xzey_1=reshape(w(:,4),mx,mz)';
% xzeydx_1=reshape(w(:,5),mx,mz)';
% xzeydz_1=reshape(w(:,6),mx,mz)';
% xzez_1=reshape(w(:,7),mx,mz)';
% xzezdx_1=reshape(w(:,8),mx,mz)';
% xzezdz_1=reshape(w(:,9),mx,mz)';
% w=load(['Efd' str_nst '002' ]);
% xzex_2=reshape(w(:,1),mx,mz)';
% xzexdx_2=reshape(w(:,2),mx,mz)';
% xzexdz_2=reshape(w(:,3),mx,mz)';
% xzey_2=reshape(w(:,4),mx,mz)';
% xzeydx_2=reshape(w(:,5),mx,mz)';
% xzeydz_2=reshape(w(:,6),mx,mz)';
% xzez_2=reshape(w(:,7),mx,mz)';
% xzezdx_2=reshape(w(:,8),mx,mz)';
% xzezdz_2=reshape(w(:,9),mx,mz)';
% w=load(['Efd' str_nst '003' ]);
% xzex_3=reshape(w(:,1),mx,mz)';
% xzexdx_3=reshape(w(:,2),mx,mz)';
% xzexdz_3=reshape(w(:,3),mx,mz)';
% xzey_3=reshape(w(:,4),mx,mz)';
% xzeydx_3=reshape(w(:,5),mx,mz)';
% xzeydz_3=reshape(w(:,6),mx,mz)';
% xzez_3=reshape(w(:,7),mx,mz)';
% xzezdx_3=reshape(w(:,8),mx,mz)';
% xzezdz_3=reshape(w(:,9),mx,mz)';
% w=load(['Efd' str_nst '004' ]);
% xzex_4=reshape(w(:,1),mx,mz)';
% xzexdx_4=reshape(w(:,2),mx,mz)';
% xzexdz_4=reshape(w(:,3),mx,mz)';
% xzey_4=reshape(w(:,4),mx,mz)';
% xzeydx_4=reshape(w(:,5),mx,mz)';
% xzeydz_4=reshape(w(:,6),mx,mz)';
% xzez_4=reshape(w(:,7),mx,mz)';
% xzezdx_4=reshape(w(:,8),mx,mz)';
% xzezdz_4=reshape(w(:,9),mx,mz)';


%end

%  if(nst>0) 
%      dvb=load(['dvb' str_nst ]);
%  end   
% %x1ps_zm   
% w=load(['x1ps_zm' str_nst '001' ]);
% pszm_rho_1=reshape(w(:,1),ms,mx)';
% pszm_p_1=reshape(w(:,2),ms,mx)';
% pszm_vx_1=reshape(w(:,3),ms,mx)';
% pszm_vy_1=reshape(w(:,4),ms,mx)';
% pszm_vz_1=reshape(w(:,5),ms,mx)';
% pszm_bx_1=reshape(w(:,6),ms,mx)';
% pszm_by_1=reshape(w(:,7),ms,mx)';
% pszm_bz_1=reshape(w(:,8),ms,mx)';
% 
% w=load(['x1ps_zm' str_nst '002' ]);
% pszm_rho_2=reshape(w(:,1),ms,mx)';
% pszm_p_2=reshape(w(:,2),ms,mx)';
% pszm_vx_2=reshape(w(:,3),ms,mx)';
% pszm_vy_2=reshape(w(:,4),ms,mx)';
% pszm_vz_2=reshape(w(:,5),ms,mx)';
% pszm_bx_2=reshape(w(:,6),ms,mx)';
% pszm_by_2=reshape(w(:,7),ms,mx)';
% pszm_bz_2=reshape(w(:,8),ms,mx)';
% 
% w=load(['x1ps_zm' str_nst '003' ]);
% pszm_rho_3=reshape(w(:,1),ms,mx)';
% pszm_p_3=reshape(w(:,2),ms,mx)';
% pszm_vx_3=reshape(w(:,3),ms,mx)';
% pszm_vy_3=reshape(w(:,4),ms,mx)';
% pszm_vz_3=reshape(w(:,5),ms,mx)';
% pszm_bx_3=reshape(w(:,6),ms,mx)';
% pszm_by_3=reshape(w(:,7),ms,mx)';
% pszm_bz_3=reshape(w(:,8),ms,mx)';
% 
% w=load(['x1ps_zm' str_nst '004' ]);
% pszm_rho_4=reshape(w(:,1),ms,mx)';
% pszm_p_4=reshape(w(:,2),ms,mx)';
% pszm_vx_4=reshape(w(:,3),ms,mx)';
% pszm_vy_4=reshape(w(:,4),ms,mx)';
% pszm_vz_4=reshape(w(:,5),ms,mx)';
% pszm_bx_4=reshape(w(:,6),ms,mx)';
% pszm_by_4=reshape(w(:,7),ms,mx)';
% pszm_bz_4=reshape(w(:,8),ms,mx)';
% 
% %x1ps_zp
% w=load(['x1ps_zp' str_nst '001' ]);
% pszp_rho_1=reshape(w(:,1),ms,mx)';
% pszp_p_1=reshape(w(:,2),ms,mx)';
% pszp_vx_1=reshape(w(:,3),ms,mx)';
% pszp_vy_1=reshape(w(:,4),ms,mx)';
% pszp_vz_1=reshape(w(:,5),ms,mx)';
% pszp_bx_1=reshape(w(:,6),ms,mx)';
% pszp_by_1=reshape(w(:,7),ms,mx)';
% pszp_bz_1=reshape(w(:,8),ms,mx)';
% 
% w=load(['x1ps_zp' str_nst '002' ]);
% pszp_rho_2=reshape(w(:,1),ms,mx)';
% pszp_p_2=reshape(w(:,2),ms,mx)';
% pszp_vx_2=reshape(w(:,3),ms,mx)';
% pszp_vy_2=reshape(w(:,4),ms,mx)';
% pszp_vz_2=reshape(w(:,5),ms,mx)';
% pszp_bx_2=reshape(w(:,6),ms,mx)';
% pszp_by_2=reshape(w(:,7),ms,mx)';
% pszp_bz_2=reshape(w(:,8),ms,mx)';
% 
% w=load(['x1ps_zp' str_nst '003' ]);
% pszp_rho_3=reshape(w(:,1),ms,mx)';
% pszp_p_3=reshape(w(:,2),ms,mx)';
% pszp_vx_3=reshape(w(:,3),ms,mx)';
% pszp_vy_3=reshape(w(:,4),ms,mx)';
% pszp_vz_3=reshape(w(:,5),ms,mx)';
% pszp_bx_3=reshape(w(:,6),ms,mx)';
% pszp_by_3=reshape(w(:,7),ms,mx)';
% pszp_bz_3=reshape(w(:,8),ms,mx)';
% 
% w=load(['x1ps_zp' str_nst '004' ]);
% pszp_rho_4=reshape(w(:,1),ms,mx)';
% pszp_p_4=reshape(w(:,2),ms,mx)';
% pszp_vx_4=reshape(w(:,3),ms,mx)';
% pszp_vy_4=reshape(w(:,4),ms,mx)';
% pszp_vz_4=reshape(w(:,5),ms,mx)';
% pszp_bx_4=reshape(w(:,6),ms,mx)';
% pszp_by_4=reshape(w(:,7),ms,mx)';
% pszp_bz_4=reshape(w(:,8),ms,mx)';
% 
% %x1ps_xm
% w=load(['x1ps_xm' str_nst '001' ]);
% psxm_rho_1=reshape(w(:,1),ms,mz)';
% psxm_p_1=reshape(w(:,2),ms,mz)';
% psxm_vx_1=reshape(w(:,3),ms,mz)';
% psxm_vy_1=reshape(w(:,4),ms,mz)';
% psxm_vz_1=reshape(w(:,5),ms,mz)';
% psxm_bx_1=reshape(w(:,6),ms,mz)';
% psxm_by_1=reshape(w(:,7),ms,mz)';
% psxm_bz_1=reshape(w(:,8),ms,mz)';
% 
% w=load(['x1ps_xm' str_nst '002' ]);
% psxm_rho_2=reshape(w(:,1),ms,mz)';
% psxm_p_2=reshape(w(:,2),ms,mz)';
% psxm_vx_2=reshape(w(:,3),ms,mz)';
% psxm_vy_2=reshape(w(:,4),ms,mz)';
% psxm_vz_2=reshape(w(:,5),ms,mz)';
% psxm_bx_2=reshape(w(:,6),ms,mz)';
% psxm_by_2=reshape(w(:,7),ms,mz)';
% psxm_bz_2=reshape(w(:,8),ms,mz)';
% 
% w=load(['x1ps_xm' str_nst '003' ]);
% psxm_rho_3=reshape(w(:,1),ms,mz)';
% psxm_p_3=reshape(w(:,2),ms,mz)';
% psxm_vx_3=reshape(w(:,3),ms,mz)';
% psxm_vy_3=reshape(w(:,4),ms,mz)';
% psxm_vz_3=reshape(w(:,5),ms,mz)';
% psxm_bx_3=reshape(w(:,6),ms,mz)';
% psxm_by_3=reshape(w(:,7),ms,mz)';
% psxm_bz_3=reshape(w(:,8),ms,mz)';
% 
% w=load(['x1ps_xm' str_nst '004' ]);
% psxm_rho_4=reshape(w(:,1),ms,mz)';
% psxm_p_4=reshape(w(:,2),ms,mz)';
% psxm_vx_4=reshape(w(:,3),ms,mz)';
% psxm_vy_4=reshape(w(:,4),ms,mz)';
% psxm_vz_4=reshape(w(:,5),ms,mz)';
% psxm_bx_4=reshape(w(:,6),ms,mz)';
% psxm_by_4=reshape(w(:,7),ms,mz)';
% psxm_bz_4=reshape(w(:,8),ms,mz)';
% 
% %x1ps_xp
% w=load(['x1ps_xp' str_nst '001' ]);
% psxp_rho_1=reshape(w(:,1),ms,mz)';
% psxp_p_1=reshape(w(:,2),ms,mz)';
% psxp_vx_1=reshape(w(:,3),ms,mz)';
% psxp_vy_1=reshape(w(:,4),ms,mz)';
% psxp_vz_1=reshape(w(:,5),ms,mz)';
% psxp_bx_1=reshape(w(:,6),ms,mz)';
% psxp_by_1=reshape(w(:,7),ms,mz)';
% psxp_bz_1=reshape(w(:,8),ms,mz)';
% 
% w=load(['x1ps_xp' str_nst '002' ]);
% psxp_rho_2=reshape(w(:,1),ms,mz)';
% psxp_p_2=reshape(w(:,2),ms,mz)';
% psxp_vx_2=reshape(w(:,3),ms,mz)';
% psxp_vy_2=reshape(w(:,4),ms,mz)';
% psxp_vz_2=reshape(w(:,5),ms,mz)';
% psxp_bx_2=reshape(w(:,6),ms,mz)';
% psxp_by_2=reshape(w(:,7),ms,mz)';
% psxp_bz_2=reshape(w(:,8),ms,mz)';
% 
% w=load(['x1ps_xp' str_nst '003' ]);
% psxp_rho_3=reshape(w(:,1),ms,mz)';
% psxp_p_3=reshape(w(:,2),ms,mz)';
% psxp_vx_3=reshape(w(:,3),ms,mz)';
% psxp_vy_3=reshape(w(:,4),ms,mz)';
% psxp_vz_3=reshape(w(:,5),ms,mz)';
% psxp_bx_3=reshape(w(:,6),ms,mz)';
% psxp_by_3=reshape(w(:,7),ms,mz)';
% psxp_bz_3=reshape(w(:,8),ms,mz)';
% 
% w=load(['x1ps_xp' str_nst '004' ]);
% psxp_rho_4=reshape(w(:,1),ms,mz)';
% psxp_p_4=reshape(w(:,2),ms,mz)';
% psxp_vx_4=reshape(w(:,3),ms,mz)';
% psxp_vy_4=reshape(w(:,4),ms,mz)';
% psxp_vz_4=reshape(w(:,5),ms,mz)';
% psxp_bx_4=reshape(w(:,6),ms,mz)';
% psxp_by_4=reshape(w(:,7),ms,mz)';
% psxp_bz_4=reshape(w(:,8),ms,mz)';


dat_mat=[ 'png1\dat' str_nst '.mat' ];
save(dat_mat);
% plot(xz_rho(128,:));fname=['rho1' str_nst ];title([ fname '   t=' str_time ]);print('-dpng', ['png1\' fname]);
% plot(xz_p(128,:));fname=['p1' str_nst ];title([ fname '   t=' str_time ]);print('-dpng', ['png1\' fname]);
% plot(xz_vx(128,:));fname=['vx1' str_nst ];title([ fname '   t=' str_time ]);print('-dpng', ['png1\' fname]);
% plot(xz_vz(128,:));fname=['vz1' str_nst ];title([ fname '   t=' str_time ]);print('-dpng', ['png1\' fname]);
% plot(xz_vy(128,:));fname=['vy1' str_nst ];title([ fname '   t=' str_time ]);print('-dpng', ['png1\' fname]);
% plot(xz_by(128,:));fname=['by1' str_nst ];title([ fname '   t=' str_time ]);print('-dpng', ['png1\' fname]);
% plot(xz_bx(128,:));fname=['bx1' str_nst ];title([ fname '   t=' str_time ]);print('-dpng', ['png1\' fname]);
% plot(xz_bz(128,:));fname=['bz1' str_nst ];title([ fname '   t=' str_time ]);print('-dpng', ['png1\' fname]);
 [C,h]=contourf(xz_xx,xz_zz,xz_rho,50);figure(gcf); set(h, 'LineStyle', 'none');fname=['rho' str_nst ];title([ fname '   t=' str_time ]);xlabel('R');ylabel('Z');colorbar;axis('equal',[xmin xmax zmin zmax ]);hold on;plot(xb,zb,'--');hold off;print('-dpng', ['png1\' fname]);
 [C,h]=contourf(xz_xx,xz_zz,xz_p,50);figure(gcf); set(h, 'LineStyle', 'none');fname=['p' str_nst ];title([ fname '   t=' str_time ]);xlabel('R');ylabel('Z');colorbar;axis('equal',[xmin xmax zmin zmax ]);hold on;plot(xb,zb,'--');hold off;print('-dpng', ['png1\' fname]);
 [C,h]=contourf(xz_xx,xz_zz,xz_vy,50);figure(gcf); set(h, 'LineStyle', 'none');fname=['vy' str_nst ];title([ fname '   t=' str_time ]);xlabel('R');ylabel('Z');colorbar;axis('equal',[xmin xmax zmin zmax ]);hold on;plot(xb,zb,'--');hold off;print('-dpng', ['png1\' fname]);
 [C,h]=contourf(xz_xx,xz_zz,xz_vx,50);figure(gcf); set(h, 'LineStyle', 'none');fname=['vx' str_nst ];title([ fname '   t=' str_time ]);xlabel('R');ylabel('Z');colorbar;axis('equal',[xmin xmax zmin zmax ]);hold on;plot(xb,zb,'--');hold off;print('-dpng', ['png1\' fname]);
 [C,h]=contourf(xz_xx,xz_zz,xz_vz,50);figure(gcf); set(h, 'LineStyle', 'none');fname=['vz' str_nst ];title([ fname '   t=' str_time ]);xlabel('R');ylabel('Z');colorbar;axis('equal',[xmin xmax zmin zmax ]);hold on;plot(xb,zb,'--');hold off;print('-dpng', ['png1\' fname]);
 [C,h]=contourf(xz_xx,xz_zz,xz_vr,50);figure(gcf); set(h, 'LineStyle', 'none');fname=['vr' str_nst ];title([ fname '   t=' str_time ]);xlabel('R');ylabel('Z');colorbar;axis('equal',[xmin xmax zmin zmax ]);hold on;plot(xb,zb,'--');hold off;print('-dpng', ['png1\' fname]);
[C,h]=contourf(xz_xx,xz_zz,xz_vp,50);figure(gcf); set(h, 'LineStyle', 'none');fname=['vp' str_nst ];title([ fname '   t=' str_time ]);xlabel('R');ylabel('Z');colorbar;axis('equal',[xmin xmax zmin zmax ]);hold on;plot(xb,zb,'--');hold off;print('-dpng', ['png1\' fname]);
 [C,h]=contourf(xz_xx,xz_zz,xz_by,50);figure(gcf); set(h, 'LineStyle', 'none');fname=['by' str_nst ];title([ fname '   t=' str_time ]);xlabel('R');ylabel('Z');colorbar;axis('equal',[xmin xmax zmin zmax ]);hold on;plot(xb,zb,'--');hold off;print('-dpng', ['png1\' fname]);
 [C,h]=contourf(xz_xx,xz_zz,xz_bx,50);figure(gcf); set(h, 'LineStyle', 'none');fname=['bx' str_nst ];title([ fname '   t=' str_time ]);xlabel('R');ylabel('Z');colorbar;axis('equal',[xmin xmax zmin zmax ]);hold on;plot(xb,zb,'--');hold off;print('-dpng', ['png1\' fname]);
 [C,h]=contourf(xz_xx,xz_zz,xz_bz,50);figure(gcf); set(h, 'LineStyle', 'none');fname=['bz' str_nst ];title([ fname '   t=' str_time ]);xlabel('R');ylabel('Z');colorbar;axis('equal',[xmin xmax zmin zmax ]);hold on;plot(xb,zb,'--');hold off;print('-dpng', ['png1\' fname]);
 [C,h]=contourf(xz_xx,xz_zz,xz_br,50);figure(gcf); set(h, 'LineStyle', 'none');fname=['br' str_nst ];title([ fname '   t=' str_time ]);xlabel('R');ylabel('Z');colorbar;axis('equal',[xmin xmax zmin zmax ]);hold on;plot(xb,zb,'--');hold off;print('-dpng', ['png1\' fname]);
 [C,h]=contourf(xz_xx,xz_zz,xz_bp,50);figure(gcf); set(h, 'LineStyle', 'none');fname=['bp' str_nst ];title([ fname '   t=' str_time ]);xlabel('R');ylabel('Z');colorbar;axis('equal',[xmin xmax zmin zmax ]);hold on;plot(xb,zb,'--');hold off;print('-dpng', ['png1\' fname]);
%  [C,h]=contourf(xz_xx,xz_zz,xz_cx,50);figure(gcf); set(h, 'LineStyle', 'none');fname=['cx' str_nst ];title(fname);colorbar;print('-dpng', ['png1\' fname]);
%  [C,h]=contourf(xz_xx,xz_zz,xz_cy,50);figure(gcf); set(h, 'LineStyle', 'none');fname=['cy' str_nst ];title(fname);colorbar;print('-dpng', ['png1\' fname]);
% [C,h]=contourf(xz_xx,xz_zz,xz_cz,50);figure(gcf); set(h, 'LineStyle', 'none');fname=['cz' str_nst ];title(fname);colorbar;print('-dpng', ['png1\' fname]);
% [C,h]=contourf(xz_xx,xz_zz,xz_ex,50);figure(gcf); set(h, 'LineStyle', 'none');fname=['ex' str_nst ];title(fname);colorbar;print('-dpng', ['png1\' fname]);
% [C,h]=contourf(xz_xx,xz_zz,xz_ey,50);figure(gcf); set(h, 'LineStyle', 'none');fname=['ey' str_nst ];title(fname);colorbar;print('-dpng', ['png1\' fname]);
% [C,h]=contourf(xz_xx,xz_zz,xz_ez,50);figure(gcf); set(h, 'LineStyle', 'none');fname=['ez' str_nst ];title(fname);colorbar;print('-dpng', ['png1\' fname]);
% 
  [C,h]=contourf(xz_xx,xz_zz,xz_vs,50);figure(gcf); set(h, 'LineStyle', 'none');fname=['vs' str_nst ];title([ fname '   t=' str_time ]);xlabel('R');ylabel('Z');colorbar;axis('equal',[xmin xmax zmin zmax ]);hold on;plot(xb,zb,'--');hold off;print('-dpng', ['png1\' fname]);

end
% 
% for jr=1:mr
%     for jt=1:mt
%         xl(jr,jt)=xzero+ar(jr)*cos(at(jt));
%         zl(jr,jt)=ar(jr)*sin(at(jt));
%     end
% end




