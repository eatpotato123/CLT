clear;
nr=256;
n=200;
np=83;
aa=1;
dir='png_line';
mkdir(dir);
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
nstart=0;
nstartime=0;
ndst=1;
phi=[0,45,90,180];
for nst=nstart:ndst:100
    str_nst=num2str(nst,'%3.3d');
    str_time=num2str(time((nst-nstartime)/ndst+1,2));
    pointr1=load(['outputdatar1' str_nst])*aa;
    pointz1=load(['outputdataz1' str_nst])*aa;
    save([dir '\line' str_nst '.mat']);
     k=1
    str_phi=num2str(phi(k),'%3.3d');   
    % 将读入的数据按照进行分组
    pointr=reshape(pointr1(:,k),np,n);
    pointz=reshape(pointz1(:,k),np,n);

    pointr=pointr';
    pointz=pointz';

    for i=1:2:np
        scatter(pointr(:,i),pointz(:,i),25,'.');
        hold on
    end
    axis('equal',[xmin xmax zmin zmax ])
    fname=['l' str_nst str_phi ];
    title(['\phi=' str_phi '   t=' str_time  ]);
    xlabel('R');
    ylabel('Z');
    print('-dpng', [dir '\' fname]);
    hold off
    
end
mkdir png_line000;
copyfile( ['./' dir '/l*000.png' ],'./png_line000');