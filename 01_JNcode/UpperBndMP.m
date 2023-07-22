function [lbs,bad]=UpperBndMP(A,s,sini,target,printlevel,reached);
n=size(A,2);
m=size(A,1);
[U,D,V]=svd(A);
dm=D(1,1);
mm=m;
for i=1:m,
    if D(i,i)<=1.e-8*dm,
        mm=i-1;
        break;
    end;
end;
B=zeros(mm,n);
for i=1:mm,
    B(i,:)=V(:,i)';
end;
lb=0;
lbs.ss=[sini:1:m]';
lbs.gs=zeros(m-sini+1,1);
bad.s=m+1;
sc=s;
disp(sprintf('Solving series of %1dx%1d saddle point reformulations of %1dx%1d LPs. Be patient...',mm,2*n,mm,3*n+1));
pause(0.01);
irepeat=0;
for itemp=1:24,
    [x,ind]=sort(rand(n,1));
    u=zeros(n,1);
    u(ind(1:s))=sign(randn(s,1));
    lbc=0;
    flag=0;
    while(1==1)
        fff=fopen('B_LP.dat','w');
        fprintf(fff,'%d %d %d \n',mm,n,mm*n);
        for i=1:n,
            fprintf(fff,'%e\n',u(i));
        end;
        for i=1:mm,
            for j=1:n,
                fprintf(fff,'%d %d %24.15e\n',i-1,j-1,B(i,j));
            end;
        end;
        fclose(fff);
        !rm Out_LP.dat;
        !./MirrorProxLP
        irepeat=irepeat+1;
        if mod(irepeat,8)==0,
            disp(sprintf('... %3d saddle point problems solved, ub=%d, more to go...',irepeat,min(bad.s,reached+1)-1));
        end;
        fff=fopen('Out_LP.dat','r');
        if fff<=0,
            disp('Something wrong with LPMirrorProx');
            return;
        end;
        z=zeros(n,1);
        for i=1:n,
            z(i)=fscanf(fff,'%e',1);
        end;
        fclose(fff);
        tmp=B*z;
        ntmp=norm(tmp);
        if ntmp>1.e-4,
            disp(sprintf('residual: %.3e',norm(tmp)));
        end;
        x=z-B'*tmp;
        if sum(abs(x))>1,
            x=x/sum(abs(x));
        end;
        [sx,ind]=sort(abs(x));
        lbn=sum(sx(n-sc+1:n));
        for i=1:m-sini+1,
            ss=lbs.ss(i);
            tmp=sum(sx(n-ss+1:n));
            if (tmp>=0.5)&(ss<bad.s),
                if printlevel>0,
                    disp(sprintf('upper bound set to %1d',ss-1));
                    pause(0.001);
                end;
                if ss-1<=target,
                    disp('target is reached');
                    pause(0.001);
                    flag=1;
                end;
                bad.s=ss;
                bad.ind=ind(n-ss+1:n);
                bad.sgn=sign(x(ind(n-ss+1:n)));
            end;
            lbs.gs(i)=max(lbs.gs(i),tmp);
        end;
        if flag==1,
            break;
        end;
        if printlevel>1,
            disp(sprintf('CPU: %5.1f %5.4f -> %5.4f',tused,lbc,lbn));
            pause(0.001);
        end;
        if lbn>=0.5,
            break;
        end;
        if lbn-lbc<0.01,
            break;
        end;
        lbc=lbn;
        lb=max(lb,lbc);
        u=zeros(n,1);
        u(ind(n-sc+1:n))=sign(x(ind(n-sc+1:n)));
    end;
    if flag==1,
        break;
    end;
    sc=floor(s+rand(1,1)*(bad.s-s));
end;
