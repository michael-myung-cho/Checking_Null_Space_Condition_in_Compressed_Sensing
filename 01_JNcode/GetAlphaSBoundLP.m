function [Ralpha,Vs,Valpha,Y]=GetAlphaSBoundLP(A,sprs, half, prn)
%Ralpha: the value of alpha found in sprs verification
%Vs: the maximum sparsity that is verified using the solution found
%Valpha: the value of alpha corresponding to Vs verification
Vs=0;
Valpha=0;
[m,n]=size(A);
cls=m*n+n+1;
ctot=m*n+n*(n+1)+1;
%
prb.a=[];
prb.blc=-inf*ones(n*(2*n+1),1);
prb.buc=zeros(n*(2*n+1),1);
prb.blx=[-inf*ones(m*n,1);zeros(n*(n+1),1);0];
prb.bux=[inf*ones(m*n,1);inf*ones(n*(n+1),1);inf];
base=0;
NNZ=n*(n*(2*(m+2))+n+2);
Nii=zeros(1,NNZ);
Njj=zeros(1,NNZ);
Nvv=zeros(1,NNZ);
rc=1;
imem=1;
for i=1:n,
    a=A(:,i)';
% block=[];
    aft=ctot-base-cls;
    for irw=1:n,
        lsp=m*(irw-1);
        rsp=m*(n-irw)+irw-1;
%
        o=lsp+1;
        oo=o+m;
        ooo=oo+rsp+base;
        oooo=ooo+n-irw+1;
        Nii(imem:imem+m+1)=rc*ones(1,m+2);
        Njj(imem:imem+m+1)=[o:oo-1, ooo, oooo];
        Nvv(imem:imem+m+1)=[a, -1, -1];
        imem=imem+m+2;
        rc=rc+1;
        Nii(imem:imem+m+1)=rc*ones(1,m+2);
        Njj(imem:imem+m+1)=[o:oo-1, ooo, oooo];        
        Nvv(imem:imem+m+1)=[-a, -1, -1];
        imem=imem+m+2;
        rc=rc+1;        
    end;
%
        o=m*n+base+1;
        oo=o+n;
        ooo=ooo+aft-1+2;
        Nii(imem:imem+n+1)=rc*ones(1,n+2);
        Njj(imem:imem+n+1)=[o:oo-1, oo, ooo];
        Nvv(imem:imem+n+1)=[ones(1,n), sprs, -1];
        imem=imem+n+2;
        rc=rc+1;    
        prb.buc((i-1)*(2*n+1)+2*i-1)=1;
        prb.buc((i-1)*(2*n+1)+2*i)=-1;
        prb.blx(m*n+(i-1)*(n+1)+n+1)=-inf;
        
    base=base+n+1;
end
prb.c=zeros(ctot,1);
prb.c(ctot)=1;
prb.a=sparse(Nii,Njj,Nvv,rc-1,ctot);
clear Nii;
clear Njj;
clear Nvv;
[pm,pn]=size(prb.a);
%
if prn, disp(['Using MOSEK to solve a ',num2str(pm),'x',num2str(pn),' LP']); end
if prn>1,
    [rr,res]=mosekopt('minimize',prb);
else
    [rr,res]=mosekopt('minimize echo(0)',prb);
end
if rr>0, 
    if isfield(res,'rcode')&&(res.rcode>0), error(res.rmsg);
    else error(['mosekopt rr: ',num2str(rr)]);
    end;
end
% Get the objective value from LP and store it as Ralpha
Ralpha=prb.c'*res.sol.itr.xx;
Y=zeros(m,n);
base=0;
for i=1:n,
    Y(:,i)=res.sol.itr.xx(base+1:base+m);
    base=base+m;
end;
Tmp=sort(abs(eye(n)-Y'*A));
for i=1:size(Tmp,1);
    tmp=ones(1,i);
    gamma=max(tmp*Tmp(n-i+1:n,:));
    if gamma<half
        Vs=max(Vs,i);
        Valpha=max(Valpha,gamma);
    else
        break;
    end;
end;
end