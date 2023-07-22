function [snew,F,tused]=GetAlphaSBound(A,s,sold,ctrl)
%GetAlphaSBound computes a lower bound for the supported by A sparsity
%
%   [s_new,Y]=GetAlphaSBound(A,target_s,s_old,control);
%
%   Inputs:
%   A; (double array) matrix to be tested
%   target_s: (scalar integer) the sparsity to test 
%   s_old: the structure of results. Should contain the fields:
%       s_mi: sparsity bound obtained using the mutual incoherence
%       s_bst: the best sparsity bound obtained
%       s_alpha: structure containing sparsity bounds 
%       Fields of s_alpha:
%           al1: an "exact" upper bound on alpha_1 (as in the J.-N.'s paper),
%           al1s: an upper bound on alpha_s obtained using the matrix eye(n)-Y'*A
%           computed in course of alpha_1 computation
%           s_alpha_1: corresponding sparsity bound
%           als: alpha_s bound
%           s_alpha_s: corresponding sparsity bound
%       control: structure with parameters of the procedure
%       Fields of control:
%           solver: either 'LP' (or 'lp', etc.) or 'MP' (or 'mp', or
%           'Mirror',etc.) -- the solver used. LP reformulation is soved
%           using MOSEK; saddle-point reformulation is solved using
%           Mirror-Prox algorithm
%       
%   Outputs: 
%   s_new: new structure (same as s_old) containing the updated sparsity bounds
%   Y: matrix certificate (matrix of the same dimension as A) as in the J.-N.'s paper
%

tstart=clock;
ni=nargin;

% if ~ni, error('At least one input argument is required!'); end
[m,n]=size(A);
printlevel=1;
half=0.499999;
cond=1e-8;
solver=0; % Using Mosek by default 

if ni>=4, 
    if isfield(ctrl,'printlevel'), printlevel=ctrl.printlevel; end
    if isfield(ctrl,'half'), half=ctrl.half; end
    if isfield(ctrl,'cond'), cond=ctrl.cond; end
    if isfield(ctrl,'solver'), 
        if ~isempty(strfind(lower(ctrl.solver),'m')), solver=1; end% MP-solver
    end
end
if printlevel, fprintf('%1dx%1d sensing matrix; s_targ=%1d \n',m,n,s); end
drawnow; 
%
if ni<3||isempty(sold)
    sold.s_alpha.al1=inf;
    sold.s_alpha.al1s=inf;
    sold.s_alpha.als=inf;
    sold.s_alpha.s_alpha_1=0;
    sold.s_alpha.s_alpha_s=0;
    sold.s_mi=0; 
end
if ni<2, s=m; end
snew=sold;
% compute mutual incoherehce
if isempty(snew.s_mi)|| snew.s_mi==inf 
    gmx=0;
    for i=1:n,
    y=A(:,i);
    y=y/norm(y)^2;
    z=A'*y;
    z(i)=z(i)-1;
    g=max(abs(z));
    g=g/(1+g);
    gmx=max(gmx,g);
    end;
    snew.s_mi=max(min(floor(half/gmx),m),snew.s_mi);
    if printlevel,
        fprintf('s.s_mi=%1d',snew.s_mi);
    	drawnow
    end;
end
% "normalize" sensing matrix
[U,D,V]=svd(A);
dm=max(diag(D));
I=find(diag(D)>cond*dm);
B=V(:,I)';
if solver,
    [~,Vs,Valpha,Y]=GetAlphaSBoundMP(B,s,half,printlevel);
else
    [~,Vs,Valpha,Y]=GetAlphaSBoundLP(A,s,half,printlevel);
end

snew.s_alpha.als=min(snew.s_alpha.als,Valpha);
snew.s_alpha.s_alpha_s=max(snew.s_alpha.s_alpha_s, Vs);
U1=U(:,I); D1=diag(1./diag(D(I,I)));
F=U1*D1*Y;
snew.bst=max([snew.s_alpha.s_alpha_s,snew.s_alpha.s_alpha_1,snew.s_mi]);

tused=etime(clock,tstart);
if printlevel,
    fprintf('Certified sparsity: %d (Alpha_%d=%2.7f); ',Vs,Vs,Valpha);
    fprintf('CPU=%5.1f  \n',tused);
    drawnow
end
end
%
%
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
%
%
%
function [Ralpha,Vs,Valpha,Y]=GetAlphaSBoundMP(B,s, half, prn)
%Ralpha: the value of alpha found in sparsity s verification
%Vs: the maximum sparsity that is verified using the solution found
%Valpha: the value of alpha corresponding to Vs verification
Ralpha=0;
Vs=0;
Valpha=0;
[m,n]=size(B);
nx=m*n+n^2+n;
ny=2*n*n;
lpnc=2*n^2+n;
lpnv=n*(m+n+1);
if prn
    fprintf('Solving %1dx%1d saddle point reformulation of %1dx%1d LP... \n',nx,ny,lpnc,lpnv);
    drawnow;
end
[sol,Y]=MPAlphaS(B,s);
lwb=sol.lb;
upb=sol.ub;
%if prn, fprintf('target s=%d; lb=%8.7f; ub=%8.7f \n',s,lwb,upb);drawnow; end
%
if upb<0,
    Ralpha=upb;
else
    if lwb>0,
        Ralpha=lwb;
    end;
end;
Z=Y'*B;
sY=sort(abs(eye(n)-Z));
for sc=1:n,
    tmp=ones(1,sc)*sY(n-sc+1:n,:);
    if max(tmp)<half,
        Vs=max(Vs,sc);
        Valpha=max(max(tmp),Valpha);
    else
        break;
    end;
end;
end