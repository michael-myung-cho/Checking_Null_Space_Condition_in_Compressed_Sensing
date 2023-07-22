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