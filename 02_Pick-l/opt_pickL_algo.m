function [ opt_UB_ak ] = opt_pickL_algo( n, kv, InList )
%% Optimized Pick-L-element algorithm
% To calculate upper bound on alpha_k from alpha_{l,L}'s via optimized
% coefficients (Section 3 of NSP paper)
%
% Sept. 12, 2017
% By Myuug (Michael) Cho
% Email: michael.myung.cho@gmail.com
%--------------------------

%% Variable assginment
f=[];
f=InList(:,end);
L=size(InList,2)-1;
[nChooseL,~]=size(InList);

%% Constraint setting in a matrix A
% Initialization
A=[];
y=[];

% 2nd constraint
A=[A;ones(1,nChooseL)];
y=[y;kv/L];

% 3rd constraint
iConst=2;
for b=1:L
    IdxB=nchoosek([1:n],b);
   for ii=1:size(IdxB,1)
       indSubset=find(sum(ismember(InList(:,1:end-1),IdxB(ii,:)),2)==b);
       A(iConst,indSubset)=1;
       y(iConst,1)=nchoosek(kv-b,L-b)/nchoosek(kv-1,L-1);
       iConst=iConst+1;
   end
end
z=zeros(nChooseL,1);
o=ones(nChooseL,1);
%[x,fval,exitflag,output,lambda]=linprog(-f,A,y,[],[],0,[],[],options)
%addpath('C:\Program Files\MATLAB\R2013b\toolbox\optim\optim');
[x,fval,exitflag,output,lambda]=linprog(-f,A,y,[],[],z,o);

opt_UB_ak=-fval;
%save('temp_opt.mat')
end
