%% Computable performance guarantees for compressed sensing matrices
% Aug. 27, 2017
% By Myung (Michael) Cho
% Email: michael.myung.cho@gmail.com
%
% Code for generating Table 4 in the paper (comparision between basic
% pick-l and optimized pick-l)
%-------------------------------------

%addpath('/Users/myucho/Documents/RESEARCH/00_NSP/mosek/8/toolbox/r2014a'); % need to add mosek path
addpath('./00_Rawdata/01_Low_dimension/03_OptPick_VS_Pick');
addpath('./02_Pick-l');

%% matrix load
% 28x40 Gaussian matrix
load('01_28x40_gaus_optpick_vs_basic.mat');
% 40x50 Gaussian matrix
%load('02_40x50_gaus_optpick_vs_basic.mat');

%% Comparison starts
% initialization
InList1=[];InList2=[];InList3=[];
T1sdp=0; T1lp=0; T1pick1=0; T1pick2=0; T1pick3=0; T1pick3opt=0; T1pick2opt=0; T1pick1opt=0;
alphaPICK1=0; alphaPICK2=0; alphaPICK3=0; alphaPICK1OPT=0; alphaPICK2OPT=0; alphaPICK3OPT=0;alphaTSA_LB=0;alphaTSA_UB=0
resbuf=[];
for kv=1:10
    %% Pick-1
    T0pick1=tic;
    if kv==1
        InList1 = pickL_element_algo_Lin(A,1);
        alphaPICK1 = sum(InList1(1:kv,end));
    else
        alphaPICK1 = sum(InList1(1:kv,end));
    end
    T1pick1=toc(T0pick1);
    
    %% opt pick-1
    T0pick1opt=tic;
    alphaPICK1OPT=opt_pickL_algo(n,kv,InList1);
    T1pick1opt=toc(T0pick1opt);
    
    %% Pick-2
    T0pick2=tic;
    if kv==1
        alphaPICK2 = alphaPICK1;
    elseif kv==2
        InList2 = pickL_element_algo_Lin(A,2);
        alphaPICK2 = sum(InList2(1:nchoosek(kv,2),end))/(nchoosek(kv-1,2-1));
    else
        alphaPICK2 = sum(InList2(1:nchoosek(kv,2),end))/(nchoosek(kv-1,2-1));
    end
    T1pick2=toc(T0pick2);
    
    %% opt pick-2
    T0pick2opt=tic;
    if kv>=2
        alphaPICK2OPT=opt_pickL_algo(n,kv,InList2);
    end
    T1pick2opt=toc(T0pick2opt);
    
    %% Pick-3
    T0pick3=tic;
    if kv==1
        alphaPICK3 = alphaPICK1;
    elseif kv==2
        alphaPICK3 = alphaPICK2;
    elseif kv==3
        InList3 = pickL_element_algo_Lin(A,3);
        alphaPICK3 = sum(InList3(1:nchoosek(kv,3),end))/(nchoosek(kv-1,3-1));
    else
        alphaPICK3 = sum(InList3(1:nchoosek(kv,3),end))/(nchoosek(kv-1,3-1));
    end
    T1pick3=toc(T0pick3);
    
    %% opt pick-3
    T0pick3opt=tic;
    if kv>=3
        alphaPICK3OPT=opt_pickL_algo(n,kv,InList3);
    end
    T1pick3opt=toc(T0pick3opt);
    
    %% TSA pick3
    T0tsa=tic;
    if kv>=3
        cd ./04_TSA_StandAlone
        T0tsaP1=tic;      [alphaTSA_LB,alphaTSA_UB]=BoundTSA(A,kv,3);       T1tsaP1=toc(T0tsaP1);
        cd ..
    end
    T1tsa=toc(T0tsa);

    %% Display results
    resbuf=[resbuf;rho,kv,alphaPICK1,alphaPICK1OPT,alphaPICK2,alphaPICK2OPT,alphaPICK3,alphaPICK3OPT,alphaTSA_LB,alphaTSA_UB,T1pick1,T1pick2,T1pick3,T1pick1opt,T1pick2opt,T1pick3opt,T1tsa];
    disp('**************************************************************************************************');
    disp(num2str(resbuf));
    disp('**************************************************************************************************');
end













