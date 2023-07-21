%% Computable performance guarantees for compressed sensing matrices 
% Aug. 27, 2017
% By Myung (Michael) Cho
% Email: michael.myung.cho@gmail.com
%
% Code for using Sandwiching Algorithm (SWA)
%-------------------------------------
clear
clc
addpath('./03_SWA_StandAlone');

% Types:
% 1    strF='Fourier';
% 3    str1='+1/-1';
% 4    strG='Gauss';
% 9    strV='DeVore';

randn('state',100);
n=40;nsamp=1;krange=1:5;rhorange=[0.5];typelist=[4*ones(1,nsamp)];

%% save
count=0;    
% save data to xls file
clk=clock;
datestmp=strcat([date,'-',num2str(clk(4)),'h',num2str(clk(5))]);
filename=strcat('T0_NSP_SWA-',datestmp,'-Gauss');
filename=strcat(filename,'.xlsx');
col_name={'type','rho','kv','alphaSWA','Tswa';};
xlswrite(filename,col_name);


T1swa=0; alphaSWA=0; stepSWA=0;
T1tsa=0; alphaTSA_LB=0; alphaTSA_UB=0;
%% Iterate over many samples
for rho=rhorange
    for tp=typelist    
    %% Generate matrix
    m = round(rho*n);
    [rss,pnt]=AAGenerateSensingMatrix(tp,n,m);
    A=rss.A;n=rss.n;m=rss.m;
    P=null(A);nk=size(P,2);
    InList1=[];InList2=[];InList3=[];
        for kv=krange
        %% SWA (based on pick-1)
        T0swa=tic;      [alphaSWA, stepSWA] = BoundSW(A,kv,1); T1swa=toc(T0swa);
        %% TSA (based on pick-1)
        cd ./04_TSA_StandAlone
            T0tsa=tic;      [alphaTSA_LB,alphaTSA_UB]=BoundTSA(A, kv,1);   T1tsa=toc(T0tsa);
        cd ..
        %% Display results
        resbuf=[kv, m, n, T1tsa, T1swa , alphaTSA_LB, alphaTSA_UB ,alphaSWA ];
        disp('**************************************************************************************************');
        count=count+1;
        
        %% save data to xls file
        xlRange = sprintf('A%d',count);
        xlsInsert = resbuf;
        xlswrite(filename,xlsInsert,'Sheet1',xlRange);        
        end
    end
end




