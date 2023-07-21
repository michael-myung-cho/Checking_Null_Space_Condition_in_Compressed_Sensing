%% Computable performance guarantees for compressed sensing matrices
% Aug. 27, 2017
% By Myung (Michael) Cho
% Email: michael.myung.cho@gmail.com
%
% Code for generating Figure 5 (Gaussian matrix) in the paper (TSA complexity)
%-------------------------------------

addpath('./00_Rawdata/01_Low_dimension/01_Gaussian/01_rho_0.5');

%% save
count=0;
% save data to xls file
clk=clock;
datestmp=strcat([date,'-',num2str(clk(4)),'h',num2str(clk(5))]);
filename=strcat('F5_TSA_complexity_Gaussian-',datestmp);
filename=strcat(filename,'.xlsx');
col_name={'rho','n','kv','alphaTSA_P1','alphaTSA_P2','alphaTSA_P3','Ttsa_p1','Ttsa_p2','Ttsa_p3';};
xlswrite(filename,col_name);

tp=4; % Gaussian
%% Generate matrix

load('Fig5_25x50_Gauss.mat');
%% Comparison starts
% initialization
T1pick3=0; T1tsaP1=0; T1tsaP2=0; T1tsaP3=0;
alphaTSA_P1=0; alphaTSA_P2=0; alphaTSA_P3=0; 
for kv=1:5
    %% TSA pick1
    cd ./04_1_TSA_StandAlone
    T0tsaP1=tic;      alphaTSA_P1=BoundTSA(A,kv,1);       T1tsaP1=toc(T0tsaP1);
    cd ..
    %% TSA pick2
    if kv>=2
        cd ./04_1_TSA_StandAlone
        T0tsaP2=tic;      alphaTSA_P2=BoundTSA(A,kv,2);       T1tsaP2=toc(T0tsaP2);
        cd ..
    end
    if kv>=3
        cd ./04_1_TSA_StandAlone
        T0tsaP3=tic;      alphaTSA_P3=BoundTSA(A,kv,3);       T1tsaP3=toc(T0tsaP3);
        cd ..
    end
    
    % Display results
    resbuf=[rho,n,kv,alphaTSA_P1,alphaTSA_P2,alphaTSA_P3,T1tsaP1,T1tsaP2,T1tsaP3];
    disp('**************************************************************************************************');
    disp(num2str(resbuf));
    disp('**************************************************************************************************');
    count=count+1;
    %% save data to xls file
    xlRange = sprintf('A%d',count+1);
    xlsInsert = resbuf;
    xlswrite(filename,xlsInsert,'Sheet1',xlRange);
end   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
