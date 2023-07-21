%% Computable performance guarantees for compressed sensing matrices
% Aug. 27, 2017
% By Myung (Michael) Cho
% Email: michael.myung.cho@gmail.com
%
% Code for generating Table 6 (Network Topology) in the paper
%-------------------------------------

addpath('./00_Rawdata/03_Network_topology');

%% save
count=0;
% save data to xls file
clk=clock;
datestmp=strcat([date,'-',num2str(clk(4)),'h',num2str(clk(5))]);
filename=strcat('T6_NSP_TSA_Large_Network-',datestmp);
filename=strcat(filename,'.xlsx');
col_name={'case','kv','alphaTSA_LB_P1','alphaTSA_UB_P1','alphaTSA_LB_P2','alphaTSA_UB_P2','Ttsa_p1','Ttsa_p2';};
xlswrite(filename,col_name);

%% Load matrix
for icase=1:1
    load('03_320x400_Network_model.mat');
    
    %% Calculation starts
    % initialization
    T1tsaP1=0; T1tsaP2=0;
    alphaTSA_LB_P1=0;  alphaTSA_UB_P1=0;
    alphaTSA_LB_P2=0;  alphaTSA_UB_P2=0;
    alphaTSA_LB_P3=0;  alphaTSA_UB_P3=0; 
    for kv=1:6
        %% TSA pick1
        cd ./04_TSA_StandAlone
            T0tsaP1=tic;      [alphaTSA_LB_P1,alphaTSA_UB_P1]=BoundTSA(A,kv,1);       T1tsaP1=toc(T0tsaP1);
        cd ..
        %% TSA pick2
        if kv>=2
            cd ./04_TSA_StandAlone
            T0tsaP2=tic;      [alphaTSA_LB_P2,alphaTSA_UB_P2]=BoundTSA(A,kv,2);       T1tsaP2=toc(T0tsaP2);
            cd ..
        end
        
        %% Display results
        resbuf=[icase,kv,alphaTSA_LB_P1,alphaTSA_UB_P1,alphaTSA_LB_P2,alphaTSA_UB_P2,T1tsaP1,T1tsaP2];
        disp('**************************************************************************************************');
        disp(num2str(resbuf));
        disp('**************************************************************************************************');
        count=count+1;
        %% save data to xls file
        xlRange = sprintf('A%d',count+1);
        xlsInsert = resbuf;
        xlswrite(filename,xlsInsert,'Sheet1',xlRange);
    end
end














