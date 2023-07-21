%% Computable performance guarantees for compressed sensing matrices 
% Aug. 27, 2017
% By Myung (Michael) Cho
% Email: michael.myung.cho@gmail.com
%
% Code for generating Table 8 (Fourier matrix) in the paper (upper bound on
% alpha_k)
%-------------------------------------

addpath('./00_Rawdata/01_Low_dimension/02_Fourier/01_rho_0.5');
addpath('./00_Rawdata/01_Low_dimension/02_Fourier/02_rho_0.6');
addpath('./00_Rawdata/01_Low_dimension/02_Fourier/03_rho_0.7');
addpath('./00_Rawdata/01_Low_dimension/02_Fourier/04_rho_0.8');
addpath('./01_JNcode');
addpath('./02_Pick-l');

%% save
count=0;
% save data to xls file
clk=clock;
datestmp=strcat([date,'-',num2str(clk(4)),'h',num2str(clk(5))]);
filename=strcat('T8_NSP_LD_comp_Fourier-',datestmp);
filename=strcat(filename,'.xlsx');
col_name={'rho','case','kv','alphaPICK1','alphaPICK1','alphaPICK3','alphaTSA_LB_P1','alphaTSA_UB_P1','alphaTSA_LB_P2','alphaTSA_UB_P2','alphaTSA_LB_P3','alphaTSA_UB_P3','alphaLP','alphaSDP','Tpick1','Tpick2','Tpick3','Ttsa_p1','Ttsa_p2','Ttsa_p3','Tlp','Tsdp';};
xlswrite(filename,col_name);

tp=1;
%% Generate matrix
for irho=1:4
    switch irho
        case 1
            msize='_20x40_Four.mat';
        case 2
            msize='_24x40_Four.mat';
        case 3
            msize='_28x40_Four.mat';
        case 4
            msize='_32x40_Four.mat';
    end
    for icase=1:10
        scase=num2str(icase);
        if icase ~= 10
            matfile=strcat('0',scase,msize);
        else
            matfile=strcat(scase,msize);
        end
        load(matfile);
        %% Comparison starts
        % initialization
        InList1=[];InList2=[];InList3=[];
        T1sdp=0; T1lp=0; T1pick1=0; T1pick2=0; T1pick3=0; T1tsaP1=0; T1tsaP2=0; T1tsaP3=0;
        alphaPICK1=0; alphaPICK2=0; alphaPICK3=0; 
        alphaTSA_LB_P1=0; alphaTSA_UB_P1=0; alphaTSA_LB_P2=0; alphaTSA_UB_P2=0; alphaTSA_LB_P3=0; alphaTSA_UB_P3=0; alphaLP=0; alphaSDP=0;
        for kv=1:5
            %% Test SDP relaxation
            T0sdp=tic; [alphaSDP,~,~,~]=BoundCVX(P,kv); T1sdp=toc(T0sdp);% Using CVX
            %% Test LP relaxation
            T0lp=tic; alphaLP=BoundJNmosek(A,kv); T1lp=toc(T0lp);
            %% Pick-1
            T0pick1=tic;
            if kv==1
                InList1 = pickL_element_algo_Lin(A,1);
                alphaPICK1 = sum(InList1(1:kv,end));
            else
                alphaPICK1 = sum(InList1(1:kv,end));
            end
            T1pick1=toc(T0pick1);
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
            %% TSA pick2            
            if kv>=3
                cd ./04_TSA_StandAlone
                    T0tsaP3=tic;      [alphaTSA_LB_P3,alphaTSA_UB_P3]=BoundTSA(A,kv,3);       T1tsaP3=toc(T0tsaP3);
                cd ..                
            end
            
            % Display results
            resbuf=[irho,icase,kv,alphaPICK1,alphaPICK2,alphaPICK3,alphaTSA_P1,alphaTSA_P2,alphaTSA_P3,alphaLP,alphaSDP,T1pick1,T1pick2,T1pick3,T1tsaP1,T1tsaP2,T1tsaP3,T1lp,T1sdp];
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
end













