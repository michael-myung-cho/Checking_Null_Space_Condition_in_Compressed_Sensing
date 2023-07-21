%% Computable performance guarantees for compressed sensing matrices
% Aug. 21, 2017
% By Myung (Michael) Cho
% Email: michael.myung.cho@gmail.com
%
% Code for generating Table 2 (n = 1024) in the paper (Lower bound on sparsity k)
%----------------------------------


clear
addpath('./00_Rawdata/02_High_dimension');
addpath('./01_JNcode');
addpath('./02_Pick-l');

filename = 'T2_NSP_HD_Comp_results.xlsx';
col_name = {'file','m,','n','k_a1','k_a1s','T1lpa1','T1lpa1s';};
xlswrite(filename,col_name);
countN = 1;

for icase=5:5
    %% Load sensing matrix
    switch icase
        case 1
            matfile = '102x1024_NSP_HD_Gaussian.mat';
        case 2
            matfile = '205x1024_NSP_HD_Gaussian.mat';
        case 3
            matfile = '307x1024_NSP_HD_Gaussian.mat';
        case 4
            matfile = '410x1024_NSP_HD_Gaussian.mat';
        case 5
            matfile = '512x1024_NSP_HD_Gaussian.mat';
        case 6
            matfile = '614x1024_NSP_HD_Gaussian.mat';
        case 7
            matfile = '717x1024_NSP_HD_Gaussia.mat';
        case 8
            matfile = '819x1024_NSP_HD_Gaussian.mat';
        case 9
            matfile = '922x1024_NSP_HD_Gaussian.mat';
    end
    load(matfile);
    [m,n]=size(A);
    H=null(A);
    rho=m/n;
    
    %% Initialization
    s1best=0; LB_kv_pick1=0; LB_kv_pick2=0; LB_kv_pick3=0; T1lpa1=0; T0pick1=0; T0pick2=0; T0pick3=0; 
    %% LP relaxation (Code from JN)
    T0lpa1=tic;
    disp(['Computing sparsity certificates for ',num2str(m),'x',num2str(n),' sensing matrix...'])
    disp('Computing lower bound for sparsity by alpha_1 bounding...')
    [s,F]=GetAlpha1Bound(A);
    disp(['Sparsity, certified using MI: ',num2str(s.s_mi)])
    disp(['Sparsity, certified using alpha_1 bound: ',num2str(s.s_alpha.s_alpha_1)])
    s1best = num2str(s.s_alpha.s_alpha_1);
    scont=1;
    s1=s;
    T1lpa1=toc(T0lpa1);
    
    %% Pick-1
    T1pick1=0;
    if rho <= 0.5
        T0pick1=tic;    SInList1 = pickL_element_algo_Lin(A,1);  T1pick1=toc(T0pick1);
    else
        T0pick1=tic;    SInList1 = pickL_element_algoH_Lin(H,1);  T1pick1=toc(T0pick1);
    end
    % lower bound on k via Pick-1
    alphaK_pre=0;
    LB_kv_pick1=0;
    for kv=1:n
        alphaK=sum(SInList1(1:kv,end));
        if (alphaK_pre < 0.5) && (alphaK >= 0.5)
            LB_kv_pick1=kv-1;
            break;
        end
        alphaK_pre=alphaK;
    end
    fprintf('(%d x %d) LB_k_pick1: %d \n',m,n,LB_kv_pick1);    
    %% Pick-2
    if rho <= 0.5
        T0pick2=tic;    SInList2 = pickL_element_algo_Lin(A,2);  T1pick2=toc(T0pick2);
    else
        T0pick2=tic;    SInList2 = pickL_element_algoH_Lin(H,2);  T1pick2=toc(T0pick2);
    end
    % lower bound on k via Pick-2
    alphaK_pre=0;
    LB_kv_pick2=0;
    for kv=2:nchoosek(n,2)
        alphaK=sum(SInList2(1:nchoosek(kv,2),end))/(nchoosek(kv-1,2-1));
        if (alphaK_pre < 0.5) && (alphaK >= 0.5)
            LB_kv_pick2=kv-1;
            break;
        end
        alphaK_pre=alphaK;
    end
    fprintf('(%d x %d) LB_k_pick2: %d \n',m,n,LB_kv_pick2);
    %% Pick-3
%     if rho <= 0.5
%         T0pick3=tic;    SInList3 = pickL_element_algo_Lin(A,3);  T1pick3=toc(T0pick3);
%     else
%         T0pick3=tic;    SInList3 = pickL_element_algoH_Lin(H,3);  T1pick3=toc(T0pick3);
%     end
%     % lower bound on k via Pick-3
%     alphaK_pre=0;
%     LB_kv_pick3=0;
%     for kv=3:nchoosek(n,3)
%         alphaK=sum(SInList3(1:nchoosek(kv,3),end))/(nchoosek(kv-1,3-1));
%         if (alphaK_pre < 0.5) && (alphaK >= 0.5)
%             LB_kv_pick3=kv-1;
%             break;
%         end
%         alphaK_pre=alphaK;
%     end
    %% Display results
    resbufTSA=[{matfile},s1best,LB_kv_pick1,LB_kv_pick2,LB_kv_pick3,T1lpa1,T0pick1,T0pick2,T0pick3 ];
    disp('**************************************************************************************************');
    countN=countN+1;
    % save data to xls file
    xlRangeTSA = sprintf('A%d',countN);
    xlsInsertTSA = resbufTSA;
    xlswrite(filename,xlsInsertTSA,'Sheet1',xlRangeTSA);
end





