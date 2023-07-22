function [alphaTSA_LB,alphaTSA_UB] = Tree_search(A, kv, L)
% Tree Search Algorithm - Verify the null space condition
% 
% by Michael Myung Cho
% On Dec. 09, 2013
%

% clear all
% profile -memory on

%% size
% n = 40;
% m = 32; % A: m*n
% kv = 5; % sparsity
% L = 2; % L-step Tree search algorithm
% tp = 4; % Gaussian

%% A matrix
% [rss,pnt]=AAGenerateSensingMatrix(tp,n,m);
% A=rss.A;n=rss.n;m=rss.m;
% A_InList = load('..\tree_search_algo.mat');
% A = A_InList.A;

% Tree search operation starts
% operation_time = tic;

%% H matrix permuate and percalcuation
[permutated_A, sorted_alphaL_list, sorted_alphaR_list, sorted_alphaL_list_for_layerL] = Matrix_Permutation(A, kv, L);
sorted_alphaL_list;
sorted_alphaR_list;
%% tree search starts

[alphaTSA_LB,alphaTSA_UB] = make_tree(permutated_A, sorted_alphaL_list, sorted_alphaR_list, sorted_alphaL_list_for_layerL, kv, L);

%% for saving date
%  operation_time_for_exact = toc(operation_time);
%  operation_time_for_exact_min = operation_time_for_exact/60;
%  save('operation_time.txt','operation_time_for_exact_min','-ASCII');
%  save('tree_search_algo.mat');
% profile viewer
% p = profile('info');
% profsave(p,'profile_results')
% fprintf('one instance memory allocation\n');
% myStruct=dlnode(0);
% whos myStruct

%end

end

