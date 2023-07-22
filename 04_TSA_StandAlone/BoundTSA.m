function [ alphaTSA_LB,alphaTSA_UB ] = BoundTSA(A,kv,L)
%GAUSSIAN Summary of this function goes here
%   Detailed explanation goes here

% if (kv == 1 || kv == 2)
%     L = 1;
% else
%     L = 1;
% end
% A matrix dimension: m*n
[alphaTSA_LB,alphaTSA_UB] = Tree_search(A, kv, L);

end 


