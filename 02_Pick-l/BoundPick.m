function [ alphaPICK1,alphaPICK2,alphaPICK3,InList1,InList2,InList3 ] = BoundPick( A,kv,InList1,InList2,InList3)
%GAUSSIAN Summary of this function goes here
%   Detailed explanation goes here

% A matrix dimension: m*n
[m,n] = size(A);
alphaPICK1 = 0;
alphaPICK2 = 0;
alphaPICK3 = 0;

if kv == 1
    InList1 = pickL_element_algo(A,1); 
    alphaPICK1 = sum(InList1(1:kv,end));
    alphaPICK2 = alphaPICK1;
    alphaPICK3 = alphaPICK1;
elseif kv == 2
    % alpha_kv from Pick-1-element algo
    alphaPICK1 = sum(InList1(1:kv,end));

    % alpha_kv from Pick-2-element algo
    InList2 = pickL_element_algo(A,2); 
    alphaPICK2 = sum(InList2(1:nchoosek(kv,2),end))/(nchoosek(kv-1,2-1));
    alphaPICK3 = alphaPICK2;
else
    % alpha_kv from Pick-1-element algo
    alphaPICK1 = sum(InList1(1:kv,end));
    
    % alpha_kv from Pick-2-element algo
    alphaPICK2 = sum(InList2(1:nchoosek(kv,2),end))/(nchoosek(kv-1,2-1));
    
    % alpha_kv from Pick-2-element algo    
    InList3 = pickL_element_algo(A,3);
    alphaPICK3 = sum(InList3(1:nchoosek(kv,3),end))/(nchoosek(kv-1,3-1));
end
 


