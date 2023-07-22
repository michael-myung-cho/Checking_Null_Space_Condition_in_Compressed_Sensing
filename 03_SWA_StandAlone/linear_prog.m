function [ min_k_sum ] = linear_prog( kv, IN, InList_Portion, set, L )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

pick2_portion(1,:) = zeros(1, size(IN,1));
cvx_begin
    variable x(1,1:kv);
    min_k_sum = sum( x(1,1:kv) );
    maximize( min_k_sum );
subject to 
    for in_index = 1:size(IN,1)
        pick2_portion(1,in_index) = find_portion(InList_Portion, set(1, IN(in_index,:)),L);
        sum( x(1,IN(in_index,:)) ) <= pick2_portion(1,in_index); 
    end
    for x_index = 1:kv
        x(1,x_index) >= 0;
    end
cvx_end
end
