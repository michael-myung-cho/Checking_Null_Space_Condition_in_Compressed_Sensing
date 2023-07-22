function [ sum_remaining_alpha ] = remaining_upper_bound( num_layer, kv, L, best_candidate_node, sorted_alphaL_list, sorted_alphaR_list )
%REMAINING_UPPER_BOUND Summary of this function goes here
%   Detailed explanation goes here

[m2, n2] = size(sorted_alphaL_list);
[m1, n1] = size(sorted_alphaR_list);

% remaining values        
sum_remaining_alpha = 0;
remaining_Lsteps = fix((kv-num_layer)/L); % divider (L step)
remaining_Rsteps = rem(kv-num_layer,L); % remaining step (R step)
if num_layer == 0
    current_max_index = 0;
else
    current_max_index = abs(best_candidate_node.Property7(1,end));
end

% Remaining step (remaining part)
% ex) k=5, num_layer=2: (1,10) -> remaining part: alpha(11,13) + alpha(13)
if remaining_Rsteps == 0
    % nothing work
else
    i = 1;
    while(m1 - i)
        % among all indexs in a branch, best_candidate_node.Property7(1,2) has maximum index
        if  (current_max_index + remaining_Lsteps*L < abs(sorted_alphaR_list(i,1)))
            sum_remaining_alpha = sum_remaining_alpha + sorted_alphaR_list(i,end);
            break;
        end
        i = i + 1;
    end    
end

% 2 step remaining part
if remaining_Lsteps == 0
    % nothing work
else
    for q = 1:remaining_Lsteps
        i = 1;
        while(m2 - i)
            % among all indexs in a branch, best_candidate_node.Property7(1,2) has maximum index
            if  current_max_index < abs(sorted_alphaL_list(i,1))
                sum_remaining_alpha = sum_remaining_alpha + sorted_alphaL_list(i,end);
                current_max_index = current_max_index + L;
                break;
            end
            i = i + 1;
        end
    end
end
           
end

