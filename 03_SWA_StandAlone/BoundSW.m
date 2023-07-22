function [ alphaSW, number_of_step ] = BoundSW( A,kv,L)
% Sandwiching algo. for exact alpha_Kv
% initialize

[m,n] = size(A);
alphaSW = 0;
number_of_step = 0; 

%% pick-l-element
if m <= 0.5*n
    InListL = pickL_element_algo_Lin_for_SW(A,L); 
else
    H = null(A);
    InListL = pickL_element_algoH_Lin_for_SW(H,L);
end

if (kv >= L) && (L >= 1)
% make table for tracking exact alpha_kv
M = nchoosek([1:n],kv); %all the (n choose k) k-sets
keeptrack = zeros(size(M,1),kv+1); % kv+1 th column of keeptrack matrix is for saving portion
keeptrack(:,1:kv) = M; % 1~k th columns of keeptrack matrix are for indexs
IN = nchoosek([1:kv],L);

% cheap upper bound
for i = 1:size(M,1)
   set = M(i,:);
   k_set_sum = 0;
   for j = 1:size(IN,1)
        indices = set(1, IN(j,:));
        pick2_portion_for_weak_bound = find_portion(InListL, indices, L); 
        k_set_sum = k_set_sum + pick2_portion_for_weak_bound;
   end
   k_set_portion = k_set_sum / (nchoosek(kv-1,L-1));
   keeptrack(i,end) = k_set_portion;
end

% sort table 
sort_keeptrack_by_portion = -sortrows(-keeptrack, kv+1);

% sandwiching operation
limit_portion = 0;
iMax = size(M,1);
for p = 1:iMax
    upper_bound = sort_keeptrack_by_portion(p,end);
    if limit_portion < upper_bound
%         k_set_portion_tighter_bound = linear_prog( kv, IN, InListL, sort_keeptrack_by_portion(p,1:kv), L);
%         if limit_portion < k_set_portion_tighter_bound
%             pick_k_index(1,1:kv) = sort_keeptrack_by_portion(p,1:kv);
%             real_portion = pickK_element_once(A, pick_k_index');
%             sort_keeptrack_by_portion(p,kv+1) = real_portion;
%         else
%             sort_keeptrack_by_portion(p,kv+1) = k_set_portion_tighter_bound;            
%         end
        %% Since it is fast enough, linear programming bound is removed
        pick_k_index(1,1:kv) = sort_keeptrack_by_portion(p,1:kv);
        real_portion = pickK_element_once(A, pick_k_index');
        sort_keeptrack_by_portion(p,kv+1) = real_portion;
    else
        number_of_step = p-1;
        break;
    end
    if real_portion > limit_portion
        limit_portion = real_portion;
    end
    fprintf('(SWA Percent: %d/%d) LB = %f, UB = %f \n',p,iMax,limit_portion,upper_bound);
end
alphaSW = limit_portion;
end
end
 


