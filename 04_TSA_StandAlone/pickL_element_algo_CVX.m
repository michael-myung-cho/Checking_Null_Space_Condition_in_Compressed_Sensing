function [ InList_Portion_pickL ] = pickL_element_algo( A, L )
% Pick-L-element algorithm

% H matrix dimension: n*m
[m,n] = size(A);

% pickL_element_algo
k = L;

% Make list having sign pattern
InList=nchoosek([1:n],k);
CtA=[];
for Cti=1:k
   CtA=[-1*ones(2^(Cti-1),1),CtA;1*ones(2^(Cti-1),1),CtA];     
end
InList_sign_pattern = [];
for In=2^(k-1)+1:2^k
    InList_row = [];
    for In2=1:k
        InList_row = [InList_row, CtA(In,In2)*InList(:,In2)]; % sign pattern indexes
    end
    InList_sign_pattern = [InList_sign_pattern; InList_row];
end
number_of_combinations = nchoosek(n,k)*2^(k-1); % sign pattern implemented
a2 = zeros(number_of_combinations,1);

% obtain each portion of Hx and get max portion 
for numerating=1:floor(number_of_combinations/100)
    fprintf('Pick Percent: %d/%d\n',numerating,number_of_combinations);
    index = InList_sign_pattern(numerating,:)';
    v_target_pickL = 0;
    % ADMM Solver
%    v_target_pickL = ADMM(A, index);
    
    % CVX Solver
    cvx_quiet true
    cvx_begin
        variable g_pickL(n);
        expression v_pickL(m);
        expression v_target_pickL;            
        v_pickL = A*g_pickL;
        for z = 1:k
           v_target_pickL = v_target_pickL + sign(index(z,1))*g_pickL(abs(index(z,1)),1);
        end             
        maximize( v_target_pickL ); 
    subject to 
        norm(g_pickL,1)  <= 1
        0 <= v_pickL <= 0
        0 <= sign(index(:,1)).*g_pickL(abs(index(:,1)),1)
    cvx_end        
    if(cvx_optval == Inf)
        v_target_pickL = inf; 
    end
    a2(numerating,1) = v_target_pickL;  
end
InList_Portion_pickL = zeros(size(InList_sign_pattern,1), size(InList_sign_pattern,2)+1);
InList_Portion_pickL(:,1:end-1) = InList_sign_pattern;
InList_Portion_pickL(:,end) = a2;
end
