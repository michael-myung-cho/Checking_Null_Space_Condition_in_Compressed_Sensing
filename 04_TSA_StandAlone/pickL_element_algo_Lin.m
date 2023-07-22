function [ InList_Portion_pickL ] = pickL_element_algo_Lin( A, L )
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
for numerating=1:number_of_combinations
    fprintf('Pick Percent: %d/%d\n',numerating,number_of_combinations);
    index = InList_sign_pattern(numerating,:)';
    v = 0;
    % ADMM Solver
%    v_target_pickL = ADMM(A, index);

    % CVX Solver
%     cvx_quiet true
%     cvx_begin
%         variable x(2*n);
% 
%         maximize( f'*x ); 
%     
%     subject to 
%         Aprime*x <= b
%         B*x == c
%     cvx_end
%     
%     v = f'*x;
%     if(cvx_optval == Inf)
%         v = inf; 
%     end

    % linear programming
    f = zeros(2*n,1);
    f(abs(index(:,1)),1) = sign(index(:,1));
    b = [zeros(2*n+1,1);1];
    B = [A, zeros(m,n)];
    Aprime = [eye(n), -eye(n);-eye(n), -eye(n);-f';zeros(1,n),ones(1,n)];
    c = zeros(m,1);    
     [x,~,exitflag] = linprog(-f,Aprime,b,B,c); 
     if exitflag == -3 
         v = inf;
     end
     v = f'*x;
     a2(numerating,1) = v;  
     
end
InList_Portion_pickL = zeros(size(InList_sign_pattern,1), size(InList_sign_pattern,2)+1);
InList_Portion_pickL(:,1:end-1) = InList_sign_pattern;
InList_Portion_pickL(:,end) = a2;
end
