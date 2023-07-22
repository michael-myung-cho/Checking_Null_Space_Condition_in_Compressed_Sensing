function [ InList_Portion_pickL ] = pickL_element_algo_Lin_for_SW( A, L )
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
number_of_combinations = nchoosek(n,k); % sign pattern NOT implemented
a2 = zeros(number_of_combinations,1);

% obtain each portion of Hx and get max portion 
for numerating=1:number_of_combinations
    fprintf('Pick Percent: %d/%d\n',numerating*2^(k-1),number_of_combinations*2^(k-1));
    v_target = zeros(1,2^(k-1));
    for v_index=1:2^(k-1)
        index = InList(numerating,:)';
        v = 0;
        % linear programming
        indexS = CtA(v_index,:)'.*index;
        f = zeros(2*n,1);
        f(abs(index(:,1)),1) = sign(indexS(:,1));
        b = [zeros(2*n+1,1);1];
        B = [A, zeros(m,n)];
        Aprime = [eye(n), -eye(n);-eye(n), -eye(n);-f';zeros(1,n),ones(1,n)];
        c = zeros(m,1);
        [x,~,exitflag] = linprog(-f,Aprime,b,B,c);
        if exitflag == -3
            v = inf;
        end
        v = f'*x;
        v_target(1,v_index) = v;
    end
    a2(numerating,1) = max(v_target);
end
InList_Portion_pickL = zeros(size(InList,1), size(InList,2)+1);
InList_Portion_pickL(:,1:end-1) = InList;
InList_Portion_pickL(:,end) = a2;
end
