function [ SortInList ] = pickL_element_algoH_Lin( H, L )
% Pick-L-element algorithm

% H matrix dimension: n*m
[n,m] = size(H);
% pickL_element_algo
k = L;
% Make list having sign pattern
InList=nchoosek([1:n],k);
CtA=[];
for Cti=1:k
   CtA=[-1*ones(2^(Cti-1),1),CtA;1*ones(2^(Cti-1),1),CtA];     
end
number_of_combinations = nchoosek(n,k); % sign pattern implemented
a2 = zeros(number_of_combinations,1);
% obtain each portion of Hx and get max portion 
for numerating=1:number_of_combinations
    fprintf('Pick Percent: %d/%d\n',numerating*2^(k-1),number_of_combinations*2^(k-1));
    v_target = zeros(1,2^k/2);
    for v_index=1:2^(k-1)
        index = InList(numerating,:);
        v = 0;
        f = zeros(n,1);
        f(abs(index(:,1)),1) = sign(index(:,1));
        b = [zeros(n,1);1;zeros(n+1,1)];
        I = eye(n);
        I(abs(index(:,1)),abs(index(:,1))) = 0;
        O = ones(1,n);
        O(1,abs(index(:,1))) = 0;
        IH = I*H; % n x m
        fTH = f'*H; % 1 x m
        fTHz = [fTH, zeros(1,n)];
        % linear programming
        Aprime = [IH, -I;zeros(1,m), O;-IH, -I;-fTHz];
        [y,~,exitflag] = linprog(-fTHz,Aprime,b);
        if exitflag == -3
            v = inf;
        end
        v = fTHz*y;
        v_target(1,v_index) = v/(v+1);
    end
    a2(numerating,1) = max(v_target);
end
InList_Portion_pickL = zeros(size(InList,1), size(InList,2)+1);
InList_Portion_pickL(:,1:end-1) = InList;
InList_Portion_pickL(:,end) = a2;
SortInList = -sortrows(-InList_Portion_pickL, L+1);
end
