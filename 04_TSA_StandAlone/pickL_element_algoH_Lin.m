function [ InList_Portion_pickL ] = pickL_element_algoH_Lin( H, L )
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
    
    % CVX Solver
%     cvx_quiet true
%     cvx_begin
%         variable g_pickL(m);
%         expression v_pickL(n);
%         expression v_target_pickL;            
%         v_pickL = H*g_pickL;
%         vC_pickL = v_pickL;
%         for z = 1:k
%            v_target_pickL = v_target_pickL + sign(index(z,1))*v_pickL(abs(index(z,1)),1);
%            vC_pickL(abs(index(z,1)),1) = 0;
%         end
%         
%         maximize( v_target_pickL ); 
%         
%     subject to 
%         norm(vC_pickL,1)  <= 1
%         0 <= sign(index(:,1)).*v_pickL(abs(index(:,1)),1)
%     cvx_end        
%     if(cvx_optval == Inf)
%         v_target_pickL = inf; 
%     end
    

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
   
%     cvx_quiet true
%     cvx_begin
%         variable y(m+n); % [y(m) z(n)]
%                         
%         maximize( fTHz*y ); 
%          
%     subject to 
%         [IH, -I;zeros(1,m), O;-IH, -I;-fTHz]*y <= b
%     cvx_end
    % linear programming
 
    Aprime = [IH, -I;zeros(1,m), O;-IH, -I;-fTHz];
    [y,~,exitflag] = linprog(-fTHz,Aprime,b);
    if exitflag == -3
        v = inf;
    end
    v = fTHz*y;
    
    a2(numerating,1) = v/(v+1);
end
InList_Portion_pickL = zeros(size(InList_sign_pattern,1), size(InList_sign_pattern,2)+1);
InList_Portion_pickL(:,1:end-1) = InList_sign_pattern;
InList_Portion_pickL(:,end) = a2;
end
