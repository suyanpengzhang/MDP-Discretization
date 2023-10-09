%% prepare data
Gs = importdata('covid-data/COVIDGs_Greedy_50.mat');
Gi = importdata('covid-data/COVIDGi_Greedy_50.mat');
Gr = importdata('covid-data/COVIDGr_Greedy_50.mat');
lgs = length(Gs)-1;
lgi = length(Gi)-1;
transitions  = sparse(lgs*lgi,lgs*lgi);
file_path = 'covid-data/beta.csv';
% Load the CSV file into a MATLAB array
beta = readmatrix(file_path);
gamma = 0.7048;
numsample=10000;
s0 = zeros(numsample,1);
i0 = zeros(numsample,1);
r0 = zeros(numsample,1);
for i =1:numsample
    s0(i) = rand(); %draw from 0.7-0.999
    i0(i) = rand(); %draw from 0.001 - 0.1
    if s0(i)+i0(i)<=1
        r0(i) = 1-s0(i)-i0(i);
    else
        s0(i) = s0(i)/(s0(i)+i0(i));
        i0(i) = i0(i)/(s0(i)+i0(i));
    end
end
%%
code = zeros(lgs*lgi,3);
for i =1:lgi*lgs
    [ss,ii,rr] = reverse_find_index(i,Gs,Gi);
    code(i,1)=ss;
    code(i,2)=ii;
    code(i,3)=rr;
end
y = zeros(numsample,3);
b0 = zeros(numsample,lgs*lgi);
for i=1:numsample
    idx0 = find_index(s0(i),i0(i),Gs,Gi);
    b0(i,idx0)=1;
    [S,I,R] = SEIR(ones(26,1).*s0(i),ones(26,1).*i0(i),ones(26,1).*r0(i),beta, gamma,0);
    [S,I,R] = compute_total_SIR(S,I,R);
    y(i,1)=S;
    y(i,2)=I;
    y(i,3)=R;
end

%%
fun = @(x)sum((b0*x*code-y).^2,'all');
fun(eye(lgs*lgi))
%%
lb = zeros(lgs*lgi,lgs*lgi);
ub = ones(lgs*lgi,lgs*lgi);
A = [];
b = [];
Aeq = [];
beq = [];
x0 = sparse(eye(lgs*lgi));
nonlcon = @matrix;
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon);

%%
function [c,ceq] = matrix(x)
ceq = x*ones(2240,1)-ones(2240,1);
c = [];
end

%% functions
function [s,i,r] = reverse_find_index(idx,Gs,Gi)
    idx = idx - 1;
    lgs = length(Gs)-1;
    lgi = length(Gi)-1;
    idx_s = floorDiv(idx,(lgi))+1;
    idx = idx - lgi*(idx_s-1);
    idx_i = idx + 1;
    %disp(idx_s)
    %disp(idx_i)
    if idx_s == lgs+1
        s = Gs(idx_s);
    else
        s = (Gs(idx_s)+Gs(idx_s+1))/2;
    end
    if idx_i == lgi+1
        i = Gi(idx_i);
    else
        i = (Gi(idx_i)+Gi(idx_i+1))/2;
    end
    if s+i<=1
        r = 1-s-i;
    else
        s = s/(s+i);
        i = i/(s+i);
        r=0;
    end
end

function idx = find_index(s,i,Gs,Gi)
    lgs = length(Gs)-1;
    lgi = length(Gi)-1;
    ls = max(find(Gs<s));
    li = max(find(Gi<i));
    idx = (ls-1)*(lgi)+li;
end

%% functions: Discretized Discrete Time SEIR Model
%G is the discretization matrix
function [S,I,R] = SEIR_trj_dis(s0,i0,r0,beta, gamma,T,pol,Gs,Gi,Gr) 
    S = zeros(T,26);
    I = zeros(T,26);
    R = zeros(T,26);
    for t = 1:T
        for i = 1:26
            if s0(i,1) <= 0
                es0 = (Gs(1));
            elseif s0(i,1) >= 1
                es0 = (Gs(end));
            else
                es0 = (max(Gs(Gs<s0(i,1)))+min(Gs(Gs>=s0(i,1))))/2;
            end
            if i0(i,1) <= 0
                ei0 = (Gi(1));
            elseif i0(i,1) >= 1
                ei0 = (Gi(end));
            else
                ei0 = (max(Gi(Gi<i0(i,1)))+min(Gi(Gi>=i0(i,1))))/2;
            end
            if r0(i,1) <= 0
                er0 = (Gr(1));
            elseif r0(i,1) >= 1
                er0 = (Gr(end));
            else
                er0 = (max(Gr(Gr<r0(i,1)))+min(Gr(Gr>=r0(i,1))))/2;
            end
            S(t,i) = es0;
            I(t,i) = ei0;
            R(t,i) = er0;
        end
        action = pol(t);
        [s0,i0,r0] = SEIR(transpose(S(t,:)),transpose(I(t,:)),transpose(R(t,:)),beta, gamma,action);
    end
end
%% functions: Discrete Time SIR Model
function [S,I,R] = SEIR_trj(s0,i0,r0,beta, gamma,T,pol)
    S = zeros(T,26);
    I = zeros(T,26);
    R = zeros(T,26);
    for t = 1:T
        S(t,:) = s0;
        I(t,:) = i0;
        R(t,:) = r0;
        action = pol(t);
        [s0,i0,r0] = SEIR(s0,i0,r0,beta, gamma,action);
    end
end

function [S,I,R] = SEIR(s0,i0,r0,beta, gamma,action)
    delta_t = 1;
    S = zeros(26,1);
    I = zeros(26,1);
    R = zeros(26,1);
    if action == 1
        beta = beta.*0.2;% lockdown is effective at reducing 30% contacts
    end
    for i = 1:26
        S(i,1) = s0(i,1);
        for j = 1:26
            S(i,1) = S(i,1)- beta(j,i)*(s0(i,1))*i0(j,1)*delta_t;
        end          
        I(i,1) = i0(i,1)- gamma*i0(i,1)*delta_t;
        for j = 1:26
            I(i,1) = I(i,1)+ beta(j,i)*(s0(i,1))*i0(j,1)*delta_t;
        end 
        R(i,1) = r0(i,1) + gamma*i0(i,1)*delta_t;
        S(i,1) = round(S(i,1),7);
        I(i,1) = round(I(i,1),7);
        R(i,1) = round(R(i,1),7);
    end
end

function [S_,I_,R_] = compute_total_SIR(S,I,R)
    totalpop = [420697,443569,344450,526877,899111,339399,419797,308499,140361,547523,354750,479505,287613,666399,278815,193899,166374,379199,356465,195082,407864,321720,201739,411617,469439,465691];
    S_ = sum(transpose(totalpop).*S)/sum(totalpop);
    I_ = sum(transpose(totalpop).*I)/sum(totalpop);
    R_ = sum(transpose(totalpop).*R)/sum(totalpop);
end
function [S_,I_,R_] = compute_total_SIR_trj(S,I,R,T)
    totalpop = [420697,443569,344450,526877,899111,339399,419797,308499,140361,547523,354750,479505,287613,666399,278815,193899,166374,379199,356465,195082,407864,321720,201739,411617,469439,465691];
    S_ = zeros(T,1);
    I_ = zeros(T,1);
    R_ = zeros(T,1);
    for t = 1:T
        S_(t,1) = sum(totalpop.*S(t,:))/sum(totalpop);
        I_(t,1) = sum(totalpop.*I(t,:))/sum(totalpop);
        R_(t,1) = sum(totalpop.*R(t,:))/sum(totalpop);
    end
end