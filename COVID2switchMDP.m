p0 = importdata('covid-data/COVIDtransitions0_Greedy_50.mat');
p1 = importdata('covid-data/COVIDtransitions1_Greedy_50.mat');
Gs = importdata('covid-data/COVIDGs_Greedy_50.mat');
Gi = importdata('covid-data/COVIDGi_Greedy_50.mat');
Gr = importdata('covid-data/COVIDGr_Greedy_50.mat');
%p0 = importdata('covid-data/COVIDtransitions0_uniform_50.mat');
%p1 = importdata('covid-data/COVIDtransitions1_uniform_50.mat');
%Gs = 0:0.02:1;
%Gi = 0:0.02:1;
%Gr = 0:0.02:1;
lgs = length(Gs)-1;
lgi = length(Gi)-1;
%in the order of historical action 00,01,10,11
P0 = sparse(lgs*lgi*4,lgs*lgi*4);
P1 = sparse(lgs*lgi*4,lgs*lgi*4);
P0(1:lgs*lgi,1:lgs*lgi) = p0;%given action 0 00->00
P0(lgs*lgi+1:2*lgs*lgi,2*lgs*lgi+1:3*lgs*lgi) = p0;%given action 0 01->10
P0(2*lgs*lgi+1:3*lgs*lgi,2*lgs*lgi+1:3*lgs*lgi) = p0;%given action 0 10->10
P0(3*lgs*lgi+1:4*lgs*lgi,2*lgs*lgi+1:3*lgs*lgi) = p0;%given action 0 11->10
P1(1:lgs*lgi,lgs*lgi+1:2*lgs*lgi) = p1;%given action 1 00->01
P1(lgs*lgi+1:2*lgs*lgi,3*lgs*lgi+1:4*lgs*lgi) = p1;%given action 1 01->11
P1(2*lgs*lgi+1:3*lgs*lgi,2*lgs*lgi+1:3*lgs*lgi) = p1;%given action 0 10->10
P1(3*lgs*lgi+1:4*lgs*lgi,3*lgs*lgi+1:4*lgs*lgi) = p1;%given action 1 11->11
T=60;
file_path = 'covid-data/beta.csv';
% Load the CSV file into a MATLAB array
beta = readmatrix(file_path);
gamma = 0.7048;
P{1} = P0;
P{2} = P1;
R = ones(length(P0),2);
costr = 0.005;
for bs = 1:lgs
    for bi = 1:lgi
        idx1 = (bs-1)*lgi+bi;
        R(idx1,1) = -(Gi(bi)+Gi(bi+1))/2;%00
        R(idx1,2) = -(Gi(bi)+Gi(bi+1))/2-costr;
        R(lgs*lgi+idx1,1) = -(Gi(bi)+Gi(bi+1))/2;%01
        R(lgs*lgi+idx1,2) = -(Gi(bi)+Gi(bi+1))/2-costr;
        R(2*lgs*lgi+idx1,1) = -(Gi(bi)+Gi(bi+1))/2;%10
        R(2*lgs*lgi+idx1,2) = -(Gi(bi)+Gi(bi+1))/2-99999;% a very large number so this action is infesaible
        R(3*lgs*lgi+idx1,1) = -(Gi(bi)+Gi(bi+1))/2;%11
        R(3*lgs*lgi+idx1,2) = -(Gi(bi)+Gi(bi+1))/2-costr;
    end
end
skip = 0;
[V, policy, cpu_time] = mdp_finite_horizon(P, R, 1, T);
%% recovery the policies
s0 = 0.9999*ones(26,1);
i0 = 0.0001*ones(26,1);
r0 = 0*ones(26,1);
SS = zeros(T,1);
II = zeros(T,1);
RR = zeros(T,1);
%idx = find_index(s0,i0,Gs,Gi);
%disp(V(idx,1));
obj = 0;
actions = zeros(T,1);
stateflag = 1; %1-00 2-01 3-10 4-11
for t =1:T
    obj = obj + i0;
    [ss,ii,rr] = compute_total_SIR(s0,i0,r0);
    SS(t,1) = ss;
    II(t,1) = ii;
    RR(t,1) = rr;
    if stateflag == 1
        idx = find_index(ss,ii,Gs,Gi);
    elseif stateflag == 2
        idx = find_index(ss,ii,Gs,Gi)+lgs*lgi;
    elseif stateflag == 3
        idx = find_index(ss,ii,Gs,Gi)+2*lgs*lgi;
    else 
        idx = find_index(ss,ii,Gs,Gi)+3*lgs*lgi;
    end
    pol = policy(idx,t)-1;
    if stateflag == 1
        if pol ==0
            stateflag = 1;
        else
            stateflag = 2;
        end
    elseif stateflag == 2
        if pol ==0
            stateflag = 3;
        else
            stateflag = 4;
        end
    elseif stateflag == 3
        stateflag =3;
    else
        if pol ==0
            stateflag = 3;
        else
            stateflag = 4;
        end
    end
    actions(t,1) = pol;
    [s0,i0,r0] = SEIR(s0,i0,r0,beta, gamma,pol);
    %[s0,i0,r0] = compute_total_SIR(S,I,R);
end
figure
trj = zeros(T,4);
trj(:,1) = SS;
trj(:,2) = II;
trj(:,3) = RR;
trj(:,4) = actions;
plot(trj)
disp(obj)
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
    r = 1-s-i;
end

function idx = find_index(s,i,Gs,Gi)
    lgs = length(Gs)-1;
    lgi = length(Gi)-1;
    ls = max(find(Gs<s));
    li = max(find(Gi<=i));
    idx = (ls-1)*(lgi)+li;
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