Gs = importdata('covid-data/Gs_density_50_COVID.mat');
Gi = importdata('covid-data/Gi_density_50_COVID.mat');
Gr = importdata('covid-data/Gr_density_50_COVID.mat');
Gs = 0:0.02:1;
Gi = 0:0.008:0.40;
Gi(51)=1;
Gr = 0:0.02:1;
%Gi = 0:0.02:1;
%Gi(51)=1;
%Gr = 0:0.02:1;
%Gs = 0:0.03333:1;
%Gi = 0:0.03333:1;
%Gi = 0:0.02:0.40;
%Gi(201)=1;
%Gr = 0:0.03333:1;
%%
lgs = length(Gs)-1;
lgi = length(Gi)-1;
transitions  = sparse(lgs*lgi,lgs*lgi);
file_path = 'covid-data/beta.csv';
% Load the CSV file into a MATLAB array
beta = readmatrix(file_path);
gamma = 0.7048;

s0 = 0.9999;
i0 = 0.0001;
r0 = 1-s0-i0;
%%
[S,I,R] = SEIR(ones(26,1).*s0,ones(26,1).*i0,ones(26,1).*r0,beta, gamma,0);
[S,I,R] = compute_total_SIR(S,I,R)
%%
tic
count = 1;
for bs = 1:lgs
    disp(bs)
    for bi = 1:lgi
        cs = (Gs(bs)+Gs(bs+1))/2;
        ci = (Gi(bi)+Gi(bi+1))/2;
        idx1 = (bs-1)*lgi+bi;
        %counter = 0;
        if cs+ci <= 1
            cr = 1-cs-ci;
            [S,I,R] = SEIR(ones(26,1).*cs,ones(26,1).*ci,ones(26,1).*cr,beta, gamma,0);
            [S,I,R] = compute_total_SIR(S,I,R);
            idx2 = find_index(S,I,Gs,Gi);
            transitions(idx1,idx2) = transitions(idx1,idx2)+1;
            count = count+1;
        end
        for iter = 1:1000
            rs = Gs(bs)+(Gs(bs+1)-Gs(bs))*rand();
            ri = Gi(bi)+(Gi(bi+1)-Gi(bi))*rand();
            if rs+ri<=1
                %counter = counter+1;
                rr = 1 - rs - ri;
                [S,I,R] = SEIR(ones(26,1).*rs,ones(26,1).*ri,ones(26,1).*rr,beta, gamma,0);
                [S,I,R] = compute_total_SIR(S,I,R);
                idx2 = find_index(S,I,Gs,Gi);
                transitions(idx1,idx2) = transitions(idx1,idx2)+1;
            end
        end
        sums = sum(transitions(idx1,:));
        if sums>0
            transitions(idx1,:) = transitions(idx1,:)./sums;
        else
            transitions(idx1,idx1) = 1;
        end
        count = count +1;
    end
end
transitions0 = transitions;

transitions  = sparse(lgs*lgi,lgs*lgi);
count = 1;
for bs = 1:lgs
    disp(bs)
    for bi = 1:lgi
        cs = (Gs(bs)+Gs(bs+1))/2;
        ci = (Gi(bi)+Gi(bi+1))/2;
        idx1 = (bs-1)*lgi+bi;
        %counter = 0;
        if cs+ci <= 1
            cr = 1-cs-ci;
            [S,I,R] = SEIR(ones(26,1).*cs,ones(26,1).*ci,ones(26,1).*cr,beta, gamma,1);
            [S,I,R] = compute_total_SIR(S,I,R);
            idx2 = find_index(S,I,Gs,Gi);
            transitions(idx1,idx2) = transitions(idx1,idx2)+1;
            count = count+1;
        end
        for iter = 1:1000
            rs = Gs(bs)+(Gs(bs+1)-Gs(bs))*rand();
            ri = Gi(bi)+(Gi(bi+1)-Gi(bi))*rand();
            if rs+ri<=1
                %counter = counter+1;
                rr = 1 - rs - ri;
                [S,I,R] = SEIR(ones(26,1).*rs,ones(26,1).*ri,ones(26,1).*rr,beta, gamma,1);
                [S,I,R] = compute_total_SIR(S,I,R);
                idx2 = find_index(S,I,Gs,Gi);
                transitions(idx1,idx2) = transitions(idx1,idx2)+1;
            end
        end
        sums = sum(transitions(idx1,:));
        if sums>0
            transitions(idx1,:) = transitions(idx1,:)./sums;
        else
            transitions(idx1,idx1) = 1;
        end
        count = count +1;
    end
end
transitions1 = transitions;
toc
%% Compute error over samples
%{
T=40;
samples = importdata('samples_for_compare.mat');
%pol_samples = importdata('policy_for_compare.mat');
pol_samples = randi(2,150,T);
num_samples = 5;
errors = zeros(num_samples,1);
for samp =1:num_samples
    pol = pol_samples(samp,1:T);
    trj_s = zeros(1,T);
    trj_i = zeros(1,T);
    trj_r = zeros(1,T);
    idx0 = find_index(samples(samp,1),samples(samp,2),Gs,Gi);
    b0 = sparse(1,lgs*lgi);
    b0(1,idx0) = 1;
    T = 10;
    for t = 1:T
        disp(t)
        if pol(t)==0
            b0 = b0*transitions0;
        else
            b0 = b0*transitions1;
        end
        nonzeros_idx = find(b0);
        values = nonzeros(b0);
        ss = 0;
        ii = 0;
        rr = 0;
        for k = 1:length(nonzeros_idx)
            [s,i,r] = reverse_find_index(nonzeros_idx(k),Gs,Gi);
            if r<0
                s = s/(s+i);
                i = i/(s+i);
                r = 0;
                %disp(nonzeros_idx(k))
                %disp(r)
                %disp(reverse_find_index(nonzeros_idx(k),Gs,Gi))
            end
            ss = ss + b0(1,nonzeros_idx(k))*(s);
            ii = ii + b0(1,nonzeros_idx(k))*(i);
            rr = rr + b0(1,nonzeros_idx(k))*(r);
        end
        trj_s(1,t) = ss;
        trj_i(1,t) = ii;
        trj_r(1,t) = rr;
    end
    [S,I,R] = SEIR_trj(ones(26,1).*samples(samp,1),ones(26,1).*samples(samp,2),ones(26,1).*samples(samp,3),beta,gamma,T,pol);
    [S,I,R] = compute_total_SIR_trj(S,I,R,T);
    err = 0;
    for t =1:T
        err = err + (S(t)-trj_s(t))^2;
        err = err + (I(t)-trj_i(t))^2;
        err = err + (R(t)-trj_r(t))^2;
    end
    errors(samp,1) = err;
end
disp('finished')
%}

%% Draw trajectory
%transitions0 = importdata('transitions0_uniform_100.mat');
%transitions1 = importdata('transitions1_uniform_100.mat');
T=40;
pol = zeros(T,1);
trj_s = zeros(1,T);
trj_e = zeros(1,T);
trj_i = zeros(1,T);
trj_r = zeros(1,T);
idx0 = find_index(s0,i0,Gs,Gi);
b0 = sparse(1,lgs*lgi);
b0(1,idx0) = 1;
for t = 1:T
    disp(t)
    if pol(t)==0
        b0 = b0*transitions0;
    else
        b0 = b0*transitions1;
    end
    nonzeros_idx = find(b0);
    values = nonzeros(b0);
    ss = 0;
    ii = 0;
    rr = 0;
    for k = 1:length(nonzeros_idx)
        [s,i,r] = reverse_find_index(nonzeros_idx(k),Gs,Gi);
        if r<0
            s = s/(s+i);
            i = i/(s+i);
            r = 0;
            %disp(nonzeros_idx(k))
            %disp(r)
            %disp(reverse_find_index(nonzeros_idx(k),Gs,Gi))
        end
        ss = ss + b0(1,nonzeros_idx(k))*(s);
        ii = ii + b0(1,nonzeros_idx(k))*(i);
        rr = rr + b0(1,nonzeros_idx(k))*(r);
    end
    trj_s(1,t) = ss;
    trj_i(1,t) = ii;
    trj_r(1,t) = rr;
end

trj1 = zeros(T,6);
%pol = zeros(T,1);
[S,I,R] = SEIR_trj(ones(26,1).*s0,ones(26,1).*i0,ones(26,1).*r0,beta,gamma,T,pol);
[S,I,R] = compute_total_SIR_trj(S,I,R,T);
%[S,I,R] = SEIR_trj_dis(s0,i0,r0,beta,gamma,T,pol,Gs,Gi,Gr);

trj1(:,1) = trj_s(1:T);
trj1(:,2) = trj_i(1:T);
trj1(:,3) = trj_r(1:T);

trj1(:,4) = S;
trj1(:,5) = I;
trj1(:,6) = R;
%%
s0 = 0.99;
i0 = 0.001;
r0 = 1-s0-i0;
time = zeros(T,1);
for i = 1:T
    time(i)=i;    
end
figure
plot(time,trj1(:,1),'r', ...
    time,trj1(:,2),'k', ...
    time,trj1(:,3),'b', ...
    time,trj1(:,4),'r--', ...
    time,trj1(:,5),'k--', ...
    time,trj1(:,6),'b--','LineWidth',3);
ax = gca;
ax.FontSize = 16; 
title({'Uniform Discretization'},'Fontsize',18)
xlabel('Time: week','FontSize',18)
ylabel('Proportions','FontSize',18)
legend('Uniform S', 'Uniform I','Uniform R','Actual S','Actual I','Actual R','Fontsize',14)
%%
%{
errors = zeros(10,1);
for iiii = 1:10
    disp(ii)
    pol = randi(2,T,1)-1;
    %pol = zeros(T,1);
    trj_s = zeros(1,T);
    trj_e = zeros(1,T);
    trj_i = zeros(1,T);
    trj_r = zeros(1,T);
    idx0 = find_index(s0,i0,Gs,Gi);
    b0 = sparse(1,lgs*lgi);
    b0(1,idx0) = 1;
    T = 10;
    for t = 1:T
        if pol(t)==0
            b0 = b0*transitions0;
        else
            b0 = b0*transitions1;
        end
        nonzeros_idx = find(b0);
        values = nonzeros(b0);
        ss = 0;
        ii = 0;
        rr = 0;
        for k = 1:length(nonzeros_idx)
            [s,i,r] = reverse_find_index(nonzeros_idx(k),Gs,Gi);
            if r<0
                s = s/(s+i);
                i = i/(s+i);
                r = 0;
                %disp(nonzeros_idx(k))
                %disp(r)
                %disp(reverse_find_index(nonzeros_idx(k),Gs,Gi))
            end
            ss = ss + b0(1,nonzeros_idx(k))*(s);
            ii = ii + b0(1,nonzeros_idx(k))*(i);
            rr = rr + b0(1,nonzeros_idx(k))*(r);
        end
        trj_s(1,t) = ss;
        trj_i(1,t) = ii;
        trj_r(1,t) = rr;
    end
    [S,I,R] = SEIR_trj_dis(s0,i0,r0,beta,gamma,T,pol,Gs,Gi,Gr);
    err = 0;
    for t =1:T
        err = err + (S(t)-trj_s(t))^2;
        err = err + (I(t)-trj_i(t))^2;
        err = err + (R(t)-trj_r(t))^2;
    end
    errors(iiii,1) = err;
end
%}
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