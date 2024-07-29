Gs = importdata('Gs_greedy_30_MAE.mat');
Gi = importdata('Gi_greedy_30_MAE.mat');
Gr = importdata('Gr_greedy_30_MAE.mat');
%Gs = 0:0.0025:1;
%Gi = 0:0.0025:1;
%Gi = 0:0.001:0.40;
%Gi(401)=1;
%Gr = 0:0.0025:1;
%%
lgs = length(Gs)-1;
lgi = length(Gi)-1;
transitions  = sparse(lgs*lgi,lgs*lgi);
count = 0;
T = 10;
beta = 0.2*7;
theta = 0.25*7;
gamma = 0.07*7;
s0 = rand();
i0 = rand();
r0 = rand();
sums = s0+r0+i0;
s0 = s0/sums;
i0 = i0/sums;
r0 = r0/sums;
s0 = 0.99;
i0 = 0.001;
r0 = 1-s0-i0;
%%

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
            [S,I,R] = SEIR(cs,ci,cr,beta, gamma,0);
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
                [S,I,R] = SEIR(rs,ri,rr,beta, gamma,0);
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
            [S,I,R] = SEIR(cs,ci,cr,beta, gamma,1);
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
                [S,I,R] = SEIR(rs,ri,rr,beta, gamma,1);
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
%% Compute error over samples

samples = importdata('samples_for_compare.mat');
pol_samples = importdata('policy_for_compare.mat');
num_samples = 9;
T=10;
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
    [S,I,R] = SEIR_trj(samples(samp,1),samples(samp,2),samples(samp,3),beta,gamma,T,pol);
    err = 0;
    for t =1:T
        err = err + (S(t)-trj_s(t))^2;
        err = err + (I(t)-trj_i(t))^2;
        err = err + (R(t)-trj_r(t))^2;
    end
    errors(samp,1) = err;
end
disp('finished')


%% Draw trajectory eval algo 1
Gs = importdata('Gs_density_100_new.mat');
Gi = importdata('Gi_density_100_new.mat');
Gr = importdata('Gr_density_100_new.mat');
Gs = importdata('Gs_greedy_100.mat');
Gi = importdata('Gi_greedy_100.mat');
Gr = importdata('Gr_greedy_100.mat');
Gs = 0:0.01:1;
Gi = 0:0.004:0.40;
Gi(101)=1;
Gr = 0:0.01:1;

lgs = length(Gs)-1;
lgi = length(Gi)-1;
beta = 0.2*7;
theta = 0.25*7;
gamma = 0.07*7;

s0 = 0.9;
i0 = 0.001;
r0 = 1-s0-i0;
transitions0 = importdata('transitions0_density_100_new.mat');
transitions1 = importdata('transitions1_density_100_new.mat');
transitions0 = importdata('transitions0_smart_100.mat');
transitions1 = importdata('transitions1_smart_100.mat');
%transitions0 = importdata('transitions0_uniform_100.mat');
%transitions1 = importdata('transitions1_uniform_100.mat');
T=10;
pol = randi(2,T,1)-1;
pol = zeros(T,1);
trj_s = zeros(1,T);
trj_e = zeros(1,T);
trj_i = zeros(1,T);
trj_r = zeros(1,T);
idx0 = find_index(s0,i0,Gs,Gi);
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

trj1 = zeros(T,6);
%pol = zeros(T,1);
[S,I,R] = SEIR_trj(s0,i0,r0,beta,gamma,T,pol);
%[S,I,R] = SEIR_trj_dis(s0,i0,r0,beta,gamma,T,pol,Gs,Gi,Gr);

trj1(:,1) = trj_s(1:T);
trj1(:,2) = trj_i(1:T);
trj1(:,3) = trj_r(1:T);

trj1(:,4) = S;
trj1(:,5) = I;
trj1(:,6) = R;
time = zeros(10,1);
for i = 1:10
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
title({'Expert Discretization'},'Fontsize',18)
xlabel('Time: week','FontSize',18)
ylabel('Proportions','FontSize',18)
ylim([0 0.9])
legend('Expert S', 'Expert I','Expert R','Actual S','Actual I','Actual R','Fontsize',14)
%% draw eval algo 2
transitions0 = importdata('transitions0_greedy_100.mat');
transitions1 = importdata('transitions1_greedy_100.mat');
T=10;
pol = randi(2,T,1)-1;
pol = zeros(T,1);
%pol(8:end)=1;
trj_s = zeros(1,T);
trj_e = zeros(1,T);
trj_i = zeros(1,T);
trj_r = zeros(1,T);
idx0 = find_index(s0,i0,Gs,Gi);
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

trj1 = zeros(T,6);
%pol = zeros(T,1);
%[S,I,R] = SEIR_trj(s0,i0,r0,beta,gamma,T,pol);
[S,I,R] = SEIR_trj_dis(s0,i0,r0,beta,gamma,T,pol,Gs,Gi,Gr);

trj1(:,1) = trj_s(1:T);
trj1(:,2) = trj_i(1:T);
trj1(:,3) = trj_r(1:T);

trj1(:,4) = S;
trj1(:,5) = I;
trj1(:,6) = R;
time = zeros(10,1);
for i = 1:10
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
title({'Markovian Trajectories V.S. Discretized Trajectories'},'Fontsize',18)
xlabel('Time: week','FontSize',18)
ylabel('Proportions','FontSize',18)
legend('Markovian S', 'Markovian I','Markovian R','Discretized S','Discretized I','Discretized R','Fontsize',14)
%%

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

%%

pol = zeros(T,1);
[S,I,R] = SEIR_trj(s0,i0,r0,beta,gamma,T,pol);
figure
trj_d = zeros(T,3);
trj_d(:,1) = S;
trj_d(:,2) = I;
trj_d(:,3) = R;
plot(trj_d)
%plot(trj1)
%%
err = 0;
for t =1:T
    err = err + (S(t)-trj_s(t))^2;
    err = err + (I(t)-trj_i(t))^2;
    err = err + (R(t)-trj_r(t))^2;
end
disp(err)
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
%%
function [S,I,R] = SEIR(s0,i0,r0,beta, gamma,action)
    delta_t = 1;
    if action == 1
        beta = beta*0.5;% lockdown is effective at reducing 70% contacts
    end
    S = s0 - beta*(s0)*i0*delta_t;
    I = i0 + beta*(s0)*i0*delta_t - gamma*i0*delta_t;
    R = r0 + gamma*i0*delta_t;
end
function [S,I,R] = SEIR_trj(s0,i0,r0,beta, gamma,T,pol)
    S = zeros(T,1);
    I = zeros(T,1);
    R = zeros(T,1);
    for t = 1:T
        S(t,1) = s0;
        I(t,1) = i0;
        R(t,1) = r0;
        action = pol(t);
        [s0,i0,r0] = SEIR(s0,i0,r0,beta, gamma,action);
    end
end
function [S,I,R] = SEIR_trj_dis(s0,i0,r0,beta, gamma,T,pol,Gs,Gi,Gr) 
    S = zeros(T,1);
    I = zeros(T,1);
    R = zeros(T,1);
    for t = 1:T
        if s0 <= 0
            es0 = (Gs(1));
        elseif s0 >= 1
            es0 = (Gs(end));
        else
            es0 = (max(Gs(Gs<s0))+min(Gs(Gs>=s0)))/2;
        end
        if i0 <= 0
            ei0 = (Gi(1));
        elseif i0 >= 1
            ei0 = (Gi(end));
        else
            ei0 = (max(Gi(Gi<i0))+min(Gi(Gi>=i0)))/2;
        end
        if r0 <= 0
            %disp(length(Gr))
            %disp(Gr(1))
            er0 = (Gr(1));
        elseif r0 >= 1
            er0 = (Gr(end));
        else
            er0 = (max(Gr(Gr<r0))+min(Gr(Gr>=r0)))/2;
        end
        %disp(es0)
        %disp(ei0)
        S(t,1) = es0;
        I(t,1) = ei0;
        R(t,1) = er0;
        action = pol(t);
        [s0,i0,r0] = SEIR(es0,ei0,er0,beta, gamma,action);
    end
end