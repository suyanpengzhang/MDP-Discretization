p0 = importdata('transitions0_greedy_400.mat');
p1 = importdata('transitions1_greedy_400.mat');
%p0 = p0*p0*p0*p0*p0*p0*p0;
%p1 = p1*p1*p1*p1*p1*p1*p1;
%p0u = importdata('transitions0_uniform.mat');
%p1u = importdata('transitions1_uniform.mat');
%p0u = p0u*p0u*p0u*p0u*p0u*p0p1u = p1u*p1u*p1u*p1u*p1u*p1u*p1u;
Gs = importdata('Gs_greedy_400.mat');
Gi = importdata('Gi_greedy_400.mat');
Gr = importdata('Gr_greedy_400.mat');
%Gs = 0:0.03333:1;
%Gi = 0:0.0025:1;
%Gi = 0:0.013333:0.40;
%Gi(31)=1;
%Gr = 0:0.03333:1;
%Gsu = 0:0.005:1;
%Giu = 0:0.005:1;
%Gru = 0:0.005:1;
beta = 0.2*7;
theta = 0.25*7;
gamma = 0.07*7;
P{1} = p0;
P{2} = p1;
lgs = length(Gs)-1;
lgi = length(Gi)-1;
R = ones(length(p0),2);
costr = 0.03;
for bs = 1:lgs
    for bi = 1:lgi
        idx1 = (bs-1)*lgi+bi;
        R(idx1,1) = -(Gi(bi)+Gi(bi+1))/2;
        R(idx1,2) = -(Gi(bi)+Gi(bi+1))/2-costr;
    end
end
skip = 0;
[V, policy, cpu_time] = mdp_finite_horizon(P, R, 1, 10);
%%

sv = 0.7:0.01:0.99;
iv = 0.001:0.001:0.01;
cdata = zeros(length(iv),length(sv));
t = 5;
for sidx = 1:length(sv)
    for iidx = 1:length(iv)
        %disp(sv(sidx))
        %disp(iv(iidx))
        vidx = find_index(sv(sidx),iv(iidx),Gs,Gi);
        polidx = policy(vidx,t);
        cdata(iidx,sidx) = polidx-1;
    end
end
h = heatmap(sv,iv,cdata);
sv_greedy = sv;
iv_greedy = iv;
cdata_greedy = cdata;
h.Title = 'Week 5 -- Greedy';
h.XLabel = 'Susceptible Proportions';
h.YLabel = 'Infectious Proportions';
h.NodeChildren(3).YDir='normal'; 
%%
%{
samples = importdata('samples_for_compare.mat');
pol_samples = importdata('policy_for_compare.mat');
idx = find_index(samples(1,1),samples(1,2),Gs,Gi);
Queue = [idx];
for samp = 1:length(samples)
    idx = find_index(samples(samp,1),samples(samp,2),Gs,Gi);
    if ismember(idx,Queue)
        continue
    else
        Queue = [Queue,idx];
    end
end
visited = [];

iter = 0;
while length(Queue)>0
    data = Queue(1);
    visited = unique([Queue(1),visited]);
    Queue(1) = [];
    prob_data = [find(p0(data,:)),find(p1(data,:))];
    for i=1:length(prob_data)
        if ismember(prob_data(i),visited)
            continue
        elseif ismember(prob_data(i),Queue)
            continue
        else
            if iter <100
                Queue = [Queue,prob_data(i)];
            end
        end
    end
    iter = iter + 1;
end
%}
sv = 0.7:0.005:0.99;
iv = 0.001:0.0005:0.01;
count = 1;
sampless = zeros(length(sv)*length(iv),2);
for is = 1:length(sv)
    for idv = 1:length(iv)
        sampless(count,1) = sv(is);
        sampless(count,2) = iv(idv);
        count = count + 1;
    end
end
 
%%
num_samples = length(sampless);
optimal_policy_brute_force = zeros(num_samples,10);
samples = zeros(num_samples,10);
%samples(1,1) = 0.99;
%samples(1,2) = 0.005;
%samples(1,3) = 0.005;
for s =1:num_samples
    %{
    rs = 0.7+0.3*rand();
    ri = 0.1*rand();
    if rs+ri<=1
        rr = 1 - rs - ri;
    else
        rs = rs/(rs+ri);
        ri = ri/(rs+ri);
        rr = 0;
    end
    %}
    %[rs,ri,rr] = reverse_find_index(visited(s),Gs,Gi);
    rs = sampless(s,1);
    ri = sampless(s,2);
    if rs+ri>1
        rs = rs/(rs+ri);
        ri = ri/(rs+ri);
        rr=0;
    else
        rr = 1-rs-ri;
    end
    %{
    if s == 1
        rs = 0.99;
        ri = 0.005;
        rr = 0.005;
    end
    %}
    samples(s,1) = rs;
    samples(s,2) = ri;
    samples(s,3) = rr;
    idx = find_index(rs,ri,Gs,Gi);
    samples(s,4) = policy(idx,1);
    samples(s,5) = V(idx,1); %v1
    mdp_action = policy(idx,:);
    %optimal_policy_brute_force(s,:) = optimal_action;
    %disp(optimal_action)
    %disp(mdp_action)
    mis_match = 0;
    for id_find = 1:length(mdp_action)-skip
        metric = evaluate_brute_force(id_find,rs,ri,beta,gamma,costr);
        optimal_action = metric(find(metric(:,11) == max(metric(:,11))),1);
        if mdp_action(id_find) ~= optimal_action(1)
            mis_match = mis_match + 1;
        end
        if id_find == 1
            samples(s,7) = max(metric(:,11));%optimal value day 1 brute force
        end
    end
    samples(s,10) = mis_match;
    %samples(s,6) = metric(find(metric(:,11) == max(metric(:,11))),1);
    %samples(s,7) = max(metric(:,11));
    samples(s,8) = samples(s,6)-samples(s,4);
    s0 = samples(s,1);
    i0 = samples(s,2);
    r0 = samples(s,3);
    cost = 0;
    for time = 1:10
        cost = cost -i0;
        if policy(find_index(s0,i0,Gs,Gi),time)==2
            cost = cost -costr;
        end
        [s0,i0,r0] = SEIR(s0,i0,r0,beta, gamma,policy(find_index(s0,i0,Gs,Gi),time)-1);
    end
    samples(s,9) = abs(cost-samples(s,7))/abs(samples(s,7));
    %disp(cost)
end
%disp(length(find(samples(:,8)~=0))/10000) %accuarcy at first day
%disp(median(abs(samples(:,5)-samples(:,7))./(abs(samples(:,7))))) %approxmiation error
disp('Greedy')
disp('Accuracy All')
disp(1-sum(samples(:,10))/(num_samples*10))
disp('MSE')
disp(mean((samples(:,5)-samples(:,7)).^2))
disp('Mean Approx Error')
disp(mean(abs(samples(:,5)-samples(:,7))./(abs(samples(:,7))))) %approxmiation error
disp('Median Approx Error')
disp(median(abs(samples(:,5)-samples(:,7))./(abs(samples(:,7))))) %approxmiation error
disp('Mean Opt Gap')
disp(mean(samples(:,9)))
disp('Median Opt Gap')
disp(median(samples(:,9)))


%%
function metric = evaluate_brute_force(T,rs,ri,beta,gamma,costr)
    T = 11-T; 
    %disp(T)
    metric = zeros(2^10,11);
    count = 1;
    for i1 = 1:2
        for i2 = 1:2
            for i3 = 1:2
                for i4 = 1:2
                    for i5 = 1:2
                        for i6 = 1:2
                            for i7 = 1:2
                                for i8 = 1:2
                                    for i9 = 1:2
                                        for i10 = 1:2
                                            s0 = rs;
                                            i0 = ri;
                                            r0 = 1-s0-i0;
                                            metric(count,1) = i1;
                                            metric(count,2) = i2;
                                            metric(count,3) = i3;
                                            metric(count,4) = i4;
                                            metric(count,5) = i5;
                                            metric(count,6) = i6;
                                            metric(count,7) = i7;
                                            metric(count,8) = i8;
                                            metric(count,9) = i9;
                                            metric(count,10) = i10;         
                                            cost = 0;
                                            if i1 == 1
                                                cost = cost -i0;
                                                [s0,i0,r0] = SEIR(s0,i0,r0,beta, gamma,0);
                                            else
                                                cost = cost -i0-costr;
                                                [s0,i0,r0] = SEIR(s0,i0,r0,beta, gamma,1);
                                            end
                                            if i2 == 1 && T>1
                                                cost = cost -i0;
                                                [s0,i0,r0] = SEIR(s0,i0,r0,beta, gamma,0);
                                            elseif i2 == 2 && T>1
                                                cost = cost -i0-costr;
                                                [s0,i0,r0] = SEIR(s0,i0,r0,beta, gamma,1);
                                            end
                                            if i3 == 1 && T>2
                                                cost = cost -i0;
                                                [s0,i0,r0] = SEIR(s0,i0,r0,beta, gamma,0);
                                            elseif i3 == 2 && T>2
                                                cost = cost -i0-costr;
                                                [s0,i0,r0] = SEIR(s0,i0,r0,beta, gamma,1);
                                            end
                                            if i4 == 1 && T>3
                                                cost = cost -i0;
                                                [s0,i0,r0] = SEIR(s0,i0,r0,beta, gamma,0);
                                            elseif i4 == 2 && T>3
                                                cost = cost -i0-costr;
                                                [s0,i0,r0] = SEIR(s0,i0,r0,beta, gamma,1);
                                            end
                                            if i5 == 1 && T>4
                                                cost = cost -i0;
                                                [s0,i0,r0] = SEIR(s0,i0,r0,beta, gamma,0);
                                            elseif i5 == 2 && T>4
                                                cost = cost -i0-costr;
                                                [s0,i0,r0] = SEIR(s0,i0,r0,beta, gamma,1);
                                            end
                                            if i6 == 1 && T>5
                                                cost = cost -i0;
                                                [s0,i0,r0] = SEIR(s0,i0,r0,beta, gamma,0);
                                            elseif i6 == 2 && T>5
                                                cost = cost -i0-costr;
                                                [s0,i0,r0] = SEIR(s0,i0,r0,beta, gamma,1);
                                            end
                                            if i7 == 1 && T>6
                                                cost = cost -i0;
                                                [s0,i0,r0] = SEIR(s0,i0,r0,beta, gamma,0);
                                            elseif i7 == 2 && T>6
                                                cost = cost -i0-costr;
                                                [s0,i0,r0] = SEIR(s0,i0,r0,beta, gamma,1);
                                            end
                                            if i8 == 1 && T>7
                                                cost = cost -i0;
                                                [s0,i0,r0] = SEIR(s0,i0,r0,beta, gamma,0);
                                            elseif i8 == 2 && T>7
                                                cost = cost -i0-costr;
                                                [s0,i0,r0] = SEIR(s0,i0,r0,beta, gamma,1);
                                            end
                                            if i9 == 1 && T>8
                                                cost = cost -i0;
                                                [s0,i0,r0] = SEIR(s0,i0,r0,beta, gamma,0);
                                            elseif i9 == 2 && T>8
                                                cost = cost -i0-costr;
                                                [s0,i0,r0] = SEIR(s0,i0,r0,beta, gamma,1);
                                            end
                                            if i10 == 1 && T>9
                                                cost = cost -i0;
                                                [s0,i0,r0] = SEIR(s0,i0,r0,beta, gamma,0);
                                            elseif i10 == 2 && T>9
                                                cost = cost -i0-costr;
                                                [s0,i0,r0] = SEIR(s0,i0,r0,beta, gamma,1);
                                            end                                       
                                            metric(count,11) = cost;
                                            count = count +1;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
      
%%
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