%This script generate discretization for SIR model
T = 20; %Time epoch
s0 = 0.99; %initial state
i0 = 0.001;
r0 = 1-s0-i0;
beta = 0.2*7; %initial parameter
theta = 0.25*7;
gamma = 0.07*7;
pol = randi(2,T,1)-1; %policy
%pol = zeros(T,1);
%%
% initialize discretization vector G
Gs = [0,0.5,1];
Gi = [0,0.5,1];
Gr = [0,0.5,1];
%Gs = importdata('Gs_greedy_100.mat');
%Gi = importdata('Gi_greedy_100.mat');
%Gr = importdata('Gr_greedy_100.mat');
budget = 50*3;
%this is sample state and policy file 
samples = importdata('samples_for_compare.mat');
pol_samples = importdata('policy_for_compare.mat');
%samples = samples(31:end,:);
%pol_samples = pol_samples(31:end,:);
tic
[Gs,Gi,Gr,costs] = greedy(budget,Gs,Gi,Gr,samples(1,1),samples(1,2),samples(1,3),beta,gamma,T,samples,pol_samples,pol_samples(1,1:T),S,I,R);
toc
%%
figure
plot(Gs)
figure
plot(Gi)
figure
plot(Gr)
%%      manual
%{
Gs = 0:0.001:1;
Ge = 0:0.001:1;
Gi = 0:0.001:1;
Gr = 0:0.001:1;
%}
pol = zeros(T,1);
[S,I,R] = SEIR_trj(s0,i0,r0,beta,gamma,T,pol);
[dS,dI,dR] = SEIR_trj_dis(s0,i0,r0,beta, gamma,T,pol,Gs,Gi,Gr);

%%
figure
trj = zeros(T,3);
trj(:,1) = S;
trj(:,2) = I;
trj(:,3) = R;
plot(trj)
figure
trj_d = zeros(T,3);
trj_d(:,1) = dS;
trj_d(:,2) = dI;
trj_d(:,3) = dR;
plot(trj_d)
%figure
%plot(log(costs))
%ylabel('Log Error')
%xlabel('Iterations')
%% function: greedy 
function [Gs,Gi,Gr,costs] = greedy(budget,Gs,Gi,Gr,s0,i0,r0,beta,gamma,T,samples,pol_samples,pol,S,I,R)
    cc_count = 1;    
    iter = 1;
    cutoff  = 10;
    old_cost = 100;
    costs = [];
    resample = 0; %reamsple s0,e0,i0,r0 every 200 iterations
    pol = zeros(T,1);
    while length(Gs)+length(Gi)+length(Gr)<budget
        if resample == 10
            iter = 1;
            resample = 0;
            s0= samples(cc_count+1,1);
            %e0 = 0.001 + rand()*0.099;
            i0 = samples(cc_count+1,2);
            if s0+i0>1
                s0 = s0/(s0+i0);
                i0 = i0/(s0+i0);
                r0 = 0;
            else
                r0 = 1-s0-i0;
            end
            pol = pol_samples(cc_count,1:T);
            cc_count = cc_count +1;
            [S,I,R] = SEIR_trj(s0,i0,r0,beta,gamma,T,pol);
        end
        [dS,dI,dR] = SEIR_trj_dis(s0,i0,r0,beta,gamma,T,pol,Gs,Gi,Gr);
        cost = cost_function(S,I,R,dS,dI,dR);
        costs = [costs,cost];
        worst_ = [0,1,1];%worst cut 1. value 2. dim 3 index
        best_ = [100,1,1];
        if iter >= cutoff
            disp(length(Gs)+length(Gi)+length(Gr))
            worst_ = [0,1,1];%worst cut 1. value 2. dim 3 index
            best_ = [100,1,1];
            for dim = 1:4
                if dim == 1
                    max_idx = length(Gs)-1;
                elseif dim == 2
                    max_idx = length(Gi)-1;
                else
                    max_idx = length(Gr)-1;
                end
                for index = 1:max_idx
                    [tGs,tGi,tGr] = cut(dim,index,Gs,Gi,Gr);
                    [dS,dI,dR] = SEIR_trj_dis(s0,i0,r0,beta,gamma,T,pol,tGs,tGi,tGr);
                    cost = cost_function(S,I,R,dS,dI,dR);
                    if cost<best_(1)
                        best_(1) = cost;
                        best_(2) = dim;
                        best_(3) = index;
                    end
                    if cost>worst_(1)
                        worst_(1) = cost;
                        worst_(2) = dim;
                        worst_(3) = index;
                    end
                end
            end
        end
        if iter >= cutoff
            if worst_(1)>best_(1)
                disp(worst_(1))
                disp(best_(1))
                [Gs,Gi,Gr] = cut(best_(2),best_(3),Gs,Gi,Gr);
                if abs(old_cost - best_(1))<=0.000001
                    iter = 1;
                    continue
                end
                old_cost = best_(1);
            else
                dim = randi(3);
                %randomly select a data point from the sample
                if dim == 1
                    num = S(randi(length(S)));
                    index = max(find(Gs<num));
                elseif dim == 2
                    num = I(randi(length(I)));
                    index = max(find(Gi<num));
                else
                    num = R(randi(length(R)));
                    index = max(find(Gr<num));
                end
                %index = randi(max_idx);
                [Gs,Gi,Gr] = cut(dim,index,Gs,Gi,Gr);
            end
        else
            dim = randi(3);
            %randomly select a data point from the sample
            if dim == 1
                num = S(randi(length(S)));
                index = max(find(Gs<=num));
            elseif dim == 2
                num = I(randi(length(I)));
                index = max(find(Gi<=num));
            else
                num = R(randi(length(R)));
                index = max(find(Gr<=num));
            end
            %index = randi(max_idx);
            
            [Gs,Gi,Gr] = cut(dim,index,Gs,Gi,Gr);
        end
        resample = resample + 1;
        iter = iter+1;
    end
    disp('xxxxx')
    disp(cc_count)
end
%% function: cut 
function [tGs,tGi,tGr] = cut(dim,index,Gs,Gi,Gr)
    if dim == 1
        tGs = [Gs(1:index),(Gs(index)+Gs(index+1))/2, Gs(index+1:end)];
        tGi = Gi;
        tGr = Gr;
    elseif dim == 2
        tGi = [Gi(1:index),(Gi(index)+Gi(index+1))/2, Gi(index+1:end)];
        tGs = Gs;
        tGr = Gr;
    else
        %disp('bug')
        %disp(index)
        %disp('above')
        tGr = [Gr(1:index),(Gr(index)+Gr(index+1))/2, Gr(index+1:end)];
        tGi = Gi;
        tGs = Gs;
    end
end

%% function: cost
function cost = cost_function(S,I,R,S1,I1,R1)
    cost = sum((S1-S).^2)+sum((I1-I).^2)+sum((R1-R).^2);
end
%% functions: Discretized Discrete Time SEIR Model
%G is the discretization matrix
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
        S(t,1) = es0;
        I(t,1) = ei0;
        R(t,1) = er0;
        action = pol(t);
        [s0,i0,r0] = SEIR(es0,ei0,er0,beta, gamma,action);
    end
end
%% functions: Discrete Time SIR Model
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

function [S,I,R] = SEIR(s0,i0,r0,beta, gamma,action)
    delta_t = 1;
    if action == 1
        beta = beta*0.5;% lockdown is effective at reducing 30% contacts
    end
    S = s0 - beta*(s0)*i0*delta_t;
    I = i0 + beta*(s0)*i0*delta_t - gamma*i0*delta_t;
    R = r0 + gamma*i0*delta_t;
end