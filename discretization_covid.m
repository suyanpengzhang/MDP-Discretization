%This script generate discretization for SIR model
T = 60; %Time epoch
s0 = ones(26,1).*0.9999; %initial state
i0 = ones(26,1).*0.0001;
r0 = ones(26,1)-s0-i0;
% Specify the path to the CSV file
file_path = 'covid-data/beta.csv';
% Load the CSV file into a MATLAB array
beta = readmatrix(file_path);
%beta = beta.*0.8;
gamma = 0.7048;
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
pol_samples = randi(2,150,T)-1;
pol_samples(1,:) = zeros(1,T);
%for i=1:150
%    pol_samples(i,:) = zeros(1,T);
%end
%%
[s,i,r] = SEIR_trj(s0,i0,r0,beta, gamma,T,zeros(1,T));
%%

tic
[Gs,Gi,Gr,costs] = greedy(budget,Gs,Gi,Gr,samples(1,1),samples(1,2),samples(1,3),beta,gamma,T,samples,pol_samples,pol_samples(1,1:T));
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
%pol(1:8,1)=1;
%pol(1:20,1)=1;
[S,I,R] = SEIR_trj(s0,i0,r0,beta,gamma,T,pol);
[dS,dI,dR] = SEIR_trj_dis(s0,i0,r0,beta, gamma,T,pol,Gs,Gi,Gr);
[S,I,R] = compute_total_SIR(S,I,R,T);
[dS,dI,dR] = compute_total_SIR(dS,dI,dR,T);
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
disp(sum(I))
%figure
%plot(log(costs))
%ylabel('Log Error')
%xlabel('Iterations')
%% function: greedy 
function [Gs,Gi,Gr,costs] = greedy(budget,Gs,Gi,Gr,s0,i0,r0,beta,gamma,T,samples,pol_samples,pol)
    cc_count = 1;    
    iter = 1;
    cutoff  = 2;
    old_cost = 100;
    costs = [];
    resample = 0; %reamsple s0,e0,i0,r0 every 10 iterations
    pol = zeros(T,1);
    [S,I,R] = SEIR_trj(ones(26,1).*s0,ones(26,1).*i0,ones(26,1).*r0,beta,gamma,T,pol);
    while length(Gs)+length(Gi)+length(Gr)<budget
        if resample == 10
            iter = 1;
            resample = 0;
            s0= samples(cc_count+1,1);
            i0 = samples(cc_count+1,2);
            %s0 = 0.999;
            %i0 = 0.0001;
            if s0+i0>1
                s0 = s0/(s0+i0);
                i0 = i0/(s0+i0);
                r0 = 0;
            else
                r0 = 1-s0-i0;
            end
            pol = pol_samples(cc_count,1:T);
            cc_count = cc_count +1;
            [S,I,R] = SEIR_trj(ones(26,1).*s0,ones(26,1).*i0,ones(26,1).*r0,beta,gamma,T,pol);
        end
        [dS,dI,dR] = SEIR_trj_dis(ones(26,1).*s0,ones(26,1).*i0,ones(26,1).*r0,beta,gamma,T,pol,Gs,Gi,Gr);
        [S,I,R] = compute_total_SIR(S,I,R,T);
        [dS,dI,dR] = compute_total_SIR(dS,dI,dR,T);
        cost = cost_function(S,I,R,dS,dI,dR);
        costs = [costs,cost];
        worst_ = [0,1,1];%worst cut 1. value 2. dim 3 index
        best_ = [100,1,1];
        if iter >= cutoff
            disp(length(Gs)+length(Gi)+length(Gr))
            worst_ = [0,1,1];%worst cut 1. value 2. dim 3 index
            best_ = [100,1,1];
            for dim = 1:3
                if dim == 1
                    max_idx = length(Gs)-1;
                elseif dim == 2
                    max_idx = length(Gi)-1;
                else
                    max_idx = length(Gr)-1;
                end
                for index = 1:max_idx
                    [tGs,tGi,tGr] = cut(dim,index,Gs,Gi,Gr);
                    [dS,dI,dR] = SEIR_trj_dis(ones(26,1).*s0,ones(26,1).*i0,ones(26,1).*r0,beta,gamma,T,pol,tGs,tGi,tGr);
                    [dS,dI,dR] = compute_total_SIR(dS,dI,dR,T);
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
        if iter >= cutoff && abs(old_cost - best_(1))>0
            if worst_(1)>best_(1)
                disp(worst_(1))
                disp(best_(1))
                [Gs,Gi,Gr] = cut(best_(2),best_(3),Gs,Gi,Gr);
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
    for t = 2:size(S,1)
        if (S(t,1)-S(t+1,1))/S(t,1) <0.001
            break
        end
    end
    %cost = sum((S1(1:t,1)-S(1:t,1)).^2)+sum((I1(1:t,1)-I(1:t,1)).^2)+sum((R1(1:t,1)-R(1:t,1)).^2);
    %T = 20;
    cost = sum((S1(1:t,1)-S(1:t,1)).^2)+sum((I1(1:t,1)-I(1:t,1)).^2)+sum((R1(1:t,1)-R(1:t,1)).^2);

end
%% functions: Discretized Discrete Time SEIR Model
%G is the discretization matrix
function [S,I,R] = SEIR_trj_dis(s0,i0,r0,beta, gamma,T,pol,Gs,Gi,Gr) 
    S = zeros(T,26);
    I = zeros(T,26);
    R = zeros(T,26);
    for t = 1:T
        %begin modify
        [s0,i0,r0] = compute_total_SIR_single(s0,i0,r0);
        s0 = s0.*ones(26,1);
        i0 = i0.*ones(26,1);
        r0 = r0.*ones(26,1);
        %end modify
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
        %begin modify
        %[s0,i0,r0] = compute_total_SIR_single(s0,i0,r0);
        %s0 = s0.*ones(26,1);
        %i0 = i0.*ones(26,1);
        %r0 = r0.*ones(26,1);
        %end modify
    end
end
function [S_,I_,R_] = compute_total_SIR_single(S,I,R)
    totalpop = [420697,443569,344450,526877,899111,339399,419797,308499,140361,547523,354750,479505,287613,666399,278815,193899,166374,379199,356465,195082,407864,321720,201739,411617,469439,465691];
    S_ = sum(transpose(totalpop).*S)/sum(totalpop);
    I_ = sum(transpose(totalpop).*I)/sum(totalpop);
    R_ = sum(transpose(totalpop).*R)/sum(totalpop);
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

function [S_,I_,R_] = compute_total_SIR(S,I,R,T)
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