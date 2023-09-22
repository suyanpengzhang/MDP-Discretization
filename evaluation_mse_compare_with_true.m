%% generate evaluation samples
%{
samples_eval =zeros(1000,3);
policy_eval = zeros(1000,20);
for idx = 1:1000
    pol = randi(2,20,1)-1;
    policy_eval(idx,:) = pol;
    s0 = 0.8+0.199*rand(); %draw from 0.7-0.999
    i0 = 0.001+0.05*rand(); %draw from 0.001 - 0.1
    if s0+i0<1
        samples_eval(idx,1) = s0;
        samples_eval(idx,2) = i0;
        samples_eval(idx,3) = 1 - s0 - i0;
    else
        samples_eval(idx,1) = s0/(s0+i0);
        samples_eval(idx,2) = i0/(s0+i0);
        samples_eval(idx,3) = 0;
    end
end
%}
%% compute error (compare with true)
% Define the file names
t_score = tinv(0.975, 99);
T = 20;
beta = 0.2*7;
theta = 0.25*7;
gamma = 0.07*7;
%policy_eval = importdata('policy_for_compare.mat');
%samples_eval = importdata('samples_for_compare.mat');
policy_eval = importdata('policy_eval.mat');
samples_eval = importdata('samples_eval.mat');
files_Gs = {'Gs_greedy_30.mat', 'Gs_greedy_50.mat', 'Gs_greedy_100.mat','Gs_greedy_400_new.mat'};
files_Gi = {'Gi_greedy_30.mat', 'Gi_greedy_50.mat', 'Gi_greedy_100.mat','Gi_greedy_400_new.mat'};
files_Gr = {'Gr_greedy_30.mat', 'Gr_greedy_50.mat', 'Gr_greedy_100.mat','Gr_greedy_400_new.mat'};
files_trans0_g = {'transitions0_greedy_30.mat', 'transitions0_greedy_50.mat', 'transitions0_greedy_100.mat','transition0_greedy_400_new.mat'};
files_trans1_g = {'transitions1_greedy_30.mat', 'transitions1_greedy_50.mat', 'transitions1_greedy_100.mat','transition1_greedy_400_new.mat'};
files_trans0_u = {'transitions0_uniform_30.mat', 'transitions0_uniform_50.mat', 'transitions0_uniform_100.mat','transitions0_uniform_400.mat'};
files_trans1_u = {'transitions1_uniform_30.mat', 'transitions1_uniform_50.mat', 'transitions1_uniform_100.mat','transitions1_uniform_400.mat'};
% Loop to load the files
Greedy_error = zeros(1,4);
Greedy_error_dis = zeros(1,4);
Uniform_error = zeros(1,4);
Uniform_error_dis = zeros(1,4);
num_samples = 100;
for ifl = 2:numel(files_Gs)
    file_name = files_Gs{ifl};
    Gs = importdata(file_name);
    file_name = files_Gi{ifl};
    Gi = importdata(file_name);
    file_name = files_Gr{ifl};
    Gr = importdata(file_name);
    file_name = files_trans0_g{ifl};
    transitions0 = importdata(file_name);
    file_name = files_trans1_g{ifl};
    transitions1 = importdata(file_name);
    disp('**********')
    disp(ifl)
    lgs = length(Gs)-1;
    lgi = length(Gi)-1;
    errors = zeros(num_samples,1);
    errors_dis = zeros(num_samples,1);
    for iiii = 1:num_samples
        %disp(iiii)
        pol = policy_eval(iiii,:);
        s0 = samples_eval(iiii,1);
        i0 = samples_eval(iiii,2);
        r0 = samples_eval(iiii,3);
        trj_s = zeros(1,T);
        trj_i = zeros(1,T);
        trj_r = zeros(1,T);
        idx0 = find_index(s0,i0,Gs,Gi);
        b0 = sparse(1,lgs*lgi);
        b0(1,idx0) = 1;
        T = 20;
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
        errors_dis(iiii,1) = err;
        [S,I,R] = SEIR_trj(s0,i0,r0,beta,gamma,T,pol);
        err = 0;
        for t =1:T
            err = err + (S(t)-trj_s(t))^2;
            err = err + (I(t)-trj_i(t))^2;
            err = err + (R(t)-trj_r(t))^2;
        end
        errors(iiii,1) = err;
        %{
        if ifl>2
            figure
            trj_d = zeros(T,6);
            trj_d(:,1) = S;
            trj_d(:,2) = I;
            trj_d(:,3) = R;
            trj_d(:,4) = trj_s;
            trj_d(:,5) = trj_i;
            trj_d(:,6) = trj_r;
            plot(trj_d)
        end
        %}
    end
    Greedy_error_dis(1,ifl) = mean(errors_dis);
    Greedy_error(1,ifl) = mean(errors);
    ci_dis_low = mean(errors_dis)-(std(errors_dis) / sqrt(100))*t_score;
    ci_dis_high = mean(errors_dis)+(std(errors_dis) / sqrt(100))*t_score;
    disp('Greedy_dis CI')
    disp(ci_dis_low)
    disp(ci_dis_high)
    ci_low = mean(errors)-(std(errors) / sqrt(100))*t_score;
    ci_high = mean(errors)+(std(errors) / sqrt(100))*t_score;
    disp('Greedy CI')
    disp(ci_low)
    disp(ci_high)
    if ifl ==1 
        Gs = 0:0.03333:1;
        Gi= 0:0.03333:1;
        Gr = 0:0.03333:1;
    elseif ifl ==2
        Gs = 0:0.02:1;
        Gi= 0:0.02:1;
        Gr = 0:0.02:1;
    elseif ifl==3
        Gs = 0:0.01:1;
        Gi= 0:0.01:1;
        Gr = 0:0.01:1;
    else
        Gs = 0:0.0025:1;
        Gi= 0:0.0025:1;
        Gr = 0:0.0025:1;
    end
    file_name = files_trans0_u{ifl};
    transitions0 = importdata(file_name);
    file_name = files_trans1_u{ifl};
    transitions1 = importdata(file_name);
    disp('**********')
    disp(ifl)
    lgs = length(Gs)-1;
    lgi = length(Gi)-1;
    errors = zeros(num_samples,1);
    errors_dis = zeros(num_samples,1);
    for iiii = 1:num_samples
        %disp(iiii)
        policy_eval(iiii,:);
        s0 = samples_eval(iiii,1);
        i0 = samples_eval(iiii,2);
        r0 = samples_eval(iiii,3);
        trj_s = zeros(1,T);
        trj_i = zeros(1,T);
        trj_r = zeros(1,T);
        idx0 = find_index(s0,i0,Gs,Gi);
        b0 = sparse(1,lgs*lgi);
        b0(1,idx0) = 1;
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
        errors_dis(iiii,1) = err;
        [S,I,R] = SEIR_trj(s0,i0,r0,beta,gamma,T,pol);
        err = 0;
        for t =1:T
            err = err + (S(t)-trj_s(t))^2;
            err = err + (I(t)-trj_i(t))^2;
            err = err + (R(t)-trj_r(t))^2;
        end
        errors(iiii,1) = err;
    end
    Uniform_error_dis(1,ifl) = mean(errors_dis);
    Uniform_error(1,ifl) = mean(errors);
    ci_dis_low = mean(errors_dis)-(std(errors_dis) / sqrt(100))*t_score;
    ci_dis_high = mean(errors_dis)+(std(errors_dis) / sqrt(100))*t_score;
    disp('Uniform_dis CI')
    disp(ci_dis_low)
    disp(ci_dis_high)
    ci_low = mean(errors)-(std(errors) / sqrt(100))*t_score;
    ci_high = mean(errors)+(std(errors) / sqrt(100))*t_score;
    disp('Uniform CI')
    disp(ci_low)
    disp(ci_high)
end



%% functions

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

function idx = find_index(s,i,Gs,Gi)
    lgs = length(Gs)-1;
    lgi = length(Gi)-1;
    ls = max(find(Gs<s));
    li = max(find(Gi<i));
    idx = (ls-1)*(lgi)+li;
end

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