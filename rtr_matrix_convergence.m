Gs = importdata('Gs_greedy_30.mat');
Gi = importdata('Gi_greedy_30.mat');
Gr = importdata('Gr_greedy_30.mat');
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
xi = [0;10;30;70;100;300;500;800;1000;1500;2000;5000;10000];
convergence = zeros(length(xi),1);
for num_ = 1:length(xi)
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
            for iter = 1:xi(num_)
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
    if num_>0
        convergence(num_,1) = sum((transitions00-transitions).^2,'all');
    end
    transitions0 = transitions;
end
%transitions00 is transitions with 10000 observations


%%

figure
plot(xi,log(convergence),'r', 'LineWidth',3);
ax = gca;
ax.FontSize = 16; 
xlabel('Iterations','FontSize',18)
ylabel({'log(\Sigma_{i,j} (P^{iter}(j|i)-P^{10000}(j|i))^2)'},'FontSize',18)

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