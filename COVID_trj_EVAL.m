p0 = importdata('covid-data/COVIDtransitions0_Greedy_50.mat');
p1 = importdata('covid-data/COVIDtransitions1_Greedy_50.mat');
Gs = importdata('covid-data/COVIDGs_Greedy_50.mat');
Gi = importdata('covid-data/COVIDGi_Greedy_50.mat');
Gr = importdata('covid-data/COVIDGr_Greedy_50.mat');

lgs = length(Gs)-1;
lgi = length(Gi)-1;


T=60;
file_path = 'covid-data/beta.csv';
% Load the CSV file into a MATLAB array
beta = readmatrix(file_path);
gamma = 0.7048;
P{1} = p0;
P{2} = p1;
time = zeros(T,1);
for i = 1:T
    time(i)=i;    
end
pol = zeros(T,1);
pol = Greedyres(:,4);
%pol(7:48,1)=1;
%pol=randi(2,T,1)-1;
s0 = 0.9999*ones(26,1);
i0 = 0.0001*ones(26,1);
r0 = 0*ones(26,1);
S = s0;
I = i0;
R = r0;
trj = zeros(T,9);


b0 = sparse(1,lgs*lgi);
idx0 = find_index(S(1,1),I(1,1),Gs,Gi);
b0(1,idx0) = 1;
for t = 1:T
    disp(t)
    [S,I,R] = SEIR(S,I,R,beta,gamma,pol(t));
    [x,y,z] = compute_total_SIR(S,I,R);
    trj(t,1) = x;
    trj(t,2) = y;
    trj(t,3) = z;
    [Sg,Ig,Rg,b0] = SEIR_markovian(b0,pol(t),p0,p1,Gs,Gi);
    trj(t,4) = Sg;
    trj(t,5) = Ig;
    trj(t,6) = Rg;
end
p0 = importdata('covid-data/COVIDtransitions0_uniform_50.mat');
p1 = importdata('covid-data/COVIDtransitions1_uniform_50.mat');
Gs = 0:0.02:1;
Gi = 0:0.02:1;
Gr = 0:0.02:1;


lgs = length(Gs)-1;
lgi = length(Gi)-1;
b0 = sparse(1,lgs*lgi);
s0 = 0.9999*ones(26,1);
i0 = 0.0001*ones(26,1);
r0 = 0*ones(26,1);
S = s0;
I = i0;
R = r0;
idx0 = find_index(S(1,1),I(1,1),Gs,Gi);
b0(1,idx0) = 1;
for t=1:T
    disp(t)
    [Su,Iu,Ru,b0] = SEIR_markovian(b0,pol(t),p0,p1,Gs,Gi);
    trj(t,7) = Su;
    trj(t,8) = Iu;
    trj(t,9) = Ru;
end

figure
plot(time,trj(:,1),'r', ...
    time,trj(:,2),'b', ...
    time,trj(:,3),'g', ...
    time,trj(:,4),'r:', ...
    time,trj(:,5),'b:', ...
    time,trj(:,6),'g:', ...
    time,trj(:,7),'r--', ...
    time,trj(:,8),'b--', ...
    time,trj(:,9),'g--', ...
    'LineWidth',3);
ax = gca;
ax.FontSize = 16; 
title({'Proportions of Susceptible over Time'},'Fontsize',18)
xlabel('Time: week','FontSize',18)
ylabel('Proportions','FontSize',18)
legend('S', ...
    'I', ...
    'R', ...
    'S-Markov', ...
    'I-Markov', ...
    'R-Markov','Fontsize',14)
%%
disp(sum((trj(:,1)-trj(:,4)).^2)+sum((trj(:,2)-trj(:,5)).^2)+sum((trj(:,3)-trj(:,6)).^2))
disp(sum((trj(:,1)-trj(:,7)).^2)+sum((trj(:,2)-trj(:,8)).^2)+sum((trj(:,3)-trj(:,9)).^2))
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
    li = max(find(Gi<=i));
    idx = (ls-1)*(lgi)+li;
end

function [S,I,R,b0] = SEIR_markovian(b0,action,t0,t1,Gs,Gi)
    if action==0
                b0(1,:) = b0(1,:)*t0;
            else
                b0(1,:) = b0(1,:)*t1;
    end
    nonzeros_idx = find(b0(1,:));
    values = nonzeros(b0(1,:));
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
    S = round(ss,7);
    I = round(ii,7);
    R = round(rr,7);
end

function [S,I,R] = SEIR(s0,i0,r0,beta, gamma,action)
    delta_t = 1;
    S = zeros(26,1);
    I = zeros(26,1);
    R = zeros(26,1);
    if action == 1
        beta = beta.*0.2;% lockdown is effective at reducing 80% contacts
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