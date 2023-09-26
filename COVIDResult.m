%This script generate discretization for SIR model
T = 40; %Time epoch
s0 = ones(26,1).*0.99; %initial state
i0 = ones(26,1).*0.0001;
r0 = ones(26,1)-s0-i0;
% Specify the path to the CSV file
file_path = 'covid-data/beta.csv';
% Load the CSV file into a MATLAB array
beta = readmatrix(file_path);
beta = beta.*0.8;
gamma = 0.7048;

Greedyres = importdata('covid-data/greedy_res.mat');
Uniformres = importdata('covid-data/uniform_res.mat');
%%
time = zeros(T,1);
for i = 1:T
    time(i)=i;    
end
pol = zeros(T,1);
pol(2:10,1) = 1;
[S,I,R] = SEIR_trj(s0,i0,r0,beta,gamma,T,pol);
[S,I,R] = compute_total_SIR(S,I,R,T);
trj = zeros(T,9);
trj(:,1) = Greedyres(:,1);
trj(:,2) = Greedyres(:,2);
trj(:,3) = Greedyres(:,3);
trj(:,4) = Uniformres(:,1);
trj(:,5) = Uniformres(:,2);
trj(:,6) = Uniformres(:,3);
trj(:,7) = S;
trj(:,8) = I;
trj(:,9) = R;
figure
plot(time,trj(:,1),'r', ...
    time,trj(:,2),'k', ...
    time,trj(:,3),'b', ...
    time,trj(:,4),'r--', ...
    time,trj(:,5),'k--', ...
    time,trj(:,6),'b--', ...
    time,trj(:,7),'r:', ...
    time,trj(:,8),'k:', ...
    time,trj(:,9),'b:', ...
    'LineWidth',3);
ax = gca;
ax.FontSize = 16; 
title({'Uniform Discretization'},'Fontsize',18)
xlabel('Time: week','FontSize',18)
ylabel('Proportions','FontSize',18)
legend('Greedy S', 'Greedy I','Greedy R', ...
    'Uniform S','Uniform I','Uniform R', ...
    'Acutal S','Actual I','Actual R','Fontsize',14)
obj = zeros(T,3);
for t = 1:T
    obj(t,1)=-trj(t,2)-0.03*Greedyres(t,4);
    obj(t,2)=-trj(t,5)-0.03*Uniformres(t,4);
    obj(t,3)=-trj(t,8)-0.03*pol(t,1);
end
disp('Greedy')
disp(sum(obj(:,1)))
disp('Uniform')
disp(sum(obj(:,2)))
disp('Actual')
disp(sum(obj(:,3)))

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