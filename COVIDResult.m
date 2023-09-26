g%This script generate discretization for SIR model
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
pol0 = zeros(T,1);
[S0,I0,R0] = SEIR_trj(s0,i0,r0,beta,gamma,T,pol0);
[S0,I0,R0] = compute_total_SIR(S0,I0,R0,T);
trj = zeros(T,12);
trj(:,1) = Greedyres(:,1);
trj(:,2) = Greedyres(:,2);
trj(:,3) = Greedyres(:,3);
trj(:,4) = Uniformres(:,1);
trj(:,5) = Uniformres(:,2);
trj(:,6) = Uniformres(:,3);
trj(:,7) = S;
trj(:,8) = I;
trj(:,9) = R;
trj(:,10) = S0;
trj(:,11) = I0;
trj(:,12) = R0;
figure
plot(time,trj(:,2),'r', ...
    time,trj(:,5),'b--', ...
    time,trj(:,8),'g:', ...
    time,trj(:,11),'k-.', ...
    'LineWidth',3);
ax = gca;
ax.FontSize = 16; 
title({'Proportions of Infectious over Time'},'Fontsize',18)
xlabel('Time: week','FontSize',18)
ylabel('Proportions','FontSize',18)
legend('Greedy policy', ...
    'Uniform policy', ...
    'Acutal policy', ...
    'No intervention','Fontsize',14)
figure
plot(time,trj(:,1),'r', ...
    time,trj(:,4),'b--', ...
    time,trj(:,7),'g:', ...
    time,trj(:,10),'k-.', ...
    'LineWidth',3);
ax = gca;
ax.FontSize = 16; 
title({'Proportions of Susceptible over Time'},'Fontsize',18)
xlabel('Time: week','FontSize',18)
ylabel('Proportions','FontSize',18)
legend('Greedy policy', ...
    'Uniform policy', ...
    'Acutal policy', ...
    'No intervention','Fontsize',14)
obj = zeros(T,3);
for t = 1:T
    obj(t,1)=-trj(t,2)-0.005*Greedyres(t,4);
    obj(t,2)=-trj(t,5)-0.005*Uniformres(t,4);
    obj(t,3)=-trj(t,8)-0.005*pol(t,1);
    obj(t,4)=-trj(t,11)-0.005*pol0(t,1);
end
disp('Greedy')
disp(sum(obj(:,1)))
disp('Uniform')
disp(sum(obj(:,2)))
disp('Actual')
disp(sum(obj(:,3)))
disp('No intervention')
disp(sum(obj(:,4)))
figure
x = categorical({'Greedy' 'Uniform' 'Actual' 'Do nothing'});
x = reordercats(x,{'Greedy' 'Uniform' 'Actual' 'Do nothing'});
bar(x,sum(obj,1))
text(1:length(sum(obj,1)),sum(obj,1)-0.05,num2str(round(sum(obj,1),4)'),'vert','bottom','horiz','center','FontSize',16); 

ax = gca;
ax.FontSize = 16; 
title({'Objective Across Different Models'},'Fontsize',18)
xlabel('Policy','FontSize',18)
ylabel('Objective value','FontSize',18)
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