%This script generate discretization for SIR model
T = 60; %Time epoch
s0 = ones(26,1).*0.9999; %initial state
i0 = ones(26,1).*0.0001;
r0 = ones(26,1)-s0-i0;
% Specify the path to the CSV file
file_path = 'covid-data/beta.csv';
% Load the CSV file into a MATLAB array
beta = 1.3.*readmatrix(file_path);
gamma = 0.7048;
costr = 0.005;
Greedyres1 = importdata('covid-data/greedy_resconstrained.mat');
Densityres1 = importdata('covid-data/density_resconstrained.mat');
Uniformres1 = importdata('covid-data/uniform_resconstrained.mat');
Greedyres = importdata('covid-data/greedy_res_30more.mat');
Uniformres = importdata('covid-data/uniform_res_30more.mat');
Densityres = importdata('covid-data/density_res_30more.mat');
%%
cdata = zeros(T,7);
cdata(3:10,1) = 1;
cdata(:,3) = Densityres(:,4);
cdata(:,4) = Uniformres(:,4);
cdata(:,2) = Greedyres(:,4);
cdata(:,7) = Uniformres1(:,4);
cdata(:,6) = Densityres1(:,4);
cdata(:,5) = Greedyres1(:,4);
cdata = transpose(cdata);
xvalues = {'Empirical','GreedyCut','Simulation-Based','Uniform','GreedyCut -- 2 switch','Simulation-Based -- 2 switch','Uniform -- 2 switch'};
yvalues = cell(1,T);
yvalues(:,:)={0};
ytik = cell(1,T);
ytik(:,:)={0};
for i =1:T
    yvalues(1,i)={i};
    if mod(i,5)==0
        ytik(1,i) = {num2str(i)};
    else
        ytik(1,i) = {' '} ;
    end
end
h = heatmap(yvalues,xvalues,cdata);
mymap = [1 1 1
    1 0 0];
colormap(mymap)
h.XDisplayLabels = ytik;
h.Title = 'Lockdown policy over time';
s = struct(h);
s.XAxis.TickLabelRotation = 0;   % horizontal
%h.Position=[0.1300 0.1100 0.8179 0.159];
h.XLabel = ('Week'); 
h.FontSize =18;
h.ColorbarVisible = 'off';
%h = heatmap(cdata);
%%
figure
x= [90,150,300,450,600,900,1200,1500];
time = [0.55,1.53,6.47,16.15,30.75,80.28,161.03,282.27];
plot(x,time,'r','LineWidth',3);
ax = gca;
ax.FontSize = 16; 
xlabel('Number of discretizations','FontSize',18)
ylabel('Runtime: seconds','FontSize',18)
%%
lambda = [0,1,2,5,10,20,50];
log_lambda = [log(0.9),log(1),log(2),log(5),log(10),log(20),log(50)];
infection_averted =[386145,401215,402393,403206,405582,405384,405169];
travel_cost = [196572213,197389518,197637565,201724792,209260624,215302070,229792389];
figure
yyaxis left
plot(log_lambda,infection_averted,'b--', ...
    'LineWidth',3);
ylim([380000 420000])
ylabel('Infections Averted','FontSize',16)
xlim([-0.5 4])
yyaxis right
plot(log_lambda,travel_cost,'r:',...
    'LineWidth',3);
ax = gca;
ax.FontSize = 16; 
title({'Objectives Based on Various \lambda'},'Fontsize',18)
xlabel('log(\lambda)','FontSize',16)
ylabel('Transportation Cost (Minutes)','FontSize',16)
figure
%x = 1900:10:2000;
%y = [75 91 105 123.5 131 150 179 203 226 249 281.5];
%bar(x,y)
%legend('Infections Averted', ...
%    'Transportation Cost (Minutes)','Fontsize',14)
%%
time = zeros(T,1);
for i = 1:T
    time(i)=i;    
end
pol = zeros(T,1);
pol(2:10,1) = 1;
pol0 = zeros(T,1);
%s0 = 0.9999; %initial state
%i0 = 0.0001;
%r0 = 1-s0-i0;
S = s0;
S0 = s0;
Sg = s0;
Su = s0;
Sd = s0;
I = i0;
I0 = i0;
Ig = i0;
Iu = i0;
Id = i0;
R = r0;
Rg = r0;
Ru = r0;
R0 = r0;
Rd = r0;
trj = zeros(T,15);
for t = 1:T
    disp(t)
    %[S,I,R] = SEIR(S*ones(26,1),I*ones(26,1),R*ones(26,1),beta,gamma,pol(t));
    %[S,I,R] = compute_total_SIR(S,I,R);
    %[S0,I0,R0] = SEIR(S0*ones(26,1),I0*ones(26,1),R0*ones(26,1),beta,gamma,pol0(t));
    %[S0,I0,R0] = compute_total_SIR(S0,I0,R0);
    %[Sg,Ig,Rg] = SEIR(Sg*ones(26,1),Ig*ones(26,1),Rg*ones(26,1),beta,gamma,Greedyres(t,4));
    %[Sg,Ig,Rg] = compute_total_SIR(Sg,Ig,Rg);
    %[Su,Iu,Ru] = SEIR(Su*ones(26,1),Iu*ones(26,1),Ru*ones(26,1),beta,gamma,Uniformres(t,4));
    %[Su,Iu,Ru] = compute_total_SIR(Su,Iu,Ru);
    [S,I,R] = SEIR(S,I,R,beta,gamma,pol(t));
    [S0,I0,R0] = SEIR(S0,I0,R0,beta,gamma,pol0(t));
    [Sg,Ig,Rg] = SEIR(Sg,Ig,Rg,beta,gamma,Greedyres(t,4));
    [Su,Iu,Ru] = SEIR(Su,Iu,Ru,beta,gamma,Uniformres(t,4));
    [Sd,Id,Rd] = SEIR(Sd,Id,Rd,beta,gamma,Densityres(t,4));
    %trj(t,1) = Sg;
    %trj(t,2) = Ig;
    %trj(t,3) = Rg;
    %trj(t,4) = Su;
    %trj(t,5) = Iu;
    %trj(t,6) = Ru;
    %trj(t,7) = S;
    %trj(t,8) = I;
    %trj(t,9) = R;
    %trj(t,10) = S0;
    %trj(t,11) = I0;
    %trj(t,12) = R0;
    [x,y,z] = compute_total_SIR(Sg,Ig,Rg);
    trj(t,1) = x;
    trj(t,2) = y;
    trj(t,3) = z;
    [x,y,z] = compute_total_SIR(Su,Iu,Ru);
    trj(t,4) = x;
    trj(t,5) = y;
    trj(t,6) = z;
    [x,y,z] = compute_total_SIR(S,I,R);
    trj(t,7) = x;
    trj(t,8) = y;
    trj(t,9) = z;
    [x,y,z] = compute_total_SIR(S0,I0,R0);
    trj(t,10) = x;
    trj(t,11) = y;
    trj(t,12) = z;
    [x,y,z] = compute_total_SIR(Sd,Id,Rd);
    trj(t,13) = x;
    trj(t,14) = y;
    trj(t,15) = z;
end
figure
plot(time,trj(:,2),'r--', ...
    time,trj(:,14),'m-', ...
    time,trj(:,5),'b', ...
    time,trj(:,8),'g:', ...
    time,trj(:,11),'k-.', ...
    'LineWidth',3);
ax = gca;
ax.FontSize = 16; 
title({'Proportion of the Population Infected Over Time'},'Fontsize',18)
xlabel('Time: week','FontSize',18)
ylabel('Proportion','FontSize',18)
legend('GreedyCut policy', ...
    'Simluation-Based policy', ...
    'Uniform policy', ...
    'Empirical policy', ...
    'No intervention','Fontsize',14)
figure
plot(time,trj(:,1),'r--', ...
    time,trj(:,13),'m-', ...
    time,trj(:,4),'b', ...
    time,trj(:,7),'g:', ...
    time,trj(:,10),'k-.', ...
    'LineWidth',3);
ax = gca;
ax.FontSize = 16; 
title({'Proportion of the Population Susceptible Over Time'},'Fontsize',18)
xlabel('Time: week','FontSize',18)
ylabel('Proportion','FontSize',18)
legend('GreedyCut policy', ...
    'Simluation-Based policy', ...
    'Uniform policy', ...
    'Empirical policy', ...
    'No intervention','Fontsize',14)
num_infections = zeros(T,5);
obj = zeros(T,5);
%costr=0
for t = 1:T
    obj(t,1)=-trj(t,2)-costr*Greedyres(t,4);
    obj(t,3)=-trj(t,5)-costr*Uniformres(t,4);
    obj(t,4)=-trj(t,8)-costr*pol(t,1);
    obj(t,5)=-trj(t,11)-costr*pol0(t,1);
    obj(t,2)=-trj(t,14)-costr*Densityres(t,4);
end
obj = -obj;
disp('Greedy')
disp(sum(obj(:,1)))
disp('Simulation-Based')
disp(sum(obj(:,2)))
disp('Uniform')
disp(sum(obj(:,3)))
disp('Empirical')
disp(sum(obj(:,4)))
disp('No intervention')
disp(sum(obj(:,5)))
figure
x = categorical({'GreedyCut' 'Simulation-Based' 'Uniform' 'Empirical' 'Do nothing'});
x = reordercats(x,{'GreedyCut' 'Simulation-Based' 'Uniform' 'Empirical' 'Do nothing'});
bar(x,sum(obj,1))
text(1:length(sum(obj,1)),sum(obj,1)+0.07,num2str(round(sum(obj,1),4)'),'vert','top','horiz','center','FontSize',16); 
ylim([0 1]);
ax = gca;
ax.FontSize = 16; 
title({'Objective Across Different Models'},'Fontsize',18)
xlabel('Policy','FontSize',18)
ylabel('Objective value','FontSize',18)
%%
time = zeros(T,1);
for i = 1:T
    time(i)=i;    
end
pol = zeros(T,1);
pol(2:10,1) = 1;
pol0 = zeros(T,1);
s0 = 0.9999; %initial state
i0 = 0.0001;
r0 = 1-s0-i0;
S = s0;
S0 = s0;
Sg = s0;
Su = s0;
I = i0;
I0 = i0;
Ig = i0;
Iu = i0;
R = r0;
Rg = r0;
Ru = r0;
R0 = r0;
trj = zeros(T,12);
for t = 1:T
    disp(t)
    [S,I,R] = SEIR(S*ones(26,1),I*ones(26,1),R*ones(26,1),beta,gamma,pol(t));
    [S,I,R] = compute_total_SIR(S,I,R);
    [S0,I0,R0] = SEIR(S0*ones(26,1),I0*ones(26,1),R0*ones(26,1),beta,gamma,pol0(t));
    [S0,I0,R0] = compute_total_SIR(S0,I0,R0);
    [Sg,Ig,Rg] = SEIR(Sg*ones(26,1),Ig*ones(26,1),Rg*ones(26,1),beta,gamma,Greedyres(t,4));
    [Sg,Ig,Rg] = compute_total_SIR(Sg,Ig,Rg);
    [Su,Iu,Ru] = SEIR(Su*ones(26,1),Iu*ones(26,1),Ru*ones(26,1),beta,gamma,Uniformres(t,4));
    [Su,Iu,Ru] = compute_total_SIR(Su,Iu,Ru);
    trj(t,1) = Sg;
    trj(t,2) = Ig;
    trj(t,3) = Rg;
    trj(t,4) = Su;
    trj(t,5) = Iu;
    trj(t,6) = Ru;
    trj(t,7) = S;
    trj(t,8) = I;
    trj(t,9) = R;
    trj(t,10) = S0;
    trj(t,11) = I0;
    trj(t,12) = R0;
end
figure
plot(time,trj(:,2),'r--', ...
    time,trj(:,5),'b', ...
    time,trj(:,11),'k-.', ...
    'LineWidth',3);
ax = gca;
ax.FontSize = 16; 
title({'Proportions of Infectious over Time'},'Fontsize',18)
xlabel('Time: week','FontSize',18)
ylabel('Proportions','FontSize',18)
legend('Greedy policy', ...
    'Uniform policy', ...
    'No intervention','Fontsize',14)
figure
plot(time,trj(:,1),'r--', ...
    time,trj(:,4),'b', ...
    time,trj(:,10),'k-.', ...
    'LineWidth',3);
ax = gca;
ax.FontSize = 16; 
title({'Proportions of Susceptible over Time'},'Fontsize',18)
xlabel('Time: week','FontSize',18)
ylabel('Proportions','FontSize',18)
legend('Greedy policy', ...
    'Uniform policy', ...
    'No intervention','Fontsize',14)
num_infections = zeros(T,4);
obj = zeros(T,3);
for t = 1:T
    obj(t,1)=-trj(t,2)-costr*Greedyres(t,4);
    obj(t,2)=-trj(t,5)-costr*Uniformres(t,4);
    obj(t,3)=-trj(t,11)-costr*pol0(t,1);
end
disp('Greedy')
disp(sum(obj(:,1)))
disp('Uniform')
disp(sum(obj(:,2)))
disp('No intervention')
disp(sum(obj(:,3)))
figure
x = categorical({'Greedy' 'Uniform'  'Do nothing'});
x = reordercats(x,{'Greedy' 'Uniform'  'Do nothing'});
bar(x,sum(obj,1))
text(1:length(sum(obj,1)),sum(obj,1)-0.05,num2str(round(sum(obj,1),4)'),'vert','bottom','horiz','center','FontSize',16); 
ylim([-1 0]);
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
        action = pol(t);
        [s0,i0,r0] = SEIR(s0,i0,r0,beta, gamma,action);
        S(t,:) = s0;
        I(t,:) = i0;
        R(t,:) = r0;
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

function [S_,I_,R_] = compute_total_SIR2(S,I,R,T)
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
