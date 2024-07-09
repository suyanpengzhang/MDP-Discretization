%This script generate discretization for SIR model
%This script generate discretization for SIR model
T = 60; %Time epoch
% Specify the path to the CSV file
file_path = 'covid-data/beta.csv';
% Load the CSV file into a MATLAB array
beta = readmatrix(file_path);
%beta = beta.*0.8;
gamma = 0.7048;
pol = randi(2,T,1)-1; %policy
%pol = zeros(T,1);
Ssample =zeros(1000,T);
Isample = zeros(1000,T);
Rsample = zeros(1000,T);
for c=1:1000
    s0 = 0.7+0.3*rand();
    i0 = 0.001+0.099*rand();
    if s0+i0>=1
        s0 = s0/(s0+i0);
        i0 = i0/(s0+i0);
        r0 = 0;
    else
        r0 = 1-s0-i0;
    end
    s0 = ones(26,1).*s0; %initial state
    i0 = ones(26,1).*i0;
    r0 = ones(26,1).*r0;
    pol = randi(2,T,1)-1; %policy
    [S,I,R] = SEIR_trj(s0,i0,r0,beta, gamma,T,pol);
    [S_,I_,R_] = compute_total_SIR(S,I,R,T);
    Ssample(c,:)=S_;
    Isample(c,:)=I_;
    Rsample(c,:)=R_;
end
%%
Nc = 50;
Isample=reshape(Isample,60000,1);
Ssample=reshape(Ssample,60000,1);
Rsample=reshape(Rsample,60000,1);
int1s = Nc*(length(Ssample(Ssample<0.2))/60000);
int2s = Nc*(length(Ssample(Ssample<0.4 &Ssample>=0.2))/60000);
int3s = Nc*(length(Ssample(Ssample<0.6 &Ssample>=0.4))/60000);
int4s = Nc*(length(Ssample(Ssample<0.8 &Ssample>=0.6))/60000);
int5s = Nc*(length(Ssample(Ssample<=1 &Ssample>=0.8))/60000);
int1i = Nc*(length(Isample(Isample<0.2))/60000);
int2i = Nc*(length(Isample(Isample<0.4 &Isample>=0.2))/60000);
int3i = Nc*(length(Isample(Isample<0.6 &Isample>=0.4))/60000);
int4i = Nc*(length(Isample(Isample<0.8 &Isample>=0.6))/60000);
int5i = Nc*(length(Isample(Isample<=1 &Isample>=0.8))/60000);
int1r = Nc*(length(Rsample(Rsample<0.2))/60000);
int2r = Nc*(length(Rsample(Rsample<0.4 &Rsample>=0.2))/60000);
int3r = Nc*(length(Rsample(Rsample<0.6 &Rsample>=0.4))/60000);
int4r = Nc*(length(Rsample(Rsample<0.8 &Rsample>=0.6))/60000);
int5r = Nc*(length(Rsample(Rsample<=1 &Rsample>=0.8))/60000);

Gs1 = 0:0.2/int1s:0.2;
Gs2 = 0.2:0.2/int2s:0.4;
Gs3 = 0.4:0.2/int3s:0.6;
Gs4 = 0.6:0.2/int4s:0.8;
Gs5 = 0.8:0.2/int5s:1;
Gs = [Gs1,Gs2,Gs3,Gs4,Gs5,1];

Gi1 = 0:0.2/int1i:0.2;
Gi2 = 0.2:0.2/int2i:0.4;
Gi3 = 0.4:0.2/int3i:0.6;
Gi4 = 0.6:0.2/int4i:0.8;
Gi5 = 0.8:0.2/int5i:1;
Gi = [Gi1,Gi2,Gi3,Gi4,Gi5,1];

Gr1 = 0:0.2/int1r:0.2;
Gr2 = 0.2:0.2/int2r:0.4;
Gr3 = 0.4:0.2/int3r:0.6;
Gr4 = 0.6:0.2/int4r:0.8;
Gr5 = 0.8:0.2/int5r:1;
Gr = [Gr1,Gr2,Gr3,Gr4,Gr5,1];



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