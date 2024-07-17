%This script generate discretization for SIR model
T = 20; %Time epoch
beta = 0.2*7; %initial parameter
theta = 0.25*7;
gamma = 0.07*7;
samples = importdata('samples_for_compare.mat');
pol_samples = importdata('policy_for_compare.mat');
N = 120;
Ssample =zeros(N,20);
Isample = zeros(N,20);
Rsample = zeros(N,20);
for c=1:N
    s0= samples(c,1);
    i0 = samples(c,2);
    if s0+i0>=1
        s0 = s0/(s0+i0);
        i0 = i0/(s0+i0);
        r0 = 0;
    else
        r0 = 1-s0-i0;
    end
    pol = pol_samples(c,1:T); %policy
    [S,I,R] = SEIR_trj(s0,i0,r0,beta, gamma,T,pol);
    Ssample(c,:)=S;
    Isample(c,:)=I;
    Rsample(c,:)=R;
end
%%
Nc = 400-4;
Isample=reshape(Isample,N*T,1);
Ssample=reshape(Ssample,N*T,1);
Rsample=reshape(Rsample,N*T,1);
int1s = Nc*(length(Ssample(Ssample<0.2))/(N*T));
int2s = Nc*(length(Ssample(Ssample<0.4 &Ssample>=0.2))/(N*T));
int3s = Nc*(length(Ssample(Ssample<0.6 &Ssample>=0.4))/(N*T));
int4s = Nc*(length(Ssample(Ssample<0.8 &Ssample>=0.6))/(N*T));
int5s = Nc*(length(Ssample(Ssample<=1 &Ssample>=0.8))/(N*T));
int1i = Nc*(length(Isample(Isample<0.2))/(N*T));
int2i = Nc*(length(Isample(Isample<0.4 &Isample>=0.2))/(N*T));
int3i = Nc*(length(Isample(Isample<0.6 &Isample>=0.4))/(N*T));
int4i = Nc*(length(Isample(Isample<0.8 &Isample>=0.6))/(N*T));
int5i = Nc*(length(Isample(Isample<=1 &Isample>=0.8))/(N*T));
int1r = Nc*(length(Rsample(Rsample<0.2))/(N*T));
int2r = Nc*(length(Rsample(Rsample<0.4 &Rsample>=0.2))/(N*T));
int3r = Nc*(length(Rsample(Rsample<0.6 &Rsample>=0.4))/(N*T));
int4r = Nc*(length(Rsample(Rsample<0.8 &Rsample>=0.6))/(N*T));
int5r = Nc*(length(Rsample(Rsample<=1 &Rsample>=0.8))/(N*T));

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