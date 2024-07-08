%%
% Define the 2D grid for x1 and x2
[x1, x2] = meshgrid(0:0.0025:1, 0:0.0025:1);  % Create a 2D grid
xxx = 0:0.0025:1;
Gs = 0:0.0025:1;
Gi = 0:0.0025:1;
Gr = 0:0.0025:1;
V_N = zeros(size(x1));
for i=1:size(x1,1)-1
    for j=1:size(x1,1)-1
        if i<size(x1,1)-1
            if j<size(x1,1)-1
                V_N(x1 >= x1(i,i) & x1 < x1(i+1,i+1) & x2 >= x1(j,j) & x2 < x1(j+1,j+1)) = (x1(j,j)+x1(j+1,j+1))/2;
            else
                V_N(x1 >= x1(i,i) & x1 < x1(i+1,i+1) & x2 >= x1(j,j) & x2 <= x1(j+1,j+1)) = (x1(j,j)+x1(j+1,j+1))/2;
            end
        else
            if j<size(x1,1)-1
                V_N(x1 >= x1(i,i) & x1 <= x1(i+1,i+1) & x2 >= x1(j,j) & x2 < x1(j+1,j+1)) = (x1(j,j)+x1(j+1,j+1))/2;
            else
                V_N(x1 >= x1(i,i) & x1 <= x1(i+1,i+1) & x2 >= x1(j,j) & x2 <= x1(j+1,j+1)) = (x1(j,j)+x1(j+1,j+1))/2;
            end
        end
    end
end


% Visualization using imagesc
figure;
imagesc([0 1], [0 1], V_N); % specify the x and y range
set(gca, 'YDir', 'normal'); % set Y direction to normal
title('Piecewise Function V_N(s, i)');
xlabel('s');
ylabel('i');
colorbar;
axis equal tight;

p0 = importdata('transitions0_uniform_400.mat');
p1 = importdata('transitions1_uniform_400.mat');



%%
%[T0,T1] = compute_transition_vec(0.7,0.001,Gs,Gi,x1,x2,p0,p1);
% Visualization using imagesc
figure;
imagesc([0 1], [0 1], T0); % specify the x and y range
set(gca, 'YDir', 'normal'); % set Y direction to normal
title('Piecewise Function T0(s, i)');
xlabel('s');
ylabel('i');
colorbar;
axis equal tight;
% Visualization using imagesc
figure;
imagesc([0 1], [0 1], T1); % specify the x and y range
set(gca, 'YDir', 'normal'); % set Y direction to normal
title('Piecewise Function T1(s, i)');
xlabel('s');
ylabel('i');
colorbar;
axis equal tight;

%%
V_10 = zeros(size(x1));
A_10 = zeros(size(x1));
for t=1:10
    disp('T')
    disp(t)
    for i=1:size(x1,1)-1
        disp(i)
        for j=1:size(x1,1)-1
            if i<size(x1,1)-1
                if j<size(x1,1)-1
                    sample_s = (x1(i,i)+x1(i+1,i+1))/2;
                    sample_i = (x1(j,j)+x1(j+1,j+1))/2;
                    idx_0 = find_index(sample_s,sample_i,Gs,Gi);
                    [T0,T1] = compute_transition_vec(sample_s,sample_i,Gs,Gi,x1,x2,p0,p1);
                    temp0 = sample_i+sum(T0.*V_N,'all');
                    temp1 = sample_i+0.03 + sum(T1.*V_N,'all');
                    if temp1<temp0
                        A_10(x1 >= x1(i,i) & x1 < x1(i+1,i+1) & x2 >= x1(j,j) & x2 < x1(j+1,j+1))=1;
                        V_10(x1 >= x1(i,i) & x1 < x1(i+1,i+1) & x2 >= x1(j,j) & x2 < x1(j+1,j+1)) = temp1;
                    else
                        A_10(x1 >= x1(i,i) & x1 < x1(i+1,i+1) & x2 >= x1(j,j) & x2 < x1(j+1,j+1))=0;
                        V_10(x1 >= x1(i,i) & x1 < x1(i+1,i+1) & x2 >= x1(j,j) & x2 < x1(j+1,j+1)) = temp0;
                    end
                else
                    sample_s = (x1(i,i)+x1(i+1,i+1))/2;
                    sample_i = (x1(j,j)+x1(j+1,j+1))/2;
                    idx_0 = find_index(sample_s,sample_i,Gs,Gi);
                    [T0,T1] = compute_transition_vec(sample_s,sample_i,Gs,Gi,x1,x2,p0,p1);
                    temp0 = sample_i+sum(T0.*V_N,'all');
                    temp1 = sample_i+0.03 + sum(T1.*V_N,'all');
                    if temp1<temp0
                        A_10(x1 >= x1(i,i) & x1 < x1(i+1,i+1) & x2 >= x1(j,j) & x2 <= x1(j+1,j+1))=1;
                        V_10(x1 >= x1(i,i) & x1 < x1(i+1,i+1) & x2 >= x1(j,j) & x2 <= x1(j+1,j+1)) = temp1;
                    else
                        A_10(x1 >= x1(i,i) & x1 < x1(i+1,i+1) & x2 >= x1(j,j) & x2 <= x1(j+1,j+1))=0;
                        V_10(x1 >= x1(i,i) & x1 < x1(i+1,i+1) & x2 >= x1(j,j) & x2 <= x1(j+1,j+1)) = temp0;
                    end
                end
            else
                if j<size(x1,1)-1
                    sample_s = (x1(i,i)+x1(i+1,i+1))/2;
                    sample_i = (x1(j,j)+x1(j+1,j+1))/2;
                    idx_0 = find_index(sample_s,sample_i,Gs,Gi);
                    [T0,T1] = compute_transition_vec(sample_s,sample_i,Gs,Gi,x1,x2,p0,p1);
                    temp0 = sample_i+sum(T0.*V_N,'all');
                    temp1 = sample_i+0.03 + sum(T1.*V_N,'all');
                    if temp1<temp0
                        A_10(x1 >= x1(i,i) & x1 <= x1(i+1,i+1) & x2 >= x1(j,j) & x2 < x1(j+1,j+1))=1;
                        V_10(x1 >= x1(i,i) & x1 <= x1(i+1,i+1) & x2 >= x1(j,j) & x2 < x1(j+1,j+1)) = temp1;
                    else
                        A_10(x1 >= x1(i,i) & x1 <= x1(i+1,i+1) & x2 >= x1(j,j) & x2 < x1(j+1,j+1))=0;
                        V_10(x1 >= x1(i,i) & x1 <= x1(i+1,i+1) & x2 >= x1(j,j) & x2 < x1(j+1,j+1)) = temp0;
                    end
                else
                    sample_s = (x1(i,i)+x1(i+1,i+1))/2;
                    sample_i = (x1(j,j)+x1(j+1,j+1))/2;
                    idx_0 = find_index(sample_s,sample_i,Gs,Gi);
                    [T0,T1] = compute_transition_vec(sample_s,sample_i,Gs,Gi,x1,x2,p0,p1);
                    temp0 = sample_i+sum(T0.*V_N,'all');
                    temp1 = sample_i+0.03 + sum(T1.*V_N,'all');
                    if temp1<temp0
                        A_10(x1 >= x1(i,i) & x1 <= x1(i+1,i+1) & x2 >= x1(j,j) & x2 <= x1(j+1,j+1))=1;
                        V_10(x1 >= x1(i,i) & x1 <= x1(i+1,i+1) & x2 >= x1(j,j) & x2 <= x1(j+1,j+1)) = temp1;
                    else
                        A_10(x1 >= x1(i,i) & x1 <= x1(i+1,i+1) & x2 >= x1(j,j) & x2 <= x1(j+1,j+1))=0;
                        V_10(x1 >= x1(i,i) & x1 <= x1(i+1,i+1) & x2 >= x1(j,j) & x2 <= x1(j+1,j+1)) = temp0;
                    end
                end
            end
        end
    end
    A{11-t}=A_10;
    Vx{11-t}=V_10;
    V_N = V_10;
    % Visualization using imagesc
    figure;
    imagesc([0 1], [0 1], A_10); % specify the x and y range
    set(gca, 'YDir', 'normal'); % set Y direction to normal
    title('Piecewise Function A_10(s, i)');
    xlabel('s');
    ylabel('i');
    colorbar;
    axis equal tight;
    % Visualization using imagesc
    figure;
    imagesc([0 1], [0 1], V_10); % specify the x and y range
    set(gca, 'YDir', 'normal'); % set Y direction to normal
    title('Piecewise Function V_10(s, i)');
    xlabel('s');
    ylabel('i');
    colorbar;
    axis equal tight;

end


% Visualization using imagesc
figure;
imagesc([0 1], [0 1], V_10); % specify the x and y range
set(gca, 'YDir', 'normal'); % set Y direction to normal
title('Piecewise Function V_10(s, i)');
xlabel('s');
ylabel('i');
colorbar;
axis equal tight;

% Visualization using imagesc
figure;
imagesc([0 1], [0 1], A_10); % specify the x and y range
set(gca, 'YDir', 'normal'); % set Y direction to normal
title('Piecewise Function A_10(s, i)');
xlabel('s');
ylabel('i');
colorbar;
axis equal tight;

%% functions


function [t0,t1] = compute_transition_vec(ss,ii,Gs,Gi,x1,x2,p0,p1)
    t0 = zeros(size(x1));
    t1 = zeros(size(x1));
    idx_0 = find_index(ss,ii,Gs,Gi);
    for i=2:length(Gs)
        for j=2:length(Gi)
            if i>0
                if j>0
                    t0(x1 > Gs(i-1) & x1 <= Gs(i) & x2 > Gi(j-1) & x2 <= Gi(j)) = p0(idx_0,find_index((Gs(i)+Gs(i-1))/2,(Gi(j)+Gi(j-1))/2,Gs,Gi));
                    t1(x1 > Gs(i-1) & x1 <= Gs(i) & x2 > Gi(j-1) & x2 <= Gi(j)) = p1(idx_0,find_index((Gs(i)+Gs(i-1))/2,(Gi(j)+Gi(j-1))/2,Gs,Gi));
                else
                    t0(x1 > Gs(i-1) & x1 <= Gs(i) & x2 >= Gi(j-1) & x2 <= Gi(j)) = p0(idx_0,find_index((Gs(i)+Gs(i-1))/2,(Gi(j)+Gi(j-1))/2,Gs,Gi));
                    t1(x1 > Gs(i-1) & x1 <= Gs(i) & x2 >= Gi(j-1) & x2 <= Gi(j)) = p1(idx_0,find_index((Gs(i)+Gs(i-1))/2,(Gi(j)+Gi(j-1))/2,Gs,Gi));
                end
            else
                if j>0
                    t0(x1 >= Gs(i-1) & x1 <= Gs(i) & x2 > Gi(j-1) & x2 <= Gi(j)) = p0(idx_0,find_index((Gs(i)+Gs(i-1))/2,(Gi(j)+Gi(j-1))/2,Gs,Gi));
                    t1(x1 >= Gs(i-1) & x1 <= Gs(i) & x2 > Gi(j-1) & x2 <= Gi(j)) = p1(idx_0,find_index((Gs(i)+Gs(i-1))/2,(Gi(j)+Gi(j-1))/2,Gs,Gi));
                else
                    t0(x1 >= Gs(i-1) & x1 <= Gs(i) & x2 >= Gi(j-1) & x2 <= Gi(j)) = p0(idx_0,find_index((Gs(i)+Gs(i-1))/2,(Gi(j)+Gi(j-1))/2,Gs,Gi));
                    t1(x1 >= Gs(i-1) & x1 <= Gs(i) & x2 >= Gi(j-1) & x2 <= Gi(j)) = p1(idx_0,find_index((Gs(i)+Gs(i-1))/2,(Gi(j)+Gi(j-1))/2,Gs,Gi));
                end
            end
        end
    end
    t0 = t0./sum(t0,'all');
    t1 = t1./sum(t1,'all');

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

function idx = find_index(s,i,Gs,Gi)
    lgs = length(Gs)-1;
    lgi = length(Gi)-1;
    ls = max(find(Gs<s));
    li = max(find(Gi<i));
    idx = (ls-1)*(lgi)+li;
end