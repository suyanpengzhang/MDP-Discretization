%Draw graph
discetizations = [15,30,60,90,120,150,180,240,300,450,600,750,900,1050,1200,1350,1500]

discetization_runtime = [0.068458,0.162452,0.307201,0.588817,1.056665,1.515952,2.204892,4.096350,6.657557,16.503736,31.314643,52.977726,82.403003,118.502493,165.701297,219.595436,286.723667]

figure
set(gca,"FontSize",16)
plot(discetizations,discetization_runtime,LineWidth=2)
xlabel('Number of Discretizations',FontSize=16) 
ylabel('Runtime (seconds)',FontSize=16) 

%%
data = importdata('covid-data/greedy_res.mat');
% Assuming your data is stored in a variable called 'data'
x = data(:,1); % Extracting the first column for x-axis
y = data(:,2); % Extracting the second column for y-axis
colorValues = data(:,4); % Extracting the fourth column for color values

% Creating a scatter plot where the color of each point is determined by 'colorValues'
figure; % Opens a new figure window
scatter(x, y, 40, colorValues, 'filled'); % '40' is the marker size, you can adjust it as needed
colorbar; % Adds a color bar to the side of the plot to indicate the scale of 'colorValues'
xlabel('X-axis Label'); % Replace 'X-axis Label' with your actual label
ylabel('Y-axis Label'); % Replace 'Y-axis Label' with your actual label
title('Your Heatmap Title'); % Provide a title for your heatmap

%%
cdata = zeros(T,5);
cdata(3:10,1) = 1;
cdata(:,3) = Uniformres(:,4);
cdata(:,2) = Greedyres(:,4);
cdata(:,5) = Uniformres1(:,4);
cdata(:,4) = Greedyres1(:,4);
cdata = transpose(cdata);
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

% Define the boundaries of the regions
x1 = [0 0.6 0.6 0];
y1 = [0 0 0.2 0.2];

x2 = [0 0.6 0.6 0];
y2 = [0.2 0.2 1 1];

x3 = [0.6 1 1 0.6];
y3 = [0 0 0.2 0.2];

x4 = [0.6 1 1 0.6];
y4 = [0.2 0.2 1 1];

% Plot the regions with different colors
figure;
hold on;

fill(x1, y1, [1 0.8 0.8]); % Light red
fill(x2, y2, [0.8 1 0.8]); % Light green
fill(x3, y3, [0.8 0.8 1]); % Light blue
fill(x4, y4, [1 1 0.6]); % Light yellow
plot(0.1,0.3,'*','MarkerSize',10,'LineWidth',2,'Color','black')
plot(0.3,0.6,'*','MarkerSize',10,'LineWidth',2,'Color','black')

% Add text to each region
text(0.15, 0.3, '$$\mathbf{X}_t$$', 'HorizontalAlignment', 'center','Interpreter', 'latex', 'fontsize', 18);
text(0.35, 0.6, '$$\bar{\mathbf{X}}_t$$', 'HorizontalAlignment', 'center','Interpreter', 'latex', 'fontsize', 18);

text(0.3, 0.1, 'State 1', 'HorizontalAlignment', 'center');
text(0.3, 0.5, 'State 2', 'HorizontalAlignment', 'center');
text(0.8, 0.1, 'State 3', 'HorizontalAlignment', 'center');
text(0.8, 0.5, 'State 4', 'HorizontalAlignment', 'center');

% Set the axis limits
axis([0 1 0 1]);
axis square;

% Label the axes
xlabel('Susceptible, G_S = [0,0.6,1]','fontsize', 18);
ylabel('Infectious, G_I = [0,0.2,1]','fontsize', 18);

hold off;
%%
%%

% Define the boundaries of the regions
x11 = [0 0.3 0.3 0];
y11 = [0 0 0.2 0.2];

x12 = [0.3 0.6 0.6 0.3];
y12 = [0 0 0.2 0.2];

x21 = [0 0.3 0.3 0];
y21 = [0.2 0.2 1 1];

x22 = [0.3 0.6 0.6 0.3];
y22 = [0.2 0.2 1 1];

x3 = [0.6 1 1 0.6];
y3 = [0 0 0.2 0.2];

x4 = [0.6 1 1 0.6];
y4 = [0.2 0.2 1 1];

% Plot the regions with different colors
figure;
hold on;

fill(x11, y11, [1 0.8 0.8]); % Light red
fill(x21, y21, [0.8 1 0.8]); % Light green
fill(x12, y12, [1 0.6 0.6]); % Light red
fill(x22, y22, [0.6 1 0.6]); % Light green
fill(x3, y3, [0.8 0.8 1]); % Light blue
fill(x4, y4, [1 1 0.6]); % Light yellow
plot(0.1,0.3,'*','MarkerSize',10,'LineWidth',2,'Color','black')
plot(0.15,0.6,'*','MarkerSize',10,'LineWidth',2,'Color','black')

% Add text to each region
text(0.15, 0.3, '$$\mathbf{X}_t$$', 'HorizontalAlignment', 'center','Interpreter', 'latex', 'fontsize', 18);
text(0.2, 0.6, '$$\bar{\mathbf{X}}_t$$', 'HorizontalAlignment', 'center','Interpreter', 'latex', 'fontsize', 18);

%text(0.15, 0.1, 'State 1', 'HorizontalAlignment', 'center');
%text(0.15, 0.5, 'State 2', 'HorizontalAlignment', 'center');
%text(0.45, 0.1, 'State 3', 'HorizontalAlignment', 'center');
%text(0.45, 0.5, 'State 4', 'HorizontalAlignment', 'center');
%text(0.8, 0.1, 'State 5', 'HorizontalAlignment', 'center');
%text(0.8, 0.5, 'State 6', 'HorizontalAlignment', 'center');

% Set the axis limits
axis([0 1 0 1]);
axis square;

% Label the axes
xlabel('Susceptible, G_S'' = [0,0.3,0.6,1]','fontsize', 18);
ylabel('Infectious, G_I'' = [0,0.2,1]','fontsize', 18);

hold off;
