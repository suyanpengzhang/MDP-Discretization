
% Define the time vector
x = 0:0.0001:1;  % Discretize the range for x'

% Define the piecewise functions
V = zeros(size(x));
V(x >= 0 & x < 0.3) = 1;
V(x >= 0.3 & x < 0.7) = 2;
V(x >= 0.7 & x <= 1) = 3;

T = zeros(size(x));
T(x >= 0 & x < 0.3) = 0.1;
T(x >= 0.3 & x < 0.5) = 0.3;
T(x >= 0.5 & x < 0.8) = 0.4;
T(x >= 0.8 & x <= 1) = 0.2;

% Perform the product of T and V over their respective intervals
product = T .* V;

% Compute the integral result (for verification)
integral_result = sum(product) * (x(2) - x(1)); % Numerical integration using the trapezoidal rule

% Find unique points in the product result
threshold = 1e-5; % Small threshold to detect changes
diffs = abs(diff(product)) > threshold; % Detect changes
breakpoints = [1, find(diffs) + 1, length(product)]; % Include first and last points

% Extract the unique values
unique_values = product(breakpoints);

% Display the integral result
disp('Integral result:');
disp(integral_result);

% Plot the piecewise functions and their product
figure;
subplot(3, 1, 1);
plot(x, T, 'DisplayName', 'T(x''|xa)');
title('Piecewise Function T(x''|xa)');
xlabel('x''');
ylabel('Value');
grid on;

subplot(3, 1, 2);
plot(x, V, 'DisplayName', 'V^n(x'')');
title('Piecewise Function V^n(x'')');
xlabel('x''');
ylabel('Value');
grid on;

subplot(3, 1, 3);
plot(x, product, 'DisplayName', 'T(x''|xa) * V^n(x'')');
hold on;
for k = 1:length(breakpoints)
    plot(x(breakpoints(k)), unique_values(k), 'ro', 'MarkerFaceColor', 'r'); % Plot breakpoints
end
title('Product of Piecewise Functions T(x''|xa) and V^n(x'')');
xlabel('x''');
ylabel('Product Value');
grid on;
hold off;

legend;

%%

% Define the 2D grid for x1 and x2
[x1, x2] = meshgrid(0:0.01:1, 0:0.01:1);  % Create a 2D grid

% Define the piecewise function V(x1, x2)
x1_dis = [0,0.3,0.7,1];
x2_dis = [0,0.3,0.7,1];
V = zeros(size(x1));
for i=1:length(x1_dis)-1
    V(x1 >= x1_dis(i) & x1 < x1_dis(i+1) & x2 >= x2_dis(i) & x2 < x2_dis(i+1)) = x1_dis(i);
end


% Define the piecewise function T(x1, x2)
T = zeros(size(x1));
T(x1 >= 0 & x1 < 0.3 & x2 >= 0 & x2 < 0.3) = 0.1;
T(x1 >= 0.3 & x1 < 0.5 & x2 >= 0.3 & x2 < 0.5) = 0.3;
T(x1 >= 0.5 & x1 < 0.8 & x2 >= 0.5 & x2 < 0.8) = 0.4;
T(x1 >= 0.8 & x1 <= 1 & x2 >= 0.8 & x2 <= 1) = 0.2;

% Perform the product of T and V over their respective intervals
product = T .* V;

% Compute the integral result using the trapezoidal rule for 2D
dx = 0.01;  % Step size in x1 direction
dy = 0.01;  % Step size in x2 direction

integral_result = sum(product(:)) * dx * dy;  % Numerical integration in 2D

% Display the integral result
disp(['Integral result: ', num2str(integral_result)]);

%%
% Define the 2D grid for x1 and x2
[x1, x2] = meshgrid(0:0.01:1, 0:0.01:1);  % Create a 2D grid

% Define the piecewise function V(x1, x2)
V = zeros(size(x1));
x1_dis = [0,0.3,0.7,1];
x2_dis = [0,0.3,0.7,1];
V = zeros(size(x1));
for i=1:length(x1_dis)-1
    if i<length(x1_dis)-1
        V(x1 >= x1_dis(i) & x1 < x1_dis(i+1) & x2 >= x2_dis(i) & x2 < x2_dis(i+1)) = x1_dis(i);
    else
        V(x1 >= x1_dis(i) & x1 <= x1_dis(i+1) & x2 >= x2_dis(i) & x2 <= x2_dis(i+1)) = x1_dis(i);
    end

end

% Define the piecewise function T(x1, x2)
T = zeros(size(x1));
T(x1 >= 0 & x1 < 0.3 & x2 >= 0 & x2 < 0.3) = 0.1;
T(x1 >= 0.3 & x1 < 0.5 & x2 >= 0.3 & x2 < 0.5) = 0.3;
T(x1 >= 0.5 & x1 < 0.8 & x2 >= 0.5 & x2 < 0.8) = 0.4;
T(x1 >= 0.8 & x1 <= 1 & x2 >= 0.8 & x2 <= 1) = 0.2;

% Visualization using surf
figure;
surf(x1, x2, V);
title('Piecewise Function V(x1, x2)');
xlabel('x1');
ylabel('x2');
zlabel('V');
colorbar;

% Visualization using imagesc
figure;
imagesc([0 1], [0 1], V); % specify the x and y range
set(gca, 'YDir', 'normal'); % set Y direction to normal
title('Piecewise Function V(x1, x2)');
xlabel('x1');
ylabel('x2');
colorbar;
axis equal tight;

% Visualization of the product T .* V
product = T .* V;

% Visualization using surf
figure;
surf(x1, x2, product);
title('Product of T(x1, x2) and V(x1, x2)');
xlabel('x1');
ylabel('x2');
zlabel('Product');
colorbar;

% Visualization using imagesc
figure;
imagesc([0 1], [0 1], product); % specify the x and y range
set(gca, 'YDir', 'normal'); % set Y direction to normal
title('Product of T(x1, x2) and V(x1, x2)');
xlabel('x1');
ylabel('x2');
colorbar;
axis equal tight;



