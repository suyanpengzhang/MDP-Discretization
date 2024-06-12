% Define the time vector
x = 0:0.01:1;  % Discretize the time range

% Define the piecewise functions
v = zeros(size(x));
v(x >= 0 & x < 0.3) = 1;
v(x >= 0.3 & x < 0.7) = 2;
v(x >= 0.7 & x <= 1) = 3;

p = zeros(size(x));
p(x >= 0 & x < 0.3) = 0.1;
p(x >= 0.3 & x < 0.5) = 0.3;
p(x >= 0.5 & x < 0.8) = 0.4;
p(x >= 0.8 & x <= 1) = 0.2;

% Perform the convolution
conv_result = conv(v, p, 'full') * (x(2) - x(1)); % Normalize by the sampling interval

% Define the time vector for the convolution result
x_conv = linspace(0, 2*x(end), length(conv_result));

% Find unique points in the convolution result
threshold = 1e-5; % Small threshold to detect changes
diffs = abs(diff(diff(conv_result))) > threshold; % Detect changes
breakpoints = [1, find(diffs) + 1, length(conv_result)]; % Include first and last points

% Extract the unique values
unique_values = conv_result(breakpoints);

% Display the piecewise linear function
disp('Breakpoints (t):');
disp(t_conv(breakpoints));
disp('Values (conv_result):');
disp(unique_values);

% Plot the functions and their convolution
figure;
hold on;
plot(x, v, 'DisplayName', 'v(t)');
plot(x, p, 'DisplayName', 'p(t)');
plot(x_conv, conv_result, 'DisplayName', '(f * g)(t)');
for k = 1:length(breakpoints)
    plot(x_conv(breakpoints(k)), unique_values(k), 'ro', 'MarkerFaceColor', 'r'); % Plot breakpoints
end
title('Convolution of Piecewise Constant Functions');
xlabel('x');
ylabel('Value');
legend;
grid on;
hold off;

%%

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

