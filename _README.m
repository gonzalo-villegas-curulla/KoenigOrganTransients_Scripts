% Load or define beta and nu values
beta_values = MX(:,:,26);
nu_values = MX(:,:,27);

beta_values = beta_values(:);
nu_values = nu_values(:);
valid_indices = ~isnan(beta_values) & ~isnan(nu_values);

beta_values = beta_values(valid_indices);
nu_values = nu_values(valid_indices);


% Define the range and step for testing bases
base_range = -2 : 0.0001 : 2.2;
num_bases = length(base_range);
residuals = zeros(1, num_bases);

% Loop over each base to find the one that best linearizes the relationship
for i = 1:num_bases
    base = base_range(i);
    transformed_nu = base .^ nu_values;
    
    % Perform linear regression: fit transformed_nu against beta_values
    p = polyfit(beta_values, transformed_nu, 1);
    fitted_values = polyval(p, beta_values);
    
    % Calculate residual as a measure of fit quality
    residuals(i) = sum((fitted_values - transformed_nu).^2);
end

% Find the base that minimizes the residual
[~, best_index] = min(residuals);
best_base = base_range(best_index);

% Display the best base and fit plot
fprintf('Best base to linearize the relationship: %.4f\n', best_base);

% Replot using the best base to visualize the linear fit
transformed_nu_best = best_base .^ nu_values;
p_best = polyfit(beta_values, transformed_nu_best, 1);
fitted_values_best = polyval(p_best, beta_values);

% Plot results
figure;
scatter(beta_values, transformed_nu_best, 'b', 'filled'); hold on;
plot(beta_values, fitted_values_best, 'r', 'LineWidth', 1.5);
xlabel('\beta values');
ylabel(['Best transformed \nu values (', num2str(best_base), '^{\nu})']);
title('Linearized Relationship between \beta and \nu');
legend('Transformed Data', 'Linear Fit');
grid on;

figure();
plot(base_range, log10(residuals)); xlabel('base range'); ylabel('residuals '); grid on;

%%

% Load or define beta and nu values
beta_values = MX(:,:,26);
nu_values   = MX(:,:,27);
beta_values = beta_values(:);
nu_values   = nu_values(:);
valid_indices = ~isnan(beta_values) & ~isnan(nu_values);
beta_values = beta_values(valid_indices);
nu_values   = nu_values(valid_indices);


% Define the range and step for testing bases
expo_range = -2 : 0.001 : -0.01;
num_expos = length(expo_range);
residuals = zeros(1, num_expos);

% Loop over each base to find the one that best linearizes the relationship
for i = 1:num_expos
    expo = expo_range(i);
    transformed_beta = beta_values.^(expo);
    
    % Perform linear regression: fit transformed_nu against beta_values
    p = polyfit(nu_values, transformed_beta,  1);
    fitted_values = polyval(p, nu_values);
    
    % Calculate residual as a measure of fit quality
    residuals(i) = sum((fitted_values - transformed_beta).^2);
end

% Find the base that minimizes the residual
[~, best_index] = min(residuals);
best_expo = expo_range(best_index);

% Display the best base and fit plot
fprintf('Best base to linearize the relationship: %.4f\n', best_expo);

% Replot using the best base to visualize the linear fit
transformed_beta_best = beta_values.^best_expo;
p_best = polyfit(nu_values, transformed_beta_best, 1);
fitted_values_best = polyval(p_best, nu_values);

% Plot results
figure;
scatter(nu_values, transformed_beta_best, 'b', 'filled'); hold on;
plot(nu_values, fitted_values_best, 'r', 'LineWidth', 1.5);
xlabel('\beta values');
ylabel(['Best transformed \nu values (', num2str(best_expo), '^{\nu})']);
title('Linearized Relationship between \beta and \nu');
legend('Transformed Data', 'Linear Fit');
grid on;

figure();
plot(expo, log10(residuals)); xlabel('base range'); ylabel('residuals '); grid on;

%%

% Flatten and filter NaN values from MX(:,:,26) and MX(:,:,27)
beta_matrix = MX(:,:,26);
nu_matrix = MX(:,:,27);

% Vectorize and filter out NaN values
beta_values = beta_matrix(:);
nu_values = nu_matrix(:);
valid_indices = ~isnan(beta_values) & ~isnan(nu_values);

% Extract only valid (non-NaN) pairs
beta_values = beta_values(valid_indices);
nu_values = nu_values(valid_indices);

% Define the range of fractional exponents to test
exponent_range = -0.1:0.001:0.8;  % Adjust this range as needed for more precision
num_exponents = length(exponent_range);
residuals = zeros(1, num_exponents);

% Loop over each fractional exponent to find the best linearization
for i = 1:num_exponents
    exponent = exponent_range(i);
    transformed_beta = beta_values .^ exponent;
    
    % Perform linear regression: fit transformed_beta against nu_values
    p = polyfit(transformed_beta, nu_values, 1);
    fitted_values = polyval(p, transformed_beta);
    
    % Calculate residual as a measure of fit quality
    residuals(i) = sum((fitted_values - nu_values).^2);
end

% Find the exponent that minimizes the residual
[~, best_index] = min(residuals);
best_exponent = exponent_range(best_index);

% Display the best exponent and plot the linear fit
fprintf('Best fractional exponent to linearize the relationship: %.4f\n', best_exponent);

% Replot using the best exponent to visualize the linear fit
transformed_beta_best = beta_values .^ best_exponent;
p_best = polyfit(transformed_beta_best, nu_values, 1);
fitted_values_best = polyval(p_best, transformed_beta_best);

% Plot results
figure;
scatter(transformed_beta_best, nu_values, 'b', 'filled'); hold on;
plot(transformed_beta_best, fitted_values_best, 'r', 'LineWidth', 1.5);
xlabel(['Transformed \beta values (\beta^{', num2str(best_exponent), '})']);
ylabel('\nu values');
title('Linearized Relationship between \beta and \nu');
legend('Transformed Data', 'Linear Fit');
grid on;


figure();
plot(exponent_range,residuals);

%%

% Flatten and filter NaN values from MX(:,:,26) and MX(:,:,27)
beta_matrix = MX(:,:,26);
nu_matrix = MX(:,:,27);

% Vectorize and filter out NaN values
beta_values = beta_matrix(:);
nu_values = nu_matrix(:);
valid_indices = ~isnan(beta_values) & ~isnan(nu_values);

% Extract only valid (non-NaN) pairs
beta_values = beta_values(valid_indices);
nu_values = nu_values(valid_indices);

% Define the range of logarithmic bases to test
log_bases = linspace(0.5, 20, 100);  % Test bases from 1.1 to 5 with 100 steps
num_bases = length(log_bases);
residuals = zeros(1, num_bases);

% Loop over each logarithmic base to find the best linearization
for i = 1:num_bases
    % Apply the logarithm with current base
    transformed_beta = log(beta_values) / log(log_bases(i));
    
    % Perform a linear regression fit between transformed_beta and nu_values
    p = polyfit(transformed_beta, nu_values, 1);
    y_fit = polyval(p, transformed_beta);
    
    % Calculate the sum of squared residuals
    residuals(i) = sum((nu_values - y_fit).^2);
end

% Find the base with the minimum residual
[~, best_idx] = min(residuals);
best_base = log_bases(best_idx);

% Display the best base and plot the result
disp(['Best base for logarithmic transformation of beta: ', num2str(best_base)]);

% Plot the linearized relationship using the best base
transformed_beta = log(beta_values) / log(best_base);
p_best = polyfit(transformed_beta, nu_values, 1);
y_fit_best = polyval(p_best, transformed_beta);

figure;
scatter(transformed_beta, nu_values, 'b');
hold on;
plot(transformed_beta, y_fit_best, 'r', 'LineWidth', 1.5);
xlabel(['log_{\it{', num2str(best_base), '}}(\beta)']);
ylabel('\nu');
title('Linearized Relationship between Transformed \beta and \nu');
legend('Data', 'Best Linear Fit');
grid on;

figure();
plot(log_bases, residuals);

%%

a = 100;
b = 2^6;
c = 1;
d = 1;

t = linspace(0,2,5e3);
y = a./( (1 + d*exp(-b*(t-c))).^(1/d) );

figure(2); clf;
plot(t, y); grid on;


%%

% Assuming MX(:,:,26) and MX(:,:,27) contain beta and nu values
beta_values = MX(:,:,26);
nu_values = MX(:,:,27);

% Vectorize, ignoring NaNs
beta_vec = beta_values(~isnan(beta_values));
nu_vec = nu_values(~isnan(nu_values));

N = 1000;
% Preallocate for residuals
residuals_log = zeros(1, N);
residuals_exp = zeros(1, N);
residuals_power = zeros(1, N);
residuals_hyperbolic = zeros(1, N);

% Define ranges for transformations
log_bases = linspace(1.1, 10, N); % Range for logarithmic transformation bases
exp_bases = linspace(1.1, 2.5, N); % Range for exponential base on nu
power_exponents = linspace(0.1, 5, N); % Range for power transformation on beta
hyperbolic_scales = linspace(-1, 1, N); % Range for scaling in atan(nu) vs. log(beta)

% --- 1. Logarithmic Transform on beta ---
for i = 1:length(log_bases)
    log_beta = log(beta_vec) / log(log_bases(i)); % Log transform with varying base
    p = polyfit(log_beta, nu_vec, 1); % Linear fit
    fitted_nu = polyval(p, log_beta);
    residuals_log(i) = sum((nu_vec - fitted_nu).^2); % Store residuals
end

% --- 2. Exponential Transform on nu ---
for i = 1:length(exp_bases)
    exp_nu = exp_bases(i).^nu_vec; % Exponential transform on nu
    p = polyfit(beta_vec, exp_nu, 1); % Linear fit
    fitted_exp_nu = polyval(p, beta_vec);
    residuals_exp(i) = sum((exp_nu - fitted_exp_nu).^2); % Store residuals
end

% --- 3. Power Transform on beta ---
for i = 1:length(power_exponents)
    power_beta = beta_vec.^(1 / power_exponents(i)); % Fractional power transform on beta
    p = polyfit(power_beta, nu_vec, 1); % Linear fit
    fitted_nu_power = polyval(p, power_beta);
    residuals_power(i) = sum((nu_vec - fitted_nu_power).^2); % Store residuals
end

% --- 4. Hyperbolic Transform on nu (atan) and Log Transform on beta ---
for i = 1:length(hyperbolic_scales)
    hyperbolic_nu = atan(nu_vec * hyperbolic_scales(i)); % atan(nu * scale)
    log_beta = log(beta_vec); % Standard log on beta for comparison
    p = polyfit(log_beta, hyperbolic_nu, 1); % Linear fit
    fitted_hyperbolic_nu = polyval(p, log_beta);
    residuals_hyperbolic(i) = sum((hyperbolic_nu - fitted_hyperbolic_nu).^2); % Store residuals
end

% Find optimal parameters and plot residuals
[~, opt_log_idx] = min(residuals_log);
[~, opt_exp_idx] = min(residuals_exp);
[~, opt_power_idx] = min(residuals_power);
[~, opt_hyperbolic_idx] = min(residuals_hyperbolic);

% Plot residuals for analysis
figure(1); clf;
subplot(2,2,1);
plot(log_bases, residuals_log, '-o');
title('Logarithmic Transform on beta');
xlabel('Logarithm Base');
ylabel('Residuals');

subplot(2,2,2);
plot(exp_bases, residuals_exp, '-o');
title('Exponential Transform on nu');
xlabel('Exponential Base');
ylabel('Residuals');

subplot(2,2,3);
plot(power_exponents, residuals_power, '-o');
title('Power Transform on beta');
xlabel('Power Exponent');
ylabel('Residuals');

subplot(2,2,4);
plot(hyperbolic_scales, residuals_hyperbolic, '-o');
title('Hyperbolic Transform (atan) on nu');
xlabel('Hyperbolic Scale Factor');
ylabel('Residuals');

% Display optimal parameters
fprintf('Optimal Logarithmic Base: %.2f\n', log_bases(opt_log_idx));
fprintf('Optimal Exponential Base: %.2f\n', exp_bases(opt_exp_idx));
fprintf('Optimal Power Exponent: %.2f\n', power_exponents(opt_power_idx));
fprintf('Optimal Hyperbolic Scale Factor: %.2f\n', hyperbolic_scales(opt_hyperbolic_idx));

% Plot transformed beta and nu values using optimal parameters for each transformation
figure(2);clf;

% 1. Logarithmic Transformation on beta
opt_log_beta = log(beta_vec) / log(log_bases(opt_log_idx));
subplot(2, 2, 1);
plot(opt_log_beta, nu_vec, 'o');
title(['Log Transform on beta (Base = ', num2str(log_bases(opt_log_idx)), ')']);
xlabel('Transformed beta');
ylabel('nu');

% 2. Exponential Transformation on nu
opt_exp_nu = exp_bases(opt_exp_idx).^nu_vec;
subplot(2, 2, 2);
plot(beta_vec, opt_exp_nu, 'o');
title(['Exp Transform on nu (Base = ', num2str(exp_bases(opt_exp_idx)), ')']);
xlabel('beta');
ylabel('Transformed nu');

% 3. Power Transformation on beta
opt_power_beta = beta_vec.^(1 / power_exponents(opt_power_idx));
subplot(2, 2, 3);
plot(opt_power_beta, nu_vec, 'o');
title(['Power Transform on beta (Exponent = ', num2str(1 / power_exponents(opt_power_idx)), ')']);
xlabel('Transformed beta');
ylabel('nu');

% 4. Hyperbolic Transformation (atan) on nu
opt_hyperbolic_nu = atan(nu_vec * hyperbolic_scales(opt_hyperbolic_idx));
log_beta_standard = log(beta_vec);
subplot(2, 2, 4);
plot(log_beta_standard, opt_hyperbolic_nu, 'o');
title(['Hyperbolic Transform on nu (Scale Factor = ', num2str(hyperbolic_scales(opt_hyperbolic_idx)), ')']);
xlabel('log(beta)');
ylabel('Transformed nu (atan)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% The best relation between beta and nu seems to be a power transformation on beta that resembles taking the square root of beta (instead of a power 0.5 it's actually 0.56). Why would that be?
% ChatGPT said:
% 
% The result of a power transformation close to a square root for the relationship between ββ and νν suggests that the growth rate parameter ββ (which controls how quickly the function reaches its steady state) has a nonlinear, sublinear relationship with the onset-sharpness parameter νν. Here’s why a transformation around the square root might fit well in this context:
% 
%     Growth Dynamics of the Logistic Function: In the generalized logistic function, ββ governs the rate at which the function approaches its asymptote, while νν determines the shape of the curve’s growth phase. If νν is high, the onset is sharper, but ββ dictates the speed of growth after the initial onset. The fact that a square root-type transformation fits best suggests that these two parameters do not scale linearly with one another in the context of how the logistic function grows over time.
% 
%     Relationship to Diffusion and Growth Processes: In many natural systems with growth or decay dynamics, such as diffusion or logistic growth, scaling laws often appear where growth rates vary with the square root of time or other parameters. This sublinear relationship (exponent < 1) is common in processes where initial growth happens quickly but then tapers off more gradually as it reaches a limiting capacity. The square root-like scaling can capture this type of transition between rapid initial change and more gradual stabilization, mirroring how ββ and νν interact.
% 
%     Approximation of Complex Interactions: The parameters ββ and νν in the logistic model control different aspects of the growth curve, but their combined effect on the rate of approach to steady state may depend on a nonlinear interaction. The optimal exponent (0.56) can be interpreted as an empirical approximation of this complex interaction. This value being close to 0.5 implies a significant yet diminishing dependence of the growth rate on the sharpness factor as ββ increases.
% 
%     Interpretation of 0.56 Exponent: The exponent of 0.56, close to the square root, may reflect a nonlinear, moderated relationship between the amplitude growth rate and the sharpness of onset. This suggests that increases in νν initially result in significant changes to ββ, but as νν becomes larger, these changes in ββ are less pronounced, reflecting a diminishing sensitivity of the growth rate to onset sharpness. This type of behavior is consistent with scaling laws in systems constrained by limits or saturation points.
% 
% In summary, a transformation around the square root of ββ implies that the growth dynamics of your system are driven by an initial rapid change modulated by the parameter νν, but with a diminishing impact of further increases in νν on ββ due to the inherent nonlinear nature of logistic growth processes.