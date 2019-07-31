Xm = input('Input data (eg: [X1, X2, ..., Xn]):\n');
fm = input('Input the density function(eg: @(x, theta) f(x, theta)):\n');
theta0m = input('Input start point:\n');

% Input
%Xm = normrnd(5, 0.1, 100, 1)
%fm = @(x, theta) 1/ (sqrt(2*pi) * theta(2)) * exp(-0.5*( (x-theta(1))/theta(2) ).^2)
%theta0m = [2, 1]
% MLE
LLm = @(theta) sum(log(fm(Xm, theta))); % Likelihood Function
[thetam, fval, exitflag, output, grad, h] = maximize(LLm, theta0m)

%CI and p-value
h_inv = inv(h);
sigma_hat_m = sqrt(diag(h_inv))

conf_level = 0.05;
l_m = norminv(1-conf_level/2,0, sigma_hat_m);
mu_CI = [thetam(1) - l_m(1), thetam(1) + l_m(1)]
sigma_CI = [thetam(2) - l_m(2), thetam(2) + l_m(2)]
mu_p = 2 * (1-normcdf(abs(thetam(1)-5), 0, sigma_hat_m(1)))
sigma_p = 2 * (1-normcdf(abs(thetam(2)-0.1), 0, sigma_hat_m(2)))







