%% CIR
syms a b c 
syms h x xs
muX=b*(a-x);
sigmaX=c*sqrt(x);
CIR_Density = Density(muX, sigmaX, 4, 5);
%% Plot density function
figure(1)
cir_fun=subs(CIR_Density,{a,b,c,h,xs},{0.1,0.1,0.2,1,0.1});
subplot(3,1,1)
fplot(cir_fun,[0,2]);
%parameters for real density function
cc=2*a/((1-exp(-a*h))*c^2);
q=2*a*b/c^2-1;
density1=cc*exp(-cc*(x+exp(-a*h)*xs))*(x/(exp(-a*h)*xs))^(q/2)*besseli(q,2*cc*sqrt(x*exp(-a*h)*xs));
v3=subs(density1,{a,b,c,h,xs},{0.1,0.1,0.2,1,0.1});
subplot(3,1,2)
fplot(v3,[0,2]);
% Error term
subplot(3,1,3)
fplot(v3-cir_fun,[0,2]);

%% Simulate Data
alpha = 0.02;
beta = 0.01;
sigma = 0.05;
dt = 1/252;
size = 1000;
initial_X = 0.02;
X = simulate(alpha, beta, sigma, dt, size, initial_X);
%%
plot(X)
%%  
Xs = X(1:size-1);
Xx = X(2:size);
names = {'alpha', 'beta', 'sigma'};
theta = [alpha, beta, sigma];
func = matlabFunction(CIR_Density);
mleProb = @(theta) -sum(log(func(theta(1),theta(2),theta(3),dt, Xx, Xs)));
theta_initial = [0.01, 0.01, 0.03];
options = optimset('LargeScale','off');
[theta_hat, ~, ~, ~, ~, hessian] = fminunc(mleProb,theta_initial,options);
n = length(theta_hat);
% %estimate the sigma
theta_Sig = sqrt(diag(inv(hessian/n)));
% %p_value and confidence interval
p_value = zeros(1,n);
ci = zeros(n,2);
 for i = 1:n
    z_score = abs((theta_hat(i)-theta(i))/(theta_Sig(i)*sqrt(n)));
    %z_score = abs(theta_hat(i)/theta_Sig(i));
    p_value(i) = 2*(1-normcdf(z_score));
    ci(i,:)=[theta_hat(i)-1.96*theta_Sig(i)/sqrt(n),theta_hat(i)+1.96*theta_Sig(i)/sqrt(n)];
    %print the result
    fprintf('Estimate of %s: %f. p_value: %f. \n', names{i}, theta_hat(i),p_value(i));
    fprintf('Confidence interval: [%f,%f]\n', ci(i,1),ci(i,2));
 end

 %% Size Analysis
 Test_Times = 100;
 reject = zeros(3,1);
 for path = 1:Test_Times
    alpha = 0.01;
    beta =  0.01;
    sigma = 0.03;
    theta = [alpha, beta, sigma];
    dt = 1/252;
    size = 1000;
    initial_X = 0.01;
    X = simulate(alpha, beta, sigma, dt, size, initial_X);
    Xs = X(1:size-1);
    Xx = X(2:size);
    names = {'alpha', 'beta', 'sigma'};
    func = matlabFunction(CIR_Density);
    mleProb = @(theta) -sum(log(func(theta(1),theta(2),theta(3),dt, Xx, Xs)));
    theta_initial = [0.02, 0, 0.05];
    options = optimset('LargeScale','off');
    [theta_hat, ~, ~, ~, ~,hessian] = fminunc(mleProb,theta_initial,options);
    n = length(theta_hat);
    theta_Sig = sqrt(diag(inv(hessian/n)));
    p_value = zeros(1,n);
    for i = 1:n
        z_score = abs((theta_hat(i)-theta(i))/(theta_Sig(i)*sqrt(n)));
        p_value(i) = 2*(1-normcdf(z_score));
        if p_value(i) < 0.05
            reject(i) = reject(i) + 1;
        end
    end        
 end
 %% Power Analysis
 Test_Times = 100;
 reject_power = zeros(3,1);
 for path = 1:Test_Times
    alpha = 0.01;
    beta =  0.01;
    sigma = 0.03;
    theta_power = [0, 0, 0];
    dt = 1/252;
    size = 1000;
    initial_X = 0.01;
    X = simulate(alpha, beta, sigma, dt, size, initial_X);
    Xs = X(1:size-1);
    Xx = X(2:size);
    names = {'alpha', 'beta', 'sigma'};
    func = matlabFunction(CIR_Density);
    mleProb = @(theta) -sum(log(func(theta(1),theta(2),theta(3),dt, Xx, Xs)));
    theta_initial = [0.02, 0, 0.05];
    options = optimset('LargeScale','off');
    [theta_hat, ~, ~, ~, ~,hessian] = fminunc(mleProb,theta_initial,options);
    n = length(theta_hat);
    % %estimate the sigma
     theta_Sig = sqrt(diag(inv(hessian/n)));
    % %p_value and confidence interval
     p_value_power = zeros(1,n);
     for i = 1:n
        z_score_power = abs((theta_hat(i)-theta_power(i))/(theta_Sig(i)*sqrt(n)));
        p_value_power(i) = 2*(1-normcdf(z_score_power));
        if p_value_power(i) < 0.05
            reject_power(i) = reject_power(i) + 1;
        end
     end        
 end
 
%% Simulate data
beta = .1;
alpha = .05;
sigma = .05;
r0 = .04;
%X = simulate(alpha, beta, sigma, dt, size, initial_X);
obj = cir(beta,alpha,sigma,'StartState',r0);
nTrials = 1;
nPeriods = 100;   % Simulate future short over the next five years
dt = 1;      % time increment = 1 day
rng('default'); 
rPaths = simByTransition(obj,nPeriods,'nTrials',nTrials,'nSteps',nSteps);
% MLE
Xx = rPaths(1:nPeriods-1);
Xs = rPaths(2:nPeriods);
names = {'alpha', 'beta', 'sigma'};
theta = [alpha, beta, sigma];
func = matlabFunction(CIR_Density);
dt = 1/252;
mleProb = @(theta) -sum(log(func(theta(1),theta(2),theta(3),dt, Xx, Xs)));
theta_initial = [0.2, 0.2, 0.3];
options = optimset('LargeScale','off');
[theta_hat, ~, ~, ~, ~, hessian] = fminunc(mleProb,theta_initial,options);
n = length(theta_hat);
% %estimate the sigma
theta_Sig = sqrt(diag(inv(hessian/n)));
% %p_value and confidence interval
p_value = zeros(1,n);
confidence_int = zeros(n,2);
 for i = 1:n
    z_score = abs((theta_hat(i)-theta(i))/(theta_Sig(i)*sqrt(n)));
    %z_score = abs(theta_hat(i)/theta_Sig(i));
    p_value(i) = 2*(1-normcdf(z_score));
    confidence_int(i,:)=[theta_hat(i)-1.96*theta_Sig(i)/sqrt(n),theta_hat(i)+1.96*theta_Sig(i)/sqrt(n)];
    %print the result
    fprintf('The estimate of %s is %f with p_value at %f. \n', names{i}, theta_hat(i),p_value(i));
    fprintf('Confidence interval is: [%f,%f]\n', confidence_int(i,1),confidence_int(i,2));
 end
 
 
 
 
 
 
 
 
 
 
 %% Empirical Part with Fed Funds Rate
load('FRBH15.mat');
data = FRBH15{:,5};
X_data = data(~isnan(data))/100;
len = length(X_data);
Xs = X_data(1:len-1);
Xx = X_data(2:len);

%%
plot(data)
ylabel('Interest Rate')
title('1-year Interest Rate 19620102 - 20190221')
%%
%-------------Using Linear Regression to Find Parameters as Theta0
x = X_data(1:end-1); % Time series of interest rates observations
dx = diff(X_data)./sqrt(x); %dx/sqrt(x)
regressors = [1./sqrt(x) sqrt(x)];
[coefficients, ~, residuals] = regress(dx,regressors);
%Get the parameters
beta = - coefficients(2)/dt;
alpha = - coefficients(1)/coefficients(2);
sigma = std(residuals,'omitnan')/sqrt(dt);
InitialParams = [alpha, beta, sigma]; % Vector of initial parameters
%%
%----------------------------------------------------------------
mleProb = @(theta) -sum(log(func(theta(1),theta(2),theta(3),dt,Xx,Xs)));
%find the estimate and hessian matrix
options = optimset('LargeScale','off');
theta_initial = [0.2, 0.1, 0.4];
[theta_hat, ~, ~, ~, ~, hessian] = fminunc(mleProb, theta_initial, options);
n = length(theta_hat);
%estimate the sigma
theta_Sig = sqrt(diag(inv(hessian/n)));
%p_value and confidence interval
p_value = zeros(1,n);
confidence_int = zeros(n,2);
for i = 1:n
    z_score = abs((theta_hat(i)-InitialParams(i))/(theta_Sig(i)*sqrt(n)));
    p_value(i) = 2*(1-normcdf(z_score));
    confidence_int(i,:)=[theta_hat(i)-1.96*theta_Sig(i)/sqrt(n),theta_hat(i)+1.96*theta_Sig(i)/sqrt(n)];
    %print the result
    fprintf('The estimate of %s is %f with p_value at %f. \n',names{i}, theta_hat(i),p_value(i));
    fprintf('Confidence interval is: [%f,%f]\n', confidence_int(i,1),confidence_int(i,2));

end
%% Subsample Test 1
X_data = data(~isnan(data))/100;
len = length(X_data)/2;
Xs = X_data(1:len-1);
Xx = X_data(2:len);
%-------------Using Linear Regression to Find Parameters as Theta0
x = X_data(1:len-1); % Time series of interest rates observations
dx = diff(X_data(1:len))./sqrt(x); %dx/sqrt(x)
regressors = [1./sqrt(x) sqrt(x)];
[coefficients, ~, residuals] = ...
    regress(dx,regressors);
%Get the parameters
beta = - coefficients(2)/dt;
alpha = - coefficients(1)/coefficients(2);
sigma = std(residuals)/sqrt(dt);
InitialParams = [alpha, beta, sigma]; % Vector of initial parameters
%----------------------------------------------------------------
mleProb = @(theta) -sum(log(func(theta(1),theta(2),theta(3),dt,Xx,Xs)));
%find the estimate and hessian matrix
options = optimset('LargeScale','off');
theta_initial = [0.2, 0.1, 0.4];
[theta_hat, ~, ~, ~, ~, hessian] = fminunc(mleProb, theta_initial, options);
n = length(theta_hat);
%estimate the sigma
theta_Sig = sqrt(diag(inv(hessian/n)));
%p_value and confidence interval
p_value = zeros(1,n);
confidence_int = zeros(n,2);
for i = 1:n
    z_score = abs((theta_hat(i)-InitialParams(i))/(theta_Sig(i)*sqrt(n)));
    p_value(i) = 2*(1-normcdf(z_score));
    confidence_int(i,:)=[theta_hat(i)-1.96*theta_Sig(i)/sqrt(n),theta_hat(i)+1.96*theta_Sig(i)/sqrt(n)];
    %print the result
    fprintf('The estimate of %s is %f with p_value at %f. \n',names{i}, theta_hat(i),p_value(i));
    fprintf('Confidence interval is: [%f,%f]\n', confidence_int(i,1),confidence_int(i,2));
end
%% Subsample 2
X_data = data(~isnan(data))/100;
len = length(X_data)/2;
Xs = X_data(len:end-1);
Xx = X_data(len+1:end);
%-------------Using Linear Regression to Find Parameters as Theta0
x = X_data(len+1:end-1);%Time series of interest rates observations
dx = diff(X_data(len+1:end))./sqrt(x); %dx/sqrt(x)
regressors = [1./sqrt(x) sqrt(x)];
[coefficients, intervals, residuals] = ...
    regress(dx,regressors);
%Get the parameters
beta = - coefficients(2)/dt;
alpha = - coefficients(1)/coefficients(2);
sigma = std(residuals)/sqrt(dt);
InitialParams = [alpha, beta, sigma]; % Vector of initial parameters
%----------------------------------------------------------------
mleProb = @(theta) -sum(log(func(theta(1),theta(2),theta(3),dt,Xx,Xs)));
%find the estimate and hessian matrix
options = optimset('LargeScale','off');
theta_initial = [0.1, 0.1, 0.4];
[theta_hat, fval, exitflag, output, grad, hessian] = fminunc(mleProb, theta_initial, options);
n = length(theta_hat);
%estimate the sigma
theta_Sig = sqrt(diag(inv(hessian/n)));
%p_value and confidence interval
p_value = zeros(1,n);
confidence_int = zeros(n,2);
for i = 1:n
    z_score = abs((theta_hat(i)-InitialParams(i))/(theta_Sig(i)*sqrt(n)));
    p_value(i) = 2*(1-normcdf(z_score));
    confidence_int(i,:)=[theta_hat(i)-1.96*theta_Sig(i)/sqrt(n),theta_hat(i)+1.96*theta_Sig(i)/sqrt(n)];
    %print the result
    fprintf('The estimate of %s is %f with p_value at %f. \n',names{i}, theta_hat(i),p_value(i));
    fprintf('Confidence interval is: [%f,%f]\n', confidence_int(i,1),confidence_int(i,2));
end
 



