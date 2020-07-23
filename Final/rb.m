syms a b x
muX_rb=a*x;
sigmaX_rb=b*x;

rb_Density = Density(muX_rb, sigmaX_rb, 4, 5);

%% Simulate Data
alpha = 0.02;
beta = 0.02;
dt = 1/252;
size = 1000;
initial_X = 0.1;
X = simulate_rb(alpha, beta, dt, size, initial_X);

%%
%%  
Xs = X(1:size-1);
Xx = X(2:size);
names = {'alpha', 'beta'};
theta = [alpha, beta];
func = matlabFunction(rb_Density);
mleProb = @(theta) -sum(log(func(theta(1),theta(2),dt, Xx, Xs)));
theta_initial = [0.2, 0.2];
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
 
 %%
  %% Empirical Part with Fed Funds Rate
load('FRBH15.mat');
data = FRBH15{:,5};
X_data = data(~isnan(data))/100;
len = length(X_data);
Xs = X_data(1:len-1);
Xx = X_data(2:len);

%%
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
mleProb = @(theta) -sum(log(func(theta(1),theta(2),dt,Xx,Xs)));
%find the estimate and hessian matrix
options = optimset('LargeScale','off');
theta_initial = [0.2, 0.2];
[theta_hat, ~, ~, ~, ~, hessian] = fminunc(mleProb, theta_initial, options);
n = length(theta_hat);
%estimate the sigma
theta_Sig = sqrt(diag(inv(hessian/n)));
%p_value and confidence interval
p_value = zeros(1,n);
confidence_int = zeros(n,2);
for i = 1:n
    z_score = abs((theta_hat(i))/(theta_Sig(i)*sqrt(n)));
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
%--------------------------------------------------------------
mleProb = @(theta) -sum(log(func(theta(1),theta(2),dt,Xx,Xs)));
%find the estimate and hessian matrix
options = optimset('LargeScale','off');
theta_initial = [0.2, 0.2];
[theta_hat, ~, ~, ~, ~, hessian] = fminunc(mleProb, theta_initial, options);
n = length(theta_hat);
%estimate the sigma
theta_Sig = sqrt(diag(inv(hessian/n)));
%p_value and confidence interval
p_value = zeros(1,n);
confidence_int = zeros(n,2);
for i = 1:n
    z_score = abs((theta_hat(i))/(theta_Sig(i)*sqrt(n)));
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
mleProb = @(theta) -sum(log(func(theta(1),theta(2),dt,Xx,Xs)));
%find the estimate and hessian matrix
options = optimset('LargeScale','off');
theta_initial = [0.2, 0.2];
[theta_hat, fval, exitflag, output, grad, hessian] = fminunc(mleProb, theta_initial, options);
n = length(theta_hat);
%estimate the sigma
theta_Sig = sqrt(diag(inv(hessian/n)));
%p_value and confidence interval
p_value = zeros(1,n);
confidence_int = zeros(n,2);
for i = 1:n
    z_score = abs((theta_hat(i))/(theta_Sig(i)*sqrt(n)));
    p_value(i) = 2*(1-normcdf(z_score));
    confidence_int(i,:)=[theta_hat(i)-1.96*theta_Sig(i)/sqrt(n),theta_hat(i)+1.96*theta_Sig(i)/sqrt(n)];
    %print the result
    fprintf('The estimate of %s is %f with p_value at %f. \n',names{i}, theta_hat(i),p_value(i));
    fprintf('Confidence interval is: [%f,%f]\n', confidence_int(i,1),confidence_int(i,2));
end
 

