% Input Example
theta0 = [1, 0.6];   % True value of theta

% Simulated data
simu_num = 1000;
X = zeros(simu_num, 1);
e = normrnd(0, 1, simu_num-1, 1);

X(1) = 2;
for i = 2:1:simu_num
    X(i) =  theta0(1) + theta0(2) * X(i-1) + e(i-1);
end
Xm = X(1:simu_num-1); 
Ym = X(2:simu_num);

% Probability density function
fm = @(x, y, theta) 1/ (sqrt(2*pi)) * exp(-0.5*( (y-(theta(1)+theta(2)*x))/2 ).^2);
theta0m = [1.2, 1]; % Start point


mle(fm, Xm, Ym, theta0m, theta0)




