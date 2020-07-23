%% Function to Simulate Data
function data = simulate(a, b, c, h, n, initial_value)
    X = zeros(n,1);
    X(1) = initial_value;
    for i = 2:n
        X(i) =b * (a - X(i-1)) *h + c * sqrt(X(i-1)*h)* normrnd(0,1) + X(i-1);
    end
    data = X;
end




