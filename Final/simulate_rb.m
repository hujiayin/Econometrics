%% Function to Simulate rb Data
function data = simulate_rb(a, b, h, n, initial_value)
    X = zeros(n,1);
    X(1) = initial_value;
    
    for i = 2:n
        X(i) = a * X(i-1) *h + b * X(i-1) * sqrt(h) * normrnd(0,1) + X(i-1);
    end
    
    data = X;
end