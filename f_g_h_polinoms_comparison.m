clear;
tic;

% Parameters
p = 5003;        % Prime modulus
N = p;            % Sequence length

% Define polynomial functions
f = @(x) mod(x.^2 + 1, p);
g = @(x) mod(x.^2 + 3.*x + 1, p);
h = @(x) mod(x.^3 - 1, p);

% Create a function to compute the Legendre sequence and measures
function [L, wd_measure, corr_measure, z_well] = analyze_poly(poly_func, N, p)
    L = zeros(1, N);
    for i = 1:N
        a = poly_func(i);
        L(i) = jacobi(a, p);
    end
    L(L == 0) = 1;  % Replace 0s with 1s

    % Well-distribution measure
    z_well = zeros(50, 50);
    for a = 1:50
        for b = 1:50
            t = floor((N - a) / b);
            W = zeros(1, t + 1);
            for j = 0:t
                W(j + 1) = sum(L(a + b * (0:j)));
            end
            z_well(a, b) = max(abs(W));
        end
    end
    wd_measure = max(z_well, [], 'all');

    % Second-order correlation measure
    z_corr = zeros(100);  % We assume d1, d2 < 100
    for d1 = 0:99
        for d2 = d1+1:99
            M = N - d2;
            C = zeros(1, M);
            for n = 1:M
                C(n) = L(n + d1) * L(n + d2);
            end
            S = cumsum(C);
            z_corr(d1+1, d2+1) = max(abs(S));
        end
    end
    corr_measure = max(z_corr, [], 'all');
end

% Run analysis for each function
[L_f, wd_f, corr_f, z_well_f] = analyze_poly(f, N, p);
[L_g, wd_g, corr_g, ~] = analyze_poly(g, N, p);
[L_h, wd_h, corr_h, ~] = analyze_poly(h, N, p);

% Theoretical bound
bound = sqrt(N) * log(N);

% Display table of results
% Display table of results
fprintf('\n%-25s | %-22s | %-25s\n', 'Function', 'Well-dist. Measure', 'Corr. Measure');
fprintf('%s\n', repmat('-', 1, 80));
fprintf('%-25s | %-22d | %-25d\n', 'f(x) = x^2 + 1', wd_f, corr_f);
fprintf('%-25s | %-22d | %-25d\n', 'g(x) = x^2 + 3x + 1', wd_g, corr_g);
fprintf('%-25s | %-22d | %-25d\n', 'h(x) = x^3 - 1', wd_h, corr_h);
fprintf('%s\n', repmat('-', 1, 80));
fprintf('Theoretical upper bound: %.2f\n', bound);

% --- Optional surface plot for one of the functions ---
figure;
surf(1:50, 1:50, z_well_f);
title('Well-Distribution Measure Surface Plot for f(x)');
xlabel('b'); ylabel('a'); zlabel('max |S(a,b)|');


% Jacobi symbol function (used for Legendre symbol since p is prime)
function j = jacobi(m,n)
    if mod(n,2)==0
        error('n must be odd');
    end
    m = mod(m,n);
    if m == 0
        j = 0;
    elseif m == 1
        j = 1;
    elseif mod(m,2) == 0
        if abs(mod(n,8)) == 1
            j = jacobi(m/2,n);
        else
            j = -jacobi(m/2,n);
        end
    else
        if mod(n,4)==3 && mod(m,4)==3
            j = -jacobi(n,m);
        else
            j = jacobi(n,m);
        end
    end
end