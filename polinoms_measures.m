clear;
tic;

% Define the primes to loop over
primes_list = [2003, 3001, 4001, 5003, 6007];
N_vals = primes_list;

% Define the functions f(x), g(x), h(x)
f = @(x, p) mod(x.^2 + 1, p);
g = @(x, p) mod(x.^2 + 3*x + 1, p);
h = @(x, p) mod(x.^3 - 1, p);

% Initialize results
wd_results = zeros(length(primes_list), 4);  % Well-distribution measures
corr_results = zeros(length(primes_list), 4);  % Correlation measures

% Loop through the primes and calculate measures for each p
for i = 1:length(primes_list)
    p = primes_list(i);
    N = p;  % Set N = p

    % Calculate measures for each function
    [~, wd_f, corr_f] = analyze_poly(@(x) f(x, p), N, p);
    [~, wd_g, corr_g] = analyze_poly(@(x) g(x, p), N, p);
    [~, wd_h, corr_h] = analyze_poly(@(x) h(x, p), N, p);

    % Calculate measures for the combined sequence E_{f,g,h}
    [E_combined, wd_combined, corr_combined] = analyze_combined(N, p, f, g, h);
    
    % Store results in the corresponding rows
    wd_results(i, :) = [wd_f, wd_g, wd_h, wd_combined];
    corr_results(i, :) = [corr_f, corr_g, corr_h, corr_combined];
end

% Create a table for Well-distribution measures
wd_table = table(primes_list', wd_results(:, 1), wd_results(:, 2), wd_results(:, 3), wd_results(:, 4), ...
    'VariableNames', {'Prime (p)', 'W(E_f)', 'W(E_g)', 'W(E_h)', 'W(E_{f,g,h})'});

% Create a table for Correlation measures
corr_table = table(primes_list', corr_results(:, 1), corr_results(:, 2), corr_results(:, 3), corr_results(:, 4), ...
    'VariableNames', {'Prime (p)', 'C_2(E_f)', 'C_2(E_g)', 'C_2(E_h)', 'C_2(E_{f,g,h})'});

% Write the tables to a LaTeX file
writetable(wd_table, 'well_distribution_measures.tex', 'WriteRowNames', true);
writetable(corr_table, 'correlation_measures.tex', 'WriteRowNames', true);

disp('Tables have been written to LaTeX files: well_distribution_measures.tex and correlation_measures.tex');

% ---- Function Definitions ----

% Analyze poly function for a single polynomial
function [L, wd_measure, corr_measure] = analyze_poly(poly_func, N, p)
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

% Analyze combined sequence function for E_{f,g,h}
function [E_combined, wd_combined, corr_combined] = analyze_combined(N, p, f, g, h)
    E_combined = zeros(1, N);
    for n = 1:N
        fn = f(n, p);
        gn = g(n, p);
        hn = h(n, p);
        Jf = jacobi(fn, p);
        Jg = jacobi(gn, p);
        Jh = jacobi(hn, p);

        if mod(fn * gn, p) == 0
            E_combined(n) = 1;
        elseif Jh == -1
            E_combined(n) = Jg;
        elseif Jh == 0 || Jh == 1
            E_combined(n) = Jf;
        else
            E_combined(n) = 1;
        end
    end
    E_combined(E_combined == 0) = 1;

    % Well-distribution measure for combined sequence
    z_well = zeros(50, 50);
    for a = 1:50
        for b = 1:50
            t = floor((N - a) / b);
            W = zeros(1, t + 1);
            for j = 0:t
                W(j + 1) = sum(E_combined(a + b * (0:j)));
            end
            z_well(a, b) = max(abs(W));
        end
    end
    wd_combined = max(z_well, [], 'all');

    % Second-order correlation measure for combined sequence
    z_corr = zeros(100);
    for d1 = 0:99
        for d2 = d1+1:99
            M = N - d2;
            C = zeros(1, M);
            for n = 1:M
                C(n) = E_combined(n + d1) * E_combined(n + d2);
            end
            S = cumsum(C);
            z_corr(d1+1, d2+1) = max(abs(S));
        end
    end
    corr_combined = max(z_corr, [], 'all');
end

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