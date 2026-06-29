function [isFullRank, critical_x] = fullcolrank(A)
    % Method Oettli-Prager / LP.
    % unverified version
    % Input :interval matrix
    %
    % Output :
    %   isFullRank : Bool : true iif A full rank
    %   critical_x : vector in the kernel if not full rank

    % See J Rohn: sufficient cond for int. mat. to have fcr
    Ac = mid(A)' * mid(A);
    Delta = rad(A)' * rad(A);
    n = size(Ac, 1);
    isFullRank = true;
    critical_x = [];
    
    num_combinations = 2^n;
    G = zeros(n, num_combinations);
    for i = 0:(num_combinations - 1)
        for j = 1:n
            if bitand(bitshift(i, -(j-1)), 1)
                G(j, i+1) = 1;
            else
                G(j, i+1) = -1;
            end
        end
    end

    for k = 1:num_combinations
        g = G(:, k);
        
        % Oettli-Prager :
        % |Ac*x| <= Delta*|x|  with diag(g)*x >= 0 et g'*x = 1 :
        % Ac*x - Delta*diag(g)*x <= 0
        % -Ac*x - Delta*diag(g)*x <= 0
        % Minimize : f(z) dummy, we just want existence w/ x = diag(g)z
        f = zeros(n, 1); 
        M1 = Ac * diag(g) - Delta;
        M2 = -Ac * diag(g) - Delta;
        A_lp = [M1; M2];
        b_lp = zeros(2*n, 1);
        Aeq = ones(1, n);
        beq = 1;
        lb = zeros(n, 1);
        ub = [];
        opt = optimoptions("linprog", "Display","none");
        [z_sol, ~, exitflag] = linprog(f, A_lp, b_lp, Aeq, beq, lb, ub, opt);
        
        if exitflag == 1 % if sol found
            isFullRank = false;
            critical_x = diag(g) * z_sol;
            return;
        end
    end
end