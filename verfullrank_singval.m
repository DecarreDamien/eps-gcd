function res = verfullrank_singval(A)
% verify if A, a point or interval matrix, is full rank
% Method : singular values
%
% res = 1 implies that A is verified to have full rank
    try
        sigma = verifysymeigall(A'*A);
        res = all(sigma.inf > 0);
    catch
        res = false;
    end
end