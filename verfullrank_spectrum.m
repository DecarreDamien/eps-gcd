function res = verfullrank_spectrum(A)
% verify if A, a point or interval matrix, is full rank
% if A point matrix, it is better to call isfullrank directly
% Method : with spectral radius, see J. Rohn
%
% res = 1 implies that A is verified to have full rank
    midA = mid(A);
    if(~isfullrank(midA)) % check midA id full rank
        res = false; 
        return;
    end
    radA = rad(A);
    if(norm(abs(verifypinv(midA))*radA, 1).sup < 1)
        res = true;
    else 
        if(norm(abs(verifyinv(midA'*midA))*(radA'*radA), 1).sup < 1)
            res = true;
        else 
            res = false;
        end
    end
end