function res = verfullrank(A, param)
% verify if A, an interval matrix, is full rank
%
% res = 1 implies that A is verified to have full rank
    if nargin == 1
        param = ["qr", "singval", "spectrum", "pmod"];
    end
    if(~isstring(param) || (~ismember('qr', param) && ...
            ~ismember('singval', param) && ~ismember('spectrum', param) ...
            && ~ismember('pmod', param)))
        error('invalid parameter');
    end
    res = false;
    if ismember('qr', param)
        res = res || isfcr(A);
    end
    if ismember('singval', param)
        res = res || verfullrank_singval(A);
    end
    if ismember('spectrum', param)
        res =  res || verfullrank_spectrum(A);
    end
    if ismember('pmod', param)
        res = res || isfullrank(A);
    end
end
