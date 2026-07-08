function [tot, fc, qr, qr_t, sp, sp_t, sv, sv_t] = test_fr_cond(max, N, eps)
    tot = 0;
    sp = 0;
    sp_t = 0;
    sv = 0;
    sv_t = 0;
    qr = 0;
    qr_t = 0;
    fc = 0;
    for i=1:N
        m = 2 + floor(rand(1)*(max-1));
        n = 2 + floor(rand(1)*(m-1));
        A = 2*rand(m,n) - 1;
        delta = rand(m,n) * eps;
        while ~isfullrank(A) 
            A = 2*rand(m,n) - 1;
        end
        if(fullcolrank(midrad(A, delta)))
            fc = fc + 1;
            tic;
            r = verfullrank(midrad(A, delta), "spectrum");
            t = toc;
            if r
                sp = sp + 1;
                sp_t = sp_t + t;
            end
            tic;
            r = verfullrank(midrad(A, delta), "qr");
            t = toc;
            if r
                qr = qr + 1;
                qr_t = qr_t + t;
            end
            tic;
            r = verfullrank(midrad(A, delta), "singval");
            t = toc;
            if r
                sv = sv + 1;
                sv_t = sv_t + t;
            end
        end
        tot = tot + 1;
    end 
end

[total_test, nb_fullrank, ...
    qr_method, qr_method_time, ...
    spectrum_methhod, spectrum_method_time, ...
    singularvalues_method, singular_values_method_time ] ...
 = test_fr_cond(12, 1000, 1e-4)

function [ sv, sv_t, sp, sp_t, qr, qr_t, s] = test_method(max, N, eps)
    sp = 0;
    sp_t = 0;
    sv = 0;
    sv_t = 0;
    qr = 0;
    qr_t = 0;
    s = 0;
    for i=1:N
        m = 2 + floor(rand(1)*(max-1));
        n = 2 + floor(rand(1)*(m-1));
        A = 2*rand(m,n) - 1;
        delta = rand(m,n) * eps;
        while ~isfullrank(A) 
            A = 2*rand(m,n) - 1;
        end
        tic;
        r1 = verfullrank(midrad(A, delta), "spectrum");
        t1 = toc;
        tic;
        r2 = verfullrank(midrad(A, delta), "qr");
        t2 = toc;
        tic;
        r3 = verfullrank(midrad(A, delta), "singval");
        t3 = toc;
        if ~(r1 && r2 && r3) && (r1 || r2 || r3)
            s = s + 1;
        end
        if r1
            sp = sp + 1;
            sp_t = sp_t + t1;
        end
        if r2
            qr = qr + 1;
            qr_t = qr_t + t2;
        end
        
        if r3
            sv = sv + 1;
            sv_t = sv_t + t3;
        end
    end 
end

%[singularvalues_method, singular_values_method_time, ...
%    spectrum_methhod, spectrum_method_time, ...
%    qr_method, qr_method_time,...
%    nb_not_all_success ] = test_method(100, 100, 1e-2)

function [A,D,res] = counter_ex()
    for i=1:10000
        m = 2 + floor(rand(1)*(20-1));
        n = 2 + floor(rand(1)*(m-1));
        A = 2*rand(m,n) - 1;
        D = rand(m,n) * rand(1);
        while ~isfullrank(A) 
            A = 2*rand(m,n) - 1;
        end
        if (verfullrank(midrad(A, D), "singval") || ...
                verfullrank(midrad(A, D), "qr")) && ...
                ~verfullrank(midrad(A, D), "spectrum")
            res = true;
            return;
        end 
    end 
    res = false;
end 
%[A,D] = counter_ex()

