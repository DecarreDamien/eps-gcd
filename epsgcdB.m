function [res,E] = epsgcdB(P,Q,k,N) 
% P,Q interval polynomials, k integer, N number of iteration of random
% choice for X
% res = true implies that there exist p in P and q in Q s.t. 
% deg gcd(p,q) = k
% E give details on why the verification failed
if(nargin == 3)
  N = 1;
end
 E.failure = "no failure";
 E.code = int32(0);

  P = P(:).';
  Q = Q(:).';
  degP = length(P) - 1;
  degQ = length(Q) - 1;
  dim = max(degP,degQ);
% Bezout matrix
  B = bezout(P,Q);
  
% first step : B_{m-k} full rank
  if(~verfullrank(B(:,k+1:end)))
    if(~verfullrank(B(:,k+2:end)))
        E.failure = 'B_{m-k-1} is not full rank either';
        E.code = int32(1);
    else 
        E.failure = 'B_{m-k} is not full rank';
        E.code = int32(1);
    end
    res = false;
    return;
  end
% second step : M_i full rank
  M1 = bezout(intval(inf(P)), intval(inf(Q)));
  M1 = M1(:,k:end);
  M2 = bezout(intval(sup(P)), intval(sup(Q)));
  M2 = M2(:,k:end);
  if(~isfullrank(M1) || ~isfullrank(M2))
    %tries an other choice of M_i
    p_mid = mid(P);
    q_mid = mid(Q);
    p_rad = rad(P);
    q_rad = rad(Q);

    X = [p_mid q_mid];
    X_grad = gradientinit(X);
    Bm_k = bezout(X_grad(1:degP+1), X_grad(degP+2:end));
    Bm_k = Bm_k(:,k:end);
    
    [Q_num, ~, P_mat] = qr(Bm_k.x);
    R_grad = Q_num' * (Bm_k * P_mat);
    idx = min(size(R_grad));
    r_nn = R_grad(idx, idx);
    grad_r_nn = r_nn.dx;
    grad_P = grad_r_nn(1:degP+1);
    grad_Q = grad_r_nn(degP+2:end);

    p1 = p_mid + sign(grad_P) .* p_rad;
    q1 = q_mid + sign(grad_Q) .* q_rad;
    p2 = p_mid - sign(grad_P) .* p_rad;
    q2 = q_mid - sign(grad_Q) .* q_rad;
    M1 = bezout(intval(p1), intval(q1));
    M1 = M1(:,k:end);
    M2 = bezout(intval(p2), intval(q2));
    M2 = M2(:,k:end);
    if(~isfullrank(M1) || ~isfullrank(M2))
        res = false;
        E.failure = 'The choice of M_i could not allow to continue';
        E.code = int32(2);
        return;
    end
  end  

% third step : X'M_1(1,1) * X'M_2(1,1) <= 0
% Random choices of X
  for i = 1:N
        X = rand([dim, dim-k+1]);
        while (~isfullrank(X))
           X = rand([dim, dim-k+1]);
        end
      
        v = 1:dim-k+1;
        en = zeros(length(v),1);
        en(1) = 1;
        Rinf = verifylss(X' * M1,en);
        Rsup = verifylss(X' * M2,en);
        if Rinf(1)*Rsup(1) <= 0
          res = true;
          return;
        else
          res = false;
        end
    end
    if(verfullrank(B(:,k:end)))
        E.failure = 'B_{m-k+1} is full rank';
        E.code = int32(4);
    else 
        E.failure = 'Did not found suitable X';
        E.code = int32(3);
    end
end  % function epsgcdB
