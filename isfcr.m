function res = isfcr(M)
% Standard method based on preconditioning and Perron-Frobenius Theory
% res=1 implies matrix A is regular
  res = 0;
  setround(0)
  % store warning mode
  wng = warning;
  warning off
  [Q,~] = qr(mid(M),0);
  A = Q'*intval(M);

  R = inv(A);
  n = size(R,1);
  C = mag(eye(n)-R*intval(A));
  x = ones(n,1);
  setround(1)
  for i=1:5
    y = C*x;
    if all(y<x)
      setround(0)
      res = 1;
      break
    else
      x = y/norm(y,inf);
    end
  end
  % restore warning and exception mode
  warning(wng)
end  % function isreg