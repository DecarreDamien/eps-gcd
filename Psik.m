function Psi = Psik(x,y,m,n)
  x = x(:);
  y = y(:);
  Psi = [ toeplitz([x ; zeros(n,1)],[x(1) zeros(1,n)]) ...
          toeplitz([y ; zeros(m,1)],[y(1) zeros(1,m)]) ];
end  % function Psik