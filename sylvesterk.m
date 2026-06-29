function S = sylvesterk(p,q,k)
  p = p(:);
  q = q(:);
  np = length(p) - 1;
  nq = length(q) - 1;
  S = [ toeplitz([p ; zeros(nq-k,1)],[p(1) zeros(1,nq-k)]) ...
        toeplitz([q ; zeros(np-k,1)],[q(1) zeros(1,np-k)]) ];
end  % function sylvesterk