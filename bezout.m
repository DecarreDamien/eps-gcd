function B = bezout(p,q)
  degP = length(p) - 1;
  degQ = length(q) - 1;
  n = max(degP,degQ);
  p = p(:).';%force row vector
  q = q(:).';
  % Compute all 4 matrices for Barnett formula in one go, up to a flipud
  % for 2 of them
  A = toeplitz([p(degP+1) zeros(1,n-1)],[fliplr(p) zeros(1,2*n-degP-1)]);
  B = toeplitz([q(degQ+1) zeros(1,n-1)],[fliplr(q) zeros(1,2*n-degQ-1)]);
  % Barnett formula, flipud needed for two matrices, but putting it outside
  % gives the same result
  B = flipud( A(:,n+1:end)*B(:,1:n) - B(:,n+1:end)*A(:,1:n) );
end  % function bezout