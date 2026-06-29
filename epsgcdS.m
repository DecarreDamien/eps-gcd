function [P,Q, res, flag] = epsgcdS(p,q,k) 
% Input :
% p,q polynomials, k input degree
%
% Output :
% P,Q interval polynomial containing respectively p and q
% and f,g such that deg gcd(f,g) >= k
% if res = true, then it is an equality
% if flag != -1, then the choice of k is supposed to be sub-optimal
% flag then contain a guess of a better value (could be wrong)

  p = p(:).';
  q = q(:).';
  np = length(p) - 1;
  nq = length(q) - 1;
  r = [p q]';
  Sk = sylvesterk(p,q,k);
    
  [~,~, V] = svd(Sk);
  z = V(:,end); %z approx. in Ker Sk if Sk rank deficient
  Psi = Psik(z(1:nq-k+1),z(nq-k+2:end),nq,np);
  [~, S, V] = svd(Psi);
  svdPsi = diag(S);
  threshhold = 1e-15;
  small = find(svdPsi<threshhold);
  if any(small)
      %disp(['It seems k = ' int2str(np+nq+2-small(1)) ' is a better choice'])
      flag = np+nq+2-small(1); % dim ker Psik
  else
      flag = -1; 
  end
  C = V(:,end-k:end)'; % C' approx. basis of ker psik
% check Cr ~= 0, should be fine with the choice of C
  if(~all(C*r > 0 | C*r < 0))
      disp('The choice of C could not allow to continue')
      P = NaN(np+1);
      Q = NaN(nq+1);
      res = false;
      return;
  end 
  G = [ Psi ; C ];
% Check G regular
  if(~isfullrank(G))
      disp('G not regular')
      P = NaN(np+1);
      Q = NaN(nq+1);
      res = false;
      return;
  end 
  res = verifylss(G, [ zeros(size(Psi,1),1) ; C*r ],'illco');
  PQ = res(1:np+nq+2);
  P = PQ(1:np+1);
  Q = PQ(np+2:end);
% check if S_{k+1} if full rank
  if((k == np || k == nq) || verfullrank(sylvesterk(P,Q,k+1), "qr"))
      res = true;
  else 
      res = false;
  end
  P = P(:)';
  Q = Q(:)';
% scaling P and Q to contain p and q
  P = P + midrad(0,norm(p-P, inf).sup); 
  Q = Q + midrad(0,norm(q-Q, inf).sup);
end  % function epsgcdS
