function out=getqr(a);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Returns a modified QR decomposition of matrix a, such that the    %
%    diagonal elements of the 'R' matrix are all positive              %
%                                                                      % 
%Version is 1.00                                                       % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  [Q,R] = QR(A), where A is m-by-n, produces an m-by-n upper triangular
%  matrix R and an m-by-m unitary matrix Q (i.e. Q Q'=I) so that A = Q*R.
[q,r]=qr(a);

% If diagonal elements of R are negative then multiply the corresponding
% column of Q by -1; Note: the modified Q matrix is still unitary.
for i=1:rows(q);
       if r(i,i)<0;
           q(:,i)=-q(:,i);
       end
end

%Return the modified Q matrix
out=q;