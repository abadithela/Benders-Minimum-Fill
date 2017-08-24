function A  = fillMatrix(n, vec)
% Function to populate matrix A given a vector with all possible locations.
% A is returned is symmetric
A = 2*n*diag(ones(n,1)) + tril(ones(n),-1);
A(A == 1) = vec;
A = A + tril(A,-1)';
end