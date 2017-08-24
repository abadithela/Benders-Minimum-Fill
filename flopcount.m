% Calculating Numeric Factorization Flop Count:
load KKT
n = size(kkt,1);
K_reg = sparse(kkt + size(kkt,1)*eye(size(kkt,1)));

% Regular LDL' factorization flop count:
[L, D, Par, fl] = ldlsparse(K_reg);

% Finding best permutation from Bender's:
[P, PEO, Fill] = bmf(K_reg);

% Calculating total number of flops required for numeric factorization:
[L_bmf, D_bmf, Par_bmf, fl_bmf] = ldlsparse(K_reg, PEO);

% Numeric factorization from other heuristics:
P_amd = amd(K_reg);
P_cp = colperm(K_reg);
P_camd = colamd(K_reg);
P_rcm = symrcm(K_reg);

[L_amd, D_amd, Par_amd, fl_amd] = ldlsparse(K_reg, P_amd);
[L_cp, D_cp, Par_cp, fl_cp] = ldlsparse(K_reg, P_cp);
[L_camd, D_camd, Par_camd, fl_camd] = ldlsparse(K_reg, P_camd);
[L_rcm, D_rcm, Par_rcm, fl_rcm] = ldlsparse(K_reg, P_rcm);

% Solve flop counts:

Lnz = nnz(tril(L,-1));
Lnz_bmf = nnz(tril(L_bmf,-1));
Lnz_amd = nnz(tril(L_amd,-1));
Lnz_cp = nnz(tril(L_cp,-1));
Lnz_camd = nnz(tril(L_camd,-1));
Lnz_rcm = nnz(tril(L_rcm,-1));

nz = [Lnz, Lnz_bmf, Lnz_amd, Lnz_cp, Lnz_camd, Lnz_rcm];
flops = [fl, fl_bmf, fl_amd, fl_cp, fl_camd, fl_rcm];

solve_times = zeros(1,6); % For a 1 GHz processor
% Average solve times assuming Newton's method is called 6 times and each
% Newton step on average takes 20 iterations
for i = 1:6
    solve_times(i) = 120*(4*nz(i) + n + flops(i))/1e9*1e3; % in ms
end