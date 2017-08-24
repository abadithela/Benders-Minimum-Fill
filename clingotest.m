%filename = 'chordal_sat.lp';
filename = 'kkt_ischordal.lp';
[status, cmdout] = system(sprintf('clingo %s 0 | grep SAT', filename));

if contains(cmdout, 'UNSATISFIABLE')
    disp('unsat');
else
    disp('sat');
end