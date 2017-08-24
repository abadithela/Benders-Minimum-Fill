% ------- ASP_clingo -------%
% Apurva Badithela
% Input to ASP to be solved by Clingo

% M = load('kkt.mat');
% M = kkt;
% M(M~=0) = 1;
% M = tril(M,-1) + tril(M,-1)';
% G = graph(M);

% Writes Graph to file to be fed in by Clingo
function filename = ASP_input(G)
    edges = string(G.Edges.EndNodes);
    var_names = mat2cell(strcat('g', edges), ones(size(edges,1),1), 2);

    [nrow, ncol] = size(var_names);
    fwriteID = fopen('kkt_ischordal.lp','w');
    formatSpec = 'edge(%s , %s).\n';
    for row = 1:nrow
        fprintf(fwriteID, formatSpec, var_names{row,:});
    end
    fprintf(fwriteID,'adj(X,Y) :- edge(X,Y).\n');
    fprintf(fwriteID,'adj(X,Y) :- edge(Y,X).\n');

    fprintf(fwriteID,'4 {cycle(X,Y) : edge(X,Y)}.\n');
    fprintf(fwriteID,'adjc(X,Y) :- cycle(X,Y).\n');
    fprintf(fwriteID,'adjc(X,Y) :- cycle(Y,X).\n');
    fprintf(fwriteID,'vertex(X) :- adjc(X,_).\n');

    fprintf(fwriteID,':- vertex(X), #count{Y : adjc(X,Y)}!=2.\n');

    fprintf(fwriteID,'r(X,Y) :- adjc(X,Y).\n');
    fprintf(fwriteID,'r(X,Z) :- adjc(X,Y), r(Y,Z).\n');
    fprintf(fwriteID,':- vertex(X), vertex(Y), not r(X,Y).\n');

    fprintf(fwriteID,':- vertex(X), vertex(Y), not adjc(X,Y), adj(X,Y).\n');
    fprintf(fwriteID,'#show cycle/2.\n');

    fclose(fwriteID);
    filename = 'kkt_ischordal.lp'; % Successfully completed writing to file
end