% ---------- Bender's Minimum Fill-in Algorithm ------------ %
% Input: Sparse Symmetric n-by-n matrix M with non-zero diagonal elements
% Output: Permutation Pattern P
% ---------------------------------------------------------- %

function [P, PEO,FILL] = bmf(M)
% Algorithm Parameters
[m,n] = size(M);
assert((m==n) && all(all(M == M')), 'Input a square symmtric matrix');
M(M~=0) = 1;
ADJ = tril(M,-1) + tril(M,-1)';
G = graph(ADJ);
MAX_FILL = n*(n-1)/2;
FILL = [];

% Edge Set and Complement Edge Set
edges = G.Edges.EndNodes; 
N = size(edges,1);
universe = [];
K = 2; m = 0; h = K; iteration = 1; combination = 1:K;
while(iteration > 0)
    [combination, m, h, iteration] = GetNextCombination(N,K,combination, m, h, iteration);
    universe = [universe; combination];
end
edge_comp = setdiff(universe, edges, 'rows');

% Initial Point:
model = struct;
model.A = [];
for i = 1:N
    temp_row = zeros(1, MAX_FILL);
    temp_row(edges(i,1) + (edges(i,2) - 1)*(edges(i,2) - 2)/2) = 1;
    model.A = [model.A; temp_row];
end
z_start = sum(model.A,1);
model.rhs(1:N,:) = 1;
assert(nnz(z_start) == N, 'Initial Solve Incorrect');

% Setting up Gurobi Model:
model.obj = 1-z_start;
model.sense(1:N,:) = ['='];
model.modelsense = 'min';
model.vtype = 'B';
model.A = sparse(model.A);
params.outputflag = 1;

filename = ASP_input(G);
Ck = clingo_ischordal(filename);

while (~isempty(Ck))
    for i = 1:length(Ck)
        C = Ck{i};
        sub_cycle = [sort(cell2mat(cellfun(@(s) str2double(s(2:end)), C, 'uni', false)),2), zeros(1, n-length(C))];
        
        rhs_c = zeros(1, MAX_FILL);
        lhs_c = zeros(1, MAX_FILL);
        
        % Identifying chords, edges and chords
        v_chord = length(C);
        
        E_chord = sort([str2num(C{1}(2:end)), str2num(C{end}(2:end))]);
        if(any(ismember(edge_comp,E_chord(end,:),'rows')))
            rhs_c(E_chord(end,1) + (E_chord(end,2) - 1)*(E_chord(end,2) - 2)/2) = 1;
        end
        for z = 1:v_chord-1
            E_chord = [E_chord; sort([str2num(C{z}(2:end)), str2num(C{z+1}(2:end))])];
            if(any(ismember(edge_comp,E_chord(end,:),'rows')))
                rhs_c(E_chord(end,1) + (E_chord(end,2) - 1)*(E_chord(end,2) - 2)/2) = 1;
            end
        end 
        
        int_chord = [];
        K = 2; m = 0; h = K; iteration = 1; combination = 1:K;
        while(iteration > 0)
            [combination, m, h, iteration] = GetNextCombination(v_chord,K,combination, m, h, iteration);
            diff = abs(combination(1) - combination(2));
            if(diff ~= 1 && diff ~= (v_chord-1)) % Not edges
                int_chord = [int_chord; sort([str2num(C{combination(1)}(2:end)),str2num(C{combination(2)}(2:end))])];
                lhs_c(int_chord(end,1) + (int_chord(end,2)-1)*(int_chord(end,2)-2)/2) = 1;
            end
        end
        assert(length(rhs_c) == length (lhs_c), 'z must be of same length on both sides of the cut')
        assert(length(E_chord) + length(int_chord) == nchoosek(v_chord,2),'No. of edges and no. of chords do not add to choose(|V|,2)');
        
        % Update Constraints for each chordless cycle:
       
        model.A = [model.A; (v_chord - 3)*rhs_c - lhs_c];
        model.rhs = [model.rhs; (v_chord - 3)*(nnz(rhs_c) - 1)]; % nnz(rhs_c) = no. of edges common to E_chord and edge_comp
        % nnz(rhs_c)  = size(intersect(E_chord, edge_comp, 'rows'),1);
        model.sense = [model.sense; '<'];
    end
    
    %Solve Gurobi Model with new constraints:
    model.A = sparse(model.A);
    assert(issparse(model.A) && ~issparse(model.rhs));
    result = gurobi(model, params); % Solution z^k with all accumulated Benders cuts
    % Getting Fill edges:
    fill_loc = triu(ones(n),1);
    fill_loc(fill_loc == 1) = result.x;
    [fill_row, fill_col] = find(fill_loc);
    FILL = unique([FILL; setdiff([fill_row, fill_col], edges, 'rows')], 'rows');
    assert(length(FILL) < MAX_FILL);
    
    % Generating New Chordal Cycles:
    filename = ASP_input(addedge(G, FILL(:,1)',FILL(:,2)',ones(1,size(FILL,1))));
    Ck = clingo_ischordal(filename);
end

% Finding Perfect Elimination Ordering using Lexicographic Breadth First Search:
if(~isempty(FILL))
    chordal_graph = addedge(G, FILL(:,1)', FILL(:,2)',ones(1, size(FILL,1)));
else
    chordal_graph = G;
end

label = strings(1,n); % All vertices are unnumbered
label(:) = '0';
elim_order(1:n) = 0;
ii = n;
while(ii >= 1)
    % Finding the largest lexicographic unnumbered vertex:
    idx = find(elim_order == 0); % Unnumbered vertices
    % Finding lexicographically largest label:
    lex_label = str2double(label(idx));
    y = find(lex_label == max(lex_label));  
    elim_order(idx(y(1))) = n+1-ii; % Picking  alexicographically largest label
    
    nbrs = neighbors(chordal_graph, idx(y(1))); % Finding all neighbors
    for kk = 1:length(nbrs)
        if(elim_order(nbrs(kk)) == 0) % If un-numbered
            label(nbrs(kk)) = strcat(label(nbrs(kk)),num2str(ii));
        end
    end
    ii = ii-1;
end
PEO = zeros(1, n);
for i = 1:n
    PEO(elim_order(i)) = i;
end
PEO = flip(PEO);
ID = speye(n);
P = ID(PEO, :); 
end