
function Ck = clingo_ischordal(filename)
Ck = cell(1,1); % Default - Grpah is chordal and no chordless cycles present
% Find the chordless cycles from each graph

% filename = 'kkt_ischordal.lp';
writefile = 'kkt_ischordal.txt';
[status, ~] = system(sprintf('clingo %s 0 > %s', filename, writefile));
[status, cmdout]  = system(sprintf('clingo %s 0 | grep SAT', filename));
if contains(cmdout, 'UNSATISFIABLE')
    disp('unsat - chordless cycles do not exist');
else
    disp('sat - chordless cycles exist');
end

% Parsing ASP output to read chordless cycles:

if(~contains(cmdout, 'UNSATISFIABLE')) % While chordless cycles exist
    fid = fopen(writefile,'rt');
    this_line = fgetl(fid);
    while ischar(this_line)
        cycle = strfind(this_line, 'cycle');
        comma = strfind(this_line, ',');
        end_brace = strfind(this_line, ')');

        if(~isempty(cycle))
            assert(length(cycle) == length(comma) && length(comma) == length(end_brace));   
            edge_col1 = {}; edge_col2 = {};
            for i = 1:length(cycle)
                edge_col1 = [edge_col1, {this_line(cycle(i) + 6 : comma(i) - 1)}];
                edge_col2 = [edge_col2, {this_line(comma(i) + 1 : end_brace(i) - 1)}];
            end

            G = graph(edge_col1, edge_col2);
            H = rmedge(G, G.Edges.EndNodes(1,1), G.Edges.EndNodes(1,2));
            node_path = shortestpath(H, G.Edges.EndNodes(1,1), G.Edges.EndNodes(1,2));
            if(~isempty(node_path))
                Ck{end+1,:} = node_path;
            end
        end

          this_line = fgetl(fid);

    end
end

Ck(1,:) = [];
end