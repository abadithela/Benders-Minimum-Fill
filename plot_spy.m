% To plot spy figures of L overlapped on M
% L is plotted as red circles, M is plotted as blue circles
function chk = plot_spy(M, L)
M_reg = M + size(M,1)*eye(size(M,1));

figure(101)
subplot(121)
hold on
for i = 1:size(M,1)
    for j = 1:i
        if(M(i,j)~=0)
            plot(j, size(M,1) - i, 'ob', 'MarkerFaceColor', 'b');
            plot(i, size(M,1) - j, 'ob', 'MarkerFaceColor', 'b');
        end
    end
end
xticks([0:10:40])
set(gca,'YTickLabel',flipud(get(gca,'YTickLabel')));
xlabel(['nz = ', num2str(nnz(M_reg))]);
title('KKT sparsity');

subplot(122)
hold on
for i = 1:size(M,1)
    for j = 1:i
        if(L(i,j)~=0)
            plot(j, size(M,1) - i, 'or', 'MarkerFaceColor', 'r');           
        end
        if(M_reg(i,j)~=0)
            plot(j, size(M,1) - i, 'ob', 'MarkerFaceColor', 'b');
        end
    end
end
xticks([0:10:40])
set(gca,'YTickLabel',flipud(get(gca,'YTickLabel')));
xlabel(['nz = ', num2str(nnz(L))]);
title('L sparsity');
chk = 1;
end