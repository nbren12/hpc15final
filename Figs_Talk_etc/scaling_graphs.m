close all
Fsize = 26;
Lwide = 8;
cmap1 = colormap(jet);
clr = [6,20,36,44,64];

smbl = ['>','^','o','s','*'];

subplot(2,1,1)
set(gca,'FontSize', Fsize)
for j = 1:5
    plot(log(nt_corse(:,j))/log(2),saving(:,j),'Marker',smbl(j),'Color',cmap1(clr(j),:),'LineWidth', Lwide)
    hold on
end
set(gca, 'XLim', [3 11], 'XTick', [3 4 5 6 7 8 9 10 11],'XTickLabel', {'8';'16'; '32'; '64'; '128'; '256'; '512';'1024';'2048'});
set(gca, 'YLim', [0 100], 'YTick', [0 20 40 60 80 100])

ylabel('Percent time saved')
xlabel('Number of coarse blocks (cores)')
legend(' nt per core = 5',' nt per core = 10',' nt per core = 20',' nt per core = 40',' nt per core = 80','Location','NorthWest')

%figure()
subplot(2,1,2)

set(gca,'FontSize', Fsize)
for j = 1:5
    plot(log(nt_corse(j,:))/log(2),saving(j,:),'Marker',smbl(j),'Color',cmap1(clr(j),:),'LineWidth', Lwide)
    hold on
end
set(gca, 'XLim', [3 11], 'XTick', [3 4 5 6 7 8 9 10 11],'XTickLabel', {'8';'16'; '32'; '64'; '128'; '256'; '512';'1024';'2048'});
set(gca, 'YLim', [0 100], 'YTick', [0 20 40 60 80 100])

ylabel('Percent time saved')
xlabel('Number of coarse blocks (cores)')
legend(' nt total = 640',' nt total = 1280',' nt total = 2560',' nt total = 5120',' nt total = 10240','Location','NorthWest')

figure
set(gca,'FontSize', Fsize)
for j = 1:max(size(rel_tol))
    clear num_iter;
    num_iter = [1:1:length(rel_tol(j).tol)];
    semilogy(num_iter,rel_tol(j).tol,'Marker',smbl(j),'Color',cmap1(clr(j),:),'LineWidth', 4)
    ymax(j) = max(rel_tol(j).tol);
    xmax(j) = length(num_iter);
    hold on
end
ylabel('Relative error')
xlabel('Number of iterations')
%xlim([0 max(xmax)]);
set(gca, 'XLim', [0 max(xmax)], 'XTick', [0:50:(xmax)])
set(gca, 'YLim', [0 max(ymax)], 'YTick', [10 10^10 10^20 10^30 10^40 10^50 10^60])

ylim([0 max(ymax)]);
legend('5 tstep/core on 512 cores','80 tstep/core on 128 cores','10 tstep/core on 64 cores','20 tstep/core on 256 cores','40 tstep/core on 64 cores','Location','NorthEast')
set(gca,'FontSize', Fsize)



