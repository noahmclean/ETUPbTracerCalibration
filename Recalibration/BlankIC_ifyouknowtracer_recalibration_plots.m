%% Originally found at the end of BlankIC_ifyouknowtracer script
% Original script placed plot in new figure. This script has been modified
% to place the plots in a uitab.

%% plot data with color

alphaVal = 0.8;
xpadbuf = 0.05;
xwidthpx = 1200; yheightpx = 650;
%hfig = figure('Position', [15 15 xwidthpx yheightpx]);
tabMcLeanFigure6 = uitab(tabGroupMcLean, "Title", "Figure 6");
tabMcLeanFigure6.AutoResizeChildren = 'off';
ax(1) = subplot(1,2,1, 'FontSize', 12, 'Parent', tabMcLeanFigure6);
hold(ax(1), 'on');

stcov = cov(dataratios'); %covariance matrix for discrete blank IC data

scatter(ax(1), dataratios(1,:), dataratios(2,:),'.','k')
pts = 200;  %number of points in ellipses
pis = 0:2*pi/(pts-1):2*pi;

%set(gca, 'DataAspectRatio', [1 1 1], 'DataAspectRatioMode', 'manual')
circlepts = 2*[cos(pis') sin(pis')];   %a 2-column matrix of the points in a circle
%plot(circlepts(:,1), circlepts(:,2))

for i = 1:n %plot measred data ellipses
    s = covmatsrat(1:2,1:2,i);
    sc = chol(s);
    elpts = circlepts*sc + repmat([dataratios(1,i) dataratios(2,i)], pts, 1);
    fill(ax(1), elpts(:,1), elpts(:,2), colorV(i,:), ...
        'LineWidth', 0.2, 'FaceAlpha', alphaVal)
    
end

xsc = chol(xs(1:2, 1:2));  
bigel = circlepts*xsc + repmat([mu(1) mu(2)], pts, 1);
plot(ax(1), bigel(:,1), bigel(:,2), 'g', 'LineWidth', 2) %plot 'overdispersion' ellipse
stcovc = chol(stcov(1:2,1:2));
stdel = circlepts*stcovc + repmat([mu(1) mu(2)], pts, 1);
plot(ax(1), stdel(:,1), stdel(:,2), 'r', 'LineWidth', 2, 'LineStyle', '--') %plot standard discrete ellipse
xlabel(ax(1), '$^{206}$ Pb/$^{204}$ Pb', 'FontSize', 14, 'Interpreter', 'latex')
ylabel(ax(1), '$^{207}$ Pb/$^{204}$ Pb', 'FontSize', 14, 'Interpreter', 'latex')
xlabh = get(ax(1),'XLabel'); set(xlabh, 'Units', 'normalized')
set(xlabh,'Position',get(xlabh,'Position') - [0 xpadbuf 0])

ax(2) = subplot(1,2,2, 'FontSize', 12, 'Parent', tabMcLeanFigure6);  %for 207/204 - 208/204

hold(ax(2), 'on');
scatter(ax(2), dataratios(2,:), dataratios(3,:),'.','k')
pts = 200;
pis = 0:2*pi/(pts-1):2*pi;

set(ax(2), 'DataAspectRatioMode', 'auto')
circlepts = 2*[cos(pis') sin(pis')];
%plot(circlepts(:,1), circlepts(:,2))

for i = 1:n
    s = covmatsrat(2:3,2:3,i);
    sc = chol(s);
    elpts = circlepts*sc + repmat([dataratios(2,i) dataratios(3,i)], pts, 1);
        fill(ax(2), elpts(:,1), elpts(:,2), colorV(i,:), ...
        'LineWidth', 0.2, 'FaceAlpha', alphaVal)
end

xsc = chol(xs(2:3, 2:3));
bigel = circlepts*xsc + repmat([mu(2) mu(3)], pts, 1);
plot(ax(2), bigel(:,1), bigel(:,2), 'g', 'LineWidth', 2)
stcovc = chol(stcov(2:3,2:3));
stdel = circlepts*stcovc + repmat([mu(2) mu(3)], pts, 1);
plot(ax(2), stdel(:,1), stdel(:,2), 'r', 'LineWidth', 2, 'LineStyle', '--')
xlabel(ax(2), '$^{207}$ Pb/$^{204}$ Pb', 'FontSize', 14, 'Interpreter', 'latex')
xlabh2 = get(ax(2),'XLabel'); set(xlabh2, 'Units', 'normalized')
set(xlabh2,'Position', get(xlabh2,'Position') - [0 xpadbuf 0])
ylabel(ax(2), '$^{208}$ Pb/$^{204}$ Pb', 'FontSize', 14, 'Interpreter', 'latex')

% for colorbar
barwidth = 0.03;
hcb = colorbar('Location','SouthOutside', 'Parent', tabMcLeanFigure6);
set(hcb, 'Position', [0.10 0.10 0.815 barwidth])
for i = 1:2
    pos = get(ax(i), 'Position');
    set(ax(i), 'PlotBoxAspectRatio', [1 1 1], ...
               'Position', [0.45*(i-1) 0.3 0.6 0.6]);
end
%colormap(jet(ncolor));
set(hcb,'XTickMode','manual', 'XLim', [0 1], 'XTick', xticLoc, 'XTickLabel', xticLab, ...
        'Box', 'off', 'TickLength', barwidth/4)
set(get(hcb, 'title'), 'string', 'Pb blank (pg)', 'FontSize', 12)
cbIm = findobj(hcb,'Type','image');
alpha(cbIm,alphaVal)
%myaa(4)
%plot2svg('ETICblanks.svg')
