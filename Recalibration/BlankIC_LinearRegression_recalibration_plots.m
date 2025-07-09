%% This code originally at the end of the script BlankIC_LinearRegression_recalibration

% plot data
% close all % remove this to allow multiple plots

%fig1 = figure('Position', [0 500 1000 1500]);
%fig1 = figure;

%stcov = cov(dataratios'); %covariance matrix for discrete blank IC data
%new:
%subplot(3,1,1,'NextPlot','add')


tabMcLeanFigure5 = uitab(tabGroupMcLean, "Title", "Figure 5");
tabMcLeanFigure5.AutoResizeChildren = 'off';
ax(1) = subplot(3,1,1, 'FontSize', 12, 'Parent', tabMcLeanFigure5);
hold(ax(1), 'on');


%hold on
scatter(ax(1), dataratios(1,:), dataratios(2,:),'.','k')
pts = 200;  %number of points in ellipses
pis = 0:2*pi/(pts-1):2*pi;

%set(gca, 'DataAspectRatio', [1 1 1], 'DataAspectRatioMode', 'manual')
circlepts = 2*[cos(pis') sin(pis')];   %a 2-column matrix of the points in a circle
%plot(circlepts(:,1), circlepts(:,2))

for i = 1:n %plot measred data ellipses
    s = covmats205(1:2,1:2,i);
    sc = chol(s);
    elpts = circlepts*sc + repmat([dataratios(1,i) dataratios(2,i)], pts, 1);
    plot(ax(1), elpts(:,1), elpts(:,2), 'b', 'LineWidth', 1)
end
xlabel(ax(1), '$^{204}$Pb/$^{205}$Pb', 'FontSize', 14, 'Interpreter', 'latex')
ylabel(ax(1), '$^{206}$Pb/$^{205}$Pb', 'FontSize', 14, 'Interpreter', 'latex')

xmin = 0; xmax = ceil( max(dataratios(1,:)*10^4))/10^4;
ymin = 0; ymax = ceil( max(dataratios(2,:)*10^3))/10^3;
set(ax(1), 'XLim',[xmin xmax], 'YLim', [ymin ymax], 'DataAspectRatio', [1 r64b 1])
linxmin = trIC45; linymin = trIC65;
linxmax = min( [xmax (ymax-a(2))/v(2)]); linymax = a(2)+v(2)*linxmax;
line(ax(1), [linxmin linxmax],[linymin linymax])
plot(ax(1), trIC45,trIC65,'+','MarkerFaceColor','k')

% make uncertainty envelope
%need the slope and intercept, plus the covariance matrix
steps = 200; tstep = (xmax-xmin)/(steps-1); t = xmin:tstep:xmax;
a1 = trIC45; v1 = 1;
a2 = trIC65; v2 = r64b;
cov456 = covtrbl_ET535; cov456([3 4 6 7],:) = []; cov456(:,[3 4 6 7]) = [];
dxda1 = 1; dxda2 = 0; dxdv2 = 0;
dyda1 = 0; dyda2 = 1; %dydv2 = t;
dxdt = 1; dydt = r64b;
vperp = [-dydt dxdt];

count = 1;
delx = zeros(size(t)); dely = zeros(size(t));
for ti = t
    dydv2 = ti;
    Jxyab = [dxda1 dxda2 dxdv2;
             dyda1 dyda2 dydv2];
    s2perp = vperp*Jxyab*cov456*Jxyab'*vperp' / (vperp*vperp');
    delx(count) = 2*cos(atan(-dxdt/dydt))*sqrt(s2perp);
    dely(count) = 2*sin(atan(-dxdt/dydt))*sqrt(s2perp);
    count = count + 1;
end

plot(ax(1), a1 + v1*t + delx, a2 + v2*t + dely, '--g')
plot(ax(1), a1 + v1*t - delx, a2 + v2*t - dely, '--g')


%%

% xsc = chol(xs(1:2, 1:2));  
% bigel = circlepts*xsc + repmat([mu(1) mu(2)], pts, 1);
% plot(bigel(:,1), bigel(:,2), 'g', 'LineWidth', 2) %plot 'overdispersion' ellipse
% stcovc = chol(stcov(1:2,1:2));
% stdel = circlepts*stcovc + repmat([mu(1) mu(2)], pts, 1);
% plot(stdel(:,1), stdel(:,2), 'r', 'LineWidth', 2, 'LineStyle', '--') %plot standard discrete ellipse

%hold off

%figure  %for 207/204 - 208/204

%hold on
ax(2) = subplot(3,1,2,'NextPlot','add', 'FontSize', 12, 'Parent', tabMcLeanFigure5);

scatter(ax(2), dataratios(2,:), dataratios(3,:),'.','k')
pts = 200;
pis = 0:2*pi/(pts-1):2*pi;

hold(ax(2), 'on');
set(ax(2), 'DataAspectRatioMode', 'auto')
circlepts = 2*[cos(pis') sin(pis')];
%plot(circlepts(:,1), circlepts(:,2))

for i = 1:n
    s = covmats205(2:3,2:3,i);
    sc = chol(s);
    elpts = circlepts*sc + repmat([dataratios(2,i) dataratios(3,i)], pts, 1);
    plot(ax(2), elpts(:,1), elpts(:,2), 'b', 'LineWidth', 1)
end

xlabel(ax(2), '$^{206}$Pb/$^{205}$Pb', 'FontSize', 14, 'Interpreter', 'latex')
ylabel(ax(2), '$^{207}$Pb/$^{205}$Pb', 'FontSize', 14, 'Interpreter', 'latex')

xmin = 0; xmax = ceil( max(dataratios(2,:)*10^3))/10^3;
ymin = 0; ymax = ceil( max(dataratios(3,:)*10^3))/10^3;
set(ax(2), 'XLim',[xmin xmax], 'YLim', [ymin ymax], 'DataAspectRatio', [1 r74b/r64b 1])
linxmin = trIC65; linymin = trIC75;
linxmax = min( [xmax (ymax-a(3))*(v(2)/v(3))+a(2)]); linymax = a(3)+v(3)*(linxmax-a(2))/v(2);
line(ax(2), [linxmin linxmax],[linymin linymax])
plot(ax(2), trIC65,trIC75,'+','MarkerFaceColor','k')

% make uncertainty envelope
%need the slope and intercept, plus the covariance matrix
a2 = trIC65; v2 = r64b;
a3 = trIC75; v3 = r74b;
cov567 = covtrbl_ET535; cov567([1 4 7],:) = []; cov567(:,[1 4 7]) = [];
dxda2 = 1; dxda3 = 0;             dxdv3 = 0; %dxdv2 = t;
dyda2 = 0; dyda3 = 1; dydv2 = 0; %dydv3 = t;
dxdt = r64b; dydt = r74b;
vperp = [-dydt dxdt];

count = 1;
delx = zeros(size(t)); dely = zeros(size(t));
for ti = t
    dxdv2 = ti;
    dydv3 = ti;
    Jxyab = [dxda2 dxda3 dxdv2 dxdv3;
             dyda2 dyda3 dydv2 dydv3];
    s2perp = vperp*Jxyab*cov567*Jxyab'*vperp' / (vperp*vperp');
    delx(count) = 2*cos(atan(-dxdt/dydt))*sqrt(s2perp);
    dely(count) = 2*sin(atan(-dxdt/dydt))*sqrt(s2perp);
    count = count + 1;
end

plot(ax(2), a2 + v2*t + delx, a3 + v3*t + dely, '--g')
plot(ax(2), a2 + v2*t - delx, a3 + v3*t - dely, '--g')


%%

% xsc = chol(xs(2:3, 2:3));
% bigel = circlepts*xsc + repmat([mu(2) mu(3)], pts, 1);
% plot(bigel(:,1), bigel(:,2), 'g', 'LineWidth', 2)
% stcovc = chol(stcov(2:3,2:3));
% stdel = circlepts*stcovc + repmat([mu(2) mu(3)], pts, 1);
% plot(stdel(:,1), stdel(:,2), 'r', 'LineWidth', 2, 'LineStyle', '--')

%figure  %for 207/204 - 208/204

%hold on
ax(3) = subplot(3,1,3,'NextPlot','add', 'FontSize', 12, 'Parent', tabMcLeanFigure5);

scatter(ax(3), dataratios(3,:), dataratios(4,:),'.','k')
pts = 200;
pis = 0:2*pi/(pts-1):2*pi;

%hold on
set(ax(3), 'DataAspectRatioMode', 'auto')
circlepts = 2*[cos(pis') sin(pis')];
%plot(circlepts(:,1), circlepts(:,2))

for i = 1:n
    s = covmats205(3:4,3:4,i);
    sc = chol(s);
    elpts = circlepts*sc + repmat([dataratios(3,i) dataratios(4,i)], pts, 1);
    plot(ax(3), elpts(:,1), elpts(:,2), 'b', 'LineWidth', 1)
end

xlabel(ax(3), '$^{207}$Pb/$^{205}$Pb', 'FontSize', 14, 'Interpreter', 'latex')
ylabel(ax(3), '$^{208}$Pb/$^{205}$Pb', 'FontSize', 14, 'Interpreter', 'latex')
xmin = 0; xmax = ceil( max(dataratios(3,:)*10^3))/10^3;
ymin = 0; ymax = ceil( max(dataratios(4,:)*10^3))/10^3;
set(ax(3), 'XLim',[xmin xmax], 'YLim', [ymin ymax], 'DataAspectRatio', [1 r84b/r74b 1])
linxmin = trIC75; linymin = trIC85;
linxmax = min( [xmax (ymax-a(4))*(v(3)/v(4))+a(3)]); linymax = a(4)+v(4)*(linxmax-a(3))/v(3);
line(ax(3), [linxmin linxmax],[linymin linymax])
plot(ax(3), trIC75,trIC85,'+','MarkerFaceColor','k')

% make uncertainty envelope
%need the slope and intercept, plus the covariance matrix
a3 = trIC75; v3 = r74b;
a4 = trIC85; v4 = r84b;
cov578 = covtrbl_ET535; cov578([1 2 5],:) = []; cov578(:,[1 2 5]) = [];
dxda3 = 1; dxda4 = 0;             dxdv4 = 0; %dxdv3 = t;
dyda3 = 0; dyda4 = 1; dydv3 = 0; %dydv4 = t;
dxdt = r74b; dydt = r84b;
vperp = [-dydt dxdt];

count = 1;
delx = zeros(size(t)); dely = zeros(size(t));
for ti = t
    dxdv3 = ti;
    dydv4 = ti;
    Jxyab = [dxda3 dxda4 dxdv3 dxdv4;
             dyda3 dyda4 dydv3 dydv4];
    s2perp = vperp*Jxyab*cov578*Jxyab'*vperp' / (vperp*vperp');
    delx(count) = 2*cos(atan(-dxdt/dydt))*sqrt(s2perp);
    dely(count) = 2*sin(atan(-dxdt/dydt))*sqrt(s2perp);
    count = count + 1;
end

plot(ax(3), a3 + v3*t + delx, a4 + v4*t + dely, '--g')
plot(ax(3), a3 + v3*t - delx, a4 + v4*t - dely, '--g')
