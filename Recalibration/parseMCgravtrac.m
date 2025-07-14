data_nm1 = dlmread('MCgravtrac1.txt');
data_sb1 = dlmread('MCgravtrac1_SDB.txt');
data_me1 = dlmread('MCgravtracEddy.txt');
data_ab1 = dlmread('MCgravtracAnnie.txt');

allData = [data_nm1' data_sb1' data_me1' data_ab1']';
usedMs = allData(:, n.totalMs+1);
allData = allData(:,1:(end-1));
covsys = cov(allData);

%% now calculate measured uncertainties
%First: run AmelinTarantola_v7_All_Mean.m

covMeas = CM - CM*uGn'/(uGn*CM*uGn'+CD)*uGn*CM;  %measured covariance matrix from ML estimate of 981IC, etc.
covTot = covsys + covMeas;  %add in systematic uncertainties (variance of MC seeds)
twosigmaM = real(2*sqrt(diag(covMeas)));
twosigmaT = real(2*sqrt(diag(covTot)));

sysVars = [MC.ics.ET535Pb(usedMs,:) MC.ics.ET2535Pb(usedMs,:) MC.ics.gtU(usedMs,:) MC.ics.gravPb(usedMs,:) ...
            MC.r206238g.ET(usedMs) MC.r206238g.RP(usedMs) MC.r206238g.JMM(usedMs) allData(:,1:2)];
covAll = cov(sysVars);
twosigmaAllVars = 2*sqrt(diag(covAll));
meansAllVars = mean(sysVars)';

%%

ex.r25m = 1.007525;  ex.r25t = sysVars(:,25); 
ex.r65m = 0.382;
ex.r85m = 0.0958;
ex.r35m = 0.9972;  ex.r35t = sysVars(:,11);
ex.r55t = sysVars(:,26);

ex.alphaPb = (ex.r25t ./ ex.r25m - 1)*(-1/3);
ex.alphaU  = (ex.r35t ./ ex.r35m - 1)*(-1/2);
ex.r65fc = ex.r65m.*(1+  ex.alphaPb);
ex.r85fc = ex.r85m.*(1+3*ex.alphaU);
ex.r68spl = ex.r65fc ./ (ex.r85fc .* ex.r55t);
ex.date68 = (1/0.000000000155125)*log(ex.r68spl+1);

2*std(ex.date68)/mean(ex.date68)*10^6

%% Visualize the gammas for Amelin Data

figure
minv = 0.9985; maxv = 1.0005;
set(gca, 'DataAspectRatio', [1 1 1], 'PlotBoxAspectRatio', [1 1 1], 'NextPlot', 'add')
xlim([minv maxv]); ylim([minv maxv]);
line([minv maxv], [minv maxv], [0 0], 'LineWidth', 2, 'Color', [0.7 0.7 0.7])
xlabel(gca, '$\gamma_{205}$', 'FontSize', 24, 'Interpreter', 'latex')
ylabel(gca, '$\gamma_{207}$', 'FontSize', 24, 'Interpreter', 'latex')

xytics = [0.9985 0.9990 0.9995 1 1.0005]';
set(gca,'XTick', xytics, 'YTick', xytics, ...
    'XTickLabel', num2str(xytics, '%1.4f'), 'YTickLabel', num2str(xytics, '%1.4f'));

g5_981v = (g5start+1):(g5start+runcount981);
g5_982v = (g5start+runcount981+1):(g5start+runcount981+runcount982);
g5_Purv = (g5start+runcount981+runcount982+1):g7start;

%plot ellipse centers
plot(um(g5_981v), um(g5_981v+totruns), '.', 'MarkerEdgeColor', [1 0 0],'MarkerSize', 8)
plot(um(g5_982v), um(g5_982v+totruns), '.', 'MarkerEdgeColor', [0 1 0],'MarkerSize', 8)
plot(um(g5_Purv), um(g5_Purv+totruns), '.', 'MarkerEdgeColor', [0 0 1],'MarkerSize', 8)

%plot ellipses
pts = 100;
pis = 0:2*pi/(pts-1):2*pi;
circlepts = 2*[cos(pis') sin(pis')]; %a 2-column matrix of the points in a circle

for i = g5_981v
    j = i+totruns;
    sub_sT = [covTot(i,i) covTot(i,j);
              covTot(j,i) covTot(j,j)];
    scT = chol(sub_sT);
    elptsT = circlepts*scT + repmat([um(i) um(j)], pts, 1);
    plot(elptsT(:,1), elptsT(:,2), 'Color', 'r' , 'LineWidth', 0.2)
end

for i = g5_982v
    j = i+totruns;
    sub_sT = [covTot(i,i) covTot(i,j);
              covTot(j,i) covTot(j,j)];
    scT = chol(sub_sT);
    elptsT = circlepts*scT + repmat([um(i) um(j)], pts, 1);
    plot(elptsT(:,1), elptsT(:,2), 'Color', 'g' , 'LineWidth', 0.2)
end

for i = g5_Purv
    j = i+totruns;
    sub_sT = [covTot(i,i) covTot(i,j);
              covTot(j,i) covTot(j,j)];
    scT = chol(sub_sT);
    elptsT = circlepts*scT + repmat([um(i) um(j)], pts, 1);
    plot(elptsT(:,1), elptsT(:,2), 'Color', 'b' , 'LineWidth', 0.2)
end

%marker at (1,1)
plot(1,1,'+','MarkerSize', 20, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')


%%  do this after running TarantolaUIC_v1.m (which is MC) and then TarantolaUIC_Mean_v1.m

% covsys = cov(ummat_104');  %appended _104 to keep from overwriting
% covmes = uCM - uCM*uGn'/(uGn*uCM*uGn'+uCD)*uGn*uCM;
% covtot = covsys + covmes;
% for i = 1:4
%     for j = 1:4
%         rhotot(i,j) = covtot(i,j)/sqrt(covtot(i,i)*covtot(j,j));
%     end
% end
% twosigmIC = 2*sqrt(diag(covmes(1:2,1:2)));
% twosigtIC = 2*sqrt(diag(covtot(1:2,1:2)));
%% Visualize the gammas for Amelin Data

figure
minv = 0.9975; maxv = 1.0005;
set(gca, 'DataAspectRatio', [1 1 1], 'PlotBoxAspectRatio', [1 1 1], 'NextPlot', 'add')
xlim([minv maxv]); ylim([minv maxv]);
line([minv maxv], [minv maxv], [0 0], 'LineWidth', 2, 'Color', [0.7 0.7 0.7])
xlabel(gca, '$\gamma_{205}$', 'FontSize', 24, 'Interpreter', 'latex')
ylabel(gca, '$\gamma_{207}$', 'FontSize', 24, 'Interpreter', 'latex')

xytics = [0.9975 0.9980 0.9985 0.9990 0.9995 1 1.0005]';
set(gca,'XTick', xytics, 'YTick', xytics, ...
    'XTickLabel', num2str(xytics, '%1.4f'), 'YTickLabel', num2str(xytics, '%1.4f'));

g5_RPv  = (st.g205+1):(st.g205+19);
g5_ETv  = (st.g205+20):(st.g205+33);
g5_JMMv = (st.g205+34):(st.g207);
totruns = 46;

%plot ellipse centers
plot(um(g5_ETv), um(g5_ETv+totruns), '.', 'MarkerEdgeColor', [1 0 0],'MarkerSize', 8)
plot(um(g5_RPv), um(g5_RPv+totruns), '.', 'MarkerEdgeColor', [0 1 0],'MarkerSize', 8)
plot(um(g5_JMMv), um(g5_JMMv+totruns), '.', 'MarkerEdgeColor', [0 0 1],'MarkerSize', 8)

%plot ellipses
pts = 100;
pis = 0:2*pi/(pts-1):2*pi;
circlepts = 2*[cos(pis') sin(pis')]; %a 2-column matrix of the points in a circle

for i = g5_ETv
    j = i+totruns;
    sub_sT = [covTot(i,i) covTot(i,j);
              covTot(j,i) covTot(j,j)];
    scT = chol(sub_sT);
    elptsT = circlepts*scT + repmat([um(i) um(j)], pts, 1);
    plot(elptsT(:,1), elptsT(:,2), 'Color', 'r' , 'LineWidth', 0.2)
end

for i = g5_RPv
    j = i+totruns;
    sub_sT = [covTot(i,i) covTot(i,j);
              covTot(j,i) covTot(j,j)];
    scT = chol(sub_sT);
    elptsT = circlepts*scT + repmat([um(i) um(j)], pts, 1);
    plot(elptsT(:,1), elptsT(:,2), 'Color', 'g' , 'LineWidth', 0.2)
end

for i = g5_JMMv
    j = i+totruns;
    sub_sT = [covTot(i,i) covTot(i,j);
              covTot(j,i) covTot(j,j)];
    scT = chol(sub_sT);
    elptsT = circlepts*scT + repmat([um(i) um(j)], pts, 1);
    plot(elptsT(:,1), elptsT(:,2), 'Color', 'b' , 'LineWidth', 0.2)
end

%marker at (1,1)
plot(1,1,'+','MarkerSize', 20, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')

%% do K-S test on 202/205 and 235/205 data

r25normv = (allData(:,1)-mean(allData(:,1)))./std(allData(:,1));
[ks.h25,ks.p25,ks.ksstat25,ks.cv25] = kstest(r25normv);

r55normv = (allData(:,2)-mean(allData(:,2)))./std(allData(:,2));
[ks.h55,ks.p55,ks.ksstat55,ks.cv55] = kstest(r55normv);

