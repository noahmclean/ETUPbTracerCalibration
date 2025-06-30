data_nm1 = dlmread('MCintercal_nm1.txt');
data_nm12 = dlmread('MCintercal_nm1-4.txt');
%data_sb1 = dlmread('MCintercal_Seth1.txt');
%data_sb22 = dlmread('MCintercal_Seth2-2.txt');
data_sb23 = dlmread('MCintercal_Seth2-4.txt');

data1short = data_nm1(:,1:20)';
data2short = data_nm12(:,1:20)';
data3short = data_sb23(:,1:20)';
allDatashort = [data1short data2short data3short]';
allData = [data_nm1(:,1:(end-1))' data_nm12(:,1:(end-1))' data_sb23(:,1:(end-1))']';
%clear data*

covsys = cov(allData);
rhosys = zeros(size(covsys));
for i = 1:20
    for j = 1:20
        rhosys(i,j) = covsys(i,j)/sqrt(covsys(i,i)*covsys(j,j));
    end
end

%% now calculate measured uncertainties
%First: run AmelinTarantola_v7_All_Mean.m

covMeas = uCM - uCM*uGn'/(uGn*uCM*uGn'+diag(uCDd))*uGn*uCM;  %measured covariance matrix from ML estimate of 981IC, etc.
covTot = covsys + covMeas;  %add in systematic uncertainties (variance of MC seeds)
twosigmaM = real(2*sqrt(diag(covMeas)));
twosigmaT = real(2*sqrt(diag(covTot)));
%% 

rhotot = zeros(20,20);
for i = 1:20
    for j = 1:20
        rhotot(i,j) = covTot(i,j)/sqrt(covTot(i,i)*covTot(j,j));
    end
end

%% Visualize the gammas

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
    %measured covariance ellipse
%     sub_sM = [covMeas(i,i) covMeas(i,j);
%               covMeas(j,i) covMeas(j,j)];
%     scM = chol(sub_sM);
%     elptsM = circlepts*scM + repmat([um(i) um(j)], pts, 1);
%     plot(elptsM(:,1), elptsM(:,2), 'r', 'LineWidth', 0.2)
    %total covariance ellipse
    sub_sT = [covTot(i,i) covTot(i,j);
              covTot(j,i) covTot(j,j)];
    scT = chol(sub_sT);
    elptsT = circlepts*scT + repmat([um(i) um(j)], pts, 1);
    plot(elptsT(:,1), elptsT(:,2), 'Color', 'r' , 'LineWidth', 0.2)
end

for i = g5_982v
    j = i+totruns;
    %measured covariance ellipse
%     sub_sM = [covMeas(i,i) covMeas(i,j);
%               covMeas(j,i) covMeas(j,j)];
%     scM = chol(sub_sM);
%     elptsM = circlepts*scM + repmat([um(i) um(j)], pts, 1);
%     plot(elptsM(:,1), elptsM(:,2), 'r', 'LineWidth', 0.2)
    %total covariance ellipse
    sub_sT = [covTot(i,i) covTot(i,j);
              covTot(j,i) covTot(j,j)];
    scT = chol(sub_sT);
    elptsT = circlepts*scT + repmat([um(i) um(j)], pts, 1);
    plot(elptsT(:,1), elptsT(:,2), 'Color', 'g' , 'LineWidth', 0.2)
end

for i = g5_Purv
    j = i+totruns;
    %measured covariance ellipse
%     sub_sM = [covMeas(i,i) covMeas(i,j);
%               covMeas(j,i) covMeas(j,j)];
%     scM = chol(sub_sM);
%     elptsM = circlepts*scM + repmat([um(i) um(j)], pts, 1);
%     plot(elptsM(:,1), elptsM(:,2), 'r', 'LineWidth', 0.2)
    %total covariance ellipse
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