%% Calculate mean and variablility of blank IC, assuming you know tracer IC

%% Load in blanks data, assume measured correlation coefficients are zero

%assume ET535 v.3.0 tracer IC
conc205t = 9.884*10^-12;
r45t = 9.000000000000000e-05;
r65t = 3.887398722782950e-04;
r75t = 2.960685875948440e-04;
r85t = 7.443761622554650e-04;

% assume ET2535 v.3.0 tracer IC
% conc205t = 1.03116*10^-11;
% r45t = 0.000105;
% r65t = 0.000482509449698359;
% r75t = 0.000432369168809505;
% r85t = 0.00104222718281;

% r45t = trIC45; % to use if doing blank/tracer IC from 1st principles
% r65t = trIC65;
% r75t = trIC75;
% r85t = trIC85;

% ***BlanksData structure:
% col 1      2       3       4        5         6       7        8        9         10       11 
% alphaPb 1s(abs) r204/205m 1s(pct) r206/205m 1s(pct) r207/205m 1s(pct) r208/205m 1s(pct) tracerMass

%note: don't forget to use the transpose here:

data = ET535LoadingBlanks';

n = size(data,2);  %number of blanks measured

covmats205 = zeros(4,4,n);   %pre-allocate for uncertainty propagation
covmatsmolb = zeros(4,4,n);
covmatsrat = zeros(3,3,n);
dataratios  = zeros(3,n);    %fractionation and tracer-corrected points
surfaceArea = zeros(n,1);
blankMasses = zeros(n,1);
for i = 1:n
   alphaPb = data(1,i);         %from Excel sheet
   s_alphaPb = data(2,i);       %not corrected for beam interpolation
   r45m = data(3,i);   
   s_r45m = data(4,i)*r45m/100; %unct. entered as 1-sigma pct
   r65m = data(5,i);   
   s_r65m = data(6,i)*r65m/100;            
   r75m = data(7,i);
   s_r75m = data(8,i)*r75m/100;
   r85m = data(9,i);
   s_r85m = data(10,i)*r85m/100;
   tracerMass = data(11,i);
   
   r45fc = r45m*(1-  alphaPb);
   r65fc = r65m*(1+  alphaPb);   %with 205 in denominator
   r75fc = r75m*(1+2*alphaPb);
   r85fc = r85m*(1+3*alphaPb);
   
   % calculate covariance matrix and jacobians

   %covariance matrix for input measured ratios, assuming measured ratios 
   %to 205 are independent
   covi = diag([s_r45m^2 s_r65m^2 s_r75m^2 s_r85m^2 s_alphaPb^2]);  

   %jacobian for transformation to fractionation-corrected values
   j1 = [diag([1 - alphaPb 1+alphaPb 1+2*alphaPb 1+3*alphaPb])...
       [-r45m r65m 2*r75m 3*r85m]'];

   covmats205(:,:,i) =    j1*covi*j1';

   molPb205t = conc205t * tracerMass;   %tracermass;   %the 1 should be tracermass, but divides out!
   molPb204b = molPb205t*(r45fc - r45t);
   molPb206b = molPb205t*(r65fc - r65t);
   molPb207b = molPb205t*(r75fc - r75t);
   molPb208b = molPb205t*(r85fc - r85t);
   blankMasses(i) = (molPb204b*204 + molPb206b*206+molPb207b*207+molPb208b*208)*10^12;
   
   j2 = diag(molPb205t*ones(1,4));
   covmatsmolb(:,:,i) = j2*covmats205(:,:,i)*j2';
   
   r64b = molPb206b/molPb204b;
   r74b = molPb207b/molPb204b;
   r84b = molPb208b/molPb204b;
   
   j3 = [-molPb206b/(molPb204b^2) 1/molPb204b 0 0;...
         -molPb207b/(molPb204b^2) 0 1/molPb204b 0;...
         -molPb208b/(molPb204b^2) 0 0 1/molPb204b];
   covmatsrat(:,:,i) = j3*covmatsmolb(:,:,i)*j3';   %covariance matrix for 206/204b, 207/204b, 208/204b
   
   dataratios(:,i) = [r64b r74b r84b]';
   
   surfaceArea(i) = det(covmatsrat(:,:,i));
   
end

[~, areaSort] = sort(surfaceArea, 'descend');
dataratios = dataratios(:, areaSort);
covmatsrat = covmatsrat(:,:,areaSort);
blankMasses = blankMasses(areaSort);

%% color bar mapping
dimsVec = blankMasses;

xwidthpx = 1200; 

nTryTics = floor(xwidthpx/200);
dimsDist = max(dimsVec)-min(dimsVec);
dTarg = dimsDist/nTryTics;
duTarg = floor(log10(dTarg));
cnt = 0;
trydiff = zeros(9,3);
for i = duTarg-1:duTarg+1
    trydiff(cnt+1,1) = 1; trydiff(cnt+1,2) = i; trydiff(cnt+1,3) = floor(dimsDist/(1*10^i));
    trydiff(cnt+2,1) = 2; trydiff(cnt+2,2) = i; trydiff(cnt+2,3) = floor(dimsDist/(2*10^i));
    trydiff(cnt+3,1) = 5; trydiff(cnt+3,2) = i; trydiff(cnt+3,3) = floor(dimsDist/(5*10^i));
    cnt = cnt + 3;
end
[~, intIndx] = min(abs(trydiff(:,3)-nTryTics));
ticdiff = trydiff(intIndx,1)*10^trydiff(intIndx,2);
ticstart = min(dimsVec) - mod(min(dimsVec), ticdiff) + ticdiff; %first tic after min
ticend  = max(dimsVec) - mod(max(dimsVec), ticdiff); %last tic before max
ntics = round((ticend-ticstart)/ticdiff) + 1;
ticvec = ticstart:ticdiff:ticend;
xticLoc = (1-0)/(max(dimsVec)-min(dimsVec))*(ticvec - min(dimsVec)) + 0;

formSpec = ['%8.' num2str(max(0,-trydiff(intIndx,2)),1) 'f'];
xticLab = cell(ntics,1);
for i = 1:ntics
    xticLab{i} = num2str(ticvec(i), formSpec);
end


%% color mapping

ncolor = 256; %colorvec = 1:ncolor (step 1)
jetmat = jet(ncolor);
colormap(jet(ncolor));

xcolorOfDim = (ncolor-1)/(max(dimsVec)-min(dimsVec))*(dimsVec-min(dimsVec))+1;
colorV = ones(n,3);
for i = 1:n
    colorV(i,:) = [interp1(1:ncolor, jetmat(:,1), xcolorOfDim(i),'pchip') ...
                 interp1(1:ncolor, jetmat(:,2), xcolorOfDim(i),'pchip') ...
                 interp1(1:ncolor, jetmat(:,3), xcolorOfDim(i),'pchip')];
end

%%

xs = cov(dataratios');
% calculate weighted mean with only random errors as first guess
covmatsrat_inv = zeros(3,3,n);
wtdpts = zeros(3,n);
for i = 1:n  %weighted mean, with overdispersion factored in
    covmatsrat_inv(:,:,i) = inv(covmatsrat(:,:,i)+xs);
    %wtdpts(:,i) = covmatsrat_inv(:,:,i)*dataratios(:,i);
    wtdpts(:,i) = (covmatsrat(:,:,i)+xs)\dataratios(:,i);
end  %make inverse covariance ratios with overdispersion addded
m1 = sum(covmatsrat_inv,3);
m2 = sum(wtdpts,2);
mu = m1\m2;

smallS = covmatsrat;  %these are the measurement uncertainties to throw the solver
datameas = dataratios; %these are the fractionation and tracer-corrected blank ratios

%% vermeeshcopy3 uses 
%oldopts = optimset('fsolve');

xs = xs*0.9;

for resoloops = 1:15   %five times is good enough for convergence
options = optimset('TolX',10^-20, 'TolFun',10^-20,...
    'MaxFunEvals',10^12,'Diagnostics','on', 'MaxIter',50, 'Display', 'on');
xs = fsolve(@(S) vermeeshcopy3(S, n , mu, smallS, datameas), xs, options);
%%%%%%%%%%%%xs = fsolve(@(S) vermeeshcopy3(S, n ,mu, smallS, datameas), cov(dataratios')/2000, options);


    for i = 1:n  %weighted mean, with overdispersion factored in
        covmatsrat_inv(:,:,i) = inv(covmatsrat(:,:,i)+xs);
        wtdpts(:,i) = covmatsrat_inv(:,:,i)*dataratios(:,i);
    end  %make inverse covariance ratios with overdispersion addded
    m1 = sum(covmatsrat_inv,3);
    m2 = sum(wtdpts,2);
    mu = m1\m2;

end

%%  make data readable
r64b = mu(1); r74b = mu(2); r84b = mu(3);
s64b = sqrt(xs(1,1));
s74b = sqrt(xs(2,2));
s84b = sqrt(xs(3,3));
rho6474 = xs(1,2)/sqrt(s64b*s74b);
rho6484 = xs(1,3)/sqrt(s64b*s84b);
rho7484 = xs(2,3)/sqrt(s74b*s84b);

disp(['mean ', num2str(r64b,'%2.6f'), ' ', num2str(r74b,'%2.6f'), ' ', num2str(r84b,'%2.6f')])
disp(['1std  ', num2str(1*s64b,'%2.6f'), '  ', num2str(1*s74b,'%2.6f'), '  ', num2str(1*s84b,'%2.6f')])
disp(['rhos  ', num2str(rho6474,'%2.6f'), '  ', num2str(rho6484,'%2.6f'), '  ', num2str(rho7484,'%2.6f')])
disp(' ')

%% plot data with color

alphaVal = 0.8;
xpadbuf = 0.05;
xwidthpx = 1200; yheightpx = 650;
hfig = figure('Position', [15 15 xwidthpx yheightpx]);

ax(1) = subplot(1,2,1, 'FontSize', 12);

stcov = cov(dataratios'); %covariance matrix for discrete blank IC data

hold on
scatter(dataratios(1,:), dataratios(2,:),'.','k')
pts = 200;  %number of points in ellipses
pis = 0:2*pi/(pts-1):2*pi;

%set(gca, 'DataAspectRatio', [1 1 1], 'DataAspectRatioMode', 'manual')
circlepts = 2*[cos(pis') sin(pis')];   %a 2-column matrix of the points in a circle
%plot(circlepts(:,1), circlepts(:,2))

for i = 1:n %plot measred data ellipses
    s = covmatsrat(1:2,1:2,i);
    sc = chol(s);
    elpts = circlepts*sc + repmat([dataratios(1,i) dataratios(2,i)], pts, 1);
    fill(elpts(:,1), elpts(:,2), colorV(i,:), ...
        'LineWidth', 0.2, 'FaceAlpha', alphaVal)
    
end

xsc = chol(xs(1:2, 1:2));  
bigel = circlepts*xsc + repmat([mu(1) mu(2)], pts, 1);
plot(bigel(:,1), bigel(:,2), 'g', 'LineWidth', 2) %plot 'overdispersion' ellipse
stcovc = chol(stcov(1:2,1:2));
stdel = circlepts*stcovc + repmat([mu(1) mu(2)], pts, 1);
plot(stdel(:,1), stdel(:,2), 'r', 'LineWidth', 2, 'LineStyle', '--') %plot standard discrete ellipse
xlabel('$^{206}$ Pb/$^{204}$ Pb', 'FontSize', 14, 'Interpreter', 'latex')
ylabel('$^{207}$ Pb/$^{204}$ Pb', 'FontSize', 14, 'Interpreter', 'latex')
xlabh = get(gca,'XLabel'); set(xlabh, 'Units', 'normalized')
set(xlabh,'Position',get(xlabh,'Position') - [0 xpadbuf 0])


hold off

ax(2) = subplot(1,2,2, 'FontSize', 12);  %for 207/204 - 208/204

hold on
scatter(dataratios(2,:), dataratios(3,:),'.','k')
pts = 200;
pis = 0:2*pi/(pts-1):2*pi;

hold on
set(gca, 'DataAspectRatioMode', 'auto')
circlepts = 2*[cos(pis') sin(pis')];
%plot(circlepts(:,1), circlepts(:,2))

for i = 1:n
    s = covmatsrat(2:3,2:3,i);
    sc = chol(s);
    elpts = circlepts*sc + repmat([dataratios(2,i) dataratios(3,i)], pts, 1);
        fill(elpts(:,1), elpts(:,2), colorV(i,:), ...
        'LineWidth', 0.2, 'FaceAlpha', alphaVal)
end

xsc = chol(xs(2:3, 2:3));
bigel = circlepts*xsc + repmat([mu(2) mu(3)], pts, 1);
plot(bigel(:,1), bigel(:,2), 'g', 'LineWidth', 2)
stcovc = chol(stcov(2:3,2:3));
stdel = circlepts*stcovc + repmat([mu(2) mu(3)], pts, 1);
plot(stdel(:,1), stdel(:,2), 'r', 'LineWidth', 2, 'LineStyle', '--')
xlabel('$^{207}$ Pb/$^{204}$ Pb', 'FontSize', 14, 'Interpreter', 'latex')
xlabh2 = get(ax(2),'XLabel'); set(xlabh2, 'Units', 'normalized')
set(xlabh2,'Position', get(xlabh2,'Position') - [0 xpadbuf 0])
ylabel('$^{208}$ Pb/$^{204}$ Pb', 'FontSize', 14, 'Interpreter', 'latex')

% for colorbar
barwidth = 0.03;
hcb = colorbar('Location','SouthOutside');
set(hcb, 'Position', [0.10 0.10 0.815 barwidth])
for i = 1:2
    pos = get(ax(i), 'Position');
    set(ax(i), 'PlotBoxAspectRatio', [1 1 1], ...
               'Position', [0.45*(i-1) 0.3 0.6 0.6]);
end
colormap(jet(ncolor));
set(hcb,'XTickMode','manual', 'XLim', [0 1], 'XTick', xticLoc, 'XTickLabel', xticLab, ...
        'Box', 'off', 'TickLength', [barwidth/4 0])
set(get(hcb, 'title'), 'string', 'Pb blank (pg)', 'FontSize', 12)
cbIm = findobj(hcb,'Type','image');
alpha(cbIm,alphaVal)
%myaa(4)
%plot2svg('ETICblanks.svg')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function F = vermeeshcopy3(S, n ,mu, smallS, datameas)
% F = 0;
% for i=1:n
%     smallSi = smallS(:,:,i);
%     A = inv([smallSi(1,1) + S(1,1), smallSi(1,2) + S(1,2), smallSi(1,3) + S(1,3)
%              smallSi(2,1) + S(2,1), smallSi(2,2) + S(2,2), smallSi(2,3) + S(2,3)
%              smallSi(3,1) + S(3,1), smallSi(3,2) + S(3,2), smallSi(3,3) + S(3,3)]);
%     B = [datameas(1,i) - mu(1); datameas(2,i) - mu(2); datameas(3,i) - mu(3)];
%     F = F + A*B*B'*A - A;
% end