%% Calculate mean and variablility of blank IC with linear regression
    %with over-dispersion calculation
    
% Instructions: run BlankIC_ifyouknowtracer.m on blanks data to get
% overdispsion matrix 'xs'.  Then run this code to do linear fit.  The
% former code requires an initial guess of tracer IC.

%% Load in Loading blanks, assume correlation coefficients are zero
% (even though they're not)

%global trIC45
%trIC45 = 0.00009;   %starting dist from y-int to tracer
%trIC45 = 0.000105;

%data = ET535LoadingBlanks;  %for blanksdata structure, see for loop below:
% ET*LoadingBlanks structure:
% col 1      2       3       4        5         6       7        8        9         10
% alphaPb 2s(abs) r204/205m 2s(abs) r206/205m 2s(abs) r207/205m 2s(abs) r208/205m 2s(abs) 

n = size(data,2);  %number of blanks measured

covmats205  = zeros(4,4,n);   %pre-allocate for uncertainty propagation
covmatsmolb = zeros(4,4,n);
dataratios  = zeros(4,n);    %fractionation and tracer-corrected points

for i = 1:n
   alphaPb = data(1,i);         %from Excel sheet
   s_alphaPb = data(2,i);       %corrected for beam interpolation
   r45m = data(3,i);   
   s_r45m = data(4,i)*r45m/100; %unct. entered as 1-sigma pct
   r65m = data(5,i);   
   s_r65m = data(6,i)*r65m/100;            
   r75m = data(7,i);
   s_r75m = data(8,i)*r75m/100;
   r85m = data(9,i);
   s_r85m = data(10,i)*r85m/100;
    
%    alphaPb = data(1,i)/100;         %from Excel sheet
%    s_alphaPb = data(2,i)/200;       %corrected for beam interpolation
%    r45m = data(3,i);   
%    s_r45m = data(4,i)/2; %unct. entered as 2-sigma abs
%    r65m = data(5,i);   
%    s_r65m = data(6,i)/2;            
%    r75m = data(7,i);
%    s_r75m = data(8,i)/2;
%    r85m = data(9,i);
%    s_r85m = data(10,i)/2;
   
   r45fc = r45m*(1-  alphaPb);
   r65fc = r65m*(1+  alphaPb);   %with 205 in denominator
   r75fc = r75m*(1+2*alphaPb);
   r85fc = r85m*(1+3*alphaPb);
   dataratios(:,i) = [r45fc r65fc r75fc r85fc]';
   % calculate covariance matrix and jacobians

   %covariance matrix for input measured ratios, assuming measured ratios 
   %to 205 are independent
   covi = diag([s_r45m^2 s_r65m^2 s_r75m^2 s_r85m^2 s_alphaPb^2]);  

   %jacobian for transformation to fractionation-corrected values
   j1 = [diag([1 - alphaPb 1+alphaPb 1+2*alphaPb 1+3*alphaPb])...
       [-r45m r65m 2*r75m 3*r85m]'];

   covmats205(:,:,i) =    j1*covi*j1';

end

% calculate weighted mean with only random errors as first guess
covmats205_inv = zeros(4,4,n);
%wtdpts = zeros(4,n);
for i = 1:n
    covmats205_inv(:,:,i) = inv(covmats205(:,:,i));
%    wtdpts(:,i) = covmats205_inv(:,:,i)*dataratios(:,i);
end
%m1 = sum(covmats205_inv,3);
%m2 = sum(wtdpts,2);
%mu = m1\m2;                

%smallS = covmats205;  %these are the measurement uncertainties to throw the solver
datameas = dataratios; %these are the fractionation and tracer-corrected blank ratios

%%   solve for line using MLE (set derivative to zero)

aguess = [0.000 0.000 0.000]; vguess = [18 15 38];
avguess = [aguess vguess];
options = optimset('TolX',10^-15, 'TolFun',10^-15,...
    'MaxFunEvals',10^12,'Diagnostics','on', 'MaxIter',10^5);
%avmin = fsolve(@(av) MLmax4Dline(av, datameas, covmats205_inv), avguess, options);
[avmin fs] = fsolve(@(av) MLmax4Dline2(av, datameas, covmats205), avguess, options);

%% solve for analytical covariance matrix for line fit

a = [0 avmin(1:3)]'; v = [1 avmin(4:6)]';

d2LdaaT = 0;
d2LdvvT = 0;
d2LdavT = 0;
chisq = 0;
for i = 1:n
    si = covmats205(:,:,i);
    xi = datameas(:,i) - a;
    d2LdaaT = d2LdaaT + -(inv(si) - ( (si\v) * (v'/si) )/(v'/si*v));
    
    d2LdavT = d2LdavT + -( (v'/si*v) * ( (si\v) * (xi'/si) + (xi'/si*v) * inv(si) ) - ...
        2*(xi'/si*v) * (si\v)* (v'/si) ) /...
        (v'/si*v)^2;
    
    d2LdvvT = d2LdvvT + ((v'/si*v)^2*(2 * (v'/si*xi) * (si\xi) * (v'/si) + ...
       (v'/si*v) * (si\xi) * (xi'/si) - (v'/si*xi)^2*inv(si) - ...
       2*(v'/si*xi) * (si\v) * (xi'/si) ) -...
       4* (v'/si*v) * ( (v'/si*v) * (v'/si*xi) * (si\xi) - (v'/si*xi)^2* (si\v)) * (v'/si) )/ ...
             (v'/si*v)^4;
    chisq = chisq + (xi'/si*xi) - (v'/si*xi)^2/(v'/si*v);
    
end

fim = -[d2LdaaT d2LdavT; d2LdavT' d2LdvvT];
fimr = fim;  %the sliced and diced version
fimr(5,:) = []; fimr(1,:) = [];   %delete rows and columns that are assumed to be 1 or 0
fimr(:,5) = []; fimr(:,1) = [];
Sav = inv(fimr);
mswd = chisq / ((4-1)*(n-2));

%% now do the overdispersion calculation with a version of the same code
% 
% %aguess = [0.00 0.00 0.000]; vguess = [18 15 37];
% a = [0 avmin(1:3)]'; v = [1 avmin(4:6)]';
% %Xiguess = zeros(1,10); 
% Xiguess = [1 1 1 1 0 0 0 0 0 0];
% %avXavGuess = [Xiguess aguess vguess];
% %avXavGuess = avminod;
% options = optimset('TolX',10^-15, 'TolFun',10^-15,...
%     'MaxFunEvals',10^5,'Diagnostics','on', 'MaxIter',10^6, 'Display', 'on'); %,...
% %    'DiffMaxChange',10^-11, 'DiffMinChange',10^-16);
% [avminod mlderiv] = fsolve(@(avXav) MLmax4DlineODXiAVknown(avXav, datameas, covmats205, a, v), Xiguess, options);
% 
% Xj_reshaped =  [avminod(1) avminod(5) avminod(6)  avminod(7);
%                 avminod(5) avminod(2) avminod(8)  avminod(9); 
%                 avminod(6) avminod(8) avminod(3)  avminod(10);
%                 avminod(7) avminod(9) avminod(10) avminod(4)];
% 
% %Xiguess = Xj_reshaped;

%%  Use over-dispersion calculation for new covmats, do linear fit again

%a = [0 avminod(11:13)]';

%get xs from BlankIC_ifyouknowtracer.m
%xs = zeros(3,3);

covxs = blkdiag(0,xs);

covmats205od = zeros(size(covmats205));
tj = zeros(n,1);
for j = 1:n
    pntj = datameas(:,j);  %measured ratios i  (a blank IC run)
    %xi = pnti - a;         %used often, no physical meaning
    tj(j) = pntj(1) - trIC45;
    covmats205od(:,:,j) = covmats205(:,:,j) + tj(j)^2*covxs ; % covariance matrix for measured ratio i WITH tj SQUARED
end
%%

aguess = avmin(1:3); vguess = avmin(4:6);
avguess = [aguess vguess];

options = optimset('TolX',10^-12, 'TolFun',10^-12,...
    'MaxFunEvals',10^12,'Diagnostics','on', 'MaxIter',10^5);
avmin2 = fsolve(@(av) MLmax4Dline2(av, datameas, covmats205od), avguess, options);

% xsi = xs;
% avXguess = [xsi(1,1) xsi(2,2) xsi(3,3) xsi(1,2) xsi(1,3) xsi(2,3) avmin(1:3) avmin(4:6)]';
% avmin2a = fsolve(@(avXav) MLmax4DlineODXiAll2(avXav, datameas, covmats205), avXguess, options);

%%
a = [0 avmin2(1:3)]'; v = [1 avmin2(4:6)]';

d2LdaaT = 0;
d2LdvvT = 0;
d2LdavT = 0;
chisq = 0;
for i = 1:n
    si = covmats205od(:,:,i);
    xi = datameas(:,i) - a;
    d2LdaaT = d2LdaaT + -(inv(si) - ( (si\v) * (v'/si) )/(v'/si*v));
    
    d2LdavT = d2LdavT + -( (v'/si*v) * ( (si\v) * (xi'/si) + (xi'/si*v) * inv(si) ) - ...
        2*(xi'/si*v) * (si\v)* (v'/si) ) /...
        (v'/si*v)^2;
    
    d2LdvvT = d2LdvvT + ((v'/si*v)^2*(2 * (v'/si*xi) * (si\xi) * (v'/si) + ...
       (v'/si*v) * (si\xi) * (xi'/si) - (v'/si*xi)^2*inv(si) - ...
       2*(v'/si*xi) * (si\v) * (xi'/si) ) -...
       4* (v'/si*v) * ( (v'/si*v) * (v'/si*xi) * (si\xi) - (v'/si*xi)^2* (si\v)) * (v'/si) )/ ...
             (v'/si*v)^4;
    chisq = chisq + (xi'/si*xi) - (v'/si*xi)^2/(v'/si*v);
    
end


%%
fim2 = -[d2LdaaT d2LdavT; d2LdavT' d2LdvvT];
fimr2 = fim2;  %the sliced and diced version
fimr2(5,:) = []; fimr2(1,:) = [];   %delete rows and columns that are assumed to be 1 or 0
fimr2(:,5) = []; fimr2(:,1) = [];
Sav2 = inv(fimr2);
mswd2 = chisq / ((4-1)*(n-2));

% aguess = [0.002 0.003 0.007]; vguess = [18 15 37];
% %Xiguess = [xs(1,1) xs(1,2) xs(1,3) xs(2,2) xs(2,3) xs(3,3)] * norm(v);
% %Xiguess = reshape(blkdiag(0, xs);
% Xiguess = zeros(1,10);
% covmats205_invp = covmats205_inv;
% 
% for i = 1:100
%     
% % FIRST - Find best line fit
% 
% 
% avguess = [aguess vguess];
% options = optimset('TolX',10^-12, 'TolFun',10^-12,...
%     'MaxFunEvals',10^12,'Diagnostics','off', 'MaxIter',10^4, 'Display', 'off');
% avmin = fsolve(@(av) MLmax4Dline(av, datameas, covmats205_invp), avguess, options);
% 
% a = [0 avmin(1:3)]'; v = [1 avmin(4:6)]';
% aguess = a; vguess = v;
% 
% options = optimset('TolX',10^-6, 'TolFun',10^-6,...
%     'MaxFunEvals',20,'Diagnostics','off', 'MaxIter',20, 'Display', 'off');
% avminod = fsolve(@(avX) MLmax4DlineODXi(avX, a, v, datameas, covmats205), Xiguess, options);
%     
% Xj_reshaped =  [avminod(1) avminod(5) avminod(6)  avminod(7);
%                 avminod(5) avminod(2) avminod(8)  avminod(9); 
%                 avminod(6) avminod(8) avminod(3)  avminod(10);
%                 avminod(7) avminod(9) avminod(10) avminod(4)];
%             
% Xiguess = Xj_reshaped;
%             
% for j = 1:n
%     pntj = datameas(:,j);  %measured ratios i  (a blank IC run)
%     %xi = pnti - a;         %used often, no physical meaning
%     tj = sqrt((pntj-a)'*(pntj-a)) - c;
%     covmats205_invp(:,:,j) = inv(covmats205(:,:,j) + tj*Xj_reshaped) ; %inverse covariance matrix for measured ratio i
% end
% 
% 
% disp(num2str(i))

%%

% end

% d2LdaaT = 0;
% d2LdvvT = 0;
% d2LdavT = 0;
% chisq_xs = 0;
% for i = 1:n
%     sinvi = covmats205_invp(:,:,i);
%     xi = datameas(:,i) - a;
%     d2LdaaT = d2LdaaT + -(sinvi-(sinvi*v*v'*sinvi)/(v'*sinvi*v));
%     d2LdavT = d2LdavT + -(v'*sinvi*v*(sinvi*v*xi'*sinvi + xi'*sinvi*v*sinvi) - ...
%         2*xi'*sinvi*v*sinvi*v*v'*sinvi)/...
%         (v'*sinvi*v)^2;
%     d2LdvvT = d2LdvvT + ((v'*sinvi*v)^2*(2*v'*sinvi*xi*sinvi*xi*v'*sinvi + ...
%        v'*sinvi*v*sinvi*xi*xi'*sinvi - (v'*sinvi*xi)^2*sinvi - ...
%        2*v'*sinvi*xi*sinvi*v*xi'*sinvi) -...
%        4*v'*sinvi*v*(v'*sinvi*v*v'*sinvi*xi*sinvi*xi - (v'*sinvi*xi)^2*sinvi*v)*v'*sinvi)/ ...
%              (v'*sinvi*v)^4;
%     chisq_xs = chisq_xs + xi'*sinvi*xi - (v'*sinvi*xi)^2/(v'*sinvi*v);
%     
% end
% 
% fim_od = -[d2LdaaT d2LdavT; d2LdavT' d2LdvvT];
% fim_odr = fim_od;  %the sliced and diced version
% fim_odr(5,:) = []; fim_odr(1,:) = [];   %delete rows and columns that are assumed to be 1 or 0
% fim_odr(:,5) = []; fim_odr(:,1) = [];
% Sav_xs = inv(fim_odr);
% mswd_xs = chisq_xs / (n-3);
%%

% %% vermeeshcopy3 uses 
% %oldopts = optimset('fsolve');
% 
% for resoloops = 1:5   %five times is good enough for convergence
% options = optimset('TolX',10^-20, 'TolFun',10^-20,...
%     'MaxFunEvals',10^12,'Diagnostics','on', 'MaxIter',10^4);
% xs = fsolve(@(S) vermeeshcopy_line(S, n ,mu, smallS, datameas), [0 0 0; 0 0 0; 0 0 0], options);
% 
%     for i = 1:n  %weighted mean, with overdispersion factored in
%         covmats205_inv(:,:,i) = inv(covmatsrat(:,:,i)+xs);
%         wtdpts(:,i) = covmats205_inv(:,:,i)*dataratios(:,i);
%     end  %make inverse covariance ratios with overdispersion addded
%     m1 = sum(covmatsrat_inv,3);
%     m2 = sum(wtdpts,2);
%     mu = m1\m2;
% 
% end

%%  make data readable
r64b = v(2); r74b = v(3); r84b = v(4); blankIC_535Line = v(2:4);
s64b = sqrt(Sav2(4,4));
s74b = sqrt(Sav2(5,5));
s84b = sqrt(Sav2(6,6));
rho6474 = Sav2(4,5)/sqrt(s64b*s74b);
rho6484 = Sav2(4,6)/sqrt(s64b*s84b);
rho7484 = Sav2(5,6)/sqrt(s74b*s84b);

disp(['mean ', num2str(r64b,'%2.3f'), ' ', num2str(r74b,'%2.3f'), ' ', num2str(r84b,'%2.3f')])
disp(['2std  ', num2str(2*s64b,'%2.3f'), '  ', num2str(2*s74b,'%2.3f'), '  ', num2str(2*s84b,'%2.3f')])
disp(['rhos  ', num2str(rho6474,'%2.3f'), '  ', num2str(rho6484,'%2.3f'), '  ', num2str(rho7484,'%2.3f')])
disp(' ')

trIC45 = trIC45; %#ok<ASGSL>
trIC65 = trIC45*r64b + a(2);
trIC75 = trIC45*r74b + a(3);
trIC85 = trIC45*r84b + a(4);

trIC_ET535 = a + trIC45*v;

trIC45s = .1*trIC45;   %10% 1-sigma uncertainty.

SaaT = blkdiag(0,Sav2(1:3,1:3));
SvaT = blkdiag(0,Sav2(1:3,4:6));
SvvT = blkdiag(0,Sav2(4:6,4:6));
Sav8 = [SaaT SvaT'; SvaT SvvT];
Sav9 = blkdiag(Sav8, trIC45s^2);
Jtrbl_1 = eye(4);       Jtrbl_2 = trIC45*eye(4);        Jtrbl_3 = v;
Jtrbl_4 = zeros(3,4);   Jtrbl_5 = [zeros(3,1) eye(3)];  Jtrbl_6 = zeros(3,1);
Jtrbl = [Jtrbl_1 Jtrbl_2 Jtrbl_3;
         Jtrbl_4 Jtrbl_5 Jtrbl_6];
covtrbl_ET535 = Jtrbl*Sav9*Jtrbl';

Sav9_535 = Sav9;
v_535 = v;
trIC45s_535 = trIC45s;
%trICs = as + trIC45^2*vs + ;

