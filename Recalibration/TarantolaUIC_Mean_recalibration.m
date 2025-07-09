% Do least squares for Critical Mixture dataset

% % Do only once:
% CRM112AllDataS = CRM112AllData;  %skipped version
% used112a = [1 1 1 1 1 0 0 1 0 1 1 0 1 0 1 1 1];
% tot112 = size(CRM112AllData,2);
% for i = tot112:-1:1
%     if used112a(i) == 0
%         CRM112AllDataS(:,i) = [];
%     end
% end
% %

diffU500 = 0;
diff112a = 10000000*eps;

%% setup
totU500s = size(U500AllData,2); tot112s = size(CRM112AllDataS,2); totICs = size(U500ICData,2);
totruns = totU500s + tot112s + totICs;
totms = 4 + 1*(totU500s + tot112s) + 2*(totU500s + tot112s + totICs);
totds = 2*(totU500s + tot112s + totICs);

uCM = zeros(totms,totms);
uCD = zeros(totds,totds);
uGn = zeros(totds,totms);

umprior = zeros(totms,1);
udobs = zeros(totds,1);
ufm = zeros(totds,1);

startbt = 4;
startst = 4 + totruns;
startbeta = 4 + totruns + tot112s+totU500s;

%% piece together dobs and CD

c = 0; %initialize counter
for i = 1:totU500s
    c = c+1;
    xi = U500AllData(2,i);  %pull out critical mix measurement data
    yi = U500AllData(1,i);
    sxi = U500AllData(4,i);
    syi = U500AllData(3,i);
    ri = U500AllData(5,i);
    covxyi = ri*sxi*syi;
    covi = [sxi^2 covxyi; covxyi syi^2];

    udobs(2*c-1) = xi; %first 233/235
    udobs(2*c  ) = yi; %then  238/235
    uCD((2*c-1):(2*c),(2*c-1):(2*c)) = covi;
end
for i = 1:tot112s
    c = c+1;
    xi = CRM112AllDataS(2,i);  %pull out critical mix measurement data
    yi = CRM112AllDataS(1,i);
    sxi = CRM112AllDataS(4,i);
    syi = CRM112AllDataS(3,i);
    ri = CRM112AllDataS(5,i);
    covxyi = ri*sxi*syi;
    covi = [sxi^2 covxyi; covxyi syi^2];

    udobs(2*c-1) = xi; %first 233/235
    udobs(2*c  ) = yi; %then  238/235
    uCD((2*c-1):(2*c),(2*c-1):(2*c)) = covi;
end
for i = 1:totICs
    c = c+1;
    xi = U500ICData(2,i);  %pull out critical mix measurement data
    yi = U500ICData(1,i);
    sxi = U500ICData(4,i);
    syi = U500ICData(3,i);
    ri = U500ICData(5,i);
    covxyi = ri*sxi*syi;
    covi = [sxi^2 covxyi; covxyi syi^2];

    udobs(2*c-1) = xi; %first 233/235
    udobs(2*c  ) = yi; %then  238/235
    uCD((2*c-1):(2*c),(2*c-1):(2*c)) = covi;
end

%% set up mprior
% from U standards inter-calibration

%constants:
m233 = 233.0396352;
m235 = 235.0439299;
m238 = 238.0507882;
r85b = 137.84;

%% Start MC madness
nM = 1;
%ummat = zeros(totms,nM);
%twoSm = zeros(nM,1);
%mswds = zeros(nM,1);
for M = 1:nM


%priors
r35t = 0.995083;
r85t = 0.00308027;
covt = [1 0;
        0 (0.003)^2]; % 100% relative uncertainties: diffuse prior

r85sU500 = 0.999781492879053 + diffU500;
r85s112a = 137.841362174735 + diff112a;
covstds = [6.6947550E-09	8.8279173E-07
           8.8279173E-07	1.3935118E-04];

r85bts = (4*10^-16)/(1*10^-9) * ones(totruns,1);
cov85bts = diag(r85bts.^2);

r85sts = zeros(totU500s + tot112s,1);
for i = 1:totU500s
    r85m = U500AllData(1,i);
    r85sts(i) = -((0.0072532*(-0.423368 + 137.84*r85m))/(-0.996787 + r85m)); %from Mathematica soln using init values
end
c = i;
for i = 1:tot112s
    c = c+1; %counter
    r85m = CRM112AllDataS(1,i);
    r85sts(c) = -((1.00001*(-0.423368 + 137.84*r85m))/(-137.428 + r85m)); %from same as above
end
cov85sts = diag((0.01*r85sts).^2);   % 1% relative unct.

betas = -0.236*ones(totruns,1);  %guess at 0.1 % per amu
covbetas = diag((0.1*betas).^2);  % 10% rel. unct.

umprior = [r35t r85t r85sU500 r85s112a r85bts' r85sts' betas']';  %values of mprior
uCM = blkdiag(covt, covstds, cov85bts, cov85sts, covbetas);

%% Monte Carlo Startup
covMC = covstds;
covstds = covstds * 10^-12;

%MCreal = mvnrnd([r85sU500 r85s112a], covMC);
%r85sU500 = MCreal(1); r85s112a = MCreal(2);

%new umprior and uCM for MC trial
%umprior(3:4) = [r85sU500 r85s112a]';
uCM(3:4,3:4) = covstds;  %now much smaller

%% evaluate fm and Gn: first time
c = 0; %initialize counter
r35t = umprior(1); r85t = umprior(2); 
for i = 1:totU500s
    c = c+1;
    r85s = umprior(3);
    r85bt = umprior(startbt+c); r85st = umprior(startst+c); beta = umprior(startbeta+c);
    ufm(2*c-1) = ((m233/m235)^beta*r35t)/(1 + r85bt/r85b + r85st/r85s);
    ufm(2*c  ) = ((m238/m235)^beta*(r85bt + r85st + r85t))/(1 + r85bt/r85b + r85st/r85s);
    uGn((2*c-1):(2*c),1:3) = [((m233/m235)^beta*r85b*r85s)/(r85bt*r85s + r85b*(r85s + r85st)), 0, ((m233/m235)^beta*r35t*r85b^2*r85st)/...
        (r85bt*r85s + r85b*(r85s + r85st))^2; 0, ((m238/m235)^beta*r85b*r85s)/(r85bt*r85s + r85b*(r85s + r85st)), ...
        ((m238/m235)^beta*r85b^2*r85st*(r85bt + r85st + r85t))/(r85bt*r85s + r85b*(r85s + r85st))^2]; %r35t r85t r85s
    uGn((2*c-1):(2*c), startbt+c) = [-(((m233/m235)^beta*r35t*r85b*r85s^2)/(r85bt*r85s + r85b*(r85s + r85st))^2), 
        ((m238/m235)^beta*r85b*r85s*(r85b*(r85s + r85st) - r85s*(r85st + r85t)))/(r85bt*r85s + r85b*(r85s + r85st))^2]; %r85bt
    uGn((2*c-1):(2*c), startst+c) = [-(((m233/m235)^beta*r35t*r85b^2*r85s)/(r85bt*r85s + r85b*(r85s + r85st))^2), 
        ((m238/m235)^beta*r85b*r85s*(r85bt*r85s - r85b*(r85bt - r85s + r85t)))/(r85bt*r85s + r85b*(r85s + r85st))^2]; %r85st
    uGn((2*c-1):(2*c), startbeta+c) = [((m233/m235)^beta*r35t*r85b*r85s*log(m233/m235))/(r85bt*r85s + r85b*(r85s + r85st)), 
        ((m238/m235)^beta*r85b*r85s*(r85bt + r85st + r85t)*log(m238/m235))/(r85bt*r85s + r85b*(r85s + r85st))]; %beta
end
for i = 1:tot112s
    c = c+1;
    r85s = umprior(4);
    r85bt = umprior(startbt+c); r85st = umprior(startst+c); beta = umprior(startbeta+c);
    ufm(2*c-1) = ((m233/m235)^beta*r35t)/(1 + r85bt/r85b + r85st/r85s);
    ufm(2*c  ) = ((m238/m235)^beta*(r85bt + r85st + r85t))/(1 + r85bt/r85b + r85st/r85s);
    uGn((2*c-1):(2*c),[1 2 4]) = [((m233/m235)^beta*r85b*r85s)/(r85bt*r85s + r85b*(r85s + r85st)), 0, ((m233/m235)^beta*r35t*r85b^2*r85st)/...
        (r85bt*r85s + r85b*(r85s + r85st))^2; 0, ((m238/m235)^beta*r85b*r85s)/(r85bt*r85s + r85b*(r85s + r85st)), ...
        ((m238/m235)^beta*r85b^2*r85st*(r85bt + r85st + r85t))/(r85bt*r85s + r85b*(r85s + r85st))^2];
    uGn((2*c-1):(2*c), startbt+c) = [-(((m233/m235)^beta*r35t*r85b*r85s^2)/(r85bt*r85s + r85b*(r85s + r85st))^2), 
        ((m238/m235)^beta*r85b*r85s*(r85b*(r85s + r85st) - r85s*(r85st + r85t)))/(r85bt*r85s + r85b*(r85s + r85st))^2]; %r85bt
    uGn((2*c-1):(2*c), startst+c) = [-(((m233/m235)^beta*r35t*r85b^2*r85s)/(r85bt*r85s + r85b*(r85s + r85st))^2), 
        ((m238/m235)^beta*r85b*r85s*(r85bt*r85s - r85b*(r85bt - r85s + r85t)))/(r85bt*r85s + r85b*(r85s + r85st))^2]; %r85st
    uGn((2*c-1):(2*c), startbeta+c) = [((m233/m235)^beta*r35t*r85b*r85s*log(m233/m235))/(r85bt*r85s + r85b*(r85s + r85st)), 
        ((m238/m235)^beta*r85b*r85s*(r85bt + r85st + r85t)*log(m238/m235))/(r85bt*r85s + r85b*(r85s + r85st))]; %beta
end
for i = 1:totICs
    c = c+1;
    r85bt = umprior(startbt+c); beta = umprior(startbeta+c);
    ufm(2*c-1) = ((m233/m235)^beta*r35t)/(1 + r85bt/r85b);
    ufm(2*c  ) = ((m238/m235)^beta*(r85bt + r85t))/(1 + r85bt/r85b);
    uGn((2*c-1):(2*c),1:2) = [((m233/m235)^beta*r85b)/(r85b + r85bt), 0; 0, ((m238/m235)^beta*r85b)/(r85b + r85bt)];%r35t r85t
    uGn((2*c-1):(2*c), startbt+c) = [-(((m233/m235)^beta*r35t*r85b)/(r85b + r85bt)^2), ...
        ((m238/m235)^beta*r85b*(r85b - r85t))/(r85b + r85bt)^2];
    uGn((2*c-1):(2*c), startbeta+c) = [((m233/m235)^beta*r35t*r85b*log(m233/m235))/(r85b + r85bt), ...
        ((m238/m235)^beta*r85b*(r85bt + r85t)*log(m238/m235))/(r85b + r85bt)];
end
%%

iters = 1; totiters = 500;
%twoSm = zeros(1,totiters);
%mswds = zeros(1,totiters);
%ummat = zeros(totms, totiters);

nun = 0.9;
um = umprior;
%gamman = CM*Gn'*invCD*(gm-dobs)+(mn-mprior);
%bn = Gn*gamman;
%dmping = (gamman'/CM*gamman)/(gamman'/CM*gamman + bn'*invCD*bn);
%dmpingv(j) = dmping;
%mn  = mn - (eye(size(CM,1)) + CM*Gn'*invCD*Gn)\gamman;
%mn = mn - nun * dmping*(CM*Gn'*(CD\(gm-dobs))+(mn-mprior)); %equation 3.89 of Tarantola (Inv. Prob. Theory)
%um = um - nun*(invCM + uGn'*invCD*uGn)\(uGn'*invCD*(ufm-udobs)+invCM*(um-umprior)); %#ok<MINV> %Newton
um = um - nun*(eye(size(uCM,1))+uCM*uGn'/uCD*uGn)\(uCM*uGn'/uCD*(ufm-udobs)+um-umprior);

%% 
for iters = 2:totiters;
    
    
    %% evaluate fm and Gn: first time
c = 0; %initialize counter
r35t = um(1); r85t = um(2); 
for i = 1:totU500s
    c = c+1;
    r85s = um(3);
    r85bt = um(startbt+c); r85st = um(startst+c); beta = um(startbeta+c);
    ufm(2*c-1) = ((m233/m235)^beta*r35t)/(1 + r85bt/r85b + r85st/r85s);
    ufm(2*c  ) = ((m238/m235)^beta*(r85bt + r85st + r85t))/(1 + r85bt/r85b + r85st/r85s);
    uGn((2*c-1):(2*c),1:3) = [((m233/m235)^beta*r85b*r85s)/(r85bt*r85s + r85b*(r85s + r85st)), 0, ((m233/m235)^beta*r35t*r85b^2*r85st)/...
        (r85bt*r85s + r85b*(r85s + r85st))^2; 0, ((m238/m235)^beta*r85b*r85s)/(r85bt*r85s + r85b*(r85s + r85st)), ...
        ((m238/m235)^beta*r85b^2*r85st*(r85bt + r85st + r85t))/(r85bt*r85s + r85b*(r85s + r85st))^2]; %r35t r85t r85s
    uGn((2*c-1):(2*c), startbt+c) = [-(((m233/m235)^beta*r35t*r85b*r85s^2)/(r85bt*r85s + r85b*(r85s + r85st))^2), 
        ((m238/m235)^beta*r85b*r85s*(r85b*(r85s + r85st) - r85s*(r85st + r85t)))/(r85bt*r85s + r85b*(r85s + r85st))^2]; %r85bt
    uGn((2*c-1):(2*c), startst+c) = [-(((m233/m235)^beta*r35t*r85b^2*r85s)/(r85bt*r85s + r85b*(r85s + r85st))^2), 
        ((m238/m235)^beta*r85b*r85s*(r85bt*r85s - r85b*(r85bt - r85s + r85t)))/(r85bt*r85s + r85b*(r85s + r85st))^2]; %r85st
    uGn((2*c-1):(2*c), startbeta+c) = [((m233/m235)^beta*r35t*r85b*r85s*log(m233/m235))/(r85bt*r85s + r85b*(r85s + r85st)), 
        ((m238/m235)^beta*r85b*r85s*(r85bt + r85st + r85t)*log(m238/m235))/(r85bt*r85s + r85b*(r85s + r85st))]; %beta
end
for i = 1:tot112s
    c = c+1;
    r85s = um(4);
    r85bt = um(startbt+c); r85st = um(startst+c); beta = um(startbeta+c);
    ufm(2*c-1) = ((m233/m235)^beta*r35t)/(1 + r85bt/r85b + r85st/r85s);
    ufm(2*c  ) = ((m238/m235)^beta*(r85bt + r85st + r85t))/(1 + r85bt/r85b + r85st/r85s);
    uGn((2*c-1):(2*c),[1 2 4]) = [((m233/m235)^beta*r85b*r85s)/(r85bt*r85s + r85b*(r85s + r85st)), 0, ((m233/m235)^beta*r35t*r85b^2*r85st)/...
        (r85bt*r85s + r85b*(r85s + r85st))^2; 0, ((m238/m235)^beta*r85b*r85s)/(r85bt*r85s + r85b*(r85s + r85st)), ...
        ((m238/m235)^beta*r85b^2*r85st*(r85bt + r85st + r85t))/(r85bt*r85s + r85b*(r85s + r85st))^2];
    uGn((2*c-1):(2*c), startbt+c) = [-(((m233/m235)^beta*r35t*r85b*r85s^2)/(r85bt*r85s + r85b*(r85s + r85st))^2), 
        ((m238/m235)^beta*r85b*r85s*(r85b*(r85s + r85st) - r85s*(r85st + r85t)))/(r85bt*r85s + r85b*(r85s + r85st))^2]; %r85bt
    uGn((2*c-1):(2*c), startst+c) = [-(((m233/m235)^beta*r35t*r85b^2*r85s)/(r85bt*r85s + r85b*(r85s + r85st))^2), 
        ((m238/m235)^beta*r85b*r85s*(r85bt*r85s - r85b*(r85bt - r85s + r85t)))/(r85bt*r85s + r85b*(r85s + r85st))^2]; %r85st
    uGn((2*c-1):(2*c), startbeta+c) = [((m233/m235)^beta*r35t*r85b*r85s*log(m233/m235))/(r85bt*r85s + r85b*(r85s + r85st)), 
        ((m238/m235)^beta*r85b*r85s*(r85bt + r85st + r85t)*log(m238/m235))/(r85bt*r85s + r85b*(r85s + r85st))]; %beta
end
for i = 1:totICs
    c = c+1;
    r85bt = um(startbt+c); beta = um(startbeta+c);
    ufm(2*c-1) = ((m233/m235)^beta*r35t)/(1 + r85bt/r85b);
    ufm(2*c  ) = ((m238/m235)^beta*(r85bt + r85t))/(1 + r85bt/r85b);
    uGn((2*c-1):(2*c),1:2) = [((m233/m235)^beta*r85b)/(r85b + r85bt), 0; 0, ((m238/m235)^beta*r85b)/(r85b + r85bt)];%r35t r85t
    uGn((2*c-1):(2*c), startbt+c) = [-(((m233/m235)^beta*r35t*r85b)/(r85b + r85bt)^2), ...
        ((m238/m235)^beta*r85b*(r85b - r85t))/(r85b + r85bt)^2];
    uGn((2*c-1):(2*c), startbeta+c) = [((m233/m235)^beta*r35t*r85b*log(m233/m235))/(r85b + r85bt), ...
        ((m238/m235)^beta*r85b*(r85bt + r85t)*log(m238/m235))/(r85b + r85bt)];
end
    
    um = um - nun*(eye(size(uCM,1))+uCM*uGn'/uCD*uGn)\(uCM*uGn'/uCD*(ufm-udobs)+um-umprior); %pre-conditioned steepest ascent
    %ummat(:,iters) = um;
    %twoSm(iters) = (ufm-udobs)'/uCD*(ufm-udobs);  %eqn 5.141 of Tarantola InvPrblmTheory 
end %for i = 1:totiters

%ummat(:,M) = um;
twoSm(M) = (ufm-udobs)'/uCD*(ufm-udobs);  %eqn 5.141 of Tarantola InvPrblmTheory
%mswds(M) = twoSm(M)/(length(udobs)-1);
dlmwrite('MCICtola.txt', [um' twoSm(M)], '-append', 'precision', '%12.12g')

end %for M = 1:nM
