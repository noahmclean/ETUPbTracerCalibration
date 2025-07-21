% Run this after making the structures MC and mix with inverseGravTrac.m

n.runsTot = length(mix.mixList);            %total runs
n.skips = n.runsTot - sum(mix.skips);       %number of skipped runs
n.used = n.runsTot - n.skips;               %number of used runs
n.runs202 = sum(mix.is202 .* mix.skips);    %ET2535 runs used
n.runs205 = n.used-n.runs202;               %ET535 runs used

n.cyclesPerBlock = 20; 

n.blocksPb = 0;  n.blocksU  = 0; n.blocks202 = 0; n.blocks205 = 0;
for i = 1:n.runsTot
    n.blocksPb = n.blocksPb + floor(length(mix.ratios.Pb{i})/n.cyclesPerBlock) * mix.skips(i);
    n.blocksU  = n.blocksU  + floor(length(mix.ratios.U{i}) /n.cyclesPerBlock) * mix.skips(i);
    n.blocks202 = n.blocks202 + floor(length(mix.ratios.Pb{i})/n.cyclesPerBlock) * mix.skips(i) * mix.is202(i);
end
n.blocks205 = n.blocksPb - n.blocks202;

%start positions
% 2 global variables: [r202205t r235205t]
% 6 run-wise Pb variables: [gamma205, gamma207, r207206b, r208206b, r238205gt, r206205bt];
% 2 run-wise U variables: [r85bt r85st]
% 2 block-wise variables: betaPb, betaU

st.g205 = 2;                 %start of gamma205
st.g207 = st.g205 + n.used;  %start of gamma207
st.r76b = st.g207 + n.used;  %start of r207/206b
st.r86b = st.r76b + n.used;  %start of r208/206b
st.rUgt = st.r86b + n.used;  %start of r238g/235t
st.rPbt = st.rUgt + n.used;  %start of r206b/205t
st.rUbt = st.rPbt + n.used;  %start of r238b/235t
st.bePb = st.rUbt + n.used;  %start of betaPb
st.betU = st.bePb + n.blocksPb; %start of betaU

n.totalMs = st.betU + n.blocksU;
n.totalDs = n.blocks205*3 + n.blocks202*4 + n.blocksU*2;

umix = mix;
for i = n.runsTot:-1:1
    if ~umix.skips(i)
        umix.ratios.Pb(i) = [];
        umix.ratios.U(i)  = [];
        umix.gravName(i)  = [];
        umix.is202(i)     = [];
        umix.mixList(i)   = [];
        umix.skips(i)     = [];
    end
end

%% Set up CM and mprior

mprior = zeros(n.totalMs,1);
CM = zeros(n.totalMs, n.totalMs);

mprior(1) = 0.9995; %202/205t
CM(1,1) = 0.05^2;
mprior(2) = 100.18;
CM(2,2) = (0.5)^2;

mprior((st.g205+1):(st.r76b)) = 0.9995;  %all gammas
CM((st.g205+1:st.r76b),(st.g205+1:st.r76b)) = diag( 0.05^2*ones(2*n.used,1) );

% Note comment here to top to implement new blank IC

%          [207/206b        208/206b], from tracer-blank IC
misc.icb = [0.837045084     2.042911461]; % From MIT-dominated blank-tracer mix
%misc.icb = [0.858862933	2.089768821]; % From NIGL Tracer IC

misc.covb = [4.789190E-05	9.275855E-05
             9.275855E-05	4.594616E-04];
mprior((st.r76b+1):st.r86b) = misc.icb(1);
mprior((st.r86b+1):st.rUgt) = misc.icb(2);
CM((st.r76b+1):st.r86b, (st.r76b+1):st.r86b) = diag( misc.covb(1,1)*ones(n.used,1) );
CM((st.r86b+1):st.rUgt, (st.r86b+1):st.rUgt) = diag( misc.covb(2,2)*ones(n.used,1) );
CM((st.r86b+1):st.rUgt, (st.r76b+1):st.r86b) = diag( misc.covb(1,2)*ones(n.used,1) );
CM((st.r76b+1):st.r86b, (st.r86b+1):st.rUgt) = diag( misc.covb(2,1)*ones(n.used,1) );

%238g/235t
for i = 1:n.used
     mprior(st.rUgt+i) = umix.ratios.U{1,i}(1,1)*1.003;  %raw 238/235 as mol238g/mol235t
end
CM((st.rUgt+1):st.rPbt,(st.rUgt+1):st.rPbt) = diag( (mprior((st.rUgt+1):st.rPbt)).^2 ); % 100% prior unct.

%r206b/205t from ResultsTable_v5.xlsx using measured tracer weights, 0.3pg loading blank
misc.r65bt = [0.000237931,0.000237931,0.000237931,0.000118966,0.000118966,0.000172915,0.000172664,0.000164760,0.000164760,0.000159650,...
    0.000159650,9.80405e-05,9.83458e-05,9.84082e-05,9.88005e-05,0.000164760,0.000164760,0.000159650,0.000159650,0.000178306,0.000178591,...
    0.000223760,0.000227845,6.94920e-05,0.000211357,0.000182044,0.000246646,0.000227251,0.000159702,0.000182681,0.000182681,0.000159702,...
    0.000159702,0.000237931,0.000118966,0.000118966,0.000160115,0.000160115,0.000160994,0.000178448,0.000113398,9.87813e-05,0.000160115,...
    0.000160115,0.000160994,0.000160994];
mprior((st.rPbt+1):st.rUbt) = misc.r65bt;
CM((st.rPbt+1):st.rUbt,(st.rPbt+1):st.rUbt) = diag( (0.1*misc.r65bt).^2*ones(n.used,1) );

%r238b/235t from ResultsTable_v5.xlsx using measured tracer weights, 0.1pg loading blank
misc.r85bt = [2.691205E-06,2.691205E-06,2.691205E-06,1.345602E-06,1.345602E-06,1.955817E-06,1.952979E-06,1.863572E-06,1.863572E-06,...
    1.805774E-06,1.805774E-06,1.108922E-06,1.112375E-06,1.113080E-06,1.117517E-06,1.863572E-06,1.863572E-06,1.805774E-06,1.805774E-06,...
    2.016790E-06,2.020020E-06,2.530914E-06,2.577124E-06,7.860133E-07,2.390624E-06,2.059070E-06,2.789777E-06,2.570396E-06,1.806366E-06,...
    2.066272E-06,2.066272E-06,1.806366E-06,1.806366E-06,2.691205E-06,1.345602E-06,1.345602E-06,1.811040E-06,1.811040E-06,1.820979E-06,...
    2.018404E-06,1.282626E-06,1.117301E-06,1.811040E-06,1.811040E-06,1.820979E-06,1.820979E-06];
mprior((st.rUbt+1):st.bePb) = misc.r85bt;
CM((st.rUbt+1):st.bePb, (st.rUbt+1):st.bePb) = diag( (0.2*misc.r85bt).^2*ones(n.used,1) );

%note: mprior and CM for betas unpacked along with the measured data
misc.r85b = 137.817;

%% assemble dobs and CD
dobs = zeros(n.totalDs,1);
CD = zeros(n.totalDs,n.totalDs);
misc.r18  = 0.00207;
misc.r18s = 0.00001;

count.ratio = 1;
count.blocksPb = 1; count.blocksU = 1;
for i = 1:n.used
    
    n.blocksi = floor(length(umix.ratios.Pb{i})/n.cyclesPerBlock);
    count.blockrangej = 1:n.cyclesPerBlock;
    for jPb = 1:n.blocksi %for the Pb analysis
        dobs(count.ratio:(count.ratio+2+umix.is202(i))) = mean( umix.ratios.Pb{i}(count.blockrangej,:) );
        CD(count.ratio:(count.ratio+2+umix.is202(i)),count.ratio:(count.ratio+2+umix.is202(i))) ...
            = cov(umix.ratios.Pb{i}(count.blockrangej,:))/(n.cyclesPerBlock-1);
        count.blockrangej = count.blockrangej + n.cyclesPerBlock;
        count.ratio = count.ratio + 3 + umix.is202(i);
        mprior(st.bePb+count.blocksPb) = -log(1+3*beta0.Pb(i,jPb))/log(208/205);
        CM(st.bePb+count.blocksPb,st.bePb+count.blocksPb) = (2*mprior(st.bePb+count.blocksPb))^2;
        count.blocksPb = count.blocksPb + 1;
    end
    
    n.blocksi = floor(length(umix.ratios.U{i})/n.cyclesPerBlock);
    count.blockrangej = 1:n.cyclesPerBlock;
    for jU = 1:n.blocksi %for the Pb analysis
        misc.uic = mean( umix.ratios.U{i}(count.blockrangej,:) );
        misc.r07m = misc.uic(1); misc.r57m = misc.uic(2);
        dobs(count.ratio:(count.ratio+1)) = [misc.r07m/(1-2*misc.r18*misc.r57m) misc.r57m/(1-2*misc.r18*misc.r57m)];
        
        misc.ucov = cov(umix.ratios.U{i}(count.blockrangej,:))/(n.cyclesPerBlock-1);
        misc.jucovr18 = [1/(1 - 2*misc.r18*misc.r57m), (2*misc.r07m*misc.r18)/(1 - 2*misc.r18*misc.r57m)^2, ...
            (2*misc.r07m*misc.r57m)/(1 - 2*misc.r18*misc.r57m)^2;
            0, 1/(1 - 2*misc.r18*misc.r57m)^2, (2*misc.r57m^2)/(1 - 2*misc.r18*misc.r57m)^2];
        CD(count.ratio:(count.ratio+1),count.ratio:(count.ratio+1)) ...
           = misc.jucovr18 * blkdiag(misc.ucov, misc.r18s^2) * misc.jucovr18';  %add in unct. for oxide correction
        count.blockrangej = count.blockrangej + n.cyclesPerBlock;
        count.ratio = count.ratio + 2;
        mprior(st.betU+count.blocksU) = -log(1+3*beta0.U(i,jU))/log(238/235);
        CM(st.betU+count.blocksU, st.betU+count.blocksU) = (2*mprior(st.betU+count.blocksU))^2;
        count.blocksU = count.blocksU + 1;
    end

end