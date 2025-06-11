%% first using 981 + 202-205 data,
    % assume 
         %1. 'minor' IC of tracer, and IC of Pb blank
         %2. IC of BaPO2 and Tl
         %3. 208/206 of 981
         %4. measured data of Amelin and Davis
    %to determine 
        %1. the 202/205 of the (Amelin) tracer
        %2. the the BaPO2 and, if present, Tl interference magnitude
        %3. the gravimetric/tracer ratio
        %4. the (inter-calibrated) 204/206 and 207/206 of 981

% Use Total Inversion technique from Tarantola and Valette, but from the InverseProblemTheory book (5.5)
% Diffuse prior on r46g and r76g, set gammaEven = 1;

warning off MATLAB:nearlySingularMatrix

%% Pick through 981 data

nbs981runs = [19:24 50:52 54];  %the indices in variable measnames from parsing code
samplesizes981 = [1863,1863,1866,1865,1864,572,1898,1899,536,1839];  %in pg
nruns981 = length(nbs981runs);

% to hold the count number of the starting block of each analysis
startblock981 = cell(4,length(nbs981runs));
for i = 1:length(nbs981runs)
    startblock981{1,i} = measnames{nbs981runs(i)};
end

totalblocks981 = 0;
runcount981 = 0;
blocklims981 = zeros(nruns981,2);
for runi = nbs981runs
    runcount981 = runcount981 + 1;
    mdata = evalin('base', ['Amelin_' measnames{runi}]);
    nblocks = size(mdata,1);
    startblock981{3,runcount981} = nblocks;
    startblock981{2,runcount981} = totalblocks981 + 1;
    blocklims981(runcount981,1) = totalblocks981 + 1;
    for i = 1:nblocks
        totalblocks981 = totalblocks981 + 1;
    end
    blocklims981(runcount981,2) = totalblocks981;
end

skips981 = [2 3 4 5 6 7 15 16 27 28 38 39 40 50 51 52 53 54 55 56 57 58 59 97 98 119 120 121 138 139 174 176 177];
skipv981 = zeros(totalblocks981,1); skipv981(skips981) = 1;   %make a vector of zeros with blocknum = 1 when block is skipped.

%% Parse through 982 data

nbs982runs = [34:35 40:42 44:47];  %the indices in variable measnames from parsing code
samplesizes982 = [1714, 1721, 1964, 1557, 1798, 1800, 1868, 2273, 2099, 2271];  %in pg
nruns982 = length(nbs982runs);

% to hold the count number of the starting block of each analysis
startblock982 = cell(4,length(nbs982runs));
for i = 1:length(nbs982runs)
    startblock982{1,i} = measnames{nbs982runs(i)};
end

totalblocks982 = 0;
runcount982 = 0;
blocklims982 = zeros(nruns982,2);
for runi = nbs982runs
    runcount982 = runcount982 + 1;
    mdata = evalin('base', ['Amelin_' measnames{runi}]);
    nblocks = size(mdata,1);
    startblock982{3,runcount982} = nblocks;
    startblock982{2,runcount982} = totalblocks982 + 1;
    blocklims982(runcount982,1) = totalblocks982 + 1;
    for i = 1:nblocks
        totalblocks982 = totalblocks982 + 1;
    end
    blocklims982(runcount982,2) = totalblocks982;
end

skips982 = [14:18 26:28 38:40 42:44 54 56 57 63 65 68 77 78 81 94 95 97 99 106 108 109 123 124 126 143 149 ...
         152 155 156 157 170 173 174 177 178 180:182 184 205 207 221];

%account for the analysis that I removed:
skips982 = skips982 - 13;

skipv982 = zeros(totalblocks982,1); skipv982(skips982) = 1;   %make a vector of zeros with blocknum = 1 when block is skipped.

%% Pick through Pur data

purruns = 36:39;  %the indices in variable measnames from parsing code
samplesizesPur = [2106 2053 2007 2006];  %in pg
nrunsPur = length(purruns);

% to hold the count number of the starting block of each analysis
startblockPur = cell(4,length(purruns));
for i = 1:length(purruns)
    startblockPur{1,i} = measnames{purruns(i)};
end

totalblocksPur = 0;
runcountPur = 0;
blocklimsPur = zeros(nrunsPur,2);
for runi = purruns
    runcountPur = runcountPur + 1;
    mdata = evalin('base', ['Amelin_' measnames{runi}]);
    nblocks = size(mdata,1);
    startblockPur{3,runcountPur} = nblocks;
    startblockPur{2,runcountPur} = totalblocksPur + 1;
    blocklimsPur(runcountPur,1) = totalblocksPur + 1;
    for i = 1:nblocks
        totalblocksPur = totalblocksPur + 1;
    end
    blocklimsPur(runcountPur,2) = totalblocksPur;
end

skipsPur = [1 2 3 4 24 25 27 35 36 37 38 39];

skipvPur = zeros(totalblocksPur,1); skipvPur(skipsPur) = 1;   %make a vector of zeros with blocknum = 1 when block is skipped.

%%

m202 = 201.972159; m204 = 203.9730436; m205 = 204.9744818; m206 = 205.9744653; m207 = 206.9758969; m208 = 207.9766521;
masses = [m202 m204 m205 m206 m207 m208];


%% initialize matrices

totalskips = length(skipsPur) + length(skips981) + length(skips982);
totalblocks = totalblocks981 + totalblocks982 + totalblocksPur;
totruns = runcount981 + runcount982 + runcountPur;

usedblocks = totalblocks - totalskips;
used981s = totalblocks981 - length(skips981);
used982s = totalblocks982 - length(skips982);
usedPurs = totalblocksPur - length(skipsPur);

totms = 20 + 4*totruns + 2*usedblocks;
totds = 6*usedblocks;

uCM = zeros(totms,totms);
uCDd = zeros(totds,1);
uGn = zeros(totds,totms);

%um = zeros(totms,1);
umprior = zeros(totms,1);
udobs = zeros(totds,1);
ufm = zeros(totds,1);

dobs981 = zeros(used981s*6,1);
CDdiag981 = zeros(used981s*6,1);
dobs982 = zeros(used982s*6,1);
CDdiag982 = zeros(used982s*6,1);
dobsPur = zeros(usedPurs*6,1);
CDdiagPur = zeros(usedPurs*6,1);

tic
%% Start MC madness
nM = 1;
ummat = zeros(totms,nM);
twoSm = zeros(nM,1);
mswds = zeros(nM,1);
for M = 1:nM

%% some initial matters

%Assumption 1: 'Minor' IC of Amelin tracer
r25t0 = 1.019836994384;   %from fmincon results for all 981
r45t = 1.09*10^-4;  r42t = r45t/r25t0;
r65t = 9.91*10^-4;  r62t = r65t/r25t0;
r75t = 8.46*10^-4;  r72t = r75t/r25t0;
r85t = 2.054*10^-3; r82t = r85t/r25t0;
r25t = r25t0;

% make a covariance matrix for above assuming 1% measurement unct. and 10%
% blank correction uncertainty
molPb202t = 1;   %not really needed
 molPb204t = r42t;      molPb206t = r62t;      molPb207t = r72t;      molPb208t = r82t;
smolPb204t = r42t/100; smolPb206t = r62t/100; smolPb207t = r72t/100; smolPb208t = r82t/100; 
smolPb204b = r42t/10;  %blank for tracer measurement is +/- 10% of total 204 from tracer
r64b = 18.574;
r74b = 15.445;
r84b = 37.514;   %from email from Yuri
r46b = 1/r64b; r76b = r74b/r64b; r86b = r84b/r64b;
covmolb = diag([smolPb204t^2 smolPb206t^2 smolPb207t^2 smolPb208t^2 smolPb204b^2]);
%           %molb204t molPb206t molPb207t molPb208t molPb204b
jrat_molb = [ 1         0           0           0        1;     %r42t
              0         1           0           0       r64b;   %r62t
              0         0           1           0       r74b;   %r72t
              0         0           0           1       r84b];  %r82t
covt = jrat_molb*covmolb*jrat_molb';   %covariance matrix for tracer ratios

%make a separate covariance matrix for the blank subtraction during the runs
sr64b = 0.592;  sr74b = 0.417; sr84b = 1.188;
rho6474b = 0.755; rho7484b = 0.864; rho6484b = 0.729;  %these are MIT rhos
covb4 = [sr64b^2 rho6474b*sr64b*sr74b rho6484b*sr64b*sr84b;
         rho6474b*sr64b*sr74b sr74b^2 rho7484b*sr74b*sr84b;
         rho6484b*sr64b*sr84b rho7484b*sr74b*sr84b sr84b^2];
Jblank4to6 = [-1/r64b^2    0      0;
              -r74b/r64b^2 1/r64b 0;
              -r84b/r64b^2 0      1/r64b];
covb6 = Jblank4to6 * covb4 * Jblank4to6';  %covariance matrix for [r46b r76b r86b]

%Assumption 2: IC of BaPO2 and Tl
aBa = 0.001;
aTl = 0.001;
r15Ba = 236933 * (1+4*aBa);  %convert fractionation-corrected to measured ratios
r25Ba = 333.008 * (1+3*aBa);
r35Ba = 973.441 * (1+2*aBa);
r45Ba = 0.52739 * (1+aBa);
r35Tl = 0.41884 * (1+2*aTl);

covBa = (0.01*diag([r15Ba r25Ba r45Ba])).^2;  % 1% uncertainties on BaPO2 ratios.

%Implicit Assumptions:
m202 = 201.972159; m204 = 203.9730436; m205 = 204.9744818; m206 = 205.9744653; m207 = 206.9758969; m208 = 207.9766521;

%from run-wise total inversion
r46g_981 = 0.0590075645109563;
r76g_981 = 0.914686223556091;
r86g_981 = 2.1681;   % RECALIBRATION: CHANGE 208/206 OF NBS981 HERE
sr86g_981 = 0.0004;

r46g_982 = 0.027205;
r76g_982 = 0.46693;
r86g_982 = 1.00016;

r46g_Pur = 0.054875;
r76g_Pur = 0.85667;
r86g_Pur = 2.1024;

%% results from run-by-run TI solutions

%981:
r62bts_981 = [0.000202057 0.000202094 0.000202049 0.000202032 0.000202045 0.000657953 0.000298934 0.000298979 0.001067412 0.000149499];
r62gts_981 = [1.775190364	1.775392258	1.777587222	1.776530687	1.776299485	1.773641785	2.676274098	2.678620371	2.69852304	1.296654845];
gamma205_981 = [0.998532278,0.998933892,0.999608621,0.999279666,0.999206740,0.999617876,0.999584161,0.999421416,0.999619743,0.998753948];
gamma207_981 = [0.998014009,0.997592927,0.999487988,0.999419468,0.998965008,0.999372064,0.999494406,0.999522242,0.999557384,0.998903404]; 

r52BaPbs_981 = [2.06790335778960e-07 2.67376140294990e-07 1.04744755591200e-07 3.10173175467470e-07 4.53758075096200e-07 2.60262667073770e-07 2.97585842793020e-07 5.44130305307020e-07 3.71398802544580e-07 4.67376367147700e-07;
1.31306850370350e-06 6.55347667359980e-07 1.16871139981030e-07 3.10373139426440e-07 5.17130617124290e-07 2.29363488400680e-07 2.89173511374060e-07 4.40816526904950e-07 3.12092605276380e-07 4.32166886139460e-07;
1.25607941382860e-06 7.80910773379350e-07 1.47702978561570e-07 3.27063342738510e-07 5.80824145983220e-07 2.32899157657810e-07 3.23038555110700e-07 4.09291807357690e-07 2.84337291579870e-07 4.65242729842430e-07;
1.22277828086170e-06 8.42443535106770e-07 2.00415528226230e-07 3.59801171115130e-07 6.24878600609290e-07 2.64433492650380e-07 4.05039046616390e-07 4.26298210899070e-07 2.86601617075950e-07 5.31744595141290e-07;
1.31686346844970e-06 8.20691257747290e-07 2.81550339890210e-07 4.10198609812120e-07 6.44744551195860e-07 3.31062365803110e-07 5.62370127435430e-07 5.11229483786670e-07 3.25809340622300e-07 6.24417683347970e-07;
1.61387803529710e-06 8.01243686021230e-07 4.12312897808300e-07 4.83799420083360e-07 6.59638260183210e-07 4.56734541192720e-07 8.65540431569980e-07 7.15489962346040e-07 4.15265912818860e-07 7.55751986455380e-07;
2.11853912000960e-06 8.84909233526270e-07 6.50847644291530e-07 5.82524325675790e-07 6.91385759450590e-07 7.05939265920730e-07 1.46325260823090e-06 1.16370364629780e-06 5.93453672297890e-07 9.49942969007230e-07;
2.92355842998490e-06 1.12464461741290e-06 1.01896263709240e-06 7.24951792062470e-07 8.24094601120180e-07 1.21067564241130e-06 2.71313691749500e-06 2.10601727099740e-06 9.70709158880370e-07 1.16675416372350e-06;
8.95579674322550e-06 1.61900231301700e-06 1.47048303261970e-06 9.18756273788860e-07 1.09998961779130e-06 2.29201544372630e-06 5.51909479736500e-06 4.12610536304490e-06 1.78501817449730e-06 1.43387268445920e-06;
1.20433061119640e-05 2.46358154671030e-06 1.84431386210470e-06 1.22075228542480e-06 1.59739511773830e-06 4.74430971117140e-06 1.23510710950900e-05 8.51688608716480e-06 3.54525319480550e-06 1.71461867528620e-06;
1.44859805886960e-05 8.66616439501890e-06 1.88985922650140e-06 1.68654993182750e-06 2.47690921232730e-06 9.71529983424140e-06 2.92575670979720e-05 1.74256299293890e-05 7.37562308109030e-06 1.94496273757860e-06;
1.56457107556220e-05 9.44469207806790e-06 1.72797171921550e-06 2.43012209669040e-06 3.80421613593060e-06 1.97790166745090e-05 0 3.29184100124860e-05 0 2.15796113470270e-06;
1.36480266228610e-05 8.89990752011030e-06 1.70189472874920e-06 3.75872683431790e-06 5.38461119812470e-06 3.94322675737130e-05 0 0 0 2.37832047736690e-06;
1.02316266834590e-05 7.29926904555280e-06 1.66262539519190e-06 5.97315687754790e-06 6.97244812500780e-06 7.44336528714620e-05 0 0 0 2.53864027269470e-06;
9.60600057546000e-06 6.64340434838620e-06 1.77833729420850e-06 9.53892715774270e-06 8.60225981769900e-06 0.000123617463007510 0 0 0 3.02594483062970e-06;
1.18126785252120e-05 5.63253257768820e-06 2.25926634339060e-06 1.47087004195410e-05 9.92837650733390e-06 0.000165846271586020 0 0 0 4.16849879708400e-06;
1.65027443521870e-05 5.83287840739030e-06 0 2.10605896510290e-05 9.55379809491240e-06 0 0 0 0 5.77927469549700e-06;
0 6.72747809339430e-06 0 2.65985597917270e-05 9.27273248591440e-06 0 0 0 0 0;
0 1.02796811956780e-05 0 3.56376620383760e-05 1.03717861695290e-05 0 0 0 0 0;
0 0 0 4.95714220372870e-05 1.13682636789090e-05 0 0 0 0 0;
0 0 0 5.33301429280280e-05 0 0 0 0 0 0];

betas_981 = [-0.222041853799820 -0.157471050590690 -0.300734194334310 -0.301677788268410 -0.306946516214370 -0.347180665783520 -0.337203775795200 -0.352512765102230 -0.353461937963190 -0.296330377418830;
-0.271305230446900 -0.228346405968370 -0.286405055563870 -0.302757689337810 -0.316563891559240 -0.334666556440150 -0.320637472165090 -0.351577482957230 -0.343045638098250 -0.329672916267860;
-0.272333706547180 -0.248148258860030 -0.272589881389450 -0.310868542161330 -0.323633297072540 -0.320189476391820 -0.303949689759170 -0.333343247471600 -0.325852148321060 -0.338209810540750;
-0.275583906595450 -0.259987122928890 -0.247128270329910 -0.313453068622230 -0.321338092585630 -0.304393753101900 -0.285389122286610 -0.315894212351820 -0.304408383351690 -0.334683938915500;
-0.266860395352220 -0.261637647842810 -0.220763032584570 -0.311322704060890 -0.319360400358890 -0.285596332949280 -0.264844352042020 -0.297083440183570 -0.285480272965600 -0.336183171824490;
-0.255772037065630 -0.251237975195270 -0.191910995786750 -0.312406589528380 -0.306503398251370 -0.263009772217990 -0.235707626746260 -0.272078107714830 -0.263346538487560 -0.323802030531080;
-0.233974898360390 -0.241779546483070 -0.170398880364600 -0.302510948411470 -0.292480138075730 -0.234442027864640 -0.197158097560390 -0.237436136727140 -0.240786971344050 -0.316133095045900;
-0.215656426431430 -0.226528171209350 -0.140364833291070 -0.293402982758770 -0.280420453158490 -0.192283181975960 -0.145602002866530 -0.190940458215180 -0.214968118975610 -0.303045913492380;
-0.183872296731760 -0.208964679529390 -0.117644656984050 -0.283902039219090 -0.262928656929430 -0.139157495782640 -0.0849698563568200 -0.138269090313510 -0.165884177376350 -0.285928058781150;
-0.181700771700120 -0.195781382717910 -0.101053870369960 -0.263085535719410 -0.247405786746540 -0.0825628319971200 -0.0253762486244600 -0.0816550242063300 -0.112736544957730 -0.261037169353260;
-0.182207755895110 -0.217967619227410 -0.0867853922552600 -0.236659715125440 -0.225212769427700 -0.0373164722546600 0.0433256011456200 -0.0427817809432400 -0.0527175346933700 -0.234711364384730;
-0.176853989758020 -0.212840477797620 -0.0782726585004300 -0.199577637815590 -0.202063396334980 0.00736649547481000 0 -0.0232736924017400 0 -0.205042588808080;
-0.162148884058110 -0.198086776285780 -0.0769229272634700 -0.159291379223770 -0.175918033246730 0.0615890770807200 0 0 0 -0.173005070625100;
-0.148797851761140 -0.175816433262200 -0.0798807096413300 -0.116064198939440 -0.160097976391220 0.132202244708830 0 0 0 -0.156203951992520;
-0.145156265937350 -0.156175027184920 -0.0907759920791700 -0.0802485658108300 -0.141294278784140 0.193147763079620 0 0 0 -0.137036592033280;
-0.133025791614860 -0.143162611427070 -0.112126675366450 -0.0341933799145300 -0.125281978303360 0.276038931750730 0 0 0 -0.124946195562690;
-0.139582481373110 -0.137005167958650 0 0.00567722592076000 -0.110010924284790 0 0 0 0 -0.123350536449690;
0 -0.138298915678870 0 0.0446585567541000 -0.101186884123100 0 0 0 0 0;
0 -0.143421391660420 0 0.0841080576062300 -0.0934590554403000 0 0 0 0 0;
0 0 0 0.101168769586760 -0.102354756419550 0 0 0 0 0;
0 0 0 0.0382176605457300 0 0 0 0 0 0];

%982:

r62bts_982 = [0.00034354 0.000340747 0.000350033 0.00035319 0.000350871 0.000364034 0.000411678 0.000318802 0.000364837];
r62gts_982 = [3.916726659	3.898989811	4.572096438	3.657922059	4.198432099	4.358654529	5.117788261	4.820270959	5.095301614];
gamma205_982 = [1.00013126 1.000957 1.000175632 0.999721633 0.998666101 0.999586058 0.99917781 0.999904674 0.998902807];
gamma207_982 = [0.998964111 0.999020637 0.998873087 0.999063628 0.999797299 0.998804496 0.998742335 0.998517688 0.99923092];

r52BaPbs_982 = 2*10^-9*ones(24,9);   %zeros(24,9);  %no 201 measured for NBS982.

betas_982 = [-0.327190184 -0.285598901 -0.220363248 -0.229472315 -0.265631685 -0.307167914 -0.318557531 -0.298345739 -0.288010334;
-0.326668294 -0.276621098 -0.309742942 -0.265092869 -0.348603426 -0.306173451 -0.342504168 -0.321000626 -0.299557386;
-0.323908061 -0.272586357 -0.299430458 -0.310520116 -0.341536016 -0.300840498 -0.354521296 -0.337641350 -0.319832850;
-0.303126604 -0.265421601 -0.298224970 -0.319066181 -0.346298603 -0.302152345 -0.341037915 -0.326590392 -0.323204329;
-0.288761242 -0.253102011 -0.295416795 -0.319794683 -0.326196689 -0.290924952 -0.351068867 -0.320872870 -0.333503988;
-0.268524014 -0.246196355 -0.293356831 -0.319678479 -0.330091103 -0.296724800 -0.339570829 -0.315222317 -0.335784817;
-0.238796006 -0.231432052 -0.285891849 -0.313895082 -0.325154424 -0.286350038 -0.344760461 -0.312747455 -0.328813785;
0 -0.211430345 -0.281476475 -0.317546679 -0.330520940 -0.293490868 -0.338143064 -0.308557386 -0.310492125;
0 -0.185190707 -0.268368726 -0.320055580 -0.331816030 -0.284694285 -0.339564553 -0.299385932 -0.306352360;
0 0 -0.267200313 -0.306619450 -0.320572658 -0.281610422 -0.332777756 -0.293139118 -0.308670803;
0 0 -0.253916270 -0.310089268 -0.323202247 -0.271438610 -0.321541263 -0.286472859 -0.290089816;
0 0 -0.224722093 -0.294678060 -0.307484038 -0.268177962 -0.321520868 -0.275213878 -0.297172237;
0 0 -0.220384973 -0.285678470 -0.318518259 -0.265417578 -0.313248894 -0.267574740 -0.280781282;
0 0 -0.212541524 -0.277312961 -0.314043363 -0.262947615 -0.304878731 -0.261623202 -0.286324569;
0 0 -0.200032497 -0.269535913 -0.311235583 -0.257114445 -0.296883421 -0.244180677 -0.276565881;
0 0 -0.197687520 -0.265092701 -0.303109017 -0.242758505 -0.288878975 -0.232419085 -0.269606956;
0 0 0 -0.264294116 -0.297101458 -0.244767425 -0.278824061 -0.217930756 -0.255422881;
0 0 0 -0.259161275 -0.286202271 -0.234389842 -0.271178981 -0.203976573 -0.249280408;
0 0 0 -0.257343690 -0.294406075 -0.223121685 0 0 -0.242574812;
0 0 0 -0.255868661 -0.284232364 -0.204184888 0 0 -0.226079772;
0 0 0 -0.249804634 -0.272756602 -0.198738968 0 0 -0.203817497;
0 0 0 -0.243037717 -0.259745619 -0.191523610 0 0 -0.189902771;
0 0 0 -0.240242977 -0.251920442 0 0 0 -0.177645835;
0 0 0 -0.230353360 0 0 0 0 0];

% Puratronic Pb

r62bts_Pur = [0.000211752 0.000214023 0.000212703 0.000211457];
r62gts_Pur = [2.964798423 2.924940621 2.839375982 2.822423477];
gamma205_Pur = [0.999519163 0.997793624 0.999556094 0.998945155];
gamma207_Pur = [0.999564945 0.998819854 0.999076807 0.99893405];

r52BaPbs_Pur = [4.33345e-09 4.00225e-09 4.68263e-09 3.58440e-09;
5.11541e-09 4.06224e-09 5.46193e-09 2.62389e-09;
6.42332e-09 4.15584e-09 5.63844e-09 2.83700e-09;
8.48288e-09 4.48006e-09 6.33304e-09 3.09913e-09;
1.24890e-08 4.82191e-09 6.82955e-09 3.29666e-09;
2.03021e-08 6.44621e-09 8.39333e-09 4.44294e-09;
3.56596e-08 9.12833e-09 1.03421e-08 6.22072e-09;
6.45656e-08 1.35561e-08 1.37473e-08 9.58972e-09;
1.19435e-07 2.01205e-08 0 1.51908e-08;
0 2.90812e-08 0 0];

betas_Pur = [-0.273495904 -0.364502017 -0.304414880 -0.365614681;
-0.266946761 -0.354953358 -0.291839056 -0.347864899;
-0.251718651 -0.339934705 -0.292271646 -0.343957770;
-0.232761224 -0.333122551 -0.286029738 -0.335942074;
-0.199722055 -0.319658937 -0.275397106 -0.320295241;
-0.158950973 -0.311832191 -0.264914179 -0.304246209;
-0.124811018 -0.300972084 -0.248688672 -0.280597010;
-0.0824069990 -0.297885062 -0.226806685 -0.251180591;
-0.0416005890 -0.292633137 0 -0.216681106;
0 -0.277456528 0 0];

% now agglomerate

r62bts = [r62bts_981 r62bts_982 r62bts_Pur];
r62gts = [r62gts_981 r62gts_982 r62gts_Pur];
gamma205s = [gamma205_981 gamma205_982 gamma205_Pur];
gamma207s = [gamma207_981 gamma207_982 gamma207_Pur];

maxblocks = max([size(betas_981,1) size(betas_982,1) size(betas_Pur,1)]);

r52BaPbs_981pad = padarray(r52BaPbs_981, [maxblocks-size(r52BaPbs_981,1) 0], 'post');
r52BaPbs_982pad = padarray(r52BaPbs_982, [maxblocks-size(r52BaPbs_982,1) 0], 'post');
r52BaPbs_Purpad = padarray(r52BaPbs_Pur, [maxblocks-size(r52BaPbs_Pur,1) 0], 'post');
r52BaPbs = [r52BaPbs_981pad r52BaPbs_982pad r52BaPbs_Purpad];

betas_981pad = padarray(betas_981, [maxblocks-size(betas_981,1) 0], 'post');
betas_982pad = padarray(betas_982, [maxblocks-size(betas_982,1) 0], 'post');
betas_Purpad = padarray(betas_Pur, [maxblocks-size(betas_Pur,1) 0], 'post');
betas = [betas_981pad betas_982pad betas_Purpad];

% define some start points
btstart = 20;
gtstart = 20 + totruns;
g5start = 20 + 2*totruns;
g7start = 20 + 3*totruns;
bpstart = 20 + 4*totruns;
bestart = 20 + 4*totruns + usedblocks;

%% Monte Carlo seeding

% give MC-valued parameters near-delta function priors
covMC = blkdiag(covt, covBa, covb6, sr86g_981^2);
covt = covt*10^-6;
covb6 = covb6*10^-6;
covBa = covBa*10^-6;
sr86g_981 = sr86g_981*10^-3;

    %NOTE: for the mean value, only the ML estimates are used a single time.
    MCreal = [r42t r62t r72t r82t r15Ba r25Ba r45Ba r46b r76b r86b r86g_981];
    r42t = MCreal(1); r62t = MCreal(2); r72t = MCreal(3); r82t = MCreal(4);
    r15Ba = MCreal(5); r25Ba = MCreal(6); r45Ba = MCreal(7); 
    r46b = MCreal(8); r76b = MCreal(9); r86b = MCreal(10); r86g_981 = MCreal(11);

umprior(1:20) = [r25t0 r42t r62t r72t r82t r15Ba r25Ba r45Ba r46b r76b r86b r46g_981 r76g_981 r86g_981...
    r46g_982 r76g_982 r86g_982 r46g_Pur r76g_Pur r86g_Pur];

CM20 = blkdiag(1, covt, covBa, covb6, 0.05^2, 1, sr86g_981^2, 0.03^2, 0.5^2, 1, 0.05^2, 1, 2);
uCM(1:20,1:20) = CM20;
%%

rct = 0;   %runcount total
bct = 0;   %blockcount total
runcount981 = 0; runcount982 = 0; runcountPur = 0;   %for skipblocks only.
for runi = [nbs981runs nbs982runs purruns]
 
rct = rct + 1;

if sum(runi == nbs981runs)
    runcount981 = runcount981 + 1;
    skipvi = skipv981(blocklims981(runcount981,1):blocklims981(runcount981,2));
    curstd = 'nbs981';
    gnGravRange = 12:14;
    r46g = umprior(12); r76g = umprior(13); r86g = umprior(14);
elseif sum(runi == nbs982runs)
    runcount982 = runcount982 + 1;
    skipvi = skipv982(blocklims982(runcount982,1):blocklims982(runcount982,2));
    curstd = 'nbs982';
    gnGravRange = 15:17;
    r46g = umprior(15); r76g = umprior(16); r86g = umprior(17);
else
    runcountPur = runcountPur + 1;
    skipvi = skipvPur(blocklimsPur(runcountPur,1):blocklimsPur(runcountPur,2));
    curstd = 'pur';
    gnGravRange = 18:20;    
    r46g = umprior(18); r76g = umprior(19); r86g = umprior(20);
end

mdata = evalin('base', ['Amelin_' measnames{runi}]);

%trim skipped blocks off of mdata
nblocks = size(mdata,1);
for i = nblocks:-1:1
    if skipvi(i) == 1;
        mdata(i,:) = [];
    end
end
nblocks = size(mdata,1);

%un-pack, re-pack run-wise prior data
r62bt = r62bts(rct);        umprior(btstart+rct) = r62bt;       uCM(btstart+rct,btstart+rct) = (0.25*r62bt)^2;
r62gt =  r62gts(rct);       umprior(gtstart+rct) = r62gt;       uCM(gtstart+rct,gtstart+rct) = r62gt^2;
gamma205 = gamma205s(rct);  umprior(g5start+rct) = gamma205;    uCM(g5start+rct,g5start+rct) = 1;
gamma207 = gamma207s(rct);  umprior(g7start+rct) = gamma207;    uCM(g7start+rct,g7start+rct) = 1;

for i = 1:nblocks  %just to calculate r62gt
    bct = bct + 1;

    %udobs and uCDd for the block
    r16m = mdata(i,4)*mdata(i,6);  udobs(6*bct-5) = r16m;  uCDd(6*bct-5) = (mdata(i,5) *r16m/100)^2;
    if mdata(i,4) == -7.092198581560280e+02  %for 982 data, no 201 measured
    r16m = 1*10^-4;                udobs(6*bct-5) = r16m;  uCDd(6*bct-5) = 10^-8;
    end
    r26m = mdata(i,6);             udobs(6*bct-4) = r26m;  uCDd(6*bct-4) = (mdata(i,7) *r26m/100)^2; 
    r46m = mdata(i,10);            udobs(6*bct-3) = r46m;  uCDd(6*bct-3) = (mdata(i,11)*r46m/100)^2; 
    r56m = mdata(i,12);            udobs(6*bct-2) = r56m;  uCDd(6*bct-2) = (mdata(i,13)*r56m/100)^2; 
    r76m = mdata(i,14);            udobs(6*bct-1) = r76m;  uCDd(6*bct-1) = (mdata(i,15)*r76m/100)^2; 
    r86m = mdata(i,17);            udobs(6*bct)   = r86m;  uCDd(6*bct)   = (mdata(i,18)*r86m/100)^2;

    %unpack, re-pack block-wise mprior data
    r52BaPb = r52BaPbs(i,rct);  umprior(bpstart+bct) = r52BaPb;  uCM(bpstart+bct,bpstart+bct) = 1;
    if r16m == 1*10^-4
        uCM(bpstart+bct,bpstart+bct) = 1;  %was 10^6... %if none there, enforce
    end
    beta = betas(i,rct);        umprior(bestart+bct) = beta;     uCM(bestart+bct,bestart+bct) = 1;
    
    %evaluate g(m) [= ufm] for the block
    gRange = (6*bct-5):(6*bct);
    ufm(gRange) = [(r15Ba*r52BaPb)/(r62bt + r62gt + r62t), (m202/m206)^beta/(r62bt + r62gt + r62t) + (r25Ba*r52BaPb)/(r62bt + r62gt + r62t), ...
        (r45Ba*r52BaPb)/(r62bt + r62gt + r62t) + ((m204/m206)^beta*(r42t + r46b*r62bt + r46g*r62gt))/(r62bt + r62gt + r62t), ((gamma205*m205)/m206)^beta/(r25t*(r62bt + r62gt + r62t)) + ...
        r52BaPb/(r62bt + r62gt + r62t), (((gamma207*m207)/m206)^beta*(r72t + r62bt*r76b + r62gt*r76g))/(r62bt + r62gt + r62t), ((m208/m206)^beta*(r82t + r62bt*r86b + r62gt*r86g))/(r62bt + r62gt + r62t)];

    %evaluate derivatives for Gm
    uGn(gRange,1:11) = [0, 0, -((r15Ba*r52BaPb)/(r62bt + r62gt + r62t)^2), 0, 0, r52BaPb/(r62bt + r62gt + r62t), 0, 0, 0, 0, 0; 0, 0, -(((m202/m206)^beta + r25Ba*r52BaPb)/(r62bt + r62gt + r62t)^2), 0, 0, 0, ...
        r52BaPb/(r62bt + r62gt + r62t), 0, 0, 0, 0; 0, (m204/m206)^beta/(r62bt + r62gt + r62t), -((r45Ba*r52BaPb + (m204/m206)^beta*(r42t + r46b*r62bt + r46g*r62gt))/(r62bt + r62gt + r62t)^2), 0, 0, 0, 0, ...
        r52BaPb/(r62bt + r62gt + r62t), ((m204/m206)^beta*r62bt)/(r62bt + r62gt + r62t), 0, 0; -(((gamma205*m205)/m206)^beta/(r25t^2*(r62bt + r62gt + r62t))), 0, ...
        -((((gamma205*m205)/m206)^beta + r25t*r52BaPb)/(r25t*(r62bt + r62gt + r62t)^2)), 0, 0, 0, 0, 0, 0, 0, 0; ...
        0, 0, -((((gamma207*m207)/m206)^beta*(r72t + r62bt*r76b + r62gt*r76g))/(r62bt + r62gt + r62t)^2), ((gamma207*m207)/m206)^beta/(r62bt + r62gt + r62t), 0, 0, 0, 0, 0, ...
        (((gamma207*m207)/m206)^beta*r62bt)/(r62bt + r62gt + r62t), 0; 0, 0, -(((m208/m206)^beta*(r82t + r62bt*r86b + r62gt*r86g))/(r62bt + r62gt + r62t)^2), 0, (m208/m206)^beta/(r62bt + r62gt + r62t), 0, ...
        0, 0, 0, 0, ((m208/m206)^beta*r62bt)/(r62bt + r62gt + r62t)]; %derivatives wrt global parameters
    uGn(gRange,gnGravRange) = [0, 0, 0; 0, 0, 0; ((m204/m206)^beta*r62gt)/(r62bt + r62gt + r62t), 0, 0; 0, 0, 0; ...
        0, (((gamma207*m207)/m206)^beta*r62gt)/(r62bt + r62gt + r62t), 0; ...
        0, 0, ((m208/m206)^beta*r62gt)/(r62bt + r62gt + r62t)]; %derivatives wrt gravimetric IC
    uGn(gRange,btstart+rct) = [-((r15Ba*r52BaPb)/(r62bt + r62gt + r62t)^2), -(((m202/m206)^beta + r25Ba*r52BaPb)/(r62bt + r62gt + r62t)^2), ((-r45Ba)*r52BaPb + (m204/m206)^beta*(-r42t - r46g*r62gt + r46b*(r62gt + r62t)))/...
        (r62bt + r62gt + r62t)^2, -((((gamma205*m205)/m206)^beta + r25t*r52BaPb)/(r25t*(r62bt + r62gt + r62t)^2)), (((gamma207*m207)/m206)^beta*(-r72t + (r62gt + r62t)*r76b - r62gt*r76g))/...
        (r62bt + r62gt + r62t)^2, ((m208/m206)^beta*(-r82t + (r62gt + r62t)*r86b - r62gt*r86g))/(r62bt + r62gt + r62t)^2]; %derivatives wrt r62bt
    uGn(gRange,gtstart+rct) = [-((r15Ba*r52BaPb)/(r62bt + r62gt + r62t)^2), -(((m202/m206)^beta + r25Ba*r52BaPb)/(r62bt + r62gt + r62t)^2), ((-r45Ba)*r52BaPb + (m204/m206)^beta*(-r42t - r46b*r62bt + r46g*(r62bt + r62t)))/...
        (r62bt + r62gt + r62t)^2, -((((gamma205*m205)/m206)^beta + r25t*r52BaPb)/(r25t*(r62bt + r62gt + r62t)^2)), (((gamma207*m207)/m206)^beta*(-r72t - r62bt*r76b + (r62bt + r62t)*r76g))/...
        (r62bt + r62gt + r62t)^2, ((m208/m206)^beta*(-r82t - r62bt*r86b + (r62bt + r62t)*r86g))/(r62bt + r62gt + r62t)^2]; %derivatives wrt r62gt
    uGn(gRange,g5start+rct) = [0, 0, 0, (beta*((gamma205*m205)/m206)^beta)/(gamma205*r25t*(r62bt + r62gt + r62t)), 0, 0]; %derivatives wrt gamma205
    uGn(gRange,g7start+rct) = [0, 0, 0, 0, (beta*((gamma207*m207)/m206)^beta*(r72t + r62bt*r76b + r62gt*r76g))/(gamma207*(r62bt + r62gt + r62t)), 0]; %derivatives wrt gamma207
    uGn(gRange,bpstart+bct) = [r15Ba/(r62bt + r62gt + r62t), r25Ba/(r62bt + r62gt + r62t), r45Ba/(r62bt + r62gt + r62t), 1/(r62bt + r62gt + r62t), 0, 0]; %derivatives wrt r52BaPb
    uGn(gRange,bestart+bct) = [0, ((m202/m206)^beta*log(m202/m206))/(r62bt + r62gt + r62t), ((m204/m206)^beta*(r42t + r46b*r62bt + r46g*r62gt)*log(m204/m206))/(r62bt + r62gt + r62t), ...
        (((gamma205*m205)/m206)^beta*log((gamma205*m205)/m206))/(r25t*(r62bt + r62gt + r62t)), (((gamma207*m207)/m206)^beta*(r72t + r62bt*r76b + r62gt*r76g)*log((gamma207*m207)/m206))/(r62bt + r62gt + r62t), ...
        ((m208/m206)^beta*(r82t + r62bt*r86b + r62gt*r86g)*log(m208/m206))/(r62bt + r62gt + r62t)]; %derivatives wrt beta

end

end

%% Set up total inversion parameters, get party started

iters = 1; totiters = 100;
%twoSm = zeros(1,totiters);
%mswds = zeros(1,totiters);
%ummat = zeros(totms, totiters);

nun = 0.9;
%invCM = inv(uCM); 
invCD = (diag(1./uCDd));
um = umprior;
%gamman = CM*Gn'*invCD*(gm-dobs)+(mn-mprior);
%bn = Gn*gamman;
%dmping = (gamman'/CM*gamman)/(gamman'/CM*gamman + bn'*invCD*bn);
%dmpingv(j) = dmping;
%mn  = mn - (eye(size(CM,1)) + CM*Gn'*invCD*Gn)\gamman;
%mn = mn - nun * dmping*(CM*Gn'*(CD\(gm-dobs))+(mn-mprior)); %equation 3.89 of Tarantola (Inv. Prob. Theory)
%um = um - nun*(invCM + uGn'*invCD*uGn)\(uGn'*invCD*(ufm-udobs)+invCM*(um-umprior)); %#ok<MINV> %Newton
um = um - nun*(eye(size(uCM,1))+uCM*uGn'*invCD*uGn)\(uCM*uGn'*((ufm-udobs)./uCDd)+um-umprior);
%ummat(:,iters) = um;
%twoSm(iters) = (ufm-udobs)'*invCD*(ufm-udobs);  %eqn 5.141 of Tarantola InvPrblmTheory
%mswds(iters) = twoSm(iters)/(length(udobs)-1);

%% do the big ol' iteration
for iters = 2:totiters  %just did number 1

rct = 0;   %runcount total
bct = 0;   %blockcount total
runcount981 = 0; runcount982 = 0; runcountPur = 0;   %for skipblocks only.
for runi = [nbs981runs nbs982runs purruns]
 
rct = rct + 1;

if sum(runi == nbs981runs)
    runcount981 = runcount981 + 1;
    skipvi = skipv981(blocklims981(runcount981,1):blocklims981(runcount981,2));
    curstd = 'nbs981';
    gnGravRange = 12:14;
    r46g = um(12); r76g = um(13); r86g = um(14);
elseif sum(runi == nbs982runs)
    runcount982 = runcount982 + 1;
    skipvi = skipv982(blocklims982(runcount982,1):blocklims982(runcount982,2));
    curstd = 'nbs982';
    gnGravRange = 15:17;
    r46g = um(15); r76g = um(16); r86g = um(17);
else
    runcountPur = runcountPur + 1;
    skipvi = skipvPur(blocklimsPur(runcountPur,1):blocklimsPur(runcountPur,2));
    curstd = 'pur';
    gnGravRange = 18:20;    
    r46g = um(18); r76g = um(19); r86g = um(20);
end

% a few other global parameters
r25t = um(1); r42t = um(2); r62t = um(3); r72t = um(4); r82t = um(5);
r15Ba = um(6); r25Ba = um(7); r45Ba = um(8);
r46b = um(9); r76b = um(10); r86b = um(11);

mdata = evalin('base', ['Amelin_' measnames{runi}]);

%trim skipped blocks off of mdata
nblocks = size(mdata,1);
for i = nblocks:-1:1
    if skipvi(i) == 1;
        mdata(i,:) = [];
    end
end
nblocks = size(mdata,1);

%un-pack new estimates
r62bt    = um(btstart+rct);
r62gt    = um(gtstart+rct);
gamma205 = um(g5start+rct);
gamma207 = um(g7start+rct);

for i = 1:nblocks  %just to calculate r62gt
    bct = bct + 1;

    %unpack, re-pack block-wise mprior data
    r52BaPb = um(bpstart+bct);
    beta    = um(bestart+bct);
    
    %evaluate g(m) [= ufm] for the block
    gRange = (6*bct-5):(6*bct);
    ufm(gRange) = [(r15Ba*r52BaPb)/(r62bt + r62gt + r62t), (m202/m206)^beta/(r62bt + r62gt + r62t) + (r25Ba*r52BaPb)/(r62bt + r62gt + r62t), ...
        (r45Ba*r52BaPb)/(r62bt + r62gt + r62t) + ((m204/m206)^beta*(r42t + r46b*r62bt + r46g*r62gt))/(r62bt + r62gt + r62t), ((gamma205*m205)/m206)^beta/(r25t*(r62bt + r62gt + r62t)) + ...
        r52BaPb/(r62bt + r62gt + r62t), (((gamma207*m207)/m206)^beta*(r72t + r62bt*r76b + r62gt*r76g))/(r62bt + r62gt + r62t), ((m208/m206)^beta*(r82t + r62bt*r86b + r62gt*r86g))/(r62bt + r62gt + r62t)];

    %evaluate derivatives for Gm
    uGn(gRange,1:11) = [0, 0, -((r15Ba*r52BaPb)/(r62bt + r62gt + r62t)^2), 0, 0, r52BaPb/(r62bt + r62gt + r62t), 0, 0, 0, 0, 0; 0, 0, -(((m202/m206)^beta + r25Ba*r52BaPb)/(r62bt + r62gt + r62t)^2), 0, 0, 0, ...
        r52BaPb/(r62bt + r62gt + r62t), 0, 0, 0, 0; 0, (m204/m206)^beta/(r62bt + r62gt + r62t), -((r45Ba*r52BaPb + (m204/m206)^beta*(r42t + r46b*r62bt + r46g*r62gt))/(r62bt + r62gt + r62t)^2), 0, 0, 0, 0, ...
        r52BaPb/(r62bt + r62gt + r62t), ((m204/m206)^beta*r62bt)/(r62bt + r62gt + r62t), 0, 0; -(((gamma205*m205)/m206)^beta/(r25t^2*(r62bt + r62gt + r62t))), 0, ...
        -((((gamma205*m205)/m206)^beta + r25t*r52BaPb)/(r25t*(r62bt + r62gt + r62t)^2)), 0, 0, 0, 0, 0, 0, 0, 0; ...
        0, 0, -((((gamma207*m207)/m206)^beta*(r72t + r62bt*r76b + r62gt*r76g))/(r62bt + r62gt + r62t)^2), ((gamma207*m207)/m206)^beta/(r62bt + r62gt + r62t), 0, 0, 0, 0, 0, ...
        (((gamma207*m207)/m206)^beta*r62bt)/(r62bt + r62gt + r62t), 0; 0, 0, -(((m208/m206)^beta*(r82t + r62bt*r86b + r62gt*r86g))/(r62bt + r62gt + r62t)^2), 0, (m208/m206)^beta/(r62bt + r62gt + r62t), 0, ...
        0, 0, 0, 0, ((m208/m206)^beta*r62bt)/(r62bt + r62gt + r62t)]; %derivatives wrt global parameters
    uGn(gRange,gnGravRange) = [0, 0, 0; 0, 0, 0; ((m204/m206)^beta*r62gt)/(r62bt + r62gt + r62t), 0, 0; 0, 0, 0; ...
        0, (((gamma207*m207)/m206)^beta*r62gt)/(r62bt + r62gt + r62t), 0; ...
        0, 0, ((m208/m206)^beta*r62gt)/(r62bt + r62gt + r62t)]; %derivatives wrt gravimetric IC
    uGn(gRange,btstart+rct) = [-((r15Ba*r52BaPb)/(r62bt + r62gt + r62t)^2), -(((m202/m206)^beta + r25Ba*r52BaPb)/(r62bt + r62gt + r62t)^2), ((-r45Ba)*r52BaPb + (m204/m206)^beta*(-r42t - r46g*r62gt + r46b*(r62gt + r62t)))/...
        (r62bt + r62gt + r62t)^2, -((((gamma205*m205)/m206)^beta + r25t*r52BaPb)/(r25t*(r62bt + r62gt + r62t)^2)), (((gamma207*m207)/m206)^beta*(-r72t + (r62gt + r62t)*r76b - r62gt*r76g))/...
        (r62bt + r62gt + r62t)^2, ((m208/m206)^beta*(-r82t + (r62gt + r62t)*r86b - r62gt*r86g))/(r62bt + r62gt + r62t)^2]; %derivatives wrt r62bt
    uGn(gRange,gtstart+rct) = [-((r15Ba*r52BaPb)/(r62bt + r62gt + r62t)^2), -(((m202/m206)^beta + r25Ba*r52BaPb)/(r62bt + r62gt + r62t)^2), ((-r45Ba)*r52BaPb + (m204/m206)^beta*(-r42t - r46b*r62bt + r46g*(r62bt + r62t)))/...
        (r62bt + r62gt + r62t)^2, -((((gamma205*m205)/m206)^beta + r25t*r52BaPb)/(r25t*(r62bt + r62gt + r62t)^2)), (((gamma207*m207)/m206)^beta*(-r72t - r62bt*r76b + (r62bt + r62t)*r76g))/...
        (r62bt + r62gt + r62t)^2, ((m208/m206)^beta*(-r82t - r62bt*r86b + (r62bt + r62t)*r86g))/(r62bt + r62gt + r62t)^2]; %derivatives wrt r62gt
    uGn(gRange,g5start+rct) = [0, 0, 0, (beta*((gamma205*m205)/m206)^beta)/(gamma205*r25t*(r62bt + r62gt + r62t)), 0, 0]; %derivatives wrt gamma205
    uGn(gRange,g7start+rct) = [0, 0, 0, 0, (beta*((gamma207*m207)/m206)^beta*(r72t + r62bt*r76b + r62gt*r76g))/(gamma207*(r62bt + r62gt + r62t)), 0]; %derivatives wrt gamma207
    uGn(gRange,bpstart+bct) = [r15Ba/(r62bt + r62gt + r62t), r25Ba/(r62bt + r62gt + r62t), r45Ba/(r62bt + r62gt + r62t), 1/(r62bt + r62gt + r62t), 0, 0]; %derivatives wrt r52BaPb
    uGn(gRange,bestart+bct) = [0, ((m202/m206)^beta*log(m202/m206))/(r62bt + r62gt + r62t), ((m204/m206)^beta*(r42t + r46b*r62bt + r46g*r62gt)*log(m204/m206))/(r62bt + r62gt + r62t), ...
        (((gamma205*m205)/m206)^beta*log((gamma205*m205)/m206))/(r25t*(r62bt + r62gt + r62t)), (((gamma207*m207)/m206)^beta*(r72t + r62bt*r76b + r62gt*r76g)*log((gamma207*m207)/m206))/(r62bt + r62gt + r62t), ...
        ((m208/m206)^beta*(r82t + r62bt*r86b + r62gt*r86g)*log(m208/m206))/(r62bt + r62gt + r62t)]; %derivatives wrt beta

end   %for i = 1:nblocks   

end   %for runi = [nbs981runs nbs982runs purruns]

%um = um - nun*(invCM + uGn'*invCD*uGn)\(uGn'*invCD*(ufm-udobs)+invCM*(um-umprior)); %#ok<MINV>
um = um - nun*(eye(size(uCM,1))+uCM*uGn'*invCD*uGn)\(uCM*uGn'*((ufm-udobs)./uCDd)+um-umprior); %pre-conditioned steepest ascent
% ummat(:,iters) = um;
% twoSm(iters) = (ufm-udobs)'*invCD*(ufm-udobs);  %eqn 5.141 of Tarantola InvPrblmTheory
% mswds(iters) = twoSm(iters)/(length(udobs)-1);

end %for i = 1:totiters
ummat(:,M) = um;
twoSm(M) = (ufm-udobs)'*((ufm-udobs)./uCDd);  %eqn 5.141 of Tarantola InvPrblmTheory
mswds(M) = twoSm(M)/(length(udobs)-1);
%dlmwrite('MCintercal.txt', [um' twoSm(M)], '-append', 'precision', '%12.12g')

end %for M = 1:nM

toc