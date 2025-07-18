%% a last go at the gravimetric-tracer data reduction, as a big inverse problem (Tarantola)

% ceate MC, containing MC realizations of systematic varialbes, and mix, containing grav-trac mix raw data.
% first, make M monte carlo simulations of the gravimetric and tracer ICs, as well as the purity and masses for Gravs
% then use these to calculate M simulations of 206g/238g.  This suite of ICs and 206/238g are used next.

nM = 10^4;
nperblock = 20;  %number of analyses per 'block'

% Monte Carlo ppb impurities from GDMS data in Condon (code not in original publications)
[concsM_981, concsM_982, concsM_Pur] = concsM_MC(nM);

%% some constants and measurements
%masses from Audi et al. AME2003, Nuclear Physics A729 p. 337-676, December 22, 2003.
masses.Pb202 = 201.972159; 
masses.Pb204 = 203.9730436; 
masses.Pb205 = 204.9744818; 
masses.Pb206 = 205.9744653; 
masses.Pb207 = 206.9758969; 
masses.Pb208 = 207.9766521;
masses.U233 = 233.039635207;
masses.U234 = 234.040952088;
masses.U235 = 235.043929918;
masses.U238 = 238.050788247;

% NOTE: NOW SET IN ETUPBTEACERRECALIBRATIONSCRIPT.M
% %means and variances of tracer Pb and U ICs...
% %from BlankIC_LinearRegressionAll3.m and BlankIC_LinearRegressionAll3_25.m
% ics.ET535Pb  = [9.00000000000000e-05 0.000388739872278295 0.000296068587594844 0.000744376162255465]';
% ics.ET2535Pb = [0.000105000000000000 0.000482509449698359 0.000432369168809505 0.00104222718281000]';
% 
% % [204/205 206/205 207/205 208/205] for ET535 and ET2535
% covs.ET535Pb  = [8.10000000000000e-11,1.49142352775740e-09,1.24860181672645e-09,3.04674573549283e-09;
%                 1.49142352775740e-09,2.85137047384809e-08,2.34779694457233e-08,5.79108827115715e-08;
%                 1.24860181672645e-09,2.34779694457233e-08,1.96403112135437e-08,4.82761791708861e-08;
%                 3.04674573549283e-09,5.79108827115715e-08,4.82761791708861e-08,1.20435254888182e-07];
% covs.ET2535Pb = [8.38768660619565e-11,1.51796787106751e-09,1.25950212394299e-09,3.06017691531573e-09;
%                 1.51796787106751e-09,2.75700242967555e-08,2.28397877479297e-08,5.55511185313650e-08;
%                 1.25950212394299e-09,2.28397877479297e-08,1.89498284903739e-08,4.60743733525839e-08;
%                 3.06017691531573e-09,5.55511185313650e-08,4.60743733525839e-08,1.12190693941688e-07];

% using MC trials from TarantolaUIC_v2.m, MLE from TarantolaUIC_Mean_v1.m, covtot from end of parseMCintercal.m
% [238/235 (CRM 115), 238/235 (CRM 112a), 233/235 (ET535), 238/235 (ET535)]
ics.gtU = [491.5480646 137.8413622 0.9950621757 0.0030799297];
covs.gtU = [1.846512295277E-03	4.185382743287E-04	-2.075050000967E-06	8.994699181763E-09;
                4.185382743287E-04	1.393511775032E-04	-5.807484183676E-07	2.517364468976E-09;
                -2.075050000967E-06	-5.807484183676E-07	2.898366694260E-09	-1.264609494813E-11;
                8.994699181763E-09	2.517364468976E-09	-1.264609615520E-11	1.564583249266E-13];
ics.r234238_112a = 0.000052841;  %from certificate, for calculating atomic weight

% NOTE: NOW SET IN ETUPBTEACERRECALIBRATIONSCRIPT.M
% % using MC trials from AmelinTarantola_v7_AllMC.m, MLE from AmelinTarantola_v8_All_Mean.m, covtot from parseMCintercal.m
% %          NBS 981                  NBS 982               Puratronic Pb
% % [204/206 207/206 208/206  204/206 207/206 208/206  204/206 207/206 208/206]
% ics.gravPb = [0.0590073543 0.9146828587 2.1681000018 0.0272057597 0.4669667737 1.0002490230 0.0548860631 0.8567195631 2.1022709594];
% covs.gravPb = [1.22890447167639e-10	-7.97132415764876e-10	-4.38543453347033e-09	5.42279273487459e-11	-3.99133391062568e-10	-2.06089730970833e-09	1.13777146683553e-10	-7.65770629215464e-10	-4.29039948170844e-09
%             -7.97540074619466e-10	5.45712931993738e-09	2.92532272393163e-08	-3.69312326446656e-10	2.45639896831742e-09	1.30696831558687e-08	-7.45935131894574e-10	5.09350671828159e-09	2.83792753955173e-08
%             -4.38543454715791e-09	2.92532272784984e-08	1.59598225190498e-07	-1.99369387122933e-09	1.38978573311883e-08	7.30662714007827e-08	-4.09863442333872e-09	2.78007101883434e-08	1.55544096654634e-07
%             5.42228552356883e-11	-3.69359175624140e-10	-1.99369387520237e-09	2.71381631979341e-11	-1.61265031347176e-10	-8.85283882620290e-10	5.12717363478250e-11	-3.48110518071626e-10	-1.94622474768884e-09
%             -3.99095870813370e-10	2.45675608913663e-09	1.38978573547492e-08	-1.61259325384347e-10	1.69151107974397e-09	7.46714149021448e-09	-3.67580970396180e-10	2.46313179939370e-09	1.37937709794046e-08
%             -2.06071206768069e-09	1.30713939640673e-08	7.30662715453605e-08	-8.85283356935586e-10	7.46734611732834e-09	3.70784157945058e-08	-1.91381737724494e-09	1.28619349674645e-08	7.21337732165744e-08
%             1.13752586447273e-10	-7.46232285706835e-10	-4.09863440855982e-09	5.12367099628891e-11	-3.67351264059843e-10	-1.91254116026332e-09	1.11072292219590e-10	-7.33633596620360e-10	-4.09778811133240e-09
%             -7.65600805962120e-10	5.09556619460290e-09	2.78007100836172e-08	-3.47865747381061e-10	2.46152441560095e-09	1.28530163611868e-08	-7.33621039919889e-10	5.39223275748153e-09	2.78105691964427e-08
%             -4.28947326679123e-09	2.83904679341531e-08	1.55544096103606e-07	-1.94491043244415e-09	1.37851568786728e-08	7.20858854140185e-08	-4.09781957786317e-09	2.78112666107912e-08	1.55397277597425e-07];
% %NOTE: this matrix is almost symmetric, but not quite.  the following averages the off-diagonal terms, but should 
% %      look into why covmeas is not perfectly symmetric: likely numerical instability of backslash operator
% covs.gravPb = (chol(covs.gravPb)'*chol(covs.gravPb) + chol(covs.gravPb')'*chol(covs.gravPb'))/2;
        
        
%% from weighing and purity measurements, corrected for buoyancy of air in BuoyancyConsiderations.xlsx

massGrav.RP_Pb.grams = 51.746291/1000;  %divide by 1000 to express result in grams
massGrav.RP_Pb.sigma  = 0.000419/2000;  %divide again by 2 because the original number is 2-sigma
massGrav.RP_U.grams = 255.636345/1000;  %ditto
massGrav.RP_U.sigma   = 0.000597/2000;  %ditto

massGrav.ET_Pb.grams = 319.729073/1000;  %divide by 1000 to express result in grams
massGrav.ET_Pb.sigma  = 0.006924/2000;  %ditto
massGrav.ET_U.grams = 5.3444762;  %this value already in grams
massGrav.ET_U.sigma = 0.0000658/2;  %ditto

massGrav.JMM_Pb.grams = 0.8946245;  %this value already in grams
massGrav.JMM_Pb.sigma = 0.0000354/2;
massGrav.JMM_U.grams  = 5.151572;
massGrav.JMM_U.sigma  = 0.000204/2;

%NOTE: Use purityCalculator_longcalc.m to generate M realizations of concsM_981/982/Pur
MC.purity.nbs981 = (1-concsM_981(1:nM)*10^-9); %convert to purity, concsM are expressed as ppb total impurities
MC.purity.nbs982 = (1-concsM_982(1:nM)*10^-9); %use as many 
MC.purity.purtPb = (1-concsM_Pur(1:nM)*10^-9);

purity.CRM112a.value = 1 - 223*10^-6; % http://www.nbl.doe.gov/docs/pdf/CRM_112A_%20Uranium_Metal_Sept_2010.pdf
purity.CRM112a.sigma = 0.00006/2;     % same ref, divided by 2 because it's reported with k=2 coverage factor
purity.CRM115.value =  0.999770;      %http://pbadupws.nrc.gov/docs/ML0512/ML051220501.pdf
purity.CRM115.sigma =  0.000046/2;    % same ref, divided by 2 because it's reported with k=2 coverage factor

%% Make MC trials for systematic parameters, calculate 206g/238g MC realizations

MC.ics.ET535Pb = mvnrnd(ics.ET535Pb, covs.ET535Pb, nM);
MC.ics.ET2535Pb = mvnrnd(ics.ET2535Pb, covs.ET2535Pb, nM);
MC.ics.gtU = mvnrnd(ics.gtU, covs.gtU, nM);
MC.ics.gravPb = mvnrnd(ics.gravPb, covs.gravPb, nM);

% scale standard normal distriubtions generated by randn to purity and mass determinations, assume all independent
MC.purity.CRM112a  = purity.CRM112a.value  + randn(nM,1)*purity.CRM112a.sigma;
MC.purity.CRM115   = purity.CRM115.value   + randn(nM,1)*purity.CRM115.sigma;

MC.massGrav.RP_Pb  = massGrav.RP_Pb.grams  + randn(nM,1)*massGrav.RP_Pb.sigma;
MC.massGrav.ET_Pb  = massGrav.ET_Pb.grams  + randn(nM,1)*massGrav.ET_Pb.sigma;
MC.massGrav.JMM_Pb = massGrav.JMM_Pb.grams + randn(nM,1)*massGrav.JMM_Pb.sigma;

MC.massGrav.RP_U  = massGrav.RP_U.grams  + randn(nM,1)*massGrav.RP_U.sigma;
MC.massGrav.ET_U  = massGrav.ET_U.grams  + randn(nM,1)*massGrav.ET_U.sigma;
MC.massGrav.JMM_U = massGrav.JMM_U.grams + randn(nM,1)*massGrav.JMM_U.sigma;

MC.r206238g.RP = ((MC.massGrav.RP_Pb .* MC.purity.nbs982) ./                                                        ...
 (MC.ics.gravPb(:,4)*masses.Pb204 + masses.Pb206 + MC.ics.gravPb(:,5)*masses.Pb207 + MC.ics.gravPb(:,6)*masses.Pb208))...
 ./                                                                                                              ...
              ((MC.massGrav.RP_U  .* MC.purity.CRM112a) ./                                                       ...
 (ics.r234238_112a*masses.U234 + masses.U235./MC.ics.gtU(:,2) + masses.U238));

MC.r206238g.ET = ((MC.massGrav.ET_Pb .* MC.purity.nbs981) ./                                                        ...
 (MC.ics.gravPb(:,1)*masses.Pb204 + masses.Pb206 + MC.ics.gravPb(:,2)*masses.Pb207 + MC.ics.gravPb(:,3)*masses.Pb208))...
 ./                                                                                                              ...
              ((MC.massGrav.ET_U  .* MC.purity.CRM112a) ./                                                       ...
 (ics.r234238_112a*masses.U234 + masses.U235 ./ MC.ics.gtU(:,2) + masses.U238));

MC.r206238g.JMM = ((MC.massGrav.JMM_Pb .* MC.purity.purtPb) ./                                                        ...
 (MC.ics.gravPb(:,7)*masses.Pb204 + masses.Pb206 + MC.ics.gravPb(:,8)*masses.Pb207 + MC.ics.gravPb(:,9)*masses.Pb208))...
 ./                                                                                                              ...
              ((MC.massGrav.JMM_U  .* MC.purity.CRM115) ./                                                       ...
 (masses.U235 ./ MC.ics.gtU(:,1) + masses.U238));  %note: look up 234 and 236 content of CRM 115.

%% set up mixlist

mixList.RP = {'ET2535_RP_Mix2', 'ET2535_RP_Mix3', 'ET2535_RP_Mix4', 'ET2535_RP_Mix5'...
    'ET2535_RP_Mix6', 'ET2535_RP_Mix7', 'ET2535_RP_Mix8', 'ET2535_RP_Mix9'...
    'ET535_RP_W1MIT', 'ET535_RP_W2MIT', 'ET535_RP_X1MIT', 'ET535_RP_X2MIT'...
    'ET535_RP_Mix3r', 'ET535_RP_Mix4r', 'ET535_RP_Mix5r', 'ET535_RP_Mix6r' ...
    'ET535_RP_WA', 'ET535_RP_WB', 'ET535_RP_XA', 'ET535_RP_XB',...
    'ET535_RP_W1UG', 'ET535_RP_W2UG', 'ET535_RP_X2UG'}; 

mixList.ET = {'ET2535_ET2_Mix1', 'ET2535_ET2_Mix3', ...
    'ET2535_ET2_Mix4', 'ET2535_ET2_Mix5', 'ET2535_ET2_MixY', 'ET2535_ET2_MixZ',...
    'ET535_ET2_Mix1','ET535_ET2_Mix2', 'ET535_ET2_Mix3', 'ET535_ET2_Mix4', ... 
    'ET535_ET2_MixW1', 'ET535_ET2_MixX2', ...
    'ET535_ET2_MixW1MIT', 'ET535_ET2_MixW2MIT', 'ET535_ET2_MixX1MIT', 'ET535_ET2_MixX2MIT',...
    'ET535_ET2_MixW2UG', 'ET535_ET2_MixX1UG', 'ET535_ET2_MixX2UG'}; 


mixList.JMM = {'ET2535_JMM_MIT_Mix1', 'ET2535_JMM_Mix3', 'ET2535_JMM_Mix4', 'ET2535_JMM_Mix5', ... 
    'ET535_JMM_MIT_W1', 'ET535_JMM_MIT_W2', 'ET535_JMM_MIT_X2', ...
    'ET535_JMM_Mix2', 'ET535_JMM_Mix3', 'ET535_JMM_Mix4', 'ET535_JMM_Mix5',...
    'ET535_JMM_W_A', 'ET535_JMM_W_B', 'ET535_JMM_X_A', 'ET535_JMM_X_B' ...
    'ET535_JMM_UNGE_W1', 'ET535_JMM_UNGE_W2', 'ET535_JMM_UNGE_X1'}; 


%% 
nMixes.RP = length(mixList.RP); nMixes.ET = length(mixList.ET); nMixes.JMM = length(mixList.JMM);
runcount = 0;
mix.ratios.Pb = []; mix.ratios.U = [];
for i = 1:nMixes.RP
runcount = runcount + 1;
runNamePb = [mixList.RP{i} '_Pb'];
mdataPb = evalin('base', runNamePb);
runNameU = [mixList.RP{i} '_U'];
mdataU = evalin('base', runNameU);

mix.ratios.Pb{i} = mdataPb;
mix.ratios.U{i} = mdataU;
mix.gravName{i} = 'RP';
if mixList.RP{i}(3) == '2';
    mix.is202(i) = 1;
else
    mix.is202(i) = 0;
end

end %for

for i = 1:nMixes.ET
runcount = runcount + 1;
runNamePb = [mixList.ET{i} '_Pb'];
mdataPb = evalin('base', runNamePb);
runNameU = [mixList.ET{i} '_U'];
mdataU = evalin('base', runNameU);

mix.ratios.Pb{runcount} = mdataPb;
mix.ratios.U{runcount} = mdataU;
mix.gravName{runcount} = 'ET';
if mixList.ET{i}(3) == '2';
    mix.is202(runcount) = 1;
else
    mix.is202(runcount) = 0;
end

end %for

for i = 1:nMixes.JMM
runcount = runcount + 1;
runNamePb = [mixList.JMM{i} '_Pb'];
mdataPb = evalin('base', runNamePb);
runNameU = [mixList.JMM{i} '_U'];
mdataU = evalin('base', runNameU);

mix.ratios.Pb{runcount} = mdataPb;
mix.ratios.U{runcount} = mdataU;
mix.gravName{runcount} = 'JMM';
if mixList.JMM{i}(3) == '2';
    mix.is202(runcount) = 1;
else
    mix.is202(runcount) = 0;
end

end %for

mix.mixList = [mixList.RP mixList.ET mixList.JMM];
mix.skips = [0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,0,0,0,1,1,0,1,1,1,1,1,1,0,1,1,1,1,1,0,0,0];

%%
ics.purity.nbs981 = 0.99999857885;
ics.purity.nbs982 = 0.99997666751;
ics.purity.purtnc = 0.99998904364;

