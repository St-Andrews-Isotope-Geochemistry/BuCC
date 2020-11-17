function [RESULTS, csysOUT] = fncsysKMgCaV2(flag1,TC1,S1,Z1,ph1,s1,hco31,co31,alk1,dic1,pco21,Mg1,Ca1)

% csys and equic were originally written by Zeebe to accompany Zeebe and
% Wolf-Gladrow 2001 textbook.
% 
% Kat Allen made some modifications to allow batch inputs.
% 
% James Rae did debugging and tweaked some of code during PhD (Feb/Mar
% 2010).  Made versions that are equivalent to CO2SYS and various other
% programs (see pdfs in Readme folder for details).
% 
% This version is what decided followed best practices most closely.  Is
% equivalent to SeaCarb.  Uses Lee [B]tot, KF determination of Perez and
% Fraga [1987], and K pressure corrections on SWS.  See pH_etc.pdf for
% further details.
% 
% csysK.m and equic2.m are JR?s modifications of Kat?s modifications of
% original Zeebe code.
% 
% csysKslim.m is a slimmed down version of above code (not really used much
% or tested).
% 
% fncsysK.m and equic3MyAMI.m are version of above that can be run as a
% function.  Gives output in standard double array (RESULTS) or table
% format (csysOUT).
%
% this version of csysKMgCa(V2) links with equic3MyAMIV2 which uses PyMyAMI
% to get MyAMI Ks as a function of Ca and Mg, rather than lookup table
% (given issue with pKb)

%% FOR TESTING
% boo
% F08_999 = readtable('Foster2008_CO2_cycle.xlsx','sheet','999');
% 
% flag1 = 8;
% TC1 = F08_999.T;
% S1 = F08_999.S;
% Z1 = F08_999.WaterDepth;
% ph1 = F08_999.pH;
% alk1 = F08_999.ALK;
% s1 = [];
% hco31 = [];
% co31 = [];
% dic1 = [];
% pco21 = [];

%% FLAG REFERENCES
% flag = 1      h and CO2 given
% flag = 101    h and pco2 given
% flag = 2      CO2 and HCO3 given
% flag = 3      CO2 and CO3 given  
% flag = 4      CO2 and ALK given
% flag = 5      CO2 and DIC given
% flag = 6      h and HCO3 given
% flag = 7      h and CO3 given
% flag = 8      h and ALK given
% flag = 9      h and DIC given
% flag = 10     HCO3 and CO3 given
% flag = 11     HCO3 and ALK given
% flag = 12     HCO3 and DIC given
% flag = 13     CO3 and ALK given
% flag = 14     CO3 and DIC given
% flag = 15     ALK and DIC given


%% MAIN SCRIPT

%Determine lengths of input vectors
veclengths=[length(TC1),length(S1),length(Z1),length(ph1),length(s1),length(hco31),length(co31),length(alk1),length(dic1),length(pco21)];

if length(unique(veclengths))>2
	disp(' '); disp('*** INPUT ERROR: Input vectors must all be of same length, or of length 1. ***'); disp(' '); return
end
numSamples = max(veclengths);

phflag = 0;                     % Total Scale (0) or Free Scale (1)
k1k2flag = 1;                   % Roy (0) or Mehrbach (1)
% T = TC;
% S = S;
P1 =  Z1./10 ;                        % Pressure      (bar) from water depth
flag = flag1.*(ones(numSamples));
%     flag = 11;                       % Given Variables (see chart)
%     TC = 25;                        % Temperature   (Celcius)
%     S =  35;                        % Salinity
%     ph1   = [8.2, 8.3, 8.4];        % pH
%     s1    = [10, 11, 12];           % [CO2]         (µmol/kg)
%     hco31 = [1700, 1750, 1800];     % [HCO3-]       (µmol/kg)
%     co31  = [200, 250, 300];        % [CO3--]       (µmol/kg)
%     alk1  = [2350, 2400, 2450];     % ALK           (µmol/kg)
%     dic1  = [1950, 2000, 2050];     % [DIC]         (µmol/kg)
%     pco21 = [355, 365, 375];        % pCO2          (µatm)


    
    % Gets user input.  If the user does not enter a selection, the
    % default (1) will be used, which corresponds to modern concs at 35 S.

            B_all = ones(numSamples);
%             Mg_all = ones(numSamples); % now input explicitly so that
%             could have variable Mg and Ca throughout a record. 
%             Ca_all = ones(numSamples);
            
 

% ---------------------------------------------------------------
% µmol/kg to mol/kg and µatm to atm conversion and other 
% declarations
% ---------------------------------------------------------------
if(flag == 1 | flag == 2 | flag == 3 | flag == 4 | flag == 5)
    s1 = s1 / 1000000;
end

if(flag == 2 | flag == 6 | flag == 10 | flag == 11 | flag == 12)
    hco31 = hco31 / 1000000;
end

if(flag == 3 | flag == 7 | flag == 10 | flag == 13 | flag == 14)
    co31 = co31 / 1000000;
end

if(flag == 4 | flag == 8 | flag == 11 | flag == 13 | flag == 15)	
    alk1 = alk1 / 1000000;	
end	

if(flag == 5 | flag == 9 | flag == 12 | flag == 14 | flag ==  15)	
    dic1 = dic1 / 1000000;	
end	

if(flag == 101)	
    pco21 = pco21 / 1000000;	
end	

if(flag == 1 | flag == 101 | flag == 6 | flag == 7 | flag == 8 | flag == 9)	
    h1 = 10.^(-ph1);	%  [H+]    
end	



% ---------------------------------------------------------------
% Define empty matrices (soon to be filled with exciting numbers!)
% ---------------------------------------------------------------
CO2_r  = zeros(size(numSamples));
HCO3_r = zeros(size(numSamples));
CO3_r  = zeros(size(numSamples));
ALK_r  = zeros(size(numSamples));
DIC_r  = zeros(size(numSamples));
PC02_r = zeros(size(numSamples));
Borate_r = zeros(size(numSamples));
K1_r = zeros(size(numSamples));
K2_r = zeros(size(numSamples));
Kb_r = zeros(size(numSamples));
Kw_r = zeros(size(numSamples));
Xco2_r = zeros(size(numSamples));

% H_r = zeros(size(numSamples));


% ---------------------------------------------------------------
% Calculation functions section
% ---------------------------------------------------------------

% ---------------------------------------------------------------    
% START LOOP (there may be a more efficient way to do this ... ?)
% ---------------------------------------------------------------

i = 1 ;
   
for i = 1:numSamples

% Case 1  (h and s given)
if(flag == 1 & i <= numSamples)
        
equic3MyAMIV2;    
    
    dic = s*(1.+K1/h+K1*K2/h/h);
    hco3 = dic/(1+h/K1+K2/h);
    co3 = dic/(1+h/K2+h*h/K1/K2);
    alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;                        
    pco2 = s/Kh;
    pco2b = pco2/(exp(1.01325e5*(BB+2*delta)/(R*T)));
    
    % ----------- change units: mumol/kg
    CO2_r(i)    = 1.e6*s;
    HCO3_r(i)   = hco3*1.e6;
    CO3_r(i)    = co3*1.e6;
    DIC_r(i)    = dic*1.e6;
    ALK_r(i)    = alk*1.e6;
    PC02_r(i)   = pco2*1.e6;
    H_r(i)      = h ;
    Borate_r(i) = bor/(1 + h/Kb)*1.e6;
    K1_r(i) = K1;
    K2_r(i) = K2;
    Kb_r(i) = Kb;
    Kw_r(i) = Kw;
    
    Kspc_r(i) = Kspc;
        Omc_r(i) = Ca1(i)/1000*co3/Kspc;

    Xco2_r(i) = (pco2b/(1-pH2O))*1e6;
    
    %alktest = (2*co31 + hco31 + bor/(1+h/Kb) + Kw/h - h)*1.e6
    
    i = i + 1;
    




% Case 101  (h and pco2 given)
elseif(flag == 101 & i <= numSamples)
    
    % note to self -- check the ppmv conversion, seems strange
        
equic3MyAMIV2;    


    h = h1(i);
    pco2 = pco21(i);
    pco2b = pco2/(exp(1.01325e5*(BB+2*delta)/(R*T)));
    
    s = Kh*pco2;
    dic = s*(1.+K1/h+K1*K2/h/h);
    hco3 = dic/(1+h/K1+K2/h);
    co3 = dic/(1+h/K2+h*h/K1/K2);
    alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
    
    % ----------- change units: mumol/kg
    CO2_r(i)    = 1.e6*s;
    HCO3_r(i)   = hco3*1.e6;
    CO3_r(i)    = co3*1.e6;
    DIC_r(i)    = dic*1.e6;
    ALK_r(i)    = alk*1.e6;
    PC02_r(i)   = pco2*1.e6;
    H_r(i)      = h;
    Borate_r(i) = bor/(1 + h/Kb)*1.e6;
    K1_r(i) = K1;
    K2_r(i) = K2;
    Kb_r(i) = Kb;
    Kw_r(i) = Kw;
    
    Kspc_r(i) = Kspc;
    Omc_r(i) = Ca1(i)/1000*co3/Kspc;
    
    Xco2_r(i) = (pco2b/(1-pH2O))*1e6;
    
    %alktest = (2*co31 + hco31 + bor/(1+h/Kb) + Kw/h - h)*1.e6
    
    i = i + 1;
 


% Case 2  (s and HCO3 given)
elseif(flag == 2 & i <= numSamples)
        
equic3MyAMIV2;    


    s = s1(i);
    hco3 = hco31(i);
    
    p3 = -hco3/K1;
    p2 = s - hco3;
    p1 = s*K1 - hco3*K2;
    p0 = s*K1*K2;
    p = [p3 p2 p1 p0];
    r = roots(p);
    h = max(real(r));
    h*1.e12;
    dic = s*(1.+K1/h+K1*K2/h/h);
    %       hco3 = dic/(1+h/K1+K2/h); 
    co3 = dic/(1+h/K2+h*h/K1/K2);
    alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
    
    pco2 = s/Kh;
    pco2b = pco2/(exp(1.01325e5*(BB+2*delta)/(R*T)));
    
    % ----------- change units: mumol/kg
    %h1 = 10^(-ph1);
    CO2_r(i)    = 1.e6*s;
    HCO3_r(i)   = hco3*1.e6;
    CO3_r(i)    = co3*1.e6;
    DIC_r(i)    = dic*1.e6;
    ALK_r(i)    = alk*1.e6;
    PC02_r(i)   = s*1.e6/Kh;
    H_r(i)      = h;
    Borate_r(i) = bor/(1 + h/Kb)*1.e6;
    K1_r(i) = K1;
    K2_r(i) = K2;
    Kb_r(i) = Kb;
    Kw_r(i) = Kw;
    
    Kspc_r(i) = Kspc;
    Omc_r(i) = Ca1(i)/1000*co3/Kspc;
    
    Xco2_r(i) = (pco2b/(1-pH2O))*1e6;
    
    i = i + 1;




% Case 3  (s and CO3 given)
elseif (flag == 3 & i <= numSamples)
        
equic3MyAMIV2;    


    s = s1(i);
    co3 = co31(i);
    
    p4 = -co3/K1/K2;
    p3 = -co3/K2;
    p2 = s-co3;
    p1 = s*K1;
    p0 = s*K1*K2;
    p = [p4 p3 p2 p1 p0];      
    r = roots(p);
    h = max(real(r));
    h*1.e12;
    dic = s*(1.+K1/h+K1*K2/h/h);
    hco3 = dic/(1+h/K1+K2/h);
    %    co3 = dic/(1+h/K2+h*h/K1/K2);
    alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
    
    pco2 = s/Kh;
    pco2b = pco2/(exp(1.01325e5*(BB+2*delta)/(R*T)));
    
    % ----------- change units: mumol/kg
    %h1 = 10^(-ph1);
    CO2_r(i)    = 1.e6*s;
    HCO3_r(i)   = hco3*1.e6;
    CO3_r(i)    = co3*1.e6;
    DIC_r(i)    = dic*1.e6;
    ALK_r(i)    = alk*1.e6;
    PC02_r(i)   = s*1.e6/Kh;
    H_r(i)      = h;
    Borate_r(i) = bor/(1 + h/Kb)*1.e6;
    K1_r(i) = K1;
    K2_r(i) = K2;
    Kb_r(i) = Kb;
    Kw_r(i) = Kw;
    
    Kspc_r(i) = Kspc;
    Omc_r(i) = Ca1(i)/1000*co3/Kspc;
    
    Xco2_r(i) = (pco2b/(1-pH2O))*1e6;
    
    i = i + 1;




% Case 4  (s and ALK given)
elseif (flag == 4 & i <= numSamples)
        
equic3MyAMIV2;    


    s = s1(i);
    alk = alk1(i);
    
    p4 = 1.;              
    p3 = Kb+alk;
    p2 = alk*Kb-s*K1-Kb*bor-Kw;
    p1 = -s*Kb*K1-s*2.*K1*K2-Kw*Kb;
    p0 = -2.*s*Kb*K1*K2;
    p = [p4 p3 p2 p1 p0];
    r = roots(p);
    h = max(real(r));
    h*1.e12;
    dic = s*(1.+K1/h+K1*K2/h/h);
    hco3 = dic/(1+h/K1+K2/h);
    co3 = dic/(1+h/K2+h*h/K1/K2);
    %    alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
    
    pco2 = s/Kh;
    pco2b = pco2/(exp(1.01325e5*(BB+2*delta)/(R*T)));
    
    % ----------- change units: mumol/kg
    %h1 = 10^(-ph1);
    CO2_r(i)    = 1.e6*s;
    HCO3_r(i)   = hco3*1.e6;
    CO3_r(i)    = co3*1.e6;
    DIC_r(i)    = dic*1.e6;
    ALK_r(i)    = alk*1.e6;
    PC02_r(i)   = s*1.e6/Kh;
    H_r(i)      = h;
    Borate_r(i) = bor/(1 + h/Kb)*1.e6;
    K1_r(i) = K1;
    K2_r(i) = K2;
    Kb_r(i) = Kb;
    Kw_r(i) = Kw;

    Kspc_r(i) = Kspc;
    Omc_r(i) = Ca1(i)/1000*co3/Kspc;
        
    Xco2_r(i) = (pco2b/(1-pH2O))*1e6;
    
    i = i + 1;




% Case 5  (s and DIC given)
elseif (flag == 5 & i <= numSamples)
        
equic3MyAMIV2;    


    s = s1(i);
    dic = dic1(i);
    
    p2 = dic - s;
    p1 = -s*K1;
    p0 = -s*K1*K2;
    p = [p2 p1 p0];
    r = roots(p);
    h = max(real(r));
    h*1.e12;
    %    dic = s*(1.+K1/h+K1*K2/h/h);
    hco3 = dic/(1+h/K1+K2/h);
    co3 = dic/(1+h/K2+h*h/K1/K2);
    alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
    
    pco2 = s/Kh;
    pco2b = pco2/(exp(1.01325e5*(BB+2*delta)/(R*T)));
    
    % ----------- change units: mumol/kg
    %h1 = 10^(-ph1);
    CO2_r(i)    = 1.e6*s;
    HCO3_r(i)   = hco3*1.e6;
    CO3_r(i)    = co3*1.e6;
    DIC_r(i)    = dic*1.e6;
    ALK_r(i)    = alk*1.e6;
    PC02_r(i)   = s*1.e6/Kh;
    H_r(i)      = h;
    Borate_r(i) = bor/(1 + h/Kb)*1.e6;
    K1_r(i) = K1;
    K2_r(i) = K2;
    Kb_r(i) = Kb;
    Kw_r(i) = Kw;
    
    Kspc_r(i) = Kspc;
    Omc_r(i) = Ca1(i)/1000*co3/Kspc;
    
    Xco2_r(i) = (pco2b/(1-pH2O))*1e6;
    
    i = i + 1;




% Case 6  (h and HCO3 given)
elseif (flag == 6 & i <= numSamples)
        
equic3MyAMIV2;    


    h = h1(i);
    hco3 = hco31(i);
    
    dic = hco3 * (1+h/K1+K2/h);
    s = dic / (1.+K1/h+K1*K2/h/h);
    h*1.e12;
    %    dic = s*(1.+K1/h+K1*K2/h/h);
    %    hco3 = dic/(1+h/K1+K2/h);
    co3 = dic/(1+h/K2+h*h/K1/K2);
    alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
    
    pco2 = s/Kh;
    pco2b = pco2/(exp(1.01325e5*(BB+2*delta)/(R*T)));
    
    % ----------- change units: mumol/kg
    %h1 = 10^(-ph1);
    CO2_r(i)    = 1.e6*s;
    HCO3_r(i)   = hco3*1.e6;
    CO3_r(i)    = co3*1.e6;
    DIC_r(i)    = dic*1.e6;
    ALK_r(i)    = alk*1.e6;
    PC02_r(i)   = s*1.e6/Kh;
    H_r(i)      = h;
    Borate_r(i) = bor/(1 + h/Kb)*1.e6;
    K1_r(i) = K1;
    K2_r(i) = K2;
    Kb_r(i) = Kb;
    Kw_r(i) = Kw;
    
    Kspc_r(i) = Kspc;
    Omc_r(i) = Ca1(i)/1000*co3/Kspc;
    
    Xco2_r(i) = (pco2b/(1-pH2O))*1e6;
    
    i = i + 1;




% Case 7  (h and CO3 given)
elseif (flag == 7 & i <= numSamples)
        
equic3MyAMIV2;    


    h = h1(i);
    co3 = co31(i);
    
    dic = co3 * (1+h/K2+h*h/K1/K2);
    s = dic / (1.+K1/h+K1*K2/h/h);
    h*1.e12;
    %    dic = s*(1.+K1/h+K1*K2/h/h);
    hco3 = dic/(1+h/K1+K2/h);
    %    co3 = dic/(1+h/K2+h*h/K1/K2);
    alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
    
    pco2 = s/Kh;
    pco2b = pco2/(exp(1.01325e5*(BB+2*delta)/(R*T)));
    
    % ----------- change units: mumol/kg
    %h1 = 10^(-ph1);
    CO2_r(i)    = 1.e6*s;
    HCO3_r(i)   = hco3*1.e6;
    CO3_r(i)    = co3*1.e6;
    DIC_r(i)    = dic*1.e6;
    ALK_r(i)    = alk*1.e6;
    PC02_r(i)   = s*1.e6/Kh;
    H_r(i)      = h;
    Borate_r(i) = bor/(1 + h/Kb)*1.e6;
    K1_r(i) = K1;
    K2_r(i) = K2;
    Kb_r(i) = Kb;
    Kw_r(i) = Kw;

    Kspc_r(i) = Kspc;
    Omc_r(i) = Ca1(i)/1000*co3/Kspc;
    
    Xco2_r(i) = (pco2b/(1-pH2O))*1e6;
    
    i = i + 1;
  



% Case 8  (h and ALK given)
elseif (flag == 8 & i <= numSamples)
        
equic3MyAMIV2;    


    h = h1(i);
    alk = alk1(i);
    
    s = (alk-Kw/h+h-Kb*bor/(Kb+h)) / (K1/h+2.*K1*K2/h/h);
    h*1.e12;
    dic = s*(1.+K1/h+K1*K2/h/h);
    hco3 = dic/(1+h/K1+K2/h);
    co3 = dic/(1+h/K2+h*h/K1/K2);
    %    alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
    
    pco2 = s/Kh;
    pco2b = pco2/(exp(1.01325e5*(BB+2*delta)/(R*T)));
    
    % ----------- change units: mumol/kg
    %h1 = 10^(-ph1);
    CO2_r(i)    = 1.e6*s;
    HCO3_r(i)   = hco3*1.e6;
    CO3_r(i)    = co3*1.e6;
    DIC_r(i)    = dic*1.e6;
    ALK_r(i)    = alk*1.e6;
    PC02_r(i)   = s*1.e6/Kh;
    H_r(i)      = h;
    Borate_r(i) = bor/(1 + h/Kb)*1.e6;
    K1_r(i) = K1;
    K2_r(i) = K2;
    Kb_r(i) = Kb;
    Kw_r(i) = Kw;
    
    Kspc_r(i) = Kspc;
    
    Omc_r(i) = Ca1(i)/1000*co3/Kspc;
    
    Xco2_r(i) = (pco2b/(1-pH2O))*1e6;
    
    i = i + 1;




% Case 9  (h and DIC given)

elseif (flag == 9 & i <= numSamples)   
        
equic3MyAMIV2;    


    h = h1(i);
    dic = dic1(i);
    
    s = dic / (1.+K1/h+K1*K2/h/h);
    h*1.e12;
    %    dic = s*(1.+K1/h+K1*K2/h/h);
    hco3 = dic/(1+h/K1+K2/h);
    co3 = dic/(1+h/K2+h*h/K1/K2);
    alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
    
    pco2 = s/Kh;
    pco2b = pco2/(exp(1.01325e5*(BB+2*delta)/(R*T)));
    
    % ----------- change units: mumol/kg
    %h1 = 10^(-ph1);
    CO2_r(i)    = 1.e6*s;
    HCO3_r(i)   = hco3*1.e6;
    CO3_r(i)    = co3*1.e6;
    DIC_r(i)    = dic*1.e6;
    ALK_r(i)    = alk*1.e6;
    PC02_r(i)   = s*1.e6/Kh;
    H_r(i)      = h;
    Borate_r(i) = bor/(1 + h/Kb)*1.e6;
    K1_r(i) = K1;
    K2_r(i) = K2;
    Kb_r(i) = Kb;
    Kw_r(i) = Kw;

    Kspc_r(i) = Kspc;
    Omc_r(i) = Ca1(i)/1000*co3/Kspc;
    
    Xco2_r(i) = (pco2b/(1-pH2O))*1e6;
    
    i = i + 1;




% Case 10  (HCO3 and CO3 given)


elseif (flag == 10 & i <= numSamples)
    
        
equic3MyAMIV2;    


    hco3 = hco31(i);  
    co3 = co31(i);
    
    p3 = -co3/K1/K2;
    p2 = -co3/K2 + hco3/K1;
    p1 = -co3 + hco3;
    p0 = hco3*K2;
    p = [p3 p2 p1 p0];
    r = roots(p);
    h = max(real(r));
    h*1.e12;
    dic = hco3 * (1+h/K1+K2/h);
    s = dic / (1.+K1/h+K1*K2/h/h);
    %    hco3 = dic/(1+h/K1+K2/h);
    %    co3 = dic/(1+h/K2+h*h/K1/K2);
    alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
    
    pco2 = s/Kh;
    pco2b = pco2/(exp(1.01325e5*(BB+2*delta)/(R*T)));
    
    % ----------- change units: mumol/kg
    %h1 = 10^(-ph1);
    CO2_r(i)    = 1.e6*s;
    HCO3_r(i)   = hco3*1.e6;
    CO3_r(i)    = co3*1.e6;
    DIC_r(i)    = dic*1.e6;
    ALK_r(i)    = alk*1.e6;
    PC02_r(i)   = s*1.e6/Kh;
    H_r(i)      = h;
    Borate_r(i) = bor/(1 + h/Kb)*1.e6;
    K1_r(i) = K1;
    K2_r(i) = K2;
    Kb_r(i) = Kb;
    Kw_r(i) = Kw;
    
    Kspc_r(i) = Kspc;
    Omc_r(i) = Ca1(i)/1000*co3/Kspc;
        
    Xco2_r(i) = (pco2b/(1-pH2O))*1e6;
    
    i = i + 1;




% Case 11  (HCO3 and ALK given)
%
%       don't split the lines,
%       matlab does not understand it.
%
elseif (flag == 11 & i <= numSamples)
            
equic3MyAMIV2;    


    hco3 = hco31(i);
    alk = alk1(i);
    
    p5  = 1.;
    p4  = alk - hco3 + K1 + Kb;
    p3  = alk*(Kb+K1)-hco3*(K1+Kb+2.*K2)-Kw+K1*Kb+K1*K2-Kb*bor;
    tmp = alk*(Kb*K1+K1*K2)-hco3*((Kb+2.*K2)*K1+2.*Kb*K2+K1*K2);
    p2  = tmp +(-K1*Kb*bor-Kw*Kb-K1*Kw+K1*K2*Kb);
    tmp = alk*Kb*K1*K2-hco3*(2.*Kb*K1*K2+K2*K1*(Kb+2.*K2));
    p1  = tmp +(-K1*K2*Kb*bor-K1*Kw*Kb-K1*K2*Kw);
    p0  = -hco3*2.*K2*Kb*K1*K2-K1*K2*Kw*Kb;
    p   = [p5 p4 p3 p2 p1 p0];
    r   = roots(p);
    h   = max(real(r));
    h*1.e12;               
    dic = hco3 * (1+h/K1+K2/h);
    s = dic / (1.+K1/h+K1*K2/h/h);
    %   dic = s*(1.+K1/h+K1*K2/h/h);
    %   hco3 = dic/(1+h/K1+K2/h);
    co3 = dic/(1+h/K2+h*h/K1/K2);
    alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
    
    pco2 = s/Kh;
    pco2b = pco2/(exp(1.01325e5*(BB+2*delta)/(R*T)));
    
    % ----------- change units: mumol/kg
    %h1 = 10^(-ph1);
    CO2_r(i)    = 1.e6*s;
    HCO3_r(i)   = hco3*1.e6;
    CO3_r(i)    = co3*1.e6;
    DIC_r(i)    = dic*1.e6;
    ALK_r(i)    = alk*1.e6;
    PC02_r(i)   = s*1.e6/Kh;
    H_r(i)      = h;
    Borate_r(i) = bor/(1 + h/Kb)*1.e6;
    K1_r(i) = K1;
    K2_r(i) = K2;
    Kb_r(i) = Kb;
    Kw_r(i) = Kw;
    
    Kspc_r(i) = Kspc;
    Omc_r(i) = Ca1(i)/1000*co3/Kspc;
        
    Xco2_r(i) = (pco2b/(1-pH2O))*1e6;
    
    i = i + 1;



% Case 12  (HCO3 and DIC given)
elseif (flag == 12 & i <= numSamples)
        
equic3MyAMIV2;    

    
    hco3 = hco31(i);
    dic = dic1(i);
    
    p2 = hco3/K1;
    p1 = hco3-dic;    
    p0 = hco3*K2;
    p = [p2 p1 p0];
    r = roots(p);
    h = min(real(r));        % min instead of max !!!!!
    h*1.e12;
    s = dic / (1.+K1/h+K1*K2/h/h);
    dic = s*(1.+K1/h+K1*K2/h/h);
    %   hco3 = dic/(1+h/K1+K2/h);
    co3 = dic/(1+h/K2+h*h/K1/K2);
    alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
    
    pco2 = s/Kh;
    pco2b = pco2/(exp(1.01325e5*(BB+2*delta)/(R*T)));
    
    % ----------- change units: mumol/kg
    %h1 = 10^(-ph1);
    CO2_r(i)    = 1.e6*s;
    HCO3_r(i)   = hco3*1.e6;
    CO3_r(i)    = co3*1.e6;
    DIC_r(i)    = dic*1.e6;
    ALK_r(i)    = alk*1.e6;
    PC02_r(i)   = s*1.e6/Kh;
    H_r(i)      = h;
    Borate_r(i) = bor/(1 + h/Kb)*1.e6;
    K1_r(i) = K1;
    K2_r(i) = K2;
    Kb_r(i) = Kb;
    Kw_r(i) = Kw;
    
    Kspc_r(i) = Kspc;
    Omc_r(i) = Ca1(i)/1000*co3/Kspc;
        
    Xco2_r(i) = (pco2b/(1-pH2O))*1e6;

    i = i + 1;




% Case 13  (CO3 and ALK given)
elseif (flag == 13 & i <= numSamples)
        
equic3MyAMIV2;    


    co3 = co31(i);
    alk = alk1(i);
    
    p5  = -co3/K2+1.;
    p4  = alk - co3*(K1/K2+(Kb+2.*K2)/K2) + Kb + K1;
    tmp = alk*(Kb+K1)-co3*(K1+K1*(Kb+2.*K2)/K2+2.*Kb);
    p3  = tmp+(-Kb*bor-Kw+K1*Kb+K1*K2);
    tmp = alk*(Kb*K1+K1*K2)-co3*(K1*(Kb+2.*K2)+2.*Kb*K1);
    p2  = tmp+(-Kw*Kb-K1*Kb*bor-K1*Kw+K1*K2*Kb);
    tmp = alk*Kb*K1*K2-co3*2.*Kb*K1*K2-K1*Kw*Kb;
    p1  = tmp+(-K1*K2*Kb*bor-K1*K2*Kw);
    p0  = -K1*K2*Kw*Kb;
    p   = [p5 p4 p3 p2 p1 p0];
    r   = roots(p);
    h   = max(real(r));
    h*1.e12;
    %
    % NOTE : calculate dic from dic = co3*(1+h/K2+h^2/K1/K2);
    %                  not from dic = hco3*(1+h/K1+K2/h);
    %
    dic = co3 * (1+h/K2+h^2/K1/K2);
    s = dic / (1.+K1/h+K1*K2/h/h);
    %   dic = s*(1.+K1/h+K1*K2/h/h);    
    hco3 = dic/(1+h/K1+K2/h);
    %   co3 = dic/(1+h/K2+h*h/K1/K2);
    %   alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
    
    pco2 = s/Kh;
    pco2b = pco2/(exp(1.01325e5*(BB+2*delta)/(R*T)));
    
    % ----------- change units: mumol/kg
    %h1 = 10^(-ph1);
    CO2_r(i)    = 1.e6*s;
    HCO3_r(i)   = hco3*1.e6;
    CO3_r(i)    = co3*1.e6;
    DIC_r(i)    = dic*1.e6;
    ALK_r(i)    = alk*1.e6;
    PC02_r(i)   = s*1.e6/Kh;
    H_r(i)      = h;
    Borate_r(i) = bor/(1 + h/Kb)*1.e6;
    K1_r(i) = K1;
    K2_r(i) = K2;
    Kb_r(i) = Kb;
    Kw_r(i) = Kw;

    Kspc_r(i) = Kspc;
    Omc_r(i) = Ca1(i)/1000*co3/Kspc;
    
    Xco2_r(i) = (pco2b/(1-pH2O))*1e6;
    
    i = i + 1;




% Case 14  (CO3 and DIC given)

elseif(flag == 14 & i <= numSamples)
    
equic3MyAMIV2;    

    
    co3 = co31(i);
    dic = dic1(i);
    
    p2 = co3/K1/K2;
    p1 = co3/K2;
    p0 = co3-dic;
    p = [p2 p1 p0];
    r = roots(p);        
    h = max(real(r));
    h*1.e12;
    s = dic / (1.+K1/h+K1*K2/h/h);
    %   dic = s*(1.+K1/h+K1*K2/h/h);
    hco3 = dic/(1+h/K1+K2/h);
    %   co3 = dic/(1+h/K2+h*h/K1/K2);
    alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
    
    pco2 = s/Kh;
    pco2b = pco2/(exp(1.01325e5*(BB+2*delta)/(R*T)));
    
    % ----------- change units: mumol/kg
    %h1 = 10^(-ph1); %        6.3096e-09
    CO2_r(i)    = 1.e6*s;
    HCO3_r(i)   = hco3*1.e6;
    CO3_r(i)    = co3*1.e6;
    DIC_r(i)    = dic*1.e6;
    ALK_r(i)    = alk*1.e6;
    PC02_r(i)   = s*1.e6/Kh;
    H_r(i)      = h;
    Borate_r(i) = bor/(1 + h/Kb)*1.e6;
    K1_r(i) = K1;
    K2_r(i) = K2;
    Kb_r(i) = Kb;
    Kw_r(i) = Kw;
    
    Xco2_r(i) = (pco2b/(1-pH2O))*1e6;
    
    Kspc_r(i) = Kspc;
    Omc_r(i) = Ca1(i)/1000*co3/Kspc;
        
    i = i + 1;
  
 


% Case 15  (ALK and DIC given)

elseif (flag == 15 & i <= numSamples)
    
equic3MyAMIV2;    
    
    
    alk = alk1(i);
    dic = dic1(i);
    
    p5  = -1.;        
    p4  = -alk-Kb-K1;
    p3  = dic*K1-alk*(Kb+K1)+Kb*bor+Kw-Kb*K1-K1*K2;
    tmp = dic*(Kb*K1+2.*K1*K2)-alk*(Kb*K1+K1*K2)+Kb*bor*K1;
    p2  = tmp+(Kw*Kb+Kw*K1-Kb*K1*K2);
    tmp = 2.*dic*Kb*K1*K2-alk*Kb*K1*K2+Kb*bor*K1*K2;
    p1  = tmp+(+Kw*Kb*K1+Kw*K1*K2);
    p0  = Kw*Kb*K1*K2;
    p   = [p5 p4 p3 p2 p1 p0];
    r   = roots(p);
    h   = max(real(r));
    
    %   test = p5*h^5+p4*h^4+p3*h^3+p2*h^2+p1*h+p0
    %h*1.e12;
    s = dic / (1.+K1/h+K1*K2/h/h);
    %   dic = s*(1.+K1/h+K1*K2/h/h);
    hco3 = dic/(1+h/K1+K2/h);
    co3 = dic/(1+h/K2+h*h/K1/K2);
    %   alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
    
    pco2 = s/Kh;
    pco2b = pco2/(exp(1.01325e5*(BB+2*delta)/(R*T)));
    
    % ----------- change units: mumol/kg
    %h1 = 10^(-ph1); %
    CO2_r(i)    = 1.e6*s;
    HCO3_r(i)   = hco3*1.e6;
    CO3_r(i)    = co3*1.e6;
    DIC_r(i)    = dic*1.e6;
    ALK_r(i)    = alk*1.e6;
    PC02_r(i)   = s*1.e6/Kh;
    H_r(i)      = h;
    Borate_r(i) = bor/(1 + h/Kb)*1.e6;
    K1_r(i) = K1;
    K2_r(i) = K2;
    Kb_r(i) = Kb;
    Kw_r(i) = Kw;
    
    Xco2_r(i) = (pco2b/(1-pH2O))*1e6;
    
    Kspc_r(i) = Kspc;
    Omc_r(i) = Ca1(i)/1000*co3/Kspc;
    
    i = i + 1;



end


end

pH_r = -log10(H_r);
      
      RESULTS = [CO2_r ; HCO3_r ; CO3_r ; Borate_r ; DIC_r ; ALK_r ; PC02_r ; Xco2_r ; H_r ; pH_r; K1_r ; K2_r ; Kb_r ; Kw_r; Kspc_r; Omc_r]' ;

% results in table form      
csysOUT = table(CO2_r', HCO3_r', CO3_r', Borate_r', DIC_r', ALK_r', PC02_r', Xco2_r', H_r', pH_r', K1_r', K2_r', Kb_r', Kw_r', Kspc_r', Omc_r',...
    'VariableNames',{'CO2_r' 'HCO3_r' 'CO3_r' 'B4_r' 'DIC_r' 'ALK_r' 'PCO2_r' 'XCO2_r' 'H_r' 'pH_r' 'K1_r' 'K2_r' 'Kb_r' 'Kw_r' 'Kspc_r' 'Omc_r'});