% This model tries to capture sex differences in the circadian rhythms of
% the HPA axis
% The model is based on Mavroudis et al. and possibly Rao et al. 
% Cytokine feedback is included in ACTH as of 04/04/16
% Further attempts to add cytokine feedback to other components of the HPA
% axis
function dy = RR_021017l_SexDifferences_CST_Nocturnal_Arnold_v2(t,ic_g,h,cp,k,k2,ksq,stress,lf,eps,tau)

% Can use parameters from Model_CALL_RA_RR_12_07_15_nopath.m
% Else, the parameters are from 'ParameterFile_JuskoOsc.xlsx and 
% ParameterFile_RA_Goodwin_10_05_15.xlsx

%% Initialize variables
CRH=ic_g(1);
ACTH=ic_g(2);
F_h=ic_g(3);
GRm_h=ic_g(4);
GR_h=ic_g(5);
FR_h=ic_g(6);
FRn_h=ic_g(7);
F_p=ic_g(8);
GRm_p=ic_g(9);
GR_p=ic_g(10);
FR_p=ic_g(11);
FRn_p=ic_g(12);
IL1Bm=ic_g(13);
IL1B_p=ic_g(14);
IL1B_h=ic_g(15);
t1=ic_g(16);
t2=ic_g(17);
t3=ic_g(18);
t4=ic_g(19);
us=ic_g(20);
t5=ic_g(21);
t6=ic_g(22);
t7=ic_g(23);
t8=ic_g(24);
t9=ic_g(25);
% CRHm=ic_g(16);
% ACTHm=ic_g(17);
% F_hm=ic_g(18);
% GRm_hm=ic_g(19);
% GR_hm=ic_g(20);
% FR_hm=ic_g(21);
% FRn_hm=ic_g(22);
% F_pm=ic_g(23);
% GRm_pm=ic_g(24);
% GR_pm=ic_g(25);
% FR_pm=ic_g(26);
% FRn_pm=ic_g(27);
% IL1Bm_m=ic_g(28);
% IL1B_pm=ic_g(29);
% IL1B_hm=ic_g(30);

%% Model Equations
% if (mod(t+kt,24.0000) >= 23)
%     f=0.4*2*n;
% else
%     f=0;
% end 

% HPA Axis
% if (mod(t,tau)>tau/2)
%     light2=eps;
% else
%     light2=0;
% end

light2=eps+eps*sin(2*pi*t/tau);

if (mod(t+3,24)>=10)
    light=1*lf;
else
    light=0;
end
%light=heaviside(mod(t,24)-7)-heaviside(mod(t,24)-21);

dt1=ksq(1)*(light-t1);
dt2=ksq(2)*(t1-t2);
dt3=ksq(2)*(t2)-ksq(2)*t3;
dt4=ksq(2)*(t3-t4);
dt5=ksq(2)*(1-light)-ksq(2)*t5;
dt6=ksq(2)*(t5-t6);
dt7=ksq(2)*1*(t6-t7);
dt8=ksq(2)*1*(t7-t8);
dt9=ksq(2)*(t8-t9);
dus=ksq(3)*t4^2/(0.01+t4^2)-ksq(2)*us*(1+25*t8^2);

% Female HPA Axis
% Corticotropin Releasing Hormone
% multiplied by 1.1

dCRH=stress*k(1)*1*k2(1)/(k2(1)+FRn_h)-k(3)*CRH*(1+us/(1+us))/(k(4)+CRH)+1*light2;% k(3)*CRH*(1+light2)/(k(4)+CRH) %% +light2; % intial 1.75 period 24 hours
%dCRH=k(1)-k(3)*CRH*(1+1*us/(1.0+us))/(k(4)+CRH)+light2; % intial 1.75 period 24 hours

% ACTH
dACTH=k(5)*k2(2)*CRH/(k2(2)+FRn_h^1)-1*k(6)*ACTH/(1*k(7)+ACTH); %10_05_11 0.3*IL1B_h/10... *(1+0.15*M^2/(10+M^2)) (1+0.2*M^1/(55^4+M^4));%0.30 works more or less for IL1B and  *****(1+3*M^2/(10+M^2))+0.35*IL1Bm/(10+FRn_h^2)***(1+3*M^2/(10+M^2))+0.35*IL1Bm/(10+FRn_h^2)**** %-0.5*IL1Bm*M^1/(10+IL1Bm^3); % 3.4*M/(10+FRn_h^2) when M is from F_h ACTH differential equation

% Corticosterone
dF_h=1*k2(3)*ACTH-k(9)*F_h/(1*k(10)+F_h);

% GR mNRA in HPA
dGRm_h=h(11)*(1.0)*(1-FRn_h/(h(12)+FRn_h))-h(13)*GRm_h; %+h(30)*1*IL1B_h^h(31)/(h(32)^h(33)+IL1B_h^h(34));   

% GR Protein in HPA
dGR_h=h(14)*GRm_h+h(15)*h(16)*FRn_h-h(17)*(F_h)*GR_h-h(18)*GR_h; % GR in HPA

% Cytoplasmic CST-Bound Receptor in HPA
dFR_h=h(17)*(F_h)*GR_h-h(19)*FR_h; % Corticosterone bound receptor in cytoplasm

% Nuclear CST-Bound Receptor
dFRn_h=h(19)*(FR_h)-h(16)*FRn_h;

% Periphery

% Corticosterone in periphery
dF_p=1/h(20)*(F_h-F_p);

%GR mRNA in periphery
dGRm_p=h(11)*(1)*(1-FRn_p/(h(12)+FRn_p))-h(13)*GRm_p; %+h(35)*1.15*IL1B_p^h(36)/(h(37)^h(38)+IL1B_p^h(39));

% GR protein in periphery
dGR_p=h(14)*GRm_p+h(15)*h(16)*FRn_p-h(17)*(F_p)*GR_p-h(18)*GR_p; % GR in periphery

% Cytoplasmic CST-bound receptor
dFR_p=h(17)*(F_p)*GR_p-h(19)*FR_p; 

% Nuclear CST-bound receptor
dFRn_p=h(19)*FR_p-h(16)*FRn_p; 

% Cytokine mRNA
%dIL1Bm=0.1*cp(9)*(1-FRn_p/(0.06*cp(11)+FRn_p))-1*h(41)*IL1Bm; %(22) Turnover of IL1-beta mRNA *(1+TC_IL1B)
dIL1Bm=6*(1-FRn_p/(1+FRn_p))-1*h(41)*IL1Bm; %(22) Turnover of IL1-beta mRNA *(1+TC_IL1B)

% Cytokine protein in periphery
dIL1B_p=h(42)*IL1Bm-1.0*IL1B_p; %Cytokine protein in the periphery

% Cytokines protein in HPA 
dIL1B_h=1/(h(43)*1*h(20))*(IL1B_p-IL1B_h); 

% Male HPA axis rhythms

% Male Corticotrophin Resleasing Hormone
% dCRHm=1*k(1)*kn(4)/(kn(4)+FRn_hm)-k(3)*CRHm*(1+1.0*light/(1+light))/(k(4)+CRHm); % intial 1.75 period 24 hours
% 
% % Male ACTH
% dACTHm=k(5)*kn(5)*CRHm/(kn(5)+FRn_hm^1)-1*k(6)*ACTHm/(1*k(7)+ACTHm)+0*k(12)*IL1B_hm/(k(13)+FRn_hm^1); %10_05_11 0.3*IL1B_h/10... *(1+0.15*M^2/(10+M^2)) (1+0.2*M^1/(55^4+M^4));%0.30 works more or less for IL1B and  *****(1+3*M^2/(10+M^2))+0.35*IL1Bm/(10+FRn_h^2)***(1+3*M^2/(10+M^2))+0.35*IL1Bm/(10+FRn_h^2)**** %-0.5*IL1Bm*M^1/(10+IL1Bm^3); % 3.4*M/(10+FRn_h^2) when M is from F_h ACTH differential equation
% 
% % Male Corticosterone
% dF_hm=1*kn(6)*ACTHm-1*k(9)*F_hm/(1*k(10)+F_hm);
% 
% % Male GR mNRA in HPA
% dGRm_hm=h(11)*(1.0)*(1-FRn_hm/(h(12)+FRn_hm))-h(13)*GRm_hm; %+h(30)*1*IL1B_h^h(31)/(h(32)^h(33)+IL1B_h^h(34));   
% 
% % Male GR Protein in HPA
% dGR_hm=h(14)*GRm_hm+h(15)*h(16)*FRn_hm-h(17)*(F_hm)*GR_h-h(18)*GR_hm; % GR in HPA
% 
% % Male Cytoplasmic CST-Bound Receptor in HPA
% dFR_hm=h(17)*(F_hm)*GR_hm-h(19)*FR_hm; % Corticosterone bound receptor in cytoplasm
% 
% % Male Nuclear CST-Bound Receptor
% dFRn_hm=h(19)*FR_hm-h(16)*FRn_hm;
% 
% % Male Periphery
% 
% % Male Corticosterone in periphery
% dF_pm=1/h(20)*(F_hm-F_pm);
% 
% % Male GR mRNA in periphery
% dGRm_pm=h(11)*(1)*(1-FRn_pm/(h(12)+FRn_pm))-h(13)*GRm_pm; %+h(35)*1.15*IL1B_p^h(36)/(h(37)^h(38)+IL1B_p^h(39));
% 
% % Male GR protein in periphery
% dGR_pm=h(14)*GRm_pm+h(15)*h(16)*FRn_pm-h(17)*(F_pm)*GR_pm-h(18)*GR_pm; % GR in periphery
% 
% % Male Cytoplasmic CST-bound receptor
% dFR_pm=h(17)*(F_pm)*GR_pm-h(19)*FR_pm; 
% 
% % Male Nuclear CST-bound receptor
% dFRn_pm=h(19)*FR_pm-h(16)*FRn_pm; 
% 
% % Male Cytokine mRNA
% %dIL1Bm=0.1*cp(9)*(1-FRn_p/(0.06*cp(11)+FRn_p))-1*h(41)*IL1Bm; %(22) Turnover of IL1-beta mRNA *(1+TC_IL1B)
% dIL1Bm_m=6*(1-FRn_pm/(1+FRn_pm))-1*h(41)*IL1Bm_m; %(22) Turnover of IL1-beta mRNA *(1+TC_IL1B)
% 
% % Male Cytokine protein in periphery
% dIL1B_pm=h(42)*IL1Bm_m-1.0*IL1B_pm; %Cytokine protein in the periphery
% 
% % Male Cytokines protein in HPA 
% dIL1B_hm=1/(h(43)*1*h(20))*(IL1B_pm-IL1B_hm); 


%% Assign Output
dy=zeros(8,1);
dy(1)=dCRH;
dy(2)=dACTH;
dy(3)=dF_h;
dy(4)=dGRm_h;
dy(5)=dGR_h;
dy(6)=dFR_h;
dy(7)=dFRn_h;
dy(8)=dF_p;
dy(9)=dGRm_p;
dy(10)=dGR_p;
dy(11)=dFR_p;
dy(12)=dFRn_p;
dy(13)=dIL1Bm;
dy(14)=dIL1B_p;
dy(15)=dIL1B_h;
dy(16)=dt1;
dy(17)=dt2;
dy(18)=dt3;
dy(19)=dt4;
dy(20)=dus;
dy(21)=dt5;
dy(22)=dt6;
dy(23)=dt7;
dy(24)=dt8;
dy(25)=dt9;
% dy(16)=dCRHm;
% dy(17)=dACTHm;
% dy(18)=dF_hm;
% dy(19)=dGRm_hm;
% dy(20)=dGR_hm;
% dy(21)=dFR_hm;
% dy(22)=dFRn_hm;
% dy(23)=dF_pm;
% dy(24)=dGRm_pm;
% dy(25)=dGR_pm;
% dy(26)=dFR_pm;
% dy(27)=dFRn_pm;
% dy(28)=dIL1Bm_m;
% dy(29)=dIL1B_pm;
% dy(30)=dIL1B_hm;
