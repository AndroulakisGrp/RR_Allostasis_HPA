% This file describes equations for a Chronic Allostatic Stress Model of the HPA axis 
function dy = RR_2019_Allostasis_HPA(t,ic_g,h,entrain,k,k2,ksq,stress,lf,eps,tau)

%% Initialize variables
CRH=ic_g(1);
ACTH=ic_g(2);
F_h=ic_g(3);
GRm_h=ic_g(4);
GR_h=ic_g(5);
FR_h=ic_g(6);
FRn_h=ic_g(7);
t1=ic_g(8);
t2=ic_g(9);
t3=ic_g(10);
t4=ic_g(11);
us=ic_g(12);
t5=ic_g(13);
t6=ic_g(14);
t7=ic_g(15);
t8=ic_g(16);
t9=ic_g(17);

%% Model Equations

light2=eps+eps*sin(2*pi*t/tau);
if entrain
    light = light2;
else
    
    if (mod(t+3,24)>=10)
        light=1*lf;
    else
        light=0;
    end
end

% Light transduction delay
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

% HPA Axis
% Corticotropin Releasing Hormone
dCRH=stress*k(1)*1*k2(1)/(k2(1)+FRn_h)-k(3)*CRH*(1+us/(1+us))/(k(4)+CRH)+1*light2;

% ACTH
dACTH=k(5)*k2(2)*CRH/(k2(2)+FRn_h^1)-1*k(6)*ACTH/(1*k(7)+ACTH); 

% Corticosterone
dF_h=1*k2(3)*ACTH-k(9)*F_h/(1*k(10)+F_h);

% GR mNRA in HPA
dGRm_h=h(11)*(1.0)*(1-FRn_h/(h(12)+FRn_h))-h(13)*GRm_h; 

% GR Protein in HPA
dGR_h=h(14)*GRm_h+h(15)*h(16)*FRn_h-h(17)*(F_h)*GR_h-h(18)*GR_h; 

% Cytoplasmic CST-Bound Receptor in HPA
dFR_h=h(17)*(F_h)*GR_h-h(19)*FR_h; 

% Nuclear CST-Bound Receptor
dFRn_h=h(19)*(FR_h)-h(16)*FRn_h;


%% Assign Output
dy=zeros(8,1);
dy(1)=dCRH;
dy(2)=dACTH;
dy(3)=dF_h;
dy(4)=dGRm_h;
dy(5)=dGR_h;
dy(6)=dFR_h;
dy(7)=dFRn_h;
dy(8)=dt1;
dy(9)=dt2;
dy(10)=dt3;
dy(11)=dt4;
dy(12)=dus;
dy(13)=dt5;
dy(14)=dt6;
dy(15)=dt7;
dy(16)=dt8;
dy(17)=dt9;
