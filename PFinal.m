t=scope13(:,1);
u=scope13(:,2);
y1=scope13(:,3); %sys ord 2
y2=scope13(:,4); %sys cu 0
plot(t-t(1),u,t-t(1),y1,'b',t-t(1),y2,'r'),
legend("u","y1","y2"),xlabel("Timp[sec]"),ylabel("Amplitudine[V]"),grid
%%
subplot(311)
plot(t-t(1),u),xlabel("Timp[sec]"),ylabel("Amplitudine[V]"),title("Semnal de intrare"),grid
subplot(312)        
plot(t-t(1),y1,'b'),xlabel("Timp[sec]"),ylabel("Amplitudine[V]"),title("Semnal de iesire sistem fara zero"),grid
subplot(313)
plot(t-t(1),y2,'r'),xlabel("Timp[sec]"),ylabel("Amplitudine[V]"),title("Semnal de iesire sistem cu zero"),grid
%%
dt=t(2)-t(1);
id_y1=iddata(y1,u,dt);
id_y2=iddata(y2,u,dt);
set(id_y2,'OutputName','y2');
%% Rezonanta Sys fara zero
plot(t,u,t,y1,'r-'),legend("u","y1"),xlabel("Timp[sec]"),ylabel("Amplitudine[V]"),grid
%%
k=mean(y1(64:994))/mean(u(64:996));

A_y=(y1(266)-y1(281))/2;
A_u=(u(260)-u(276))/2;
Mr=A_y/A_u;
zeta=sqrt((Mr-sqrt(Mr^2-1))/2/Mr); 
Tr=t(297)-t(265);
wr=2*pi/Tr;
wn=wr/(sqrt(1-2*zeta^2));

H_y1=tf(k*wn^2,[1 2*zeta*wn wn^2]);
%figure;compare(id_y1,H_y1);

A1=[0 1; -wn^2 -2*wn*zeta];
B1=[0; k*wn^2];
C1=[1 0];
D1=0;
dy1=(y1(2)-y1(1))/dt; %derivata lui y1
sys=ss(A1,B1,C1,D1);
ysim=lsim(sys,u,t,[y1(1) dy1]);
plot(t,[y1 ysim]),legend("y mas","y sim"),xlabel("Timp[sec]"),ylabel("Amplitudine[V]"),grid
%% Validare numerica

%eroare medie patratica relativa
eMPN=norm(ysim-y1)/norm(ysim-mean(y1))
%% Rezonanta y2
plot(t,u,t,y2,'r-'),legend("u","y2"),xlabel("Timp[sec]"),ylabel("Amplitudine[V]"),grid
%%
k=mean(y2(69:993))/mean(u(67:987));

%A_y=(y2(266)-y2(281))/2; %vf max-vf min/2
A_y=(y2(231)-y2(247))/2;
A_u=(u(260)-u(276))/2;
Mr=A_y/A_u;

zeta=sqrt((1-sqrt(1-1/Mr^2))/2); 
Tr=t(264)-t(231);
wr=2*pi/Tr;
wn=wr/(sqrt(1-2*zeta^2));

%ID Zero dupa Faza
fir=wr*(t(264)-t(190));
T0=tan(rad2deg(fir)-90-720)/wr;

H_y2=tf(k*wn^2*[T0 1],[1 2*zeta*wn wn^2]);
zpk(H_y2)
figure;compare(H_y2,id_y2);

A2=[0 1; -wn^2 -2*wn*zeta];
B2=[k*wn^2*T0;k*wn^2-2*zeta*k*wn^3*T0];
C2=[1 0];
D2=0;

dy2=(y2(2)-y2(1))/dt;
ysim=lsim(A2,B2,C2,D2,u,t,[y2(1) dy2+u(1)*(-k*wn^2*T0)]);

figure,plot(t,[y2,ysim]),legend("y mas","y sim"),xlabel("Timp[sec]"),ylabel("Amplitudine[V]"),grid
%% Validare numerica

%eroare medie patratica relativa
eMPN=norm(ysim-y2)/norm(ysim-mean(y2))

%-----------------------------------------------------------

%% ARMAX Y1
Mamx_y1=armax(id_y1,[2 1 2 0]);
Hz_y1=tf(Mamx_y1.B,Mamx_y1.A,dt);
resid(id_y1,Mamx_y1,'corr',5); figure; compare(id_y1,Mamx_y1);
Hs_y1=d2c(Hz_y1,'matched')
Hz_y1.Variable='z^-1';
Hz_y1
%% OE Y1
Moe_y1=oe(id_y1,[2 2 0]);
Hz_y1=tf(Moe_y1.B,Moe_y1.F,dt);
resid(Moe_y1,id_y1,'corr',5); figure; compare(id_y1,Moe_y1);
Hs_y1=d2c(Hz_y1)
Hz_y1.Variable='z^-1';
Hz_y1
%% ARMAX Y2
Mamx_y2=armax(id_y2,[2 2 2 0]);
Hz_y2=tf(Mamx_y2.B,Mamx_y2.A,dt);
figure;resid(Mamx_y2,id_y2,'corr',5); figure; compare(id_y2,Mamx_y2);
Hs_y2=d2c(Hz_y2,'matched');
Hz_y2.Variable='z^-1';
Hz_y2
%% PEM Y2
Pem_y2=pem(id_y2,[2 2 1 0 2 0]);
Hz_y2_=tf([Pem_y2.B 0]/Pem_y2.F,Pem_y2.A,dt,'Variable','z^-1')
figure;resid(Pem_y2,id_y2,'corr',5); figure; compare(id_y2,Pem_y2);
Hs_y2=d2c(Hz_y2_,'tustin');
zpk(Hs_y2)

%%
Mn_y2=n4sid(id_y2);
Mn_pem_y2=pem(id_y2,Mn_y2)
[num,den]=ss2tf(Mn_pem_y2.A,Mn_pem_y2.B,Mn_pem_y2.C,Mn_pem_y2.D);
Hz_y2_=tf(num,den,dt)
figure;resid(Mn_pem_y2,id_y2,'corr',5); figure; compare(id_y2,Mn_pem_y2);
Hs_y2=d2c(Hz_y2_)
Hz_y2_.Variable="z^-1";
Hz_y2_
%-----------------------------------------------------------
%% RF Y2
plot(t,u,t,y2,'r-'),legend("u","y2"),xlabel("Timp[sec]"),ylabel("Amplitudine[V]"),grid
%%
A_y=[(y2(99)-y2(127))/2,(y2(152)-y2(173))/2,...
    (y2(193)-y2(211))/2,(y2(230)-y2(248))/2,...
    (y2(264)-y2(280))/2,(y2(295)-y2(310))/2,...
    (y2(477)-y2(489))/2,(y2(500)-y2(511))/2,...
    (y2(707)-y2(716))/2,(y2(907)-y2(916))/2];

uind=[u(99)-u(124),u(150)-u(171),u(190)-u(210),u(227)-u(244),u(258)-u(276),...
    u(290)-u(305),u(471)-u(482),u(494)-u(506),u(702)-u(712),u(902)-u(911)]./2;
A_u=uind;
uind=(u(623)-u(612))/2;
A_u=uind*ones(size(A_y));

%Pe perioada (1)
Trf=[t(152)-t(99),t(193)-t(152),t(230)-t(193),...
    t(264)-t(230),t(295)-t(264),...
    t(295)-t(264),t(500)-t(477),t(500)-t(477),...
    t(726)-t(707),t(924)-t(907)];


dt1=[t(99)-t(98),t(152)-t(150),t(193)-t(190),t(230)-t(227),...
    t(264)-t(258),t(295)-t(291),t(477)-t(471),t(500)-t(494),...
    t(707)-t(702),t(907)-t(902)];
%%
wvect=2*pi./Trf; %(1)
Mvect=A_y./A_u;
 
figure,subplot(211);
semilogx(wvect,20*log10(Mvect),'x'),hold on;
bp1=bodeplot(H_y2);
bp1.PhaseVisible='off';


Fivect=-rad2deg(wvect.*dt1);

subplot(212),bp2=bodeplot(H_y2);hold on
semilogx(wvect,Fivect,'x'),hold on;
bp2.MagVisible='off';


