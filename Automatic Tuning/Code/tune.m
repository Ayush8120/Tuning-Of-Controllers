% tune3 Inicialización simple del integredor con retraso
% Incluye rechaso a perturbaciones por carga durante la prueba it includes
% Hay tres modelos de inicialización: 0 = IDDT, 1 = FODT, 2 = IDT
%----------------------------------------------------------------

clear all;
lw=1.5;   % Ancho de línea de los graficos

% Algunas definiciones por defecto y G(s)
stepsize = 0.5;   % Amplitud del escalón de entrada
SetPoint = 1;     % set-point del lazo a probar
load_s = -0.000;  % [medio-tiempo] carga en la simulación de sintonización
T_l = 0;          % Se hace desde el inicio...
tstep = 0;        % step despues de 0, si así es requerido
Td = 0;           % Tiempo muerto de la planta si se requiere
B  = [-1 1];      % La mayoría de las plantas no tienen ceros
load = 0.0;       % Carga de perturbación, al medio tiempo en la simulación final
sigma = 0.0;      % 0 = "No quiero ruido"
Snset = sigma^2;
beta = 1;         % Beta del VBPF
a = 0.01*SetPoint;% Amplitud de entrada (design for about 1:1 when noisy)
Kfiddle = 0.1;    % Ajustar ganancia adaptiva para casos más difíciles .... 
zeta = .5;        % Sobre tiro deseado en lazo cerrado
Y1 = max(5*sigma,0.05); Y2 = 0.8*stepsize; %0.4; % on/off  niveles de salida
test = '6';

%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% Escoge una planta
test = '2';% input('Test: eg *NU3 ','s');
m = length(test);
n = str2num(test(length(test)));
if m > 1,       % Encontrar el orden por los números 
    jk = test(length(test)-1);
if ~isletter(jk) & jk ~='*', n = n + 10*str2num(jk);
end; end;
plant = ['G_{',num2str(n),'}']; 
Tmax = 4*n; Tmax_c = 100+100*n; Tmax_s = 20*(n+Td);
%Tmax = 100
%Tmax_c = 2000
%Tmax_s = 100

% Contruir A(s) convulucionando plantas de primer orden
A = [1 0.5];
for ii = 2:n,
    A = conv(A,[1 1]);
end;


% Pon un integreador
intshow = 0;    % no integrador
if (~isempty(strfind(test,'I')))|(~isempty(strfind(test,'i'))),
    intshow = 1;    % hay uno, cambia los cálculos de diseño
    ades = 1/(1+2*zeta);
    stepsize = 1/Tmax;
    A = conv(A,[1 0]); plant = ['I',plant]; 
    Kfiddle = 1; Tmax_c = 10*Tmax_c; Tmax_s =5*Tmax_s;
    if n==1, Kfiddle = 0.25; end
end;

% Pon un polo subamortiguado
if (~isempty(strfind(test,'U')))|(~isempty(strfind(test,'u'))),
A = conv(A,[1 .5 1]); plant = ['U',plant]; Tmax = Tmax + 5; Tmax_s = Tmax_s + 40;
end;
% Pon un NMP cer si hay una 'n' 
if (~isempty(strfind(test,'N')))|(~isempty(strfind(test,'n'))),
B = conv(B,[-1.5 1]); plant = ['N',plant]; Tmax = Tmax + 5; Tmax_s = Tmax_s+40;
%Tmax_c = Tmax_c+100;
end;
% Ahora el ruido si hay un '*'
if strfind(test,'*'),
Snset = sigma^2; alpha = 0; Tmax = Tmax; Tmax_c = Tmax_c; plant = [plant,'^n'];
end
%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%%=================================================================
%   Aplicamos el escalón
%%=================================================================
T_s = Tmax;             % Tiempo de asentamiento
Npts = 200;             % número de muestras en la prueba
h = Tmax/Npts;           % intervalo de muestreo, por lo que Npts puntos de datos para el ajuste
Td = max(h*round(Td/h),h);
d = [0*(1:Td/h) 1];     % Para el tiempo muerto simularlo como z^{-n}
Sn = Snset*h;       % Si se quiere ruido da SD aprox para sigma
sim('steptest2');       % escalón G(s)

%%------------------------------------------------------------------
%   Se analizan los resultados obtenidos de diversas maneras
%%------------------------------------------------------------------
jy = find(y>=Y1);   % Punto inicuial: donde cruza Y1
ky = find(y>=Y2);   % Punto final: donde cruza Y2
k2 = ky(1);
ty2 = tout2(k2);
jstart = jy(1); t1 = tout2(jstart);
yp = max(y,0);      % clip respuesta de fase no mínima
yi  = h*cumsum(yp); 
f = jy*h - 2*yi(jy)./y(jy);    % modelo con integrador y tiempo muerto
pp = polyfit(tout2(jy),f,2);   % Ajusta f a segundo orden para rechazo al ruido
% El máximo de ax^2+bx+c es cuando x=-b/2a
% y tiene un valor de c-b^2/2a, porlo que usamos Tm y Td
if pp(1) == 0,
    Tdhat = pp(3); khat = jy(1);
else
    Tdhat = (pp(3)-pp(2)^2/4/pp(1));
    khat = ceil(-pp(2)/2/pp(1)/h);  % Indicamos el índice de mestreo para el máximo
end;
Tdhat = max(Tdhat,h);           % clip
Tdhat = h*round(Tdhat/h);       % simulado digitalmente

Tmhat = khat*h;                 %  Tm estimada
Tmhat1 = khat*h;                 % Tm estimada (antes de fiddles)
mode = 0;                       % Bandera si Tm está en la mitad
delfact = 0.2;
if Tmhat < t1 + 0.25*(tstop-t1),    % FODT si Tm está muy pronto
    mode = 1; 
    Tmhat = Tdhat+ delfact*(ty2-Tdhat);
end;
if Tmhat >= tstop,              % IDT si Tm está muy tarde
    Tmhat = tstop; mode = 2; 
end; % en el final o después
khat = round(Tmhat/h);
Tmhat = tout2(khat);

ydot0=y(khat)/(Tmhat-Tdhat); % rampa inicial 
tp = Tdhat + 2*delfact*(ty2-Tdhat); 
ktp = round(tp/h); 
Yp = y(ktp); tp=h*ktp;
yL = ydot0*(tp-Tdhat);
TThat = 0.5*yL*(tp-Tdhat)/(yL-Yp);      % no se usa con IDDT, pero...
K_hat = ydot0/stepsize;  % IDT/IDDT modelo
%Ts = max((2*Tmhat-Tdhat),1.5*(ty2-Tdhat)); 
Ts = max((2*Tmhat-Tdhat),ty2); 

if mode == 1, 
ktd = (Tdhat/h)+1;      % nuevo método: ktd = muestras justo antes de la respuesta
ksecond = k2;
kfit = jstart:ksecond;      % ajuste lineal a los datos Y1 ... Y2
tkfit = h*(kfit'-ktd);      % mueve el origen a Td
ppp = polyfit(tkfit,y(kfit)./tkfit,1);      % y ajuste lineal a y/t
TThat = -0.5*ppp(2)/ppp(1);
K_hat = TThat*ppp(2)/stepsize;
end;        % del 1er orden especial a FODT

%%------------------------------------------------------------------
%  Se muesra la respuesta al escalón en contra de la planta
%%------------------------------------------------------------------
d1 = [0*(1:round(Tdhat/h)) 1];     % para el tiempo muerto, simular como z^{-n}
Tmax = max(4*Tmhat,Tmax);
Tmax = max(Ts+Tmhat,Tmax);
if mode == 2, Ts=Tmax; end;
sim('final1');
% grafica la respuesta
tt =ky:length(tout);
clf; 
subplot(221);
plot(tout2(jy),f,tout2(jy),polyval(pp,tout2(jy)),':');
v = axis;
%axis([v(1) v(2) 0 v(4)]);
legend('\it{f}','fit');
title('\it{f} \rm{and quadratic fit}')
xlabel(['Mode = ',num2str(mode),' T_m = ',num2str(Tmhat1)]);
subplot(222); 
if mode == 0|mode==2,
    te = 'K_p exp(-sT_d)/s'; if mode == 0, te = 'IDDT'; end;
plot(tout2,y,'-',tout,ymodel,':',tout(tt),yplant(tt),'-.',tout,Y1,'--',tout,Y2,'--','LineWidth',lw); 
title(['Fitting by ', te])
xlabel(['K_p = ',num2str(K_hat),' T_d = ',num2str(Tdhat)])
end;
if mode == 1,
plot(tout2,y,'-',tout,ymodel1,':',tout(tt),yplant(tt),'-.',tout,Y1,'--',tout,Y2,'--','LineWidth',lw);
title(['Fitting by FOPDT'])
xlabel(['K = ',num2str(K_hat),' T_d = ',num2str(Tdhat),' T = ',num2str(TThat)])
end
legend([plant,' data'],'model fit','location','SouthEast')

%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%   Ahora la sintonización adaptativa
%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
phi_c = -pi/2; 
ades = 2*zeta;
if intshow == 1, ades = 1/(1+2*zeta); phi_c = -1.25; end;
Kmult = 1/(1+ades^2); 
%Kmult = Kmult*(0.5+zeta);
Kmult=1*Kmult;

% caso IDDT 
if mode == 0,
    wbar_n = 0.5*pi/Tmhat;
    Teff = Tmhat;
end;
% caso especial FODT 
if mode == 1, 
    wbar_n = 1/sqrt(TThat*Tdhat);
    Teff = max(Tdhat,TThat);
end;
% caso especial IDT
if mode == 2, 
    wbar_n =  0.1*pi/Tdhat;
    Teff = Tdhat;
%    Kfiddle=Kfiddle/8; 
    Kw = K_hat/wbar_n/sqrt(2);       % ganancia de la planta a esta frecuencia IDT
end;

Kp=K_hat;
wfinal = wbar_n;      % para ti Ti
wfinalo = wfinal;     % guardar para graficar
if mode == 0, 
    Kw = Kp*abs(exp(-j*wbar_n*Tdhat)-exp(-j*wbar_n*Ts))/wbar_n;
end;
if mode == 1, Kw = Kp/sqrt((1+(wbar_n*TThat)^2)); end;
kfinal = 1/Kw;
kfinalo = kfinal;   % guardar para después
st = 0.01;

%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%   Respuesta en lazo cerrado para sintonización del PI inicial 
%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Tmax_s=Tmax_s/2;
sim('final');
yto=yt; uto=ut;

Tmax_n  = 10*Tmax/Tmhat;  % permitir para 20X el tiempo de asentamiento nominal
Kp0 = Kw;
K0 = 1/Kp0;               % primera ganancia del controlador
wbar = wbar_n; st = 0.1/wbar; Tmax = Tmax_n;
K = Kfiddle*0.25/Teff^2;   % Ganancia adaptativa
Tf = 5/wbar;       % TC del filtro de amplitud-estimación (LAG)
K = Kfiddle*0.1*wbar^2;   % ganacial adaptativa común?

st = 0.1
% si está ruidoso, detune? adaptación, toma más tiempo para los filtro de
% banda
if Snset ~= 0, 
    Sn=Snset*st; % para un caso ruidoso, la misma salida SD que antes
    a = sqrt(2)*sigma;
    fact = 2;
    Tf = Tf*fact;
    beta = beta/fact; 
    K = K/fact; 
    Tmax_c=Tmax_c*fact; 
end;
Tsw = max(10*Tmhat,5/(beta*wbar));    %enciende
TswK = Tsw/2;
%T_l = Tmax_c/2; % medio tiempo para la carga
T_l=0;
Tmax_c = Tmax_c+T_l;
%Tf = 10*Tf

% adapta, simula el últimp PI y plotea
sim('cloop1')
sim('final')
subplot(223); plot(wh(:,1),wh(:,2)/wbar,Khat(:,1),K0./Khat(:,2),':','LineWidth',lw);
title(['Adaptation (normalised)']);
legend('\omega_o','|G_p(j\omega_o)|')
xlabel([' K_a = ',num2str(K)]);
subplot(224); plot(yto(:,1),yto(:,2),':',yt(:,1),yt(:,2),'LineWidth',lw);grid;
title(['CL step: \it{y }\rm for \zeta_D= ',num2str(zeta)]);
xlabel(['K (initial) = ',num2str(kfinalo*Kmult),' (final) = ',num2str((kfinal*Kmult))]);
legend('initial PI','final PI')
c00 = Kmult*kfinalo; c10 = Kmult*kfinalo*wfinalo*ades;
c0 = Kmult*kfinal; c1 = Kmult*kfinal*wfinal*ades;
initcont = [num2str(c00),' + ',num2str(c10),'/s']
finalcont = [num2str(c0),' + ',num2str(c1),'/s']


%subplot(224); plot(uto(:,1),uto(:,2),':',ut(:,1),ut(:,2),0,0,'LineWidth',lw);grid;
%xlabel(['Final K = ',num2str(kfinal*Kmult),' T_i = ',num2str(1/(wfinal*ades))]);
%title('Control signal \it{u}')

% guarda
%print('-depsc2',plant);
%saver_td;
