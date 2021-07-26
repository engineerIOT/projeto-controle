%Projeto 1 Elvis
clc
clear all
close all

format long; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLANTA

%Dados de entrada do projeto (delta1,delta2,tp)
delta1 = 0.144; %retirado da figura 
delta2 = 0.496; %retirado da figura 
tp = 0.0192;    %retirado da figura 

%Cálculo do sobressinal(Mp)
disp('Mp: Sobre-Sinal da planta:')
Mp_planta = delta1/delta2

%Calculo de Zeta a partir do Sobre-sinal
disp('Zeta: Fator de amortecimento da planta:')
zeta_planta = sqrt(((log(Mp_planta))^2)/(pi^2+(log(Mp_planta))^2))

% Cálculo da freqüência não-amortecida da planta(Wn)    
disp('Wn: Frequência natural não amortecida da planta:')
Wn_planta = pi/(tp*sqrt(1-zeta_planta^2))

%Cálculo do Numerador
num = Wn_planta^2;

%Cálculo do Denominador
den = [1 2*Wn_planta*zeta_planta Wn_planta^2];

%Cálculo da Função de Transferência da Planta
disp('Função de Transferência da Planta G(s)')
G = tf(num,den)   

%Cálculo do Tempo de subida da planta (tr)
disp('Tempo de subida da (tr) planta  G(s)')
tr_planta =  (2.16*zeta_planta + 0.6)/Wn_planta

%Cálculo do Tempo de acomodação (ts5%_planta)
disp('Tempo de acomodação (ts5%_planta G(s)')
ts5_planta = 3/(zeta_planta*Wn_planta)

%Cálculo do Tempo de pico da planta (tp_planta)
disp('Tempo de pico (tp) da planta G(s)')
tp_planta = pi / (Wn_planta*sqrt(1-zeta_planta.^2))

%Cálculo do sobressinal da planta (Mp_planta)
disp('Sobressinal (Mp%) da planta G(s)')
Mp_planta=exp(-pi*zeta_planta/sqrt(1-zeta_planta^2))

figure(1)
step(G)
title(' Função de Transferência da Planta G(s)')

% NT = 42 (Elvis Roberto de Jesus Avila Carvalho Fernandes
NT = 0.040; 
Ts5 = NT;

% N = 5 (Elvis)
N = 5;
disp('Mp: Sobre-Sinal desejado:')
Mp = 2*N/100

% Calculo de Zeta desejado a partir do Sobre-sinal do controlador digital 
disp('Zeta: Fator de amortecimento desejado do controlador digital')
zeta = sqrt(((log(Mp))^2)/(pi^2+(log(Mp))^2))

% Cálculo da freqüência não-amortecida do controlador digital(Wn) 
disp('Wn: Frequencia natural não amortecida do controlador digital:')
Wn = 3/(Ts5 * zeta)

% Cálculo da freqüência amortecida do controlador digital(Wd)
disp('Wd: Frequência natural amortecida do controlador digital:')
Wd = Wn * (sqrt(1-zeta^2))

% Cálculo do período de amostragem (T)
% Uma boa prática é adotar o Período de Amostragem igual a
% 10 ou 15 vezes menos que o Tempo de Acomodação Ts5%
disp('Ta: Período de Amostragem:')
T = Ts5 / 15

% Cálculo da freqüência de amostragem (Ws)
disp('Ws: Frequência de amostragem desejada:')
Ws = (2*pi) / T

% Cálculo do número de amostras (Na)
disp('Na: Número de amostras por ciclo de oscilação:')
Na = Ws/Wd

% Cálculo do Módulo  de Z do pólo dominante(Z)
disp(' Módulo de Z:')
z1_modulo = exp( ((-2 * pi * zeta) / sqrt( 1 - zeta^2 )) * (Wd/Ws) )

% Cálculo do ângulo de Z do pólo dominante em graus(Z)
disp('Ângulo de Z em Graus:')
z1_angulo = (2 * pi * Wd) / Ws
z1_angulo = rad2deg( z1_angulo )

% Valor do Pólo desejado
s1 = -zeta*Wn + Wd*i
disp('Valor do Pólo desejado:')
z1 = exp ( T * s1)

% Localização do Pólo desejado no Plano Z
figure(2)
zplane( 0, z1 )
title(' Pólo de Malha Fechada Desejado no Plano Z ')

%Cálculo da Função de Transferência discreta da Planta
disp(' Função de Transferência Discreta G(z) ')
Gz = c2d( G, T)
zpk(Gz)
figure(3)
step(Gz)
title(' Função de Transferência da Planta G(z)')
figure(4)
%rlocus(Gz);
pzmap(Gz)

% considalfaerando que o zero do controlador (alfa) cancela os polos de Gz
[numz,denz]=tfdata(Gz,'v')
alfa = denz

n2=numz;
d2=[1 -1];

% fi2 é o ângulo de G2z quando z=z1
fi2=angle(polyval(n2,z1)/polyval(d2,z1))

display('-----Contribuições angulares: -----')
anguloGz1 = angle(fi2)*(180/pi)

% fi1 é o ângulo de G1z quando z=z1
fi1=-pi-fi2
beta=(imag(z1)-real(z1)*tan(-fi1))/tan(-fi1)

num2=[1];
den2=[1 beta];
Cz1=tf(num2,den2,T)

num3=[alfa];
den3=[1 -1];  %integrador
Cz2=tf(num3,den3,T)

Cz=Cz1*Cz2
zpk(Cz)

FTMA1=minreal(Cz*Gz)
disp('FTMA sem K')
zpk(FTMA1)
figure(5)
rlocus(FTMA1);
%zplane;

[NFTMA,DFTMA]=tfdata(FTMA1,'v');
K=abs((polyval(NFTMA,z1)/polyval(DFTMA,z1)));
valork=K
K=1/real(K)

FTMA2=minreal(Cz*Gz*K)
disp('FTMA2 com K')
zpk(FTMA2)

FTMF=feedback(FTMA2,1)

figure(6)
step(FTMF) %para FTMF
%step(Gz) %para Gz

%step(Cz*K)


%Resposta ao degrau
[sys,kT]= step(FTMF); %para Gz
%[sys,kT]= step(Gz); %para Gz
%[sys,kT]= step(Gz);
%[sys,kT]= step(Cz*K);

ck=zeros(1,length(kT));

e(1)=1;
e(2)=1;
c(1) = 0;
c(2) = 0.04941;
y(1)=0;
y(2)=0.09662;
u(1)=0.51141*e(1);
u(2)=0.5114*e(2)-0.7807*e(1)+0.7807*u(1);

for k=3:length(kT)
   
   c(k)=c(k-1)*1.577-c(k-2)*0.6703+0.04941+0.04403 %FTMF
   y(k)=y(k-1)*1.527-y(k-2)*0.7093+0.09662+0.0861  %Gz
   e(k)=ck(k)-y(k); %Somador
   u(k)=u(k-1)*1.626-u(k-2)*0.6263+e(k)*0.5114-e(k-1)*0.7807+e(k-2)*0.3627 %Cz*K
 end
hold on
step(G)
plot(kT,c,'*r'); %para FTMF
%plot(kT,y,'*r'); %para Gz
hold off
legend('FTMF','G','Equação Recursiva')

%erro FTMF
[num1,den1]=tfdata(FTMF,'v');
n1 = polyval(num1,1);
d1 = polyval(den1,1);
Valor_final_FTMF = n1/d1

erroFTMF = 1/(Valor_final_FTMF+1)

%erro Gz
[num2,den2]=tfdata(Gz,'v');
n2 = polyval(num2,1);
d2 = polyval(den2,1);
Valor_final_Gz = n2/d2

erroGz= 1/(Valor_final_Gz+1)

