%Projeto 1 Controle II
%Elvis Fernandes
clear;clc;close all;
%G(s)=175,84^2/(s^2 +128,829s +175,84^2)
numG = [0 0 175.84^2]; denG = [1 128.829 175.84^2];
G = tf(numG,denG);

%Requisitos de Projeto:
%ts5%=42ms   MP=10%  %Ta=2.67ms  para entrada degrau de 1 para 1,5
Mp=0.1; Ta = 0.002666666666667; 

%Calculo Fator de Amortecimento (Zeta)
Zeta=sqrt(log(Mp)^2/(pi^2+log(Mp)^2))

% C�lculo da freq��ncia de amostragem (Ws)
Ws = (2*pi) / Ta

% C�lculo do n�mero de amostras (Na)
Na =  10

% C�lculo da freq��ncia amortecida Wd)
Wd = Ws/Na

% C�lculo da freq��ncia amortecida do controlador anal�gico(Wd)
%Calculo Frequ�ncia Natural n�o amortecida (Wn)
Wn = Wd / (sqrt(1-Zeta^2))

%Calculo Frequ�ncia Natural n�o amortecida (Wn)
Ts5=3/(Zeta*Wn)

Gz = c2d(G,Ta);
%Determina��o do polo dominante no plano s
s1 = -Zeta*Wn+1j*Wn*sqrt(1-Zeta^2)

%M�dulo de s1
Ms = abs(s1)

%Angulo de s1
thetaS = angle(s1)

% g1 � o �ngulo de G1 quando s=s1
g1 = polyval(numG,s1)/polyval(denG,s1)

%M�dulo de g1
Mg = abs(g1)

%Angulo de g1
thetaG = angle(g1)

Ki = 110
Kp = (-sin(thetaG-thetaS)/(Mg*sin(thetaS))) - (2*Ki*cos(thetaS))/Ms
Kd = (sin(thetaG)/(Ms*Mg*sin(thetaS))) + Ki/Ms^2

numC = [Kd+(Ki*Ta^2)+Kp*Ta -(2*Kd+Kp*Ta) Kd];
denC = [Ta -Ta 0];
Cz = zpk(tf(numC,denC,Ta))

%Fun��o de Transfer�ncia  em malha aberta
FTMA = zpk(minreal(Cz*Gz))

%Fun��o de Transfer�ncia  em malha fechada
FTMF = feedback(FTMA,1)

%Resposta ao degrau
figure(2);
step(FTMF)
figure(3);
rlocus(FTMF)
hold off

figure (4)
hold on
step(G)
hold on
step(FTMF)
legend('Gs','FTMF')

