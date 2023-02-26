clc
clear all
close all

A = [   -0.018223   -0.088571   -9.78   0;
        -0.003038   -1.2563     0       1;
        0           0           0       1;
        0.0617      -28.078     0       -4.5937];

B = [   0           1.1962;
        0           -0.0012;
        0           0;
        7.84        -4.05];

C = [   1           0           0       0;
        0           57.296      0       0;
        0           0           57.296  0;
        0           0           0       57.296;
        0           -57.296     57.296  0];

D = [   0           0;
        0           0;
        0           0;
        0           0;
        0           0];

states = {'v' 'alpha' 'teta' 'q'};
inputs = {'deltaC' 'aprop'};
outputs = {'v' 'alpha', 'teta', 'q', 'gamma'};

FTBO = ss(A,B,C,D,'statename',states,...
'inputname',inputs,...
'outputname',outputs);

G = tf([33.7002], [1 5.8402 33.702])
G2 = tf([0.0472], [ 1 0.0282 0.0472])
x = tf(FTBO)
step(FTBO)
grid minor
figure()
subplot(3,2,1)
bode(x(1,1))
grid minor

subplot(3,2,2)
bode(x(2,1))
grid minor

subplot(3,2,3)
bode(x(3,1))
grid minor

subplot(3,2,4)
bode(x(4,1))
grid minor

subplot(3,2,5)
bode(x(5,1))
grid minor

[num_tfva, den_tfva] = tfdata(x(1,2));
num_tfva = num_tfva{1};
den_tfva = den_tfva{1};

%% Analyse des différents mode dynamique de l'avion
clc
figure();

[zeros, poles, K] = tf2zp(num_tfva, den_tfva);
disp(["Valeur de TF pour la sortie v et l'entrée a: "])
x(1,2)
disp(["Valeur des poles pour le mode phygoide : ", poles(1)])
disp(["Valeur des poles pour le mode short-period : ", poles(3)])
disp(["Analyse du mode phygoide : "])
phi_phy = atan(imag(poles(1))/abs(real(poles(1))));
zeta_phy = cos(phi_phy);
R_phy = abs(poles(1));
wn_phy = cos(phi_phy)*R_phy/zeta_phy;
wa_phy = wn_phy*(1-zeta_phy^2)^0.5;
Mp_phy = 100*exp(-pi/tan(phi_phy));
ts_phy = 4/(zeta_phy*wn_phy);

disp(["Valeur de phi : ", phi_phy]);
disp([ "Valeur de zeta : ", zeta_phy]);
disp([ "Valeur de R : ", R_phy]);
disp(["Valeur de wn : ", wn_phy]);
disp(["Valeur de wa :", wa_phy]);
disp([ "Valeur de Mp : " Mp_phy]);
disp([ "Temps de stabilisation (2%)", ts_phy]);
[phygoid_num, phygoid_den] = zp2tf( [], poles(1:2), 1);
disp("Fonction de transfert de Phygoid: ")
G_phygoid = 33.74*tf(phygoid_num, phygoid_den)

disp("Valeur théorique de matlab:")
damp(G_phygoid)
step(G_phygoid)

figure();
disp(["Analyse du mode short-period : "])
phi_sp = atan(imag(poles(3))/abs(real(poles(3))));
zeta_sp = cos(phi_sp);
R_sp = abs(poles(1));
wn_sp = cos(phi_sp)*R_sp/zeta_sp;
wa_sp = wn_sp*(1-zeta_sp^2)^0.5;
Mp_sp = 100*exp(-pi/tan(phi_sp));
ts_sp = 4/(zeta_sp*wn_sp);

disp(["Valeur de phi : ", phi_sp]);
disp([ "Valeur de zeta : ", zeta_sp]);
disp([ "Valeur de R : ", R_sp]);
disp(["Valeur de wn : ", wn_sp]);
disp(["Valeur de wa :", wa_sp]);
disp([ "Valeur de Mp : " Mp_sp]);
disp([ "Temps de stabilisation (2%)", ts_sp]);
[sp_num, sp_den] = zp2tf( [], poles(3:4), 1);
disp("Fonction de transfert de Phygoid: ")
G_sp = 0.04719*tf(sp_num, sp_den)

disp("Valeur théorique de matlab:")
damp(G_sp)
step(G_sp)

%% Dessin du lieu des racine

rlocus(x(1,2))

n = size(poles);
n = n(1);
m = size(zeros);
m = m(1);
disp(["nombre de branche = ", n])
disp(["Nombre d'asymptote = ",n-m])
disp(["Direction des asymptotes = ", 180/(n-m)*(2*0+1)])
phi = (sum(real(poles)) - sum(real(zeros)))/(n-m)
disp(["Intersection sur l'axe des réels = ", phi])

ptA = zeros(1)
ptB = poles(1)
ptC = poles(3)
ptD = zeros(3)
ptE = poles(4)
ptF = poles(2)
ptG = zeros(2)
calc_angle_racine

disp(["Angle A = ", angle_A])
disp(["Angle B = ", angle_B])
disp(["Angle C = ", angle_C])
disp(["Voir règle 7 et 9 sur document papier"])
Kv = 1.0263;
disp(["Kv trouvé : ", Kv])
disp("Aucune intersection possible sur l'axe des imaginaire")

p = rlocus(x(1,2),1.0263);
hold on
grid minor
plot(p, 'p', 'markerSize', 15)

%[Kv,POLES] = rlocfind(FTBO(1,2))

%% n1/d1 Conception de la boucle interne
C1 = C(5, :); %Enlever C(1) car c'est une sortie qui ne sera pas utilisé
A1 = A - B(:,2)*Kv*C(1,:); 
B1 = B(:,1);
D1 = [0]'; %Une sortie seulement

states = {'v' 'alpha' 'teta' 'q'};
inputs = {'deltaC'};
outputs = {'gamma'};

boucle_int = ss(A1, B1, C1, D1,'statename',states,...
'inputname',inputs,...
'outputname',outputs);

figure()
bode(Kv*FTBO(1,2))
[Gm,Pm,Wcg,Wcp] = margin(FTBO(1,2))
disp(["Gain margin:", Gm])
disp(["Phase margin:", Pm])
disp(["Wcg : ", Wcg])
disp(["Wcp : ", Wcp])
grid on



%% Pour méthode analytique
[num, den] = tfdata(FTBO(1,2));
[R,P,K] = residue(num{1},den{1});
C = abs(R)./(abs(real(P)))
[num_red, den_red] = residue(R(3:4), P(3:4), K)
FTBO_red = tf(num_red, den_red)
disp(["Voir démarche pour nouveau Kv : ", -1.3])
disp("Fonction de transfert d'ordre 2 en boucle ouverte: ")
figure();
hold on
FTBO_red
rlocus(FTBO_red)
rlocus(x(1,2))

%% Analyse de A1 B1 C1 D1
figure();
hold on

C1 = C(5, :); %Enlever C(1) car c'est une sortie qui ne sera pas utilisé
A1 = A - B(:,2)*Kv*C(1,:); 
B1 = B(:,1);
D1 = [0]'; %Une sortie seulement

[num_1, den_1] = ss2tf(A1, B1, C1, D1);
TFBF_1 = tf(num_1,den_1);

C1_2 = C(5, :); %Enlever C(1) car c'est une sortie qui ne sera pas utilisé
A1_2 = A - B(:,2)*1.3*C(1,:); 
B1_2 = B(:,1);
D1_2 = [0]'; %Une sortie seulement

[num_2, den_2] = ss2tf(A1_2, B1_2, C1_2, D1_2);
TFBF_2 = tf(num_2,den_2);

rlocus(TFBF_1);
rlocus(TFBF_2);



%%
%Calcul simple pour déterminer kp: Prendre le point à
%-180 (fréquence) trouver le gain à cette fréquence (11dB), Trouver combien
%on doit diminuer la courbe pour avoir un GM de 6dB (-17dB ou -18dB).
%Déterminer ce gain et le multiplier à TFBF_1
figure();
hold on
bode(TFBF_1)
margin(TFBF_1)
[Gm, Pm, wcg, wcp] = margin(TFBF_1)
Kp = 10^(((20*log10(Gm)-6))/20);
bode(Kp*TFBF_1)
margin(Kp*TFBF_1)
erreur = 1/(1+10^(5.34/20));
disp(["L'erreur est de : " , erreur])


%% Afficher la réponse à l'échelon
%L'erreur en régime permanent est la même et l'erreur n'est pas nul
figure()
TFBF_1_FB = feedback(Kp*TFBF_1,1)
step(TFBF_1_FB)
xlim([0 14])
disp(["L'erreur est de : " , 1-0.65])




%% PD, PI, PID
num_1_PD = Kp.*[1 1];
den_1_PD = [1];
num_1_PI = Kp.*[1 1];
den_1_PI = [1 0];
num_1_PID = Kp.*[1 1 1];
den_1_PID = [1 0 0];


tf_pd = tf(num_1_PD, den_1_PD);
tf_pi = tf(num_1_PI, den_1_PI);
tf_pid = tf(num_1_PID, den_1_PID);

TFBO_1_p = Kp*TFBF_1;
TFBO_1_pd = tf_pd*TFBF_1;
TFBO_1_pi = tf_pi*TFBF_1;
TFBO_1_pid = tf_pid*TFBF_1;

TFBF_2_p = feedback(TFBO_1_p, 1);
TFBF_2_pd = feedback(TFBO_1_pd, 1);
TFBF_2_pi = feedback(TFBO_1_pi, 1);
TFBF_2_pid = feedback(TFBO_1_pid, 1);

figure();
hold on
step(TFBF_2_p)
step(TFBF_2_pd)
step(TFBF_2_pi)
step(TFBF_2_pid)
legend('p', 'pd', 'pi', 'pid')


