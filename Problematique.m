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
FTVA = tf(num_tfva, den_tfva);

%% Analyse des différents mode dynamique de l'avion
clc
close all
figure();

[R,P,K_res] = residue(num_tfva,den_tfva);

[zeros, poles, K] = tf2zp(num_tfva, den_tfva);
disp(["Valeur de TF pour la sortie v et l'entrée a: "])
x(1,2)
disp(["Valeur des poles pour le mode phugoide : ", poles(1)])
disp(["Valeur des poles pour le mode short-period : ", poles(3)])
disp(["Analyse du mode phugoide : "])
phi_phu = atan(imag(poles(3))/abs(real(poles(3))));
zeta_phu = cos(phi_phu);
R_phu = abs(poles(3));
wn_phu = R_phu;
wa_phu = wn_phu*(1-zeta_phu^2)^0.5;
Mp_phu = 100*exp(-pi/tan(phi_phu));
ts_phu = 4/(zeta_phu*wn_phu);

disp(["Valeur de phi : ", phi_phu]);
disp([ "Valeur de zeta : ", zeta_phu]);
disp([ "Valeur de R : ", R_phu]);
disp(["Valeur de wn : ", wn_phu]);
disp(["Valeur de wa :", wa_phu]);
disp([ "Valeur de Mp : " Mp_phu]);
disp([ "Temps de stabilisation (2%)", ts_phu]);

disp("Fonction de transfert de Phugoid: ")
G_phu_num = [R(3)+R(4) -R(3)*P(4)-R(4)*P(3)];
G_phu_den = [1 -P(3)-P(4) P(3)*P(4)];
G_phu = tf(G_phu_num, G_phu_den)

disp("Valeur théorique de matlab:")
damp(G_phu)
[y,t] = step(G_phu);
plot(t,y,'LineWidth',3)
xlabel('temps(s)', 'FontSize',20)
ylabel('Amplitude', 'FontSize',20)
ax = gca;
ax.FontSize = 15;  % Font Size of 15
title('Réponse échelon du mode dynamique phugoide', FontSize=25)
grid minor

figure();
disp(["Analyse du mode short-period : "])
phi_sp = atan(imag(poles(1))/abs(real(poles(1))));
zeta_sp = cos(phi_sp);
R_sp = abs(poles(1));
wn_sp = cos(phi_sp)*R_sp/zeta_sp;
wa_sp = wn_sp*(1-zeta_sp^2)^0.5;
Mp_sp = 100*exp(-pi/tan(phi_sp));
ts_sp = 4/(zeta_sp*wn_sp);
grid minor

disp(["Valeur de phi : ", phi_sp]);
disp([ "Valeur de zeta : ", zeta_sp]);
disp([ "Valeur de R : ", R_sp]);
disp(["Valeur de wn : ", wn_sp]);
disp(["Valeur de wa :", wa_sp]);
disp([ "Valeur de Mp : " Mp_sp]);
disp([ "Temps de stabilisation (2%)", ts_sp]);

disp("Fonction de transfert courte période: ")

G_sp_num = [R(1)+R(2) -R(1)*P(2)-R(2)*P(1)];
G_sp_den = [1 -P(1)-P(2) P(1)*P(2)];
G_sp = tf(G_sp_num, G_sp_den)

disp("Valeur théorique de matlab:")
damp(G_sp)
[y,t] = step(G_sp);
plot(t,y,'LineWidth',3)
xlabel('temps(s)', 'FontSize',20)
ylabel('Amplitude', 'FontSize',20)
ax = gca;
ax.FontSize = 15;  % Font Size of 15
title('Réponse échelon du mode dynamique à courte période', FontSize=25)
grid minor

[y,t] = step(G_sp*G_phu);
plot(t,y,'LineWidth',3)
xlabel('temps(s)', 'FontSize',20)
ylabel('Amplitude', 'FontSize',20)
ax = gca;
ax.FontSize = 15;  % Font Size of 15
title('Réponse échelon des modes(ordre 4)', FontSize=25)
grid minor

%% Dessin du lieu des racine
close all
rlocusplot(x(1,2))

hm = findobj(gca, 'Type', 'Line');          % Handle To 'Line' Objects
hm(6).MarkerSize = 12;  
hm(7).MarkerSize = 12; 
hm(2).LineWidth = 2;
hm(3).LineWidth = 2;
hm(4).LineWidth = 2;
hm(5).LineWidth = 2;
title('Lieu des racines', FontSize=25)
xlabel('Axe reel', 'FontSize',20)
ylabel('Axe imaginaire','FontSize',20)
xlim([-4 0])
grid minor

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
plot(p, 'p', 'markerSize', 15)

% [Kv,POLES] = rlocfind(FTBO(1,2))

%% n1/d1 Conception de la boucle interne
close all
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

h = bodeplot(Kv*FTBO(1,2));
options = getoptions(h);
options.Title.String = 'Lieu de Bode';
options.XLabel.FontSize = 20;
options.YLabel.FontSize = 20;
options.Title.FontSize = 25;
setoptions(h,options);

[Gm,Pm,Wcg,Wcp] = margin(Kv*FTBO(1,2)) %%Kv
disp(["Gain margin:", Gm])
disp(["Phase margin:", Pm])
disp(["Wcg : ", Wcg])
disp(["Wcp : ", Wcp])
grid on



%% Pour méthode analytique
close all
[num, den] = tfdata(FTBO(1,2));
[R,P,K] = residue(num{1},den{1});
Coef = abs(R)./(abs(real(P)))
[num_red, den_red] = residue(R(3:4), P(3:4), K)
FTBO_red = tf(num_red, den_red)
disp(["Voir démarche pour nouveau Kv : ", 1.2987])
disp("Fonction de transfert d'ordre 2 en boucle ouverte: ")
FTBO_red
figure();
subplot(2,1,1)
rlocus(FTBO_red)
hm = findobj(gca, 'Type', 'Line');          % Handle To 'Line' Objects
hm(5).MarkerSize = 12;  
hm(4).MarkerSize = 12; 
hm(1).LineWidth = 2;
hm(2).LineWidth = 2;
hm(3).LineWidth = 2;

title('Lieu des racines', FontSize=25)
xlabel('Axe reel', 'FontSize',20)
ylabel('Axe imaginaire','FontSize',20)
xlim([-4 0])
grid minor



subplot(2,1,2)
rlocus(x(1,2))
hm = findobj(gca, 'Type', 'Line');          % Handle To 'Line' Objects
hm(6).MarkerSize = 12;  
hm(7).MarkerSize = 12; 
hm(2).LineWidth = 2;
hm(3).LineWidth = 2;
hm(4).LineWidth = 2;
hm(5).LineWidth = 2;
title('Lieu des racines', FontSize=25)
xlabel('Axe reel', 'FontSize',20)
ylabel('Axe imaginaire','FontSize',20)
xlim([-4 0])
grid minor

%% Analyse de A1 B1 C1 D1
close all
figure();
hold on

C1 = C(5, :); %Enlever C(1) car c'est une sortie qui ne sera pas utilisé
% C1 = C([1,5], :);
A1 = A - B(:,2)*Kv*C(1,:); 
B1 = B(:,1);
D1 = [0]'; %Une sortie seulement

[num_1, den_1] = ss2tf(A1, B1, C1, D1);
TFBF_1 = tf(num_1,den_1);

C1_2 = C(5, :); %Enlever C(1) car c'est une sortie qui ne sera pas utilisé
A1_2 = A - B(:,2)*1.2987*C(1,:); 
B1_2 = B(:,1);
D1_2 = [0]'; %Une sortie seulement

[num_2, den_2] = ss2tf(A1_2, B1_2, C1_2, D1_2);
TFBF_2 = tf(num_2,den_2);

rlocus(TFBF_1);
rlocus(TFBF_2);

hm = findobj(gca, 'Type', 'Line');          % Handle To 'Line' Objects
hm(6).MarkerSize = 12;  
hm(7).MarkerSize = 12; 
hm(12).MarkerSize = 12;  
hm(13).MarkerSize = 12; 
hm(2).LineWidth = 2;
hm(3).LineWidth = 2;
hm(4).LineWidth = 2;
hm(5).LineWidth = 2;
hm(8).LineWidth = 2;
hm(9).LineWidth = 2;
hm(10).LineWidth = 2;
hm(11).LineWidth = 2;
title('Comparaison lieu de racine', FontSize=25)
xlabel('Axe reel', 'FontSize',20)
ylabel('Axe imaginaire','FontSize',20)

legend('Kv = 1.0263', 'Kv = 1.2987', 'FontSize', 15)
grid on
xlim([-5 0])
ylim([-8 8])


%%
%Calcul simple pour déterminer kp: Prendre le point à
%-180 (fréquence) trouver le gain à cette fréquence (11dB), Trouver combien
%on doit diminuer la courbe pour avoir un GM de 6dB (-17dB ou -18dB).
%Déterminer ce gain et le multiplier à TFBF_1
close all
figure();
hold on

h = bodeplot(TFBF_1);
options = getoptions(h);
options.XLabel.FontSize = 20;
options.YLabel.FontSize = 20;
options.Title.FontSize = 25;
setoptions(h,options);




margin(TFBF_1)
[Gm, Pm, wcg, wcp] = margin(TFBF_1)
Kp = 10^(((20*log10(Gm)-6))/20);
% bode(Kp*TFBF_1)
% margin(Kp*TFBF_1)
erreur = 1/(1+10^(5.34/20));
disp(["L'erreur est de : " , erreur])
grid on
% legend('FTBF','FTBF*Kp', 'FontSize', 15)


%% Afficher la réponse à l'échelon
%L'erreur en régime permanent est la même et l'erreur n'est pas nul
figure()
TFBF_1_FB = feedback(Kp*TFBF_1,1)
step(TFBF_1_FB)
xlim([0 14])
disp(["L'erreur est de : " , 0.3510])




%% PD, PI, PID
close all

% num_1_PD = Kp.*[1 1];
% den_1_PD = [1];
num_1_PI = Kp.*[1 1];
den_1_PI = [1 0];
num_1_PID = Kp.*[1 1 1];
den_1_PID = [0 1 0];
test = ss(A1, B1, C1, D1)

C2 = C1
A2 = A1 - B1*Kp*C2;
B2 = B1*Kp;
D2 = [0]';

[n2, d2] = ss2tf(A2, B2, C2, D2);
TFBF_2_p = tf(n2,d2)


tf_pd = tf(num_1_PD, den_1_PD);
tf_pi = tf(num_1_PI, den_1_PI);
tf_pid = tf(num_1_PID, den_1_PID);

% TFBO_1_p = Kp*TFBF_1;
TFBO_1_pd = tf_pd*TFBF_1;
TFBO_1_pi = tf_pi*TFBF_1;
TFBO_1_pid = tf_pid*TFBF_1;

% TFBF_2_p = feedback(TFBO_1_p, 1);
TFBF_2_pd = feedback(TFBO_1_pd, 1);
TFBF_2_pi = feedback(TFBO_1_pi, 1);
TFBF_2_pid = feedback(TFBO_1_pid, 1);

figure();
hold on
step(TFBF_2_p)
step(TFBF_2_pd)
step(TFBF_2_pi)
step(TFBF_2_pid)
xlabel('Temps', 'Fontsize',20)
ylabel('Amplitude', 'Fontsize',20)
legend('p', 'pd', 'pi', 'pid', 'Fontsize', 15)
title('Réponse à l''échelon', 'Fontsize', 25)
grid minor


