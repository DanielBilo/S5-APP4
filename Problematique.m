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
% step(FTBO)
% grid minor
% figure()
% subplot(3,2,1)
% bode(x(1,1))
% grid minor
% 
% subplot(3,2,2)
% bode(x(2,1))
% grid minor
% 
% subplot(3,2,3)
% bode(x(3,1))
% grid minor
% 
% subplot(3,2,4)
% bode(x(4,1))
% grid minor
% 
% subplot(3,2,5)
% bode(x(5,1))
% grid minor


%% Dessin du lieu des racine
[num, den] = tfdata(x(1,2));
num = num{1};
den = den{1};
 

[zeros, poles, K] = tf2zp(num, den);
% rlocus(x(1,2))
% [Kv,POLES] = rlocfind(FTBO(1,2))
Kv = 1.0263;
p = rlocus(x(1,2),1.0263);
hold on
grid minor
% plot(p, 'p', 'markerSize', 15)

% Regle 5
phi = (sum(real(poles)) - sum(real(zeros)))

%règle 6

%% Règle 7
ptA = zeros(1)
ptB = poles(1)
ptC = poles(3)
ptD = zeros(3)
ptE = poles(4)
ptF = poles(2)
ptG = zeros(2)

calc_angle_racine

%% n1/d1
C1 = C([1,5], :);
A1 = A - B(:,2)*Kv*C(1,:);
B1 = B(:,1);
D1 = [0 0]';

states = {'v' 'alpha' 'teta' 'q'};
inputs = {'deltaC'};
outputs = {'v', 'gamma'};

boucle_int = ss(A1, B1, C1, D1,'statename',states,...
'inputname',inputs,...
'outputname',outputs);

figure()
bode(Kv*FTBO(1,2))
[Gm,Pm,Wcg,Wcp] = margin(FTBO(1,2))
grid on

%% Pour méthode analytique
[num, den] = tfdata(FTBO(1,2));
[R,P,K] = residue(num{1},den{1});

C = abs(R)./(abs(real(P)))

[num_red, den_red] = residue(R(3:4), P(3:4), K)

FTBO_red = tf(num_red, den_red)
rlocus(FTBO_red)
