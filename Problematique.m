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

sys = ss(A,B,C,D,'statename',states,...
'inputname',inputs,...
'outputname',outputs);

G = tf([33.7002], [1 5.8402 33.702])
G2 = tf([0.0472], [ 1 0.0282 0.0472])
x = tf(sys)
% step(sys)
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
x(1,2)
[num, den] = tfdata(x(1,2))
num = num{1};
den = den{1};
 

[zeros, poles, K] = tf2zp(num, den)
rlocus(x(1,2))

% Ràgle 5
phi = (sum(real(poles)) - sum(real(zeros)))

%règle 6

%% Règle 7
A = zeros(1)
B = poles(1)
C = poles(3)
D = zeros(3)
E = poles(4)
F = poles(2)
G = zeros(2)

calc_angle_racine