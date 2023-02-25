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
step(sys)
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

subplot(3,2,6)
bode(x(6,1))
grid minor

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

%Règle 7
A = zeros(1)
B = poles(1)
C = poles(3)
D = zeros(3)
E = poles(4)
F = poles(2)
G = zeros(2)

ca

angle_AC = rad2deg(atan((imag(A) - imag(C))/(real(A) - real(C))))
angle_BC = rad2deg(atan((imag(B) - imag(C))/(real(B) - real(C))))
angle_DC = rad2deg(atan((imag(D) - imag(C))/(real(D) - real(C))))
angle_EC = 90
angle_FC = rad2deg(atan((imag(F) - imag(C))/(real(F) - real(C))))
angle_GC = rad2deg(atan((imag(G) - imag(C))/(real(G) - real(C))))

angle_C = 180-(angle_BC + angle_EC + angle_FC) + (angle_AC + angle_DC + angle_GC)

angle_AE = rad2deg(atan((imag(A) - imag(E))/(real(A) - real(E))))
angle_BE = rad2deg(atan((imag(B) - imag(E))/(real(B) - real(E))))
angle_DE = rad2deg(atan((imag(D) - imag(E))/(real(D) - real(E))))
angle_EE = -90
angle_FE = rad2deg(atan((imag(F) - imag(E))/(real(F) - real(E))))
angle_GE = rad2deg(atan((imag(G) - imag(E))/(real(G) - real(E))))

angle_E = 180-(angle_BC + angle_EC + angle_FC) + (angle_AC + angle_DC + angle_GC)

angle_BA = rad2deg(atan((imag(B) - imag(A))/(real(B) - real(A))))
angle_CA = rad2deg(pi + atan((imag(C) - imag(A))/(real(C) - real(A))))
angle_DA = rad2deg(pi + atan((imag(D) - imag(A))/(real(D) - real(A))))
angle_EA = rad2deg(pi + atan((imag(E) - imag(A))/(real(E) - real(A))))
angle_FA = rad2deg(atan((imag(F) - imag(A))/(real(F) - real(A))))
angle_GA = 90

angle_A = 180-(angle_DA + angle_GA) + (angle_BA + angle_CA + angle_EA + angle_FA)

angle_BA = rad2deg(atan((imag(B) - imag(A))/(real(B) - real(A))))
angle_CA = rad2deg(pi + atan((imag(C) - imag(A))/(real(C) - real(A))))
angle_DA = rad2deg(pi + atan((imag(D) - imag(A))/(real(D) - real(A))))
angle_EA = rad2deg(pi + atan((imag(E) - imag(A))/(real(E) - real(A))))
angle_FA = rad2deg(atan((imag(F) - imag(A))/(real(F) - real(A))))
angle_GA = 90


angle_A = 180-(angle_DA + angle_GA) + (angle_BA + angle_CA + angle_EA + angle_FA)



