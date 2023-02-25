angle_BA = rad2deg(atan((imag(B) - imag(A))/(real(B) - real(A))))
angle_CA = rad2deg(pi + atan((imag(C) - imag(A))/(real(C) - real(A))))
angle_DA = rad2deg(pi + atan((imag(D) - imag(A))/(real(D) - real(A))))
angle_EA = rad2deg(pi + atan((imag(E) - imag(A))/(real(E) - real(A))))
angle_FA = rad2deg(atan((imag(F) - imag(A))/(real(F) - real(A))))
angle_GA = 90

angle_A = 180-(angle_DA + angle_GA) + (angle_BA + angle_CA + angle_EA + angle_FA)

angle_AB = rad2deg(pi + atan((imag(A) - imag(B))/(real(A) - real(B))))
angle_CB = rad2deg(pi + atan((imag(C) - imag(B))/(real(C) - real(B))))
angle_DB = rad2deg(pi + atan((imag(D) - imag(B))/(real(D) - real(B))))
angle_EB = rad2deg(pi + atan((imag(E) - imag(B))/(real(E) - real(B))))
angle_FB = 90;
angle_GB = rad2deg(pi + atan((imag(G) - imag(B))/(real(G) - real(B))))

angle_B = 180 -(angle_CB + angle_EB + angle_FB) + (angle_AB + angle_DB + angle_GB)

angle_AC = rad2deg(atan((imag(A) - imag(C))/(real(A) - real(C))))
angle_BC = rad2deg(atan((imag(B) - imag(C))/(real(B) - real(C))))
angle_DC = rad2deg(atan((imag(D) - imag(C))/(real(D) - real(C))))
angle_EC = 90
angle_FC = rad2deg(atan((imag(F) - imag(C))/(real(F) - real(C))))
angle_GC = rad2deg(atan((imag(G) - imag(C))/(real(G) - real(C))))

angle_C = 180-(angle_BC + angle_EC + angle_FC) + (angle_AC + angle_DC + angle_GC)



