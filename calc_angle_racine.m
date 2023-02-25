angle_BA = rad2deg(atan((imag(ptB) - imag(ptA))/(real(ptB) - real(ptA))))
angle_CA = rad2deg(pi + atan((imag(ptC) - imag(ptA))/(real(ptC) - real(ptA))))
angle_DA = rad2deg(pi + atan((imag(ptD) - imag(ptA))/(real(ptD) - real(ptA))))
angle_EA = rad2deg(pi + atan((imag(ptE) - imag(ptA))/(real(ptE) - real(ptA))))
angle_FA = rad2deg(atan((imag(ptF) - imag(ptA))/(real(ptF) - real(ptA))))
angle_GA = 90

angle_A = 180-(angle_DA + angle_GA) + (angle_BA + angle_CA + angle_EA + angle_FA)

angle_AB = rad2deg(pi + atan((imag(ptA) - imag(ptB))/(real(ptA) - real(ptB))))
angle_CB = rad2deg(pi + atan((imag(ptC) - imag(ptB))/(real(ptC) - real(ptB))))
angle_DB = rad2deg(pi + atan((imag(ptD) - imag(ptB))/(real(ptD) - real(ptB))))
angle_EB = rad2deg(pi + atan((imag(ptE) - imag(ptB))/(real(ptE) - real(ptB))))
angle_FB = 90;
angle_GB = rad2deg(pi + atan((imag(ptG) - imag(ptB))/(real(ptG) - real(ptB))))

angle_B = 180 -(angle_CB + angle_EB + angle_FB) + (angle_AB + angle_DB + angle_GB)

angle_AC = rad2deg(atan((imag(ptA) - imag(ptC))/(real(ptA) - real(ptC))))
angle_BC = rad2deg(atan((imag(ptB) - imag(ptC))/(real(ptB) - real(ptC))))
angle_DC = rad2deg(atan((imag(ptD) - imag(ptC))/(real(ptD) - real(ptC))))
angle_EC = 90
angle_FC = rad2deg(atan((imag(ptF) - imag(ptC))/(real(ptF) - real(ptC))))
angle_GC = rad2deg(atan((imag(ptG) - imag(ptC))/(real(ptG) - real(ptC))))

angle_C = 180-(angle_BC + angle_EC + angle_FC) + (angle_AC + angle_DC + angle_GC)



