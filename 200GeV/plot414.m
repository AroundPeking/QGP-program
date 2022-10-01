clc;clear;close all
load Tu_Jpsi_c020.dat
load Ts_Jpsi_c020.dat
load Tc_Jpsi_c020.dat
load Sj_c.dat
load Sj_cbar.dat
figure(1)
semilogy(Tu_Jpsi_c020(:,1),Tu_Jpsi_c020(:,2),'r');hold on
semilogy(Ts_Jpsi_c020(:,1),Ts_Jpsi_c020(:,2),'g');hold on
semilogy(Tc_Jpsi_c020(:,1),Tc_Jpsi_c020(:,2),'b');
legend('T_u','T_s','T_c')
figure(2)
semilogy(Sj_c(:,1),Sj_c(:,2),'r');hold on
semilogy(Sj_cbar(:,1),Sj_cbar(:,2),'b');hold on
legend('S_c','S_{cbar}')
figure(3)
semilogy(Sj_cbar(:,1),(1-exp(-Sj_cbar(:,1)/2)).*Sj_c(:,2).*Sj_cbar(:,2),'b');hold on