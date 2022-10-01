clc;clear;close all
load int_Sj_pt2_q2_50.dat
load int_Sj_pt2_q2d3_50.dat
semilogy(int_Sj_pt2_q2_50(:,1),int_Sj_pt2_q2_50(:,2),'r');hold on
semilogy(int_Sj_pt2_q2d3_50(:,1),int_Sj_pt2_q2d3_50(:,2),'g');hold on
legend('2.0','2.3')