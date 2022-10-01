clc;clear;close all
mat=readtable('dNptdpt_D0_c020_200GeV_final.dat');
mat=table2array(mat);
%%
pt_Exper=[0.315 0.968 1.35 1.89 2.56 3.75 5.88];
y_Exper=[0.0855 0.0963 0.0732 0.0236 0.00352 0.000461 1.47*10^-5];
y_error=[0.0256 0.0198 0.0126 0.0041 0.00104 0.000117 4.00*10^-6];
pt=double(mat(:,1));

TT=mat(:,2);
TS=mat(:,3);
SS=mat(:,4);
Total=TT+TS+SS;
%%
figure;
ax=axes;
plt_0=errorbar(pt_Exper,y_Exper,y_error,'k');ax.YScale='log';
plt_0.Marker='o';
plt_0.MarkerSize=4;
plt_0.MarkerFaceColor='k';
plt_0.LineStyle='none';
hold on
%=====
plt_1=semilogy(pt,Total,'g-');
plt_1.LineWidth=1.2;
hold on
%=====
plt_2=semilogy(pt,TT,'k-.');
plt_2.LineWidth=1;
hold on
%=====
plt_3=semilogy(pt,TS,'k--');
plt_3.LineWidth=1;
hold on
%=====
plt_4=semilogy(pt,SS,'k-');
plt_4.LineWidth=1;
hold on
%=====
xlim([0,7])
 ylim([10^-8,10^0])
xlabel('p_T(GeV/c)');xticks(0:7);
ylabel('1/2\pip_Td^2N/dydp_T(GeV/c)^{-2}');yticks([10^-8,10^-6,10^-4,10^-2,10^0]);
legend('STAR D^0','sum','therml-thermal','thermal-shower','shower-shower')
ax.FontName='Times New Roman';
ax.FontSize=12;
ax.FontWeight='bold';

% hold on;
% load SS_Jpsi_c020.dat;
% semilogy(SS_Jpsi_c020(:,1),SS_Jpsi_c020(:,2));