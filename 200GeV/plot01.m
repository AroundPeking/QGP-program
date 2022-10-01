clc;clear;close all
mat=readtable('dNptdpt_Jpsi_c020.dat');
mat=table2array(mat);
%%
pt_Exper=[0.5 1.5 2.5 3.5 4.5];
y_Exper=10^-6*[24.800 12.800 4.020 0.431 0.420 ];
y_error=10^-6*[3.180 1.460 0.739 0.337 0.179];
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
xlim([0,6.5])
 ylim([10^-10,10^-4])
xlabel('p_T(GeV/c)');xticks(0:6);
ylabel('1/2\pip_TBd^2N/dydp_T(GeV/c)^{-2}');yticks([10^-10,10^-8,10^-6,10^-4]);
legend('PHENIX J/\psi','sum','therml-thermal','thermal-shower','shower-shower')
ax.FontName='Times New Roman';
ax.FontSize=12;
ax.FontWeight='bold';

% hold on;
% load SS_Jpsi_c020.dat;
% semilogy(SS_Jpsi_c020(:,1),SS_Jpsi_c020(:,2));