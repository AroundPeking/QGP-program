clc;clear;close all
mat=readtable('dNptdpt_Jpsi_c020_55_final.dat');
mat=table2array(mat);
%%
pt=double(mat(:,1));
TT=mat(:,2);
TS=mat(:,3)*10^4;
SS=mat(:,4);
SS2j=mat(:,5);
Total=TT+TS+SS+SS2j;
%%
figure;
ax=axes;
%=====
%plt_1=semilogy(pt,Total,'g-');
%plt_1.LineWidth=1.2;
%hold on
%=====
%plt_2=semilogy(pt,TT,'k-.');
%plt_2.LineWidth=1;
%hold on
%=====
plt_3=semilogy(pt,TS,'b--');
plt_3.LineWidth=1;
hold on
%=====
plt_4=semilogy(pt,SS,'k-');
plt_4.LineWidth=1;
hold on
%=====
plt_5=semilogy(pt,SS2j,'r-');
plt_5.LineWidth=1;
hold on
%=====
xlim([10,20])
ylim([10^-15,10^-4])
xlabel('p_T(GeV/c)');xticks(10:20);
ylabel('d^2N/p_Tdp_Tdy(GeV/c)^{-2}');yticks([10^-15,10^-10,10^-8,10^-6,10^-5,10^-4,10^-3,10^-2,10^-1]);
legend('thermal-shower*10^4','shower-shower','shower-shower(2j)')
title('ALICE J/\psi@5.5 TeV')
ax.FontName='Times New Roman';
ax.FontSize=12;
ax.FontWeight='bold';