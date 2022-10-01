clc;clear;close all
mat=readtable('dNptdpt_D0_c020_276_1st.dat');
mat=table2array(mat);
%%
pt_Exper=[1.5 2.5 3.5 4.5 5.5 7 10 14 20];
y_Exper=[1.766666667	0.412	0.074857143	0.016533333	0.004109091	0.000922857	0.000152	1.84286E-05	0.00000358 ]/2/pi;
y_error=[0.258	0.02728	0.004628571	0.001282222	0.0004	7.94286E-05	0.0000125	2.71429E-06	0.00000092]/2/pi;
pt=double(mat(:,1));
TT=mat(:,2);
TS=mat(:,3);
SS=mat(:,4);
SS2j=mat(:,5);
Total=TT+TS+SS+SS2j;
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
xlim([0,20])
 ylim([10^-10,10^1])
xlabel('p_T(GeV/c)');xticks(0:20);
ylabel('d^2N/p_Tdp_Tdy(GeV/c)^{-2}');yticks([10^-10,10^-8,10^-6,10^-4,10^-3,10^-2,10^-1,10^0,10^1]);
legend('ALICE D^0@2.76 TeV','sum','therml-thermal','thermal-shower','shower-shower','shower-shower(2j)')
ax.FontName='Times New Roman';
ax.FontSize=12;
ax.FontWeight='bold';