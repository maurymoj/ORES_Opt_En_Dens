load('Curvas_de_carga_brasil.mat')

Time_Winter_fr.Format = 'HH:mm:ss';
Time_Winter_sun.Format = 'HH:mm:ss';
Time_Summer_fr.Format = 'HH:mm:ss';
Time_Summer_sun.Format = 'HH:mm:ss';

figure('color',[1 1 1]);  % Cost per unit energy of CAES
plot(Time_Summer_fr,Load_Summer_fr./1000,'k-')
hold all
plot(Time_Summer_sun,Load_Summer_sun./1000,'-','Color',[0.5 0.5 0.5])

plot(Time_Winter_fr,Load_Winter_fr./1000,'k--')

plot(Time_Winter_sun,Load_Winter_sun./1000,'--','Color',[0.5 0.5 0.5])
ytickformat('%,1.0f')
datetick('x','HH:MM','keeplimits','keepticks')

%xtickformat('HH:mm:ss')
grid on
xlabel('Time')
ylabel('Hourly Load [GWh/h]')
legend('Summer weekday','Summer weekend','Winter weekday','Winter weekend','location','SouthEast')

fig = gcf;
fig.PaperOrientation = 'landscape';
fig.PaperSize = [6 4.5];