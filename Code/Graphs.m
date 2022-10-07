%% ib crit val
clear; clc;
T = readtable('crit_ib_table.xlsx');

for n = 1:length(T.observation_date)
    T.observation_date(n) = T.observation_date(n) + 693960;
    %T.observation_date(n) = datestr(T.observation_date(n))
end

T.date= datestr(T.observation_date)

close;
plot(datenum(T.date), T.DWRate,'k', datenum(T.date), T.crit_1,'-.',datenum(T.date), T.crit_2,'-.',datenum(T.date), T.crit_10,'-.',datenum(T.date), T.crit_30,'-.');
dateFormat = 10; datetick('x',dateFormat); xlabel('Time'); ylabel('Interest Rate');
legend('DW Rate','1Y T-bill','2Y T-bill','10Y T-bill','30Y T-bill');

print -djpeg fig_ibstar


%% int rate gap

clear; clc;
T = readtable('dw_effr_spread.xls');
T.date= datestr(T.observation_date)

close;
plot(datenum(T.observation_date), T.Spread);
dateFormat = 10; datetick('x',dateFormat); xlabel('Time'); ylabel('Interest Rate');
legend('DW - EFFR')

print -djpeg int_spread
