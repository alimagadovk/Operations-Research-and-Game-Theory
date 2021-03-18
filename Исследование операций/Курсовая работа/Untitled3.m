clc
clear
close all
v = 140000;
S = 4;
k = 700;
price = 20;
warehouse_vol = 10000;
vv = xlsread('D_Z_2019_all_MIET.xls',['Продажи']);
vv = vv(243:494,6);
N = length(vv);
% figure
% plot(vv)
delay = 14;
% for i = 1:N
%     mean_vv(i) = mean(vv(max([1,i - delay]):i));
%     std_vv(i) = sqrt(var(vv(max([1,i - delay]):i)));
%     %mean_vv(i) = mean(vv(1:i));
%     %std_vv(i) = sqrt(var(vv(1:i)));
% end




%warehouse_1 = 500;
warehouse = 100;
%orders1 = zeros(1,N);
orders = zeros(1,N);
mean_k(1) = 0;
order = 0;
j = 1;
for i = 1:N
    
    [mean_vv, std_vv] = mean_std(i, vv, delay);
    sales(i) = min(warehouse(i), vv(i));
    warehouse(i + 1) = warehouse(i) - sales(i) + order;

    if ((warehouse(i + 1) < (mean_vv + 3*std_vv)))
        tau = warehouse(j) / mean_vv;
        i
        j
        %tau = max([warehouse(j) - warehouse(i + 1) / mean_vv, 0]);
        T(i) = tau;
%         warehouse_2(i + 1) = min(warehouse_2(i + 1) + tau*(mean_vv + 3*std_vv), warehouse_vol);
        order = min(tau*(mean_vv + 3*std_vv), warehouse_vol - warehouse(i + 1));
        orders(i) = 1;
        j = i;
    else
        order = 0;
    end
end
figure
hold on
grid on
plot(T)

figure
hold on
grid on
plot(vv)
title('Значение спроса в каждый день')
xlabel('День')
ylabel('Спрос')

mas = find(orders == 1) + 1;
figure
hold on
grid on
plot(warehouse, 'b')
plot(mas,warehouse(mas), 'b*')
title('Состояния складов в течение ' + string(N) + ' дней')
xlabel('День')
ylabel('Кол-во единиц товара')
%legend('Без учёта статистики')

Income = cumsum(price * sales);
figure
hold on
grid on
plot(Income, 'b')
title('Графики доходов')
xlabel('День')
ylabel('Кол-во рублей')
%legend('Без учёта статистики')

Costs = cumsum([0,orders].*(-k)) + cumsum(- S/30 .* warehouse);
figure
hold on
grid on
plot(Costs, 'b')
title('Графики расходов')
xlabel('День')
ylabel('Кол-во рублей')
%legend('Без учёта статистики')

Itog = [0, Income] + Costs;

figure
hold on
grid on
plot(Itog, 'b')
title('Прибыль')
xlabel('День')
ylabel('Кол-во рублей')
legend('Без учёта статистики')