clc
clear
close all
S = 4;
k = 700;
price = 20;
warehouse_vol = 1000;
vv = xlsread('D_Z_2019_all_MIET.xls',['Продажи']);
vv = vv(243:494,6);
N = length(vv);
depth = 20;





warehouse = 100;
orders = zeros(1,N);
order = 0;
vw(1) = warehouse;
flag = 0;
warehouse_l = warehouse;
for i = 1:N
    mean_vv = mean(vv(max([1,i - depth]):i));
    std_vv = sqrt(var(vv(max([1,i - depth]):i)));
    sales(i) = min(warehouse, vv(i));
    if (flag == 1)
        warehouse_l = warehouse - sales(i) + order;
        flag = 0;
    end
    warehouse = warehouse - sales(i) + order;
    if (warehouse < mean_vv + 3*std_vv)
        tau = (warehouse_l - warehouse) / mean_vv;
        order = min(floor(tau*(mean_vv + 3*std_vv)), warehouse_vol - warehouse);
        orders(i) = 1;
        flag = 1;
    else
        order = 0;
    end
    orders_size(i) = order;
    vw(i + 1) = warehouse;
end



figure
hold on
grid on
plot(vv, '*')
title('Значение спроса в каждый день')
xlabel('День')
ylabel('Спрос')

mas = find(orders == 1) + 1;
figure
hold on
grid on
plot(vw, 'b')
plot(mas,vw(mas), 'b*')
title('Состояния складов в течение ' + string(N) + ' дней')
xlabel('День')
ylabel('Кол-во единиц товара')


Income = cumsum(price * sales);
figure
hold on
grid on
plot(Income, 'b')
title('График доходов')
xlabel('День')
ylabel('Кол-во рублей')


time = 1:30:N;
vS = zeros(1,N);
for i = 2:length(time)
    vS(time(i)) = S .* max(vw(time(i - 1):time(i)));
end
Costs = cumsum([0,orders].*(-k)) - cumsum([0, vS]);
figure
hold on
grid on
plot(Costs, 'b')
title('График расходов')
xlabel('День')
ylabel('Кол-во рублей')

Itog = [0, Income] + Costs;

figure
hold on
grid on
plot(Itog, 'b')
title('Прибыль')
xlabel('День')
ylabel('Кол-во рублей')
%%
xlswrite('result.xls',[vv,vw(2:end)',orders',orders_size',Itog(2:end)'])