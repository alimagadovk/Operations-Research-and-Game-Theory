clc
clear
close all
N = 90;
v = 140000;
S = 4;
k = 700;
price = 40;
warehouse_vol = 20000;
vv = v + 40000*randn(1,N);
for i = 1:N
    if (vv(i) < 0)
        vv(i) = -vv(i);
    end
    mean_vv(i) = mean(vv(1:i));
    std_vv(i) = sqrt(var(1:i));
end


warehouse_1 = v/30;
warehouse_2 = mean_vv(1)/30;
orders1 = zeros(1,N);
orders2 = zeros(1,N);
for i = 1:N
    sales1(i) = min(warehouse_1(i), vv(i)/30);
    sales2(i) = min(warehouse_2(i), vv(i)/30);
    warehouse_1(i + 1) = warehouse_1(i) - sales1(i);
    warehouse_2(i + 1) = warehouse_2(i) - sales2(i);
    if (warehouse_1(i + 1) < v/30)
        warehouse_1(i + 1) = min(warehouse_1(i + 1) + v/30,warehouse_vol);
        orders1(i) = 1;
    end
    if ((warehouse_2(i + 1) < (mean_vv(i) + 3*std_vv(i))/30))
        tau = warehouse_2(i) / (mean_vv(i) / 30);
        warehouse_2(i + 1) = warehouse_2(i + 1) + tau*(mean_vv(i) + 3*std_vv(i))/30;
        orders2(i) = 1;
    end
end
figure
hold on
grid on
plot(mean_vv)
title('Среднее значение спроса в каждый день')
xlabel('День')
ylabel('Среднее значение спроса')

mas1 = find(orders1 == 1) + 1;
mas2 = find(orders2 == 1) + 1;
figure
hold on
grid on
plot(warehouse_1, 'r')
plot(warehouse_2, 'b')
plot(mas1,warehouse_1(mas1), 'r*')
plot(mas2,warehouse_2(mas2), 'b*')
title('Состояния складов в течение 90 дней')
xlabel('День')
ylabel('Кол-во единиц товара')
legend('Без учёта статистики','С учётом статистики')

Income1 = cumsum(price * sales1);
Income2 = cumsum(price * sales2);
figure
hold on
grid on
plot(Income1, 'r')
plot(Income2, 'b')
title('Графики доходов')
xlabel('День')
ylabel('Кол-во рублей')
legend('Без учёта статистики','С учётом статистики')

Costs1 = cumsum([0,orders1].*(-k)) + cumsum(- S/30 .* warehouse_1);
Costs2 = cumsum([0,orders2].*(-k)) + cumsum(- S/30 .* warehouse_2);
figure
hold on
grid on
plot(Costs1, 'r')
plot(Costs2, 'b')
title('Графики расходов')
xlabel('День')
ylabel('Кол-во рублей')
legend('Без учёта статистики','С учётом статистики')

Itog1 = [0, Income1] + Costs1;
Itog2 = [0, Income2] + Costs2;

figure
hold on
grid on
plot(Itog1, 'r')
plot(Itog2, 'b')
title('Прибыль')
xlabel('День')
ylabel('Кол-во рублей')
legend('Без учёта статистики','С учётом статистики')
xlswrite('data1_1.xlsx',[warehouse_1; warehouse_2])
%%
clc
clear
close all
N = 90;
v = 140000;
S = 4;
k = 700;
price = 20;
warehouse_vol = 20000;
vv = v + 50000.*cos([1:N].*pi.*2./14);
% figure
% plot(vv)
for i = 1:N
    mean_vv(i) = mean(vv(1:i));
    std_vv(i) = sqrt(var(1:i));
end




warehouse_1 = v/30;
warehouse_2 = mean_vv(1)/30;
orders1 = zeros(1,N);
orders2 = zeros(1,N);
for i = 1:N
    sales1(i) = min(warehouse_1(i), vv(i)/30);
    sales2(i) = min(warehouse_2(i), vv(i)/30);
    warehouse_1(i + 1) = warehouse_1(i) - sales1(i);
    warehouse_2(i + 1) = warehouse_2(i) - sales2(i);
    if (warehouse_1(i + 1) < v/30)
        warehouse_1(i + 1) = min(warehouse_1(i + 1) + v/30,warehouse_vol);
        orders1(i) = 1;
    end
    if ((warehouse_2(i + 1) < (mean_vv(i) + 3*std_vv(i))/30))
        tau = warehouse_2(i) / (mean_vv(i) / 30);
        warehouse_2(i + 1) = warehouse_2(i + 1) + tau*(mean_vv(i) + 3*std_vv(i))/30;
        orders2(i) = 1;
    end
end
figure
hold on
grid on
plot(mean_vv)
title('Среднее значение спроса в каждый день')
xlabel('День')
ylabel('Среднее значение спроса')

mas1 = find(orders1 == 1) + 1;
mas2 = find(orders2 == 1) + 1;
figure
hold on
grid on
plot(warehouse_1, 'r')
plot(warehouse_2, 'b')
plot(mas1,warehouse_1(mas1), 'r*')
plot(mas2,warehouse_2(mas2), 'b*')
title('Состояния складов в течение 90 дней')
xlabel('День')
ylabel('Кол-во единиц товара')
legend('Без учёта статистики','С учётом статистики')

Income1 = cumsum(price * sales1);
Income2 = cumsum(price * sales2);
figure
hold on
grid on
plot(Income1, 'r')
plot(Income2, 'b')
title('Графики доходов')
xlabel('День')
ylabel('Кол-во рублей')
legend('Без учёта статистики','С учётом статистики')

Costs1 = cumsum([0,orders1].*(-k)) + cumsum(- S/30 .* warehouse_1);
Costs2 = cumsum([0,orders2].*(-k)) + cumsum(- S/30 .* warehouse_2);
figure
hold on
grid on
plot(Costs1, 'r')
plot(Costs2, 'b')
title('Графики расходов')
xlabel('День')
ylabel('Кол-во рублей')
legend('Без учёта статистики','С учётом статистики')

Itog1 = [0, Income1] + Costs1;
Itog2 = [0, Income2] + Costs2;

figure
hold on
grid on
plot(Itog1, 'r')
plot(Itog2, 'b')
title('Прибыль')
xlabel('День')
ylabel('Кол-во рублей')
legend('Без учёта статистики','С учётом статистики')
xlswrite('data2_1.xlsx',[warehouse_1; warehouse_2])
%%
clc
clear
close all
N = 90;
v = 140000;
S = 4;
k = 700;
price = 20;
warehouse_vol = 20000;
vv = v + 50000.*cos([1:N].*pi.*2./14) + 25000*randn(1,N);
% figure
% plot(vv)
for i = 1:N
    mean_vv(i) = mean(vv(1:i));
    std_vv(i) = sqrt(var(1:i));
end




warehouse_1 = v/30;
warehouse_2 = mean_vv(1)/30;
orders1 = zeros(1,N);
orders2 = zeros(1,N);
for i = 1:N
    sales1(i) = min(warehouse_1(i), vv(i)/30);
    sales2(i) = min(warehouse_2(i), vv(i)/30);
    warehouse_1(i + 1) = warehouse_1(i) - sales1(i);
    warehouse_2(i + 1) = warehouse_2(i) - sales2(i);
    if (warehouse_1(i + 1) < v/30)
        warehouse_1(i + 1) = min(warehouse_1(i + 1) + v/30,warehouse_vol);
        orders1(i) = 1;
    end
    if ((warehouse_2(i + 1) < (mean_vv(i) + 3*std_vv(i))/30))
        tau = warehouse_2(i) / (mean_vv(i) / 30);
        warehouse_2(i + 1) = warehouse_2(i + 1) + tau*(mean_vv(i) + 3*std_vv(i))/30;
        orders2(i) = 1;
    end
end
figure
hold on
grid on
plot(mean_vv)
title('Среднее значение спроса в каждый день')
xlabel('День')
ylabel('Среднее значение спроса')

mas1 = find(orders1 == 1) + 1;
mas2 = find(orders2 == 1) + 1;
figure
hold on
grid on
plot(warehouse_1, 'r')
plot(warehouse_2, 'b')
plot(mas1,warehouse_1(mas1), 'r*')
plot(mas2,warehouse_2(mas2), 'b*')
title('Состояния складов в течение 90 дней')
xlabel('День')
ylabel('Кол-во единиц товара')
legend('Без учёта статистики','С учётом статистики')

Income1 = cumsum(price * sales1);
Income2 = cumsum(price * sales2);
figure
hold on
grid on
plot(Income1, 'r')
plot(Income2, 'b')
title('Графики доходов')
xlabel('День')
ylabel('Кол-во рублей')
legend('Без учёта статистики','С учётом статистики')

Costs1 = cumsum([0,orders1].*(-k)) + cumsum(- S/30 .* warehouse_1);
Costs2 = cumsum([0,orders2].*(-k)) + cumsum(- S/30 .* warehouse_2);
figure
hold on
grid on
plot(Costs1, 'r')
plot(Costs2, 'b')
title('Графики расходов')
xlabel('День')
ylabel('Кол-во рублей')
legend('Без учёта статистики','С учётом статистики')

Itog1 = [0, Income1] + Costs1;
Itog2 = [0, Income2] + Costs2;

figure
hold on
grid on
plot(Itog1, 'r')
plot(Itog2, 'b')
title('Прибыль')
xlabel('День')
ylabel('Кол-во рублей')
legend('Без учёта статистики','С учётом статистики')
xlswrite('data3_1.xlsx',[warehouse_1; warehouse_2])