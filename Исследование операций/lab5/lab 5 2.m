clc
clear
close all
k = 700;
v = 140000;
S = 4;
C = @(q) k.*v./q + S.*q./2;
vq = 1000:1:20000;
plot(vq,C(vq))
Q = sqrt(2.*k.*v./S);
hold on
grid on
plot(Q,C(Q),'r*')
title('Функция расходов C(q)')
xlabel('q')
ylabel('C(q)')
tau1 = Q / v
C1 = C(Q)
%
tau2 = 0.1;
Q2 = tau2 * v;
C2 = C(Q2);
S = C2 - C1
plot(Q2, C2, 'b*')
legend('C(q)','Оптимальный размер заказа q', 'Текущий размер заказа')
%
figure
hold on
grid on
vy = [];
vx = [];
for i = 1:5
    f = @(t) f1(t,Q,v,tau1*i);
    vt = (i - 1)*tau1:0.01:i*tau1;
    for j = 1:length(vt)
        y(j) = f(vt(j) - tau1*(i - 1));
    end
    vx = [vx, vt];
    vy = [vy, y];
end
plot(vx,vy)
title('Состояние «запаса (склада)» при оптимальных параметрах')
xlabel('t')
ylabel('q(t)')
axis([0 max(vx) 0 max(vy)])
figure
hold on
grid on
vy = [];
vx = [];
for i = 1:5
    f = @(t) f1(t,Q2,v,tau2*i);
    vt = (i - 1)*tau2:0.01:i*tau2;
    for j = 1:length(vt)
        y(j) = f(vt(j) - tau2*(i - 1));
    end
    vx = [vx, vt];
    vy = [vy, y];
end
plot(vx,vy)
title('Состояние «запаса (склада)» при неоптимальных параметрах')
xlabel('t')
ylabel('q(t)')
axis([0 max(vx) 0 max(vy)])
%%
clc
clear
close all
k = 700;
v = 140000;
S = 4;
lamb = 900000;
C = @(q) k.*v./q + q.*S.*(1 - v/lamb)./2;
vq = 1000:1:70000;
plot(vq,C(vq))
Q = sqrt(2*k*v/(S*(1 - v/lamb)))
hold on
grid on
plot(Q,C(Q),'r*')
title('Функция расходов C(q)')
xlabel('q')
ylabel('C(q)')
legend('C(q)','Оптимальный размер заказа q')
tau1 = Q / lamb;
Qm = (lamb - v)*tau1;
tau2 = Qm / v;
figure
t = [0 tau1];
plot(t, (lamb - v).*t)
hold on
grid on
t = [tau1 tau1 + tau2];
plot(t, Qm - v.*(t - tau1))

figure
hold on
grid on
vy = [];
vx = [];
for i = 1:5
    vt = [0, tau1];
    y = (lamb - v).*vt;
    vx = [vx, vt + (i - 1)*(tau1 + tau2)];
    vy = [vy, y];
    vt = [tau1, tau1 + tau2];
    y = Qm - v.*(vt - tau1);
    vx = [vx, vt  + (i - 1)*(tau1 + tau2)];
    vy = [vy, y];
end
plot(vx,vy)
title('Состояние «запаса (склада)» при оптимальных параметрах')
xlabel('t')
ylabel('q(t)')
axis([0 max(vx) 0 max(vy)])
%%
clc
clear
close all
k = 21;
v = 300;
S = 14;
C = @(q) k.*v./q + S.*q./2;
vq = 10:0.001:80;
plot(vq,C(vq))
Q = sqrt(2.*k.*v./S);
hold on
grid on
plot(Q,C(Q),'r*')
title('Функция расходов C(q)')
xlabel('q')
ylabel('C(q)')
legend('C(q)','Оптимальный размер заказа q')
tau1 = Q / v
C1 = C(Q)
figure
hold on
grid on
vy = [];
vx = [];
for i = 1:5
    f = @(t) f1(t,Q,v,tau1*i);
    vt = (i - 1)*tau1:0.01:i*tau1;
    for j = 1:length(vt)
        y(j) = f(vt(j) - tau1*(i - 1));
    end
    vx = [vx, vt];
    vy = [vy, y];
end
plot(vx,vy)
title('Состояние «запаса (склада)» при оптимальных параметрах')
xlabel('t')
ylabel('q(t)')
axis([0 max(vx) 0 max(vy)])
%% 4
clc
clear
close all
k = 700;
v = 140000;
S = 4;
C = @(q) k.*v./q + S.*q./2;
vq = 1000:1:20000;
plot(vq,C(vq))
Q = sqrt(2.*k.*v./S);
hold on
grid on
plot(Q,C(Q),'r*')
tau1 = Q / v
C1 = C(Q)
tau2 = 0.1;
Q2 = tau2 * v;
C2 = C(Q2);
S = C2 - C1
plot(Q2, C2, 'b*')
title('Функция расходов C(q)')
xlabel('q')
ylabel('C(q)')
tau0 = 1/30;
N = 4;
reserv = 0:N - 1;
reserv = reserv .* tau1;
figure
hold on
grid on
vy = [];
vx = [];
for i = 1:N
    f = @(t) f1(t,Q,v,tau1*i);
    vt = (i - 1)*tau1:0.01:i*tau1;
    for j = 1:length(vt)
        y(j) = f(vt(j) - tau1*(i - 1));
    end
    vx = [vx, vt];
    vy = [vy, y];
end
plot(vx + tau0,vy)
plot(reserv,reserv.*0, 'r*')
title('Состояние «запаса (склада)» при оптимальных параметрах')
xlabel('t')
ylabel('q(t)')
legend('q(t)','Заявки на заказ')
axis([0 max(vx + tau0) 0 max(vy)])
%%
clc
clear
close all
k = 700;
v = 140000;
S = 4;
lamb = 900000;
C = @(q) k.*v./q + q.*S.*(1 - v/lamb)./2;
vq = 1000:1:70000;
plot(vq,C(vq))
Q = sqrt(2*k*v/(S*(1 - v/lamb)))
hold on
grid on
plot(Q,C(Q),'r*')
title('Функция расходов C(q)')
xlabel('q')
ylabel('C(q)')
tau1 = Q / lamb;
Qm = (lamb - v)*tau1;
tau2 = Qm / v;
figure
t = [0 tau1];
plot(t, (lamb - v).*t)
hold on
grid on
t = [tau1 tau1 + tau2];
plot(t, Qm - v.*(t - tau1))

tau0 = 1/30;
N = 5;
reserv = 0:N - 1;
reserv = reserv .* (tau1 + tau2);
figure
hold on
grid on
vy = [];
vx = [];
for i = 1:N
    vt = [0, tau1];
    y = (lamb - v).*vt;
    vx = [vx, vt + (i - 1)*(tau1 + tau2)];
    vy = [vy, y];
    vt = [tau1, tau1 + tau2];
    y = Qm - v.*(vt - tau1);
    vx = [vx, vt  + (i - 1)*(tau1 + tau2)];
    vy = [vy, y];
end
plot(vx + tau0,vy)
plot(reserv,reserv.*0, 'r*')
title('Состояние «запаса (склада)» при оптимальных параметрах')
xlabel('t')
ylabel('q(t)')
legend('q(t)','Заявки на заказ')
axis([0 max(vx + tau0) 0 max(vy)])
%%
clc
clear
close all
k = 21;
v = 300;
S = 14;
C = @(q) k.*v./q + S.*q./2;
vq = 10:0.001:80;
plot(vq,C(vq))
Q = sqrt(2.*k.*v./S);
hold on
grid on
plot(Q,C(Q),'r*')
title('Функция расходов C(q)')
xlabel('q')
ylabel('C(q)')
tau1 = Q / v
C1 = C(Q)

tau0 = 1/12;
N = 4;
reserv = 0:N - 1;
reserv = reserv .* tau1;
figure
hold on
grid on
vy = [];
vx = [];
for i = 1:N
    f = @(t) f1(t,Q,v,tau1*i);
    vt = (i - 1)*tau1:0.01:i*tau1;
    for j = 1:length(vt)
        y(j) = f(vt(j) - tau1*(i - 1));
    end
    vx = [vx, vt];
    vy = [vy, y];
end
plot(vx + tau0,vy)
plot(reserv,reserv.*0, 'r*')
title('Состояние «запаса (склада)» при оптимальных параметрах')
xlabel('t')
ylabel('q(t)')
legend('q(t)','Заявки на заказ')
axis([0 max(vx + tau0) 0 max(vy)])