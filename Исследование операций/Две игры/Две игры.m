clc
clear
close all
A = [20, 30;
     10, 25;
     35, 15]
 
 A1 = -A'
 
f = ones(1,3);
b = -ones(2,1);
D = zeros(3,3);
[A1;-diag(ones(3,1))]
[b;zeros(3,1)]
%
[y,fmin] = linprog(f,[A1;-diag(ones(3,1))],[b;zeros(3,1)]);
v = 1/fmin
p = v*y

q = [0.5 0.5]
p2 = [1/3 1/3 1/3]

N = 10000;

for i = 1:N
    j = sum(rand >= cumsum(p)) + 1;
    k = sum(rand >= cumsum(q)) + 1;
    res(i) = A(j,k);
    j2 = sum(rand >= cumsum(p2)) + 1;
    k2 = sum(rand >= cumsum(q)) + 1;
    res2(i) = A(j2,k2);
end
figure
S1 = cumsum(res);
plot(S1)
hold on
grid on
S2 = cumsum(res2);
plot(S2)
%%
clc
clear
close all
A = [1,-2;
    -0.3, 0.7]
 A1 = (A + 2) * 10
 
f = ones(1,2);
b = ones(2,1);
D = zeros(2,2);
[A1;-diag(ones(2,1))]
[b;zeros(2,1)]
[qy,fmin] = linprog(-f,[A1;D-diag(ones(2,1))],[b;zeros(2,1)]);
v = -1/fmin
q = v*qy
%
A1 = A + 2
A1 = -A1'
 
f = ones(1,2);
b = -ones(2,1);
D = zeros(2,2);
[A1;-diag(ones(2,1))]
[b;zeros(2,1)]
[py,fmin] = linprog(f,[A1;D-diag(ones(2,1))],[b;zeros(2,1)]);
v = 1/fmin
p = v*py


p2 = [0.5 0.5]
%q2 = [0.5 0.5]

N = 10000;

for i = 1:N
    j = sum(rand >= cumsum(p)) + 1;
    k = sum(rand >= cumsum(q)) + 1;
    res(i) = A(j,k);
    j2 = sum(rand >= cumsum(p2)) + 1;
    k2 = sum(rand >= cumsum(q)) + 1;
    res2(i) = A(j2,k2);
end
figure
S1 = cumsum(res);
plot(S1)
hold on
grid on
S2 = cumsum(res2);
plot(S2)
xlabel('Кол-во игр')
ylabel('Кол-во очков допрашивающиего')
title('Игра "Допрос"')
legend('Оптимальные смежные стратегии','Неоптимальные стратегии')

x = [0 0.05];

figure
hold on
grid on
plot(x.*0 + 1/30, x,'b')
plot(x, (1 - 17.*x)./27,'b')

plot(x, x.*0, 'k')
plot(x.*0, x, 'k')

plot(x, 0.01 - x, 'g')
plot(qy(1), qy(2), 'r*')
axis([0 0.04 0 0.04])
xlabel('x')
ylabel('y(x)')
%axis equal