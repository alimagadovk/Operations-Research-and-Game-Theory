clc
clear
close all

tt = 1000*(1.2 - 0.8) + 6000*(0.8 - 0.5) - 200;
th = 1000*(1.2 - 0.8) - 0.5 * 6000 + 0.8 * 1200 + 0.2 * 4800 - 200;
ht = -0.8 * 4000 + 1000 * 1.2 + 3000 * 0.3 + 1200*(0.8 - 0.5) - 200;
hh = 4000*(1.2 - 0.8) + 1200*(0.8 - 0.5) - 200;

A = [tt, th;
    ht, hh]

A1 = (A + 940) / 10
%

f = ones(1,2);
b = ones(2,1);
D = zeros(2,2);
[A1;D-diag(ones(2,1))]
[b;zeros(2,1)]
[y,fmin] = linprog(-f,[A1;D-diag(ones(2,1))],[b;zeros(2,1)]);
v = -1/fmin
q = v*y
%
A2 = -A1';

f = ones(1,2);
b = -ones(2,1);
D = zeros(2,2);
[y,fmin] = linprog(f,[A2;D-diag(ones(2,1))],[b;zeros(2,1)]);
v = 1/fmin
p = v*y

q2 = [0.4 0.6];
p2 = [2/3 1/3];

N = 1000;

for i = 1:N
    j = sum(rand >= cumsum(p)) + 1;
    k = sum(rand >= cumsum(q)) + 1;
    res(i) = A(j,k);
    j2 = sum(rand >= cumsum(p2)) + 1;
    k2 = sum(rand >= cumsum(q2)) + 1;
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

A = [2, 5, 6;
     5, 3, 1];
A1 = -A
 
f = ones(1,3);
b = -ones(2,1);
D = zeros(3,3);
[y,fmin] = linprog(f,[A1;D-diag(ones(3,1))],[b;zeros(3,1)]);
v = 1/fmin
q = v*y


N = 1000;

q2 = [0.625 0 0.375];
p = [1/2 1/2];

for i = 1:N
    j = sum(rand >= cumsum(p)) + 1;
    k = sum(rand >= cumsum(q)) + 1;
    res(i) = A(j,k);
    j2 = sum(rand >= cumsum(p)) + 1;
    k2 = sum(rand >= cumsum(q2)) + 1;
    res2(i) = A(j2,k2);
end
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

A = [[180 80 -20 -120]
     [60 360 260 160]
     [-60 240 540 440]
     [-180 120 420 720]];
 
 mas_max = max(A);
 for i = 1:size(A,1)
     for j = 1:size(A,2)
         R(i,j) = mas_max(j) - A(i,j);
     end
 end
 R
 
%  A1 = -(A' + 180) / 20;
% %
% f = ones(1,4);
% b = -ones(1,4);
% D = zeros(4,4);
% [A1;D-diag(ones(4,1))]
% [b,zeros(1,4)]
% [y,fmin] = linprog(f,[A1;D-diag(ones(4,1))],[b,zeros(1,4)]);
% v = 1/fmin
% q = v*y
%
%
R1 = R / 40;

f = ones(1,4);
b = ones(1,4);
D = zeros(4,4);
[R;D-diag(ones(4,1))]
[b,zeros(1,4)]
[y,fmin] = linprog(-f,[R1;D-diag(ones(4,1))],[b,zeros(1,4)]);
v = -1/fmin
q = v*y

for i = 1:size(R,1)
    a(i) = A(i,:)*q;
    r(i) = R(i,:)*q;
end
max_a = max(a)
min_r = min(r)

q = [0.2; 0.35; 0.25; 0.2];

for i = 1:size(R,1)
    a(i) = A(i,:)*q;
    r(i) = R(i,:)*q;
end
max_a = max(a)
min_r = min(r)