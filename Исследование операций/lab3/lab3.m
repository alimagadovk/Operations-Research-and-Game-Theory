clc
clear
close all

n=5;
C=100;

A = zeros(n,n);
for i=1:n
    for j=1:n
        if (i < j)
            A(i,j)=C*(j-i);
        end
        if (i > j)
            A(i,j)=C*(n-i+1);
        end
        if (i == j)
            A(i,j)=C*(n-j+1)/2;
        end
    end
end

A

f = ones(1,5);
b = ones(5,1);
D = zeros(5,5);
[y,fmin] = linprog(-f,[A;D-diag(ones(5,1))],[b;zeros(5,1)]);
v = -1/fmin
q = v*y
%
A1 = -A';
f = ones(1,5);
b = -ones(5,1);
D = zeros(5,5);
[y,fmin] = linprog(f,[A1;D-diag(ones(5,1))],[b;zeros(5,1)]);
v = 1/fmin
p = v*y
%%
%q = [3/7 0 2/7 2/7 0];
%p = [0 2/7 4/7 1/7 0];

N = 1000;

for i = 1:N
    j = sum(rand >= cumsum(p)) + 1;
    k = sum(rand >= cumsum(q)) + 1;
    res1(i) = A(j,k);
    if (j == k)
        res2(i) = res1(i);
    else
        res2(i) = 500 - res1(i);
    end
end
S1 = cumsum(res1);
S2 = cumsum(res2);

m1 = mean(res1)
m2 = mean(res2)

M = 0;
for i = 1:n
    for j = 1:n
        M = M + A(i,j)*p(i)*q(j);
    end
end
M
figure(1)
plot(res1)
grid on
axis([0 N 100 500])
figure(2)
plot(res2)
grid on
axis([0 N 100 500])
figure(3)
hold on
grid on
axis([0 N 0 max([S1, S2])])
plot(S1,'.b')
plot(S2,'.r')
% for i = 1:N
%     figure(3)
%     plot(i,S1(i),'.b')
%     plot(i,S2(i),'.r')
%     pause(0.001)
% end
%%
clc
clear
close all

n=5;
C=100;

A = zeros(n,n);
for i=1:n
    for j=1:n
        if (i < j)
            A(i,j)=C*(j-i);
        end
        if (i > j)
            A(i,j)=C*(n-i+1);
        end
        if (i == j)
            A(i,j)=C*(n-j+1)/2;
        end
    end
end

A = A.*0.1 + 5

f = ones(1,5);
b = ones(5,1);
D = zeros(5,5);
[y,fmin] = linprog(-f,[A;D-diag(ones(5,1))],[b;zeros(5,1)]);
v = -1/fmin
q = v*y