clc
clear
close all
q = [0 3/5 2/5 0] ;
p = [0 3/5 2/5 0];
N = 10000;

mas = [1 1; 1 2; 2 1; 2 2];
res = zeros(N,1);

for i = 2:N
    j = sum(rand >= cumsum(q)) + 1;
    k = sum(rand >= cumsum(p)) + 1;
    res(i) = res(i-1) + sud(mas(j,:),mas(k,:));
end
grid on
hold on
if (min(res) ~= max(res))
    axis([1 N min(res) max(res)])
else
    axis([1 N -5 5])
end
comet(res)
res(end)
%%
clc
clear
close all

n=5;
C=100;

matrix = zeros(n,n);
for i=1:n
    for j=1:n
        if (i < j)
            matrix(i,j)=C*(j-i);
        end
        if (i > j)
            matrix(i,j)=C*(n-i+1);
        end
        if (i == j)
            matrix(i,j)=C*(n-j+1)/2;
        end
    end
end

q = [3/7 0 2/7 2/7 0];
p = [0 2/7 4/7 1/7 0];

N = 1000;

for i = 1:N
    j = sum(rand >= cumsum(q)) + 1;
    k = sum(rand >= cumsum(p)) + 1;
    res1(i) = matrix(j,k);
    if (j == k)
        res2(i) = res1(i);
    else
        res2(i) = 500 - res1(i);
    end
end
figure
plot(res1)
grid on
figure
plot(res2)
grid on
m1 = mean(res1)
m2 = mean(res2)



M = 0;
for i = 1:n
    for j = 1:n
        M = M + matrix(i,j)*p(i)*q(j);
    end
end
M
%%
clc
clear
close all

c=7.5;
S=[8.5, 9, 9.5];
N=[6,7,8];
matrix = zeros(3,3);
for i=1:3
    for j=1:3
        if (i<j)
            matrix (i,j) = -c *N(i)-S(j)*(j-i);
        else
            matrix (i,j) = -c *N(i);
        end
    end
end

q = [1/3 1/3 1/3];
p = [1/3 1/3 1/3];

N = 10000;
for i = 1:N
    j = sum(rand >= cumsum(q)) + 1;
    k = sum(rand >= cumsum(p)) + 1;
    res(i) = matrix(j,k);
end

plot(res)
grid on
mean(res)
%%
x = [0 5];

figure
hold on
grid on
plot(x, 3.5 - x./2,'b')
plot(x, x.*0 + 3,'b')


plot(x, 8 - x.*2,'b')

plot(x, x.*0, 'k')
plot(x.*0, x, 'k')

plot(x, 6 - 1.5.*x, 'g')
axis([0 5 0 5])
axis equal


f = [3 2]';
A = [1 2; 2 1; 0 1; -1 0; 0 -1];
b = [7 8 3 0 0]';
[x, fv] = linprog(-f,A,b)
%%
close all
x = [0 5];

figure
hold on
grid on
plot(x, 3.5 - x./2,'b')
plot(x, x.*0 + 3,'b')
plot(x, 8 - x.*2,'b')
plot(x, 1 - x./2,'b')

plot(x, x.*0, 'k')
plot(x.*0, x, 'k')

plot(x, 0.9 - 0.5.*x, 'g')
axis([0 5 0 5])
axis equal


f = [-3 -2]';
A = [1 2; 2 1; 0 1; 1 2; -1 0; 0 -1];
b = [7 8 3 2 0 0]';
[x, fv] = linprog(-f,A,b)
%%
close all
x = [0 2];

figure
hold on
grid on
plot(x, 1 - x,'b')


plot(x, x.*0, 'k')
plot(x.*0, x, 'k')

plot(x, 0.9 - 0.5.*x, 'g')
axis([0 2 0 2])
axis equal


f = [3 2]';
A = [1 1; -1 0; 0 -1];
b = [1 0 0]';
[x, fv] = linprog(-f,A,b)
%%
close all
x = [0 5];

figure
hold on
grid on
plot(x, 4.25 - x./2,'b')
plot(x, x.*0 + 3,'b')


plot(x, 8 - x.*2,'b')

plot(x, x.*0, 'k')
plot(x.*0, x, 'k')

plot(x, 6 - 1.5.*x, 'g')
axis([0 5 0 5])
axis equal


f = [3 2]';
A = [1 2; 2 1; 0 1];
b = [8.5 8 3]';
[x, fv] = linprog(-f,[],[],A,b)
%%
close all
x = [0 5];

figure
hold on
grid on
plot(x, 3.5 - x./2,'b')
plot(x, x.*0 + 3,'b')
plot(x, 8 - x.*2,'b')
plot(x, 1 - x./2,'b')

plot(x, x.*0, 'k')
plot(x.*0, x, 'k')

plot(x, 0.9 - 0.5.*x, 'g')
axis([0 5 0 5])
axis equal


f = [-3 -2]';
A = [1 2; 2 1; 0 1; -1 0; 0 -1];
b = [7 8 3 0 0]';
Aeq = [1 2];
beq = 2;
[x, fv] = linprog(-f,A,b,Aeq,beq)