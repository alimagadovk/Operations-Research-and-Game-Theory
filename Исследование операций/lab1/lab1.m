clc
clear
A = [1 2];
B = [1 2];
[AV, BV] = sud(A,B,0,0)
A = [1 2];
B = [2 2];
[AV, BV] = sud(A,B,0,0)
A = [[1 1]
     [1 2]
     [2 1]
     [2 2]];
B = A;
matr = zeros(4);
for i = 1:size(matr,1)
    for j = 1:size(matr,2)
        [AV, BV] = sud(A(i,:),B(j,:),0,0);
        matr(i,j) = matr(i,j) + AV - BV;
    end
end
matr
%%
f = [-1 -1]';
A = [-2 1; 2 1; 0 -1; -1 0];
b = [2 -4 0 0];
[x, y] = linprog(f,A,b);
x
y
%%
f = [-2 -1]';
A = [0 1; 1 1];
b = [3; 5];
Aeq = [1 0];
beq = [3];
[x, y] = linprog(f,A,b,Aeq,beq);
x
y
%%
f = [-2 -1]';
A = [0 1; 0 -1];
b = [1; -1];
[x, y] = linprog(f,A,b);
x
y