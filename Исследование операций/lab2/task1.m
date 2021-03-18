function task1(q,p,N)

if (sum(q) ~= 1 || sum(p) ~= 1)
    error('Закон распределения задан неверно!')
end

mas = [1 1; 1 2; 2 1; 2 2];
res = zeros(N,1);

for i = 2:N
    res(i) = res(i-1) + sud(mas(myrand(q),:),mas(myrand(p),:));
end
grid on
hold on
if (min(res) ~= max(res))
    axis([1 N min(res) max(res)])
else
    axis([1 N -5 5])
end
comet(res)
end

function res = myrand(q)
    res = sum(rand >= cumsum(q)) + 1;
end