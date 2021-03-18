function res = f1(t,Q,v,tau1)
if (t <= tau1)
    res = Q - v.*t;
else
    res = 0;
end
end