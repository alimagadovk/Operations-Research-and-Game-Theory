clc
clear
close all
k = 700;
v = 140000;
S = 4;
C = @(q) k.*v./q + S.*q.^2./(2.*v);
vq = 5000:1:20000;
plot(vq,C(vq))
Q = (k*(v^2)/S)^(1/3);
hold on
grid on
plot(Q,C(Q),'r*')
tau1 = Q / v
C1 = C(Q)
%
tau2 = 0.1;
Q2 = tau2 * v;
C2 = C(Q2);
S = C2 - C1
plot(Q2, C2, 'b*')
%%
clc
clear
close all
k = 700;
v = 140000;
S = 4;
lamb = 900000;
C = @(q) k.*v./q + q.*S.*(lamb - v)./(2.*lamb);
vq = 1000:1:70000;
plot(vq,C(vq))
Q = sqrt(2*k*v/(S*(1 - v/lamb)))
hold on
grid on
plot(Q,C(Q),'r*')
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
%%
clc
clear
close all
k = 21;
v = 300;
S = 14;
C = @(q) k.*v./q + S.*q.^2./(2.*v);
vq = 10:0.001:80;
plot(vq,C(vq))
Q = (k*(v^2)/S)^(1/3);
hold on
grid on
plot(Q,C(Q),'r*')
tau1 = Q / v
C1 = C(Q)
%%
clc
clear
close all
lamb = 4;
v = 2;
tau0 = 2;
tau1 = 10;
tau2 = 15;
f = @(x) treug1(x,v,tau0,tau1,tau2) + treug2((x - 15),v,tau0,tau1,tau2);
vx = 0:0.01:70;
for i = 1:length(vx)
    vf(i) = f(vx(i));
end
figure
plot(vx,vf)
figure
plot(vx,(vx - (tau0 + tau1).*(lamb - v))./lamb)