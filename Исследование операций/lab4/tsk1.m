A1=[2000 -880;
    -940 1760];
A1=A1+940;

N=1000;
S=zeros(1, N+1);
for k=1:N
    x=myrand([1/2 1/2])+1;
    y=myrand([1/2 1/2])+1;
    S(k+1)=S(k)+A1(x,y);
end
figure; hold on; grid on
plot(1:N+1,S)
 
% по оптимальной
S=zeros(1, N+1);
for k=1:N
    x=myrand([0.4731 1-0.4731])+1;
    y=myrand([1/2 1/2])+1;
    S(k+1)=S(k)+A1(x,y);
end
plot(1:N+1,S,'r')
legend('optimum', 'rand')

%%
A1=[2000 -880;
    -940 1760];
A1=A1+940;
f=ones(1,2);b=ones(2,1);
D=zeros(2,2);
[y,fmin]=linprog(-f,[A1;D-diag(ones(2,1))],[b;zeros(2,1)])
v=-1/fmin
q=v*y

%%
A2=[2 5 6;
    5 3 1];
f=ones(1,3);b=ones(2,1);
D=zeros(3,3);
[y,fmin]=linprog(-f,[A2;D-diag(ones(3,1))],[b;zeros(3,1)])
v=-1/fmin
q=v*y

%%
N=1000;
S=zeros(1, N+1);
for k=1:N
    x = myrand([1/2 1/2])+1;
    y = myrand([1/3 1/3 1/3])+1;
    S(k+1)=S(k)+A2(x,y);
end
figure; hold on; grid on
plot(1:N+1,S)
 
% по оптимальнй
S=zeros(1, N+1);
for k=1:N
    S(k+1)=S(k)+A2(x,y);
    y = myrand([0.6250 0 1-0.6250])+1;
    x=myrand([1/2 1/2])+1;
    S(k+1)=S(k)+A2(x,y);
end
plot(1:N+1,S,'r')
legend('optimum','rand')

%%
A3=[ 180 80 -20 -120;
     60 360 260 160;
     -60 240 540 440;
     180 120 420 720];
q=[0.2 0.35 0.25 0.2];
for n=1:4
a(n)=sum(A3(n,:)*q(n));
end
a
max(a)

%%

for n=1:4
a(n)=sum(A3(n,:)*q(1));
end
a
max(a)



