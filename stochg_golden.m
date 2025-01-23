clear; 

format compact

% set parameters
alpha = 0.33; % capital's share
beta = 0.95;
Del = 0.1; % depreciation rate (annual)
s = 2;

nk = 20; % grid points + 1
kstar = (alpha/(1/beta - (1-Del)))^(1/(1-alpha)); % steady state k

kmin = 0.25*kstar;
kmax = 4*kstar;
gridk = linspace(kmin,kmax,nk)';

na = 11; % 125 number of transitory shock grids
mua = 0;
m = 4.2;
rhoa = 0.9;
sda = 0.01;
[grida,pdfa] = tauchen(na,mua,rhoa,sda,m);
grida = exp(grida);

v0 = zeros(nk,na);
v1 = zeros(nk,na);
kpol = zeros(nk,na);
diff = 1;
tol = 0.0001;  % tolerance level for value function iteration
its = 1;

smctime   = tic;
totaltime = 0;

while diff > tol

%for ia = 1:na
for ia = 1:na
    Ev = v0*pdfa(ia,:)';   
    for ik = 1:nk
        k1 = gridk(ik);
        income = grida(ia)*k1^alpha + (1 - Del)*k1;
        [v1(ik,ia),kpol(ik,ia)] = golden(income,Ev,gridk,beta,s,kmin);
    end
end

diff = max(max(abs(v0-v1)));
v0 = v1;

totaltime = totaltime + toc(smctime);
avgtime   = totaltime/its;
if mod(its, 20) == 0 || diff<=tol; fprintf('%8.0f ~%8.8f ~%8.5fs ~%8.5fs \n', its, diff, totaltime, avgtime); end
its = its+1;
smctime = tic; % re-start clock

end

consum = zeros(nk,na);
for ia = 1:na
  for ik =1:nk
   consum(ik,ia) = grida(ia)*gridk(ik)^(alpha) - kpol(ik) + (1-Del)*gridk(ik);
  end
end

%% draw figures
lwidth = 1.25;


figure
plot(gridk,v1(:,ceil(na/2 - 0.8*(na-na/2))),'Linewidth',lwidth)
hold on
plot(gridk,v1(:,ceil(na/2)),'Linewidth',lwidth)
plot(gridk,v1(:,ceil(na/2 + 0.8*(na-na/2))),'Linewidth',lwidth)
hold off
xlabel('k')
ylabel('V(k)')
legend('A low','A medium','A high','location','southeast')
title('Value function')

figure
plot(gridk,consum(:,ceil(na/2 - 0.8*(na-na/2))),'Linewidth',lwidth)
hold on
plot(gridk,consum(:,ceil(na/2)),'Linewidth',lwidth)
plot(gridk,consum(:,ceil(na/2 + 0.8*(na-na/2))),'Linewidth',lwidth)
hold off
xlabel('k')
ylabel('c')
legend('A low','A medium','A high','location','southeast')


%%

save sg_golden.mat kpol v1 gridk grida consum na



