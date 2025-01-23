clear;

tic
global v0 beta Del alpha kmat k0 s a0 ia pdfa

% set parameters
alpha = 0.33; % capital's share
beta = 0.95;
Del = 0.1; % depreciation rate (annual)
s = 2;

kgrid = 20; % grid points + 1
kstar = (alpha/(1/beta - (1-Del)))^(1/(1-alpha)); % steady state k

cstar = kstar^(alpha) - Del*kstar;
istar = Del*kstar;
ystar = kstar^(alpha);

kmin = 0.25*kstar;
kmax = 4*kstar;
grid = (kmax-kmin)/kgrid;

kmat = kmin:grid:kmax;
kmat = kmat';
[nk,~] = size(kmat);

na = 7; % number of transitory shock grids
mua = 0;
m = 4.2;
rhoa = 0.9;
sda = 0.01;
[amat,pdfa] = tauchen(na,mua,rhoa,sda,m);
amat = exp(amat);

v0 = ones(nk,na);
v1 = ones(nk,na);
k11 = zeros(nk,na);
diff = 1;      
tol = 0.0001;
its = 1;
maxits = 2000;

smctime   = tic;
totaltime = 0;

while diff>tol 

  for ia = 1:na
      a0 = amat(ia,1);
    for ik = 1:nk
      k0 = kmat(ik,1);
      k1 = fminbnd(@valfun_stoch,kmin,kmax);
      v1(ik,ia) = - valfun_stoch(k1);
      k11(ik,ia) = k1;
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
   consum(ik,ia) = amat(ia)*kmat(ik)^(alpha) - k11(ik) + (1-Del)*kmat(ik);
  end
end

toc
%% draw figures
lwidth = 1.25;


figure
plot(kmat,v1(:,ceil(na/2 - 0.8*(na-na/2))),'Linewidth',lwidth)
hold on
plot(kmat,v1(:,ceil(na/2)),'Linewidth',lwidth)
plot(kmat,v1(:,ceil(na/2 + 0.8*(na-na/2))),'Linewidth',lwidth)
hold off
xlabel('k')
ylabel('V(k)')
legend('A low','A medium','A high','location','southeast')
title('Value function')

figure
plot(kmat,consum(:,ceil(na/2 - 0.8*(na-na/2))),'Linewidth',lwidth)
hold on
plot(kmat,consum(:,ceil(na/2)),'Linewidth',lwidth)
plot(kmat,consum(:,ceil(na/2 + 0.8*(na-na/2))),'Linewidth',lwidth)
hold off
ylabel('c')
xlabel('k')
legend('A low','A medium','A high','location','southeast')


save sg.mat v1 kmat consum k11 amat

