function [val,k11] = golden(income,Ev,gridk,beta,s,kmin)

wgss = 0.9;
mink = 0.1*kmin;
maxk = 0.999*income;
bgss = wgss*mink + (1-wgss)*maxk;%%% min amount of assets for next period
cgss= wgss*maxk + (1-wgss)*mink; %%% max amount of assets for next period

%  Evaluate at lower bound for capital
con = income - bgss;
if con <=0
    ubgss = -8888888888888888-800*abs(c); % keeps it from going negative
else
    nextv = interp1(gridk,Ev,bgss,'spline');
    ubgss = con.^(1-s)/(1-s) + beta*nextv;
end

%  Evaluate at lower bound for capital
con = income - cgss;
if con <=0
    ucgss = -8888888888888888-800*abs(c); % keeps it from going negative
else
    nextv = interp1(gridk,Ev,cgss,'spline');
    ucgss = con.^(1-s)/(1-s) + beta*nextv;
end
%% bracket the maximum

dist = abs(cgss-bgss);

while dist>=.0002
    
    if ubgss < ucgss
        mink = bgss;
        bgss = cgss;
        ubgss = ucgss;
        cgss = wgss*bgss + (1-wgss)*maxk;
        
        con = income-cgss;
        if con <=0
            ucgss = -8888888888888888-800*abs(c); % keeps it from going negative
        else
            nextv = interp1(gridk,Ev,cgss,'spline');
            ucgss = con.^(1-s)/(1-s) + beta*nextv;
        end
        
    else % ucgss<ubgss
        maxk = cgss;
        cgss = bgss;
        ucgss = ubgss;
        bgss = wgss*cgss + (1-wgss)*mink;
        
        con = income-bgss;
        if con <=0
            ubgss = -8888888888888888-800*abs(c); % keeps it from going negative
        else
            nextv = interp1(gridk,Ev,bgss,'spline');
            ubgss = con.^(1-s)/(1-s) + beta*nextv;
        end       
    end
  
    dist=abs(cgss-bgss);
end

if (ubgss>ucgss) 
  val = ubgss;
  k11 = bgss;
else
  val = ucgss;
  k11 = cgss;

end

