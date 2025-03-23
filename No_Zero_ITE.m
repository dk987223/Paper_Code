clear 
clc 
mi=1607;
ci=0.414;
fi=0.02;
g=9.8;
W_rad=0.28;
eff = 0.85;

oumiga_hat = 0.01*ones(6,1);   
cij=0.10*[-1 -0.5 0 0.5 1 0.5;
          -1 -0.5 0 0.5 1 0.5];
bj = 5.0;
h = zeros(6,1); 
u = zeros(6,1);
h_i = 1;
gwwi = 0.0001;
QUAN_wwi = 0.001;
adp_wwi =  0.1;

gwwwi = 0.0001;
QUAN_wwwi = 0.001;
adp_wwwi =  0.1;
v0 = 2.2;
vi = [1.9 1.8 2 2.1 2 2]';
p0 = 0;
p1 = 0.4;
standstill_distance = 10;
 x0 = 60;
 xi = [50.1  40.2  30.3 20.4 10.5 0.2 ]';
deta=0.01;
q = 0.85;
Q = [q -1 0 0 0 0;0 q -1 0 0 0;0 0 q -1 0 0;0 0 0 q -1 0;0 0 0 0 q -1;0 0 0 0 0 q];
inte = 0;
langda = 2;
adp_C = 0.8;
adp_f = 0.05;
gci = 0.0001;
gfi = 0.0001;
gyi = 0.0001;
QUAN_C = 0.001;
QUAN_f = 0.001;
gou = 0.001;
Siinte = 0;
aaaa=1;
ecci = 2;
di = [x0-xi(1),xi(1)-xi(2),xi(2)-xi(3),xi(3)-xi(4),xi(4)-xi(5),xi(5)-xi(6)]'; 
eevi= [vi(1)-v0,vi(2)-vi(1),vi(3)-vi(2),vi(4)-vi(3),vi(5)-vi(4),vi(6)-vi(5)]';
ex0 = di - p0*vi.*vi- p1.*vi - standstill_distance;
period = 0.001;
xita1 = 0.3;
xita2 = 0.5;
epexilon_1 = 0.7;
epexilon_2 = 0.6;

fai = 0.4;
R0 = 0.5;
Rinf =0.08;
i_wi = 0;
wii = 0;
k3 = 20;
terminal = 30;
zihao = 14;
xiankuan = 2;

for t = 0:period:terminal
    aaa = t ;

    if    t > 30
       a0 = 0 ;
    elseif t>0 && t<= 5 
       a0 = 3;
    elseif t>15 && t< 18
        a0 = -2;
    else 
        a0 = 0;
    end 
   
 for i = 1:6
     fi_1(i) = (xi(i)^2/(1+xi(i)^2)+0.1*vi(i)^2+0.5*vi(i))*exp(-(t-1+0.5*i)^2)+ 0.1*sin(t);
 end
 
       di = [x0-xi(1),xi(1)-xi(2),xi(2)-xi(3),xi(3)-xi(4),xi(4)-xi(5),xi(5)-xi(6)]'; 
       evi= [vi(1)-v0,vi(2)-v0,vi(3)-v0,vi(4)-v0,vi(5)-v0,vi(6)-v0]';
       evi1= [vi(1)-v0,vi(2)-vi(1),vi(3)-vi(2),vi(4)-vi(3),vi(5)-vi(4),vi(6)-vi(5)]';
       exi = di - p0*vi.*vi- p1.*vi - standstill_distance;
       ro = (R0 - Rinf)*exp(-fai*t) + Rinf;
       dro = -fai*(R0 - Rinf)*exp(-fai*t);
       yyi = (ex0 + ((ecci*ex0) + eevi)*t)*exp(-ecci*t);
       dyyi = -ecci*(ex0 + (ecci*ex0 + eevi)*t)*exp(-ecci*t) + (ecci*ex0 + eevi)*exp(-ecci*t);
       eexi =  exi - yyi;   

       if abs(max(ex0)) < 10e-8 &&  abs(min(ex0)) < 10e-8
           PL = -xita1 * ro ;
           PR = xita2 * ro ;
           DPL = -xita1 * dro;
           DPR = xita2 * dro ;
       else 
           PL = (epexilon_1*sign(min(ex0)) - xita1)*ro - Rinf * sign(min(ex0));
           DPL = (epexilon_1*sign(min(ex0)) - xita1)*dro;
           PR = (epexilon_2*sign(max(ex0)) + xita2)*ro - Rinf * sign(max(ex0));
           DPR = (epexilon_2*sign(max(ex0)) - xita2)*dro;
       end
       va = (eexi - PL)./(PR - PL);
       r1 = 1./((1-va).*(PR-PL).*va);
       b1 = (PL.*DPL - DPL.*PR - eexi.*(DPR - DPL))./(PR - PL);
       QUAN_y = 10*exp(-4*t);
       wi = log(va./(1-va)); 
       i_wi = i_wi + wii*period;
       wii = wi;
       si = wi + langda * i_wi ;
       Si = Q*si;
       dci = gci* (2*p0.*vi+p1).*Si.*vi.*vi - QUAN_C.*gci.*adp_C;
       adp_C = adp_C + dci.*period;
       dfi = mi*g*gfi*(2.*p0.*vi+p1).*Si - QUAN_f.*gfi.*adp_f;
       adp_f = adp_f + dfi*period;
       
      dexi=[v0-vi(1)-p1*u(1)-dyyi(1);
      vi(1)-vi(2)-p1*u(2)-dyyi(2);
      vi(2)-vi(3)-p1*u(3)-dyyi(3);
      vi(3)-vi(4)-p1*u(4)-dyyi(4);
      vi(4)-vi(5)-p1*u(5)-dyyi(5)
      vi(5)-vi(6)-p1*u(6)--dyyi(6)];
        mmi = [exi,dexi]';
for j = 1:1:6
h(j) = exp(-norm(mmi-cij(:,j))^2/(2*bj^2));
end

       dci = r1.*gci.* (2*p0.*vi+p1).*Si.*vi.*vi - QUAN_C.*gci.*adp_C;
       adp_C = adp_C + dci.*period;
       dfi = r1.*mi.*g.*gfi.*(2.*p0.*vi+p1).*Si - QUAN_f.*gfi.*adp_f;
       adp_f = adp_f + dfi*period;
       dwwi = -r1.*mi.*gwwi.*(2.*p0.*vi+p1).*Si.* h- QUAN_wwi.*gwwi.*adp_wwi;
       adp_wwi = adp_wwi + dwwi*period;
       dwwwi = r1.*mi.*gwwwi.*(2.*p0.*vi+p1).*Si.*tanh(Si./0.1) - QUAN_wwwi.*gwwwi.*adp_wwwi;
       adp_wwwi = adp_wwwi + dwwwi*period;
       
 


         for i = 1:6
            Siinte =  abs(Si) * period + Siinte ; 
            if  i == 1
                u(6) = (W_rad/eff)*(mi*b1(6)./(2*p0*vi(6)+p1) + mi*dyyi(6)./(2*p0*vi(6)+p1) + mi*adp_wwwi(6)*tanh(Si(6)./0.1)- mi*adp_wwi(6)*h(6)+ mi*(vi(5)-vi(6))./(2*p0*vi(6)+p1) + adp_C(6).*vi(6).^2+ mi*g*adp_f(6)+ mi*(k3*Si(6)/(norm(Si(6))+deta) + q*langda*wi(6))./(q*r1(6)*(2*p0*vi(6)+p1)));
            elseif i > 1 && i < 6
                u(7-i) = (W_rad/eff)*(mi*b1(7-i)./(2*p0*vi(7-i)+p1)+ mi*adp_wwwi(7-i)*tanh(Si(7-i)./0.1)- mi*adp_wwi(7-i)*h(7-i) + mi*dyyi(7-1)./(2*p0*vi(7-i)+p1)+ mi*(vi(7-1-i)-vi(7-i))./(2*p0*vi(7-i)+p1) +  adp_C(7-i).*vi(7-i).^2+ mi*g*adp_f(7-i)...
                 + mi./(q*r1(7-i).*(2*p0*vi(7-i)+p1))*(k3*Si(7-i)/(norm(Si(7-i))+deta) - q*langda*wi(7-i) - langda*wi(7-i+1)...
                 -r1(7-i+1)* (vi(7-i)-vi(7-i+1)-(2*p0*vi(6-i+1)+p1)./mi*(eff/W_rad*u(7-i+1) - ci*vi(7-i+1)^2 - mi*g*fi + mi*fi_1(7-i+1)) + b1(7-i+1))));                  
            else  
               u(1) = (W_rad/eff)*(mi*b1(1)./(2*p0*vi(1)+p1)+ mi*adp_wwwi(1)*tanh(Si(1)./0.1)- mi*adp_wwi(6)*h(6) + mi*dyyi(1)./(2*p0*vi(1)+p1) + mi*(v0-vi(1))./(2*p0*vi(1)+p1) + adp_C(1).*vi(1).^2+ mi*g*adp_f(1)...
                + mi./(q*r1(1).*(2*p0*vi(1)+p1))*(k3*Si(1)/(norm(Si(1))+deta) - q*langda*wi(1) - langda*wi(2)...
               -r1(2)* (vi(1)-vi(2)-(2*p0*vi(2)+p1)./mi*(eff/W_rad*u(2) - ci*vi(2)^2 - mi*g*fi + mi*fi_1(2)) + b1(2))));  
            end 
         end
       
    da0 = a0;
    v0 =  v0 + da0*period ; 
    x0 =  x0 + v0*period; 
    
    dvi = 1/mi*(eff*u/W_rad - ci*vi.*vi-mi*g*fi) + fi_1';
    vi  = vi + dvi*period;
    xi  = xi + vi*period;
    ccc(aaaa) = PL;
    ddd(aaaa) = PR;
    x1(aaaa) = xi(1);
    v1(aaaa) = vi(1);
    x2(aaaa) = xi(2);
    v2(aaaa) = vi(2);
    x3(aaaa) = xi(3);
    v3(aaaa) = vi(3);
    x4(aaaa) = xi(4);
    v4(aaaa) = vi(4);
    x5(aaaa) = xi(5);
    v5(aaaa) = vi(5);
    x6(aaaa) = xi(6);
    v6(aaaa) = vi(6);
    xx0(aaaa) = x0;
    vv0(aaaa) = v0;
    ee1(aaaa) = eexi(1);
    ee2(aaaa) = eexi(2);
    ee3(aaaa) = eexi(3);
    ee4(aaaa) = eexi(4);
    ee5(aaaa) = eexi(5);
    ee6(aaaa) = eexi(6);
    SS1(aaaa) = Si(1);
    SS2(aaaa) = Si(2);
    SS3(aaaa) = Si(3);
    SS4(aaaa) = Si(4);
    SS5(aaaa) = Si(5);
    SS6(aaaa) = Si(6);
    
    
    aaaa = aaaa+1;
end



t = 0:period:terminal;



A = plot(t,xx0,t,x1,t,x2,t,x3,t,x4,t,x5,t,x6)
A(1).LineWidth = xiankuan;
A(2).LineWidth = xiankuan;
A(3).LineWidth = xiankuan;
A(4).LineWidth = xiankuan;
A(5).LineWidth = xiankuan;
A(6).LineWidth = xiankuan;
A(7).LineWidth = xiankuan;
legend('$x_0$' ,'$x_1$','$x_2$','$x_3$','$x_4$','$x_5$','$x_6$','Interpreter', 'latex')
set(gca,'FontName','Times New Roman','FontSize',zihao,'FontWeight','normal')
xlabel('time (s)','FontSize',zihao)
set(gca,'FontName','Times New Roman','FontSize',zihao,'FontWeight','normal')
ylabel('Position (m)','FontSize',zihao)
set(gca,'FontName','Times New Roman','FontSize',zihao,'FontWeight','normal')
figure;



B=plot(t,vv0,t,v1,t,v2,t,v3,t,v4,t,v5,t,v6)
B(1).LineWidth = xiankuan;
B(2).LineWidth = xiankuan;
B(3).LineWidth = xiankuan;
B(4).LineWidth = xiankuan;
B(5).LineWidth = xiankuan;
B(6).LineWidth = xiankuan;
B(7).LineWidth = xiankuan;
legend('$v_0$' ,'$v_1$','$v_2$','$v_3$','$v_4$','$v_5$','$v_6$','Interpreter', 'latex')
set(gca,'FontName','Times New Roman','FontSize',zihao,'FontWeight','normal')
xlabel('time (s)','FontSize',zihao)
set(gca,'FontName','Times New Roman','FontSize',zihao,'FontWeight','normal')
ylabel('Velocity (m/s)','FontSize',zihao)
set(gca,'FontName','Times New Roman','FontSize',zihao,'FontWeight','normal')
figure;




C = plot(t,ccc,'-.',t,ee1,t,ee2,t,ee3,t,ee4,t,ee5,t,ee6,t,ddd,'--')
C(1).LineWidth = xiankuan;
C(2).LineWidth = xiankuan;
C(3).LineWidth = xiankuan;
C(4).LineWidth = xiankuan;
C(5).LineWidth = xiankuan;
C(6).LineWidth = xiankuan;
C(7).LineWidth = xiankuan;
C(8).LineWidth = xiankuan;
legend(	'$Z_{li}$','$\bar e_1$' ,'$\bar e_2$','$\bar e_3$','$\bar e_4$','$\bar e_5$','$\bar e_6$','$Z_{ri}$','Interpreter', 'latex')
set(gca,'FontName','Times New Roman','FontSize',zihao,'FontWeight','normal')
xlabel('time (s)','FontSize',zihao)
set(gca,'FontName','Times New Roman','FontSize',zihao,'FontWeight','normal')
ylabel('Tracking error (m)','FontSize',zihao)
set(gca,'FontName','Times New Roman','FontSize',zihao,'FontWeight','normal')
figure;

C = plot(t,SS1,t,SS2,t,SS3,t,SS4,t,SS5,t,SS6)
C(1).LineWidth = xiankuan;
C(2).LineWidth = xiankuan;
C(3).LineWidth = xiankuan;
C(4).LineWidth = xiankuan;
C(5).LineWidth = xiankuan;
C(6).LineWidth = xiankuan;
legend(	'$\bar S_1$' ,'$\bar S_2$','$\bar S_3$','$\bar S_4$','$\bar S_5$','$\bar S_6$','Interpreter', 'latex')
set(gca,'FontName','Times New Roman','FontSize',zihao,'FontWeight','normal')
xlabel('time (s)','FontSize',zihao)
set(gca,'FontName','Times New Roman','FontSize',zihao,'FontWeight','normal')
ylabel('Sliding surface ','FontSize',zihao)
set(gca,'FontName','Times New Roman','FontSize',zihao,'FontWeight','normal')


(sum(abs(ee1))+sum(abs(ee2))+sum(abs(ee3))+sum(abs(ee4))+sum(abs(ee5))+sum(abs(ee6)))*0.001
max(ee1),max(ee2),max(ee3),max(ee4),max(ee5),max(ee6)



