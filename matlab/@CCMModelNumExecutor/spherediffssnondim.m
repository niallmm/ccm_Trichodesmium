function dy = spherediffssnondim(t,y,pin)
% has both HCO3- and CO2 diffusing
% the equations in x have been discritized using finite-difference method
p = pin.ccm_params;
x = pin.x;
dx = pin.dx;
xnum = pin.xnum;


dh = zeros(xnum, 1);        %derivative of h
dc = zeros(xnum, 1);        %derivative of c

% the vector y contains both species, where y(1) = h(1), y(2) = c(1),
% y(3) = h(2), ect... we want to convert this to two vectors h = HCO3-, and
% c = CO2
h = y(1:2:xnum*2)';
c = y(2:2:xnum*2)';
Rca = (h-p.gamma*c)./(1+ h + c); %non-dimensional carbonic anhydrase reaction
Rrub = p.Vmax*c./(p.Vba*(p.Km/p.Kca+c)); % non-dimensional RuBisCO reaction
% 


% figure(11)
% plot(h./(1+h+c)-Rrub, 'k')
% % 
% figure(11)
% plot(c, 'r')
% hold on
% plot(h, 'b')
% hold off
% drawnow



% boundary condition at the 

beta = p.xi/(dx^2); % grouping of parameters

dh(1) = beta*((3*(x(1)+x(2))/2)^(4/3))*(h(2)-h(1)) - Rca(1); % change in bicarbonate
dc(1) = beta*((3*(x(1)+x(2))/2)^(4/3))*(c(2)-c(1)) + p.kappa*Rca(1)-p.kappa*Rrub(1); % change in carbon dioxide

% inside carboxysome

xhp(2:xnum-1) = (3*(x(2:xnum-1)+x(3:xnum))/2).^(4/3);
xhm(2:xnum-1) = (3*(x(2:xnum-1)+x(1:xnum-2))/2).^(4/3);
dh(2:xnum-1) = beta*(xhp(2:xnum-1).*(h(3:xnum)-h(2:xnum-1))-xhm(2:xnum-1).*(h(2:xnum-1)-h(1:xnum-2)))- Rca(2:xnum-1);
dc(2:xnum-1) = beta*(xhp(2:xnum-1).*(c(3:xnum)-c(2:xnum-1))-xhm(2:xnum-1).*(c(2:xnum-1)-c(1:xnum-2)))+ p.kappa*Rca(2:xnum-1)- p.kappa*Rrub(2:xnum-1);
% 
% figure(10)
% plot(h./(1+h+c), 'b')
% hold on
% plot(p.gamma*c./(1+h+c), '--b')
% plot(Rrub, 'k')
% plot(beta*(xhp(2:xnum-1).*(h(3:xnum)-h(2:xnum-1))-xhm(2:xnum-1).*(h(2:xnum-1)-h(1:xnum-2))), 'c')
% hold off
% drawnow

% figure(12)
% plot(beta*(xhp(2:xnum-1).*(h(3:xnum)-h(2:xnum-1))-xhm(2:xnum-1).*(h(2:xnum-1)-h(1:xnum-2))))
% drawnow
%at carboxysome boundary

i = xnum;
    
xhp(i) = (3*(x(i)+x(i)+0.5*dx)/2)^(2/3); % 
xhm(i) = (3*(x(i)+x(i-1))/2)^(4/3);
dh(xnum) = beta*(xhp(i)*(p.beta_h*h(i)-p.epsilon_h+p.beta_c2*c(i))*dx - xhm(i)*(h(i)-h(i-1)));
dc(xnum) = beta*(xhp(i)*(p.beta_c*c(i)-p.epsilon_c)*dx - xhm(i)*(c(i)-c(i-1)));

%     

% now we need to convert the derivative vectors dh and dc into one vector
dy = zeros(xnum*2,1);
dy(1:2:xnum*2)= dh;
dy(2:2:xnum*2)= dc;
% 
% figure(14)
% plot(4*pi*p.Rc^2*(dh(xnum)+dc(xnum)) - p.Vcsome*sum(Rrub), 'o')
% drawnow
