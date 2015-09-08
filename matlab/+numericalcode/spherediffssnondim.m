function dy = spherediffssnondim(t,y,p)
% has both HCO3- and CO2 diffusing
% the equations in x have been discritized using finite-difference method

Rc = 1;

dh = zeros(p.xnum, 1);        %derivative of h
dc = zeros(p.xnum, 1);        %derivative of c

% the vector y contains both species, where y(1) = h(1), y(2) = c(1),
% y(3) = h(2), ect... we want to convert this to two vectors h = HCO3-, and
% c = CO2
h = y(1:2:p.xnum*2)';
c = y(2:2:p.xnum*2)';

% boundary condition at the center
beta = p.xi/(p.dx^2); % grouping of parameters
Rca = (h(1)-p.gamma*c(1))/(1 + h(1) + c(1)); %non-dimensional carbonic anhydrase reaction
Rrub = p.Vmax*c(1)/(p.Vba*(p.Km/p.Kca+c(1))); % non-dimensional RuBisCO reaction
dh(1) = beta*((3*(p.x(1)+p.x(2))/2)^(4/3))*(h(2)-h(1)) - Rca; % change in bicarbonate
dc(1) = beta*((3*(p.x(1)+p.x(2))/2)^(4/3))*(c(2)-c(1)) + p.kappa*Rca-p.kappa*Rrub; % change in carbon dioxide

% inside carboxysome

for i = 2:p.xnum-1
    xhp = (3*(p.x(i)+p.x(i+1))/2)^(4/3); % factor from finite difference scheme in spherical coordinates
    xhm = (3*(p.x(i)+p.x(i-1))/2)^(4/3);
    Rca = (h(i)-p.gamma*c(i))/(1+ h(i) + c(i));
    Rrub = p.Vmax*c(i)/(p.Vba*(p.Km/p.Kca+c(i)));
    dh(i) = beta*(xhp*(h(i+1)-h(i))-xhm*(h(i)-h(i-1)))- Rca;
    dc(i) = beta*(xhp*(c(i+1)-c(i))-xhm*(c(i)-c(i-1)))+ p.kappa*Rca- p.kappa*Rrub;
end

%at carboxysome boundary

i = p.xnum;
    
xhp = (3*(p.x(i)+p.x(i)+0.5*p.dx)/2)^(2/3); % 
xhm = (3*(p.x(i)+p.x(i-1))/2)^(4/3);
dh(i) = beta*(xhp*(p.beta_h*h(i)-p.epsilon_h+p.beta_c2*c(i))*p.dx - xhm*(h(i)-h(i-1)));
dc(i) = beta*(xhp*(p.beta_c*c(i)-p.epsilon_c)*p.dx - xhm*(c(i)-c(i-1)));
    

% now we need to convert the derivative vectors dh and dc into one vector
dy = zeros(p.xnum*2,1);
dy(1:2:p.xnum*2)= dh;
dy(2:2:p.xnum*2)= dc;