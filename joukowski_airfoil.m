% This code is used to solve the flow over an joukowski airfoil by using joukowski 
% transformation method
%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg
% 12 - 5 - 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all ; clc;
%% givens
number=100;  % length of x1 and y1 (must be even)
c=1; %chord
Cc=5/c/100; %max.camber/c
tc=9/c/100; %max.thickness/c
v_inf=100; % stream velocity
alpha=pi/180*[0,5,10,15]; % AOA
%% circle parameter
b=c/4; 
e=tc/1.3; 
beta=2*Cc; 
a=b*(1+e)/cos(beta);
x0=-b*e; 
y0=a*beta;
%% airfoil coordinates
theta=linspace(3,357,number)*pi/180; 
x1=2*b*cos(theta); 
y1=2*b*e*(1-cos(theta)).*sin(theta)+2*b*beta*(sin(theta)).^2;
r=b*(1+e*(1-cos(theta))+beta*sin(theta));
if length(y1)/2 == floor(length(y1)/2)
    Camb=0.5*(y1(1:length(y1)/2)+y1(end:-1:length(y1)/2+1));
else
    Camb=0.5*(y1(1:floor(length(y1)/2))+y1(end:-1:floor(length(y1)/2)+2));
end
x1=x1+c/2;
x1(1)=c;
x1(number/2)=0;
x1(end)=c;
y1(1)=0;
y1(end)=0;
y1(number/2)=0;
y1(number/2+1)=0;
Camb(1)=0;
Camb(end)=0;
%% Z-plane
x=r.*cos(theta);
y=r.*sin(theta);
%% Zd-plane
xd=x-x0*ones(1,length(x0)); 
yd=y-y0*ones(1,length(y0)); 
rd=sqrt(xd.^2+yd.^2); 
thetad=atan2(yd,xd);
for n=1:length(alpha)
    v_thetad{n}=-v_inf*(sin(thetad-alpha(n)*ones(1,length(thetad))).*(1+a^2*rd.^-2)+2*(a*rd.^-1)*sin(alpha(n)+beta));
end
%% Z1-plane
for m=1:length(alpha)
    v1{m}=v_thetad{m}.*(sqrt(1*ones(1,length(theta))-2*b^2*r.^-2.*cos(2*theta)+b^4*r.^-4)).^-1;
    cp{m}=1-(v1{m}.*(v_inf*ones(1,length(theta))).^-1).^2;
    cl(m)=2*pi*(1+e)*sin(alpha(m)+beta);
end
for l=1:length(alpha)
    [ geom, iner, cpmo ] = polygeom( x1, cp{l} );
    cm0(l)=-geom(1)*geom(2);
    xcp(l)=geom(2);
    n(l)=geom(1);
    m(l)=geom(2);
end
for f=2:length(n)
    xac(f-1)=m(f-1)-((m(f)-m(f-1))*n(f)/(n(f-1)-n(f)));
end
%% saving data
save('ComformalMappingSolution.mat', 'x1', 'y1','cl','cm0','xcp','xac','v1','alpha','v_inf','c');
%% airfoil plot
figure(1);
plot(x1,y1,x1(1:length(Camb)),Camb,'--','LineWidth',2);
xlim([0,c])
grid on;
title('Joukwoski Airfoil','Fontsize',12)
xlabel('Airfoil x-axis','Fontsize',12)
ylabel('Airfoil y-axis','Fontsize',12)
set(gcf,'Color','w')
axis equal
legend('Airfoil','Camber')
%% Cl alpha plot
if length(alpha) >=2
    figure(2);
    plot(alpha*180/pi,cl,'LineWidth',2);
    grid on;
    title('C_l VS \alpha','Fontsize',12)
    xlabel('\alpha ^o','Fontsize',12)
    ylabel('C_l','Fontsize',12)
    set(gcf,'Color','w')
else
    display(cl);
end
%% cmo alpha plot
if length(alpha) >=2
    figure(3);
    plot(alpha*180/pi,cm0,'LineWidth',2);
    grid on;
    title('C_m_0 VS \alpha','Fontsize',12)
    xlabel('\alpha ^o','Fontsize',12)
    ylabel('C_m_0','Fontsize',12)
    set(gcf,'Color','w')
else
    display(cmo);
end
%% xcp alpha plot
if length(alpha) >=2
    figure(4);
    plot(alpha*180/pi,xcp,'LineWidth',2); 
    grid on;
    title('X_c_p VS \alpha','Fontsize',12)
    xlabel('\alpha ^o','Fontsize',12)
    ylabel('X_c_p','Fontsize',12)
    set(gcf,'Color','w')
else
    display(xcp);
end
%% xac plotting
figure(5);
plot(alpha(2:end)*180/pi,xac,'LineWidth',2);
grid on;
title('X_a_c VS \alpha','Fontsize',12)
xlabel('\alpha ^o','Fontsize',12)
ylabel('X_a_c','Fontsize',12)
set(gcf,'Color','w')
%% v1/v-inf & cp plotting
for XxX=1:length(alpha)
    % v1/v-inf
    figure(2*XxX-1+5);
    hold all
    plot(x1,abs(v1{XxX}.*(v_inf*ones(1,length(theta))).^-1),'LineWidth',2);
    grid on;
    title(['Velocity ratio (V/V_i_n_f) \alpha = '  num2str(alpha(XxX)*180/pi) ' ^o'] ,'Fontsize',12)
    xlabel('Airfoil x-axis','Fontsize',12)
    ylabel('V/V_i_n_f','Fontsize',12)
    set(gcf,'Color','w')
    % cp
    figure(2*XxX+5);
    hold all
    h6=plot(x1,cp{XxX},'LineWidth',2); 
    grid on;
    title(['C_p distribution \alpha = '  num2str(alpha(XxX)*180/pi) ' ^o'] ,'Fontsize',12)
    xlabel('Airfoil x-axis','Fontsize',12)
    ylabel('C_p','Fontsize',12)
    set(gcf,'Color','w')
end