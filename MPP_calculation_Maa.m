%calculation of MPP absorber based on Maa's Theory
%==================================================

clear all
clc

%parameter definition
%====================

rho0=1.2;                       %density of air
c=343;                          %sound speed in air
freq=logspace(1,3.2,600);        %frequency range of sound wave  
omega=2*pi*freq;                %angular frequency
eta= 1.789e-5; % coefficient of viscosity
r0=0.15e-3;                     %radius of hole in meter
d=2*r0;                         %diameter of hole in meter
t=0.9e-3;                       %thickness of the plate in meter
k=d*sqrt(omega*rho0/eta/4);     %perforation constant



%calculate perforation area ratio (tau)
%======================================

b=2e-3; % distance between two holes; % distance between two holes

tau=(pi/4)*(d/b)^2;
D=10e-3;
%calculate relative impedance
%============================
%1998 Journal of acoustical society of america
kr=sqrt(1+k.^2/32)+ (1.4142/32)*k*d/t;  %resistance coefficient
r=32*eta*t*kr./(tau*rho0*c*d^2);        %resistance
km=1+(1./sqrt(1+0.5*(k.^2)))+(0.85*d/t);%reactance coefficient
xm=omega.*t.*km./(tau*c);               %reactance
alfa=4*r./((1+r).^2+...                 %absorption coefficient for MPP
    (xm-(cot(omega*D./c))).^2);

figure(1)
plot(freq,alfa,'r-')
hold on

