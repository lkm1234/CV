clc
close all
clear all
%% define r space 
xmin = 0;
xmax = 40;
Lx = xmax-xmin;
N = 2^7+1;
dx = Lx/(N-1);
x = xmin:dx:xmax;
[x,y] = meshgrid(x);
%% define k space
dkx = 2*pi/Lx;
%Kx = fftshift(-pi/dx:dkx:pi/dx);
Kx= -pi/dx:dkx:pi/dx
[kx,ky]=meshgrid(Kx);
dt = 1i*1e-3;
%% Potential
V = 50./(1+(x-y).^6/3^6);  
%% define T and V
vg = 10
T = -vg*kx./2;
A = exp(dt/2.*T);
B = exp(-dt.*V);
%% initial state
x0 = 10;
y0 = 20;
gammax = 1;
gammay = 1;
u=exp(-(gammax*(x-x0).^2+gammay*(y-y0).^2)/2)/sqrt(pi^(3/2));
nor1=sum(sum(abs(u).^2))*dx^2;
u=u./sqrt(nor1); % 归一化 centered gaussian
%contourf(x,y,100.*abs(u).^2)
%% time evolution
n=0;
while n<20000
    n=n+1;
    %step1
    U1 = fftshift(fft(u,[],2),2);
    u1 = ifft(ifftshift(A.*U1,2),[],2);
    %step2
    u2 = B.*u1;
    %step3
    U3 = fftshift(fft(u2,[],2),2);
    u3 = ifft(ifftshift(A.*U3,2),[],2);
    u = u3;
    %nor1=sum(sum(abs(u).^2))*dx^2;
    %u=u./sqrt(nor1);
    if mod(n,100)==0
    figure(1)
    subplot(2,2,1)
    contourf(x,y,abs(u))
    title('abs|\phi|');
    colorbar;
    subplot(2,2,2)
    contourf(x,y,abs(u).^2)
    title('abs|\phi|^2');
    colorbar;
    subplot(2,2,3)
    contourf(x,y,real(u))
    title('real(\phi)');
    colorbar;
    subplot(2,2,4)
    contourf(x,y,imag(u))
    title('imag(\phi)');
    colorbar;
    %{
    xlabel('x')
    ylabel('y')
    zlabel('EE')
    %}
    pause(0.1)
    end
end


