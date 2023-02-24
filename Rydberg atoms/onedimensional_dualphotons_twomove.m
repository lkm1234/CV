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
Kx= -pi/dx:dkx:pi/dx;
[kx,ky]=meshgrid(Kx);
dt = 1i*1e-3;
%% Potential
V = 50./(1+(x-y).^6/3^6);  
%% define T and V
vg1 = -10;
vg2 = 10;
T1 = vg1*kx;
T2 = vg2*ky;
A1 = exp(-dt/2.*T1); %x direction
A2 = exp(-dt/2.*T2); %y direction
B = exp(-dt.*V);
%% initial state
x0 = 30;
y0 = 10;
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
    %step1 x direction
    U1 = fftshift(fft(u,[],1),1);
    u1 = ifft(ifftshift(A2.*U1,1),[],1);
    %step2 y direction
    U2 = fftshift(fft(u1,[],2),2);
    u2 = ifft(ifftshift(A1.*U2,2),[],2);
    %step3
    u3 = B.*u2;
    %step4
    U4 = fftshift(fft(u3,[],2),2);
    u4 = ifft(ifftshift(A1.*U4,2),[],2);
    %step5
    U5 = fftshift(fft(u4,[],1),1);
    u5 = ifft(ifftshift(A2.*U5,1),[],1);
    u = u5;
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


