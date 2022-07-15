

close all;

N=500;

Gss=5 * 1e-3;  % 10 mT/m
g=ones(N,1)*Gss;

gamma=42.576*1e6;
slice_thick=5e-3;  % 5 mm
omega=gamma*Gss*slice_thick;

deltat=1e-5;
Trf=deltat * N;

times=([1:N]'-0.5-N/2)*deltat;
x=omega*times;

f1=sinc(x);
f3=0.5 + 0.5*cos( 2*pi*times/Trf );  % hanning window
        
b1=f1.*f3;

% add gradient ramps
Nramp=5;
gup=linspace(0,Gss,Nramp)';
gdown=linspace(Gss,0,Nramp)';
g=[gup;g;gdown];
b1=[zeros(Nramp,1);b1;zeros(Nramp,1)];

b1 = [0 ;b1; 0];	
g = [0 ;g; 0];

Gmax=40 * 1e-3 ;  % 40 mT/m
Smax=150;  % 150 T/m/s

% Gmax=Gmax / Gss;  % unit of Gss
% Smax=Smax / Gss * deltat;  % unit of Gss/deltat

[b1v gv] = mintverse(b1,g,deltat,1.0*max(b1),Gmax,Smax,deltat);

figure;
subplot(2,1,1); plot(b1);
subplot(2,1,2); plot(g);

figure;
subplot(2,1,1); plot(b1v);
subplot(2,1,2); plot(gv);

% Nv=size(b1v,1);
% slew=zeros(Nv,1);
% for i=2:Nv
%     slew(i)=( gv(i) - gv(i-1) ) /deltat;
% end
% 
% figure;
% plot(slew);


