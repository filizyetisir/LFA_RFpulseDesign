

function [mxf myf mzf]=bloch_simulation_z(spokes,system,zs,b0_off,pulse,adj)
% function [mxf myf mzf]=bloch_simulation_z(spokes,system,zs,b0_off)

mxy0=0.0; mz0=1.0;  % initial magnetization state
gamma=42.576*1e6;
q=zeros(size(zs,1),4); q(:,1)=1; q(:,4)=1;
B1combined = 0;
for i = 1:size(pulse.LFA.rf_full,2)
    B1combined = B1combined + pulse.LFA.rf_full(:,i).*adj.b1mapsWhole(i,end/2,end/2);
end
for time=1:spokes.ntimes
    
    % Cayley-Klein parameters
    grad=zs*spokes.g(time,3) + b0_off/gamma;
    
    b1s = B1combined(time)*ones(size(grad));
%     if time<=spokes.ntimesperspoke    
%         b1s=repmat( 6.4e-8*spokes.sinc(time),size(zs,1),1 );
%     else
%         b1s=zeros(size(zs,1),1);
%     end
    norm=sqrt( abs(b1s).^2+grad.^2 );
    phi=-2*pi*gamma*system.deltat*norm;
    if sum(abs(phi))<=0.0
        continue;
    end
    n=[real(b1s) imag(b1s) grad];
    n=n./repmat(norm,1,3);
    a=cos(phi/2.0) - 1j*(n(:,3)).*sin(phi/2.0);
    b=sin(phi/2.0).*n(:,2) - 1j*sin(phi/2.0).*n(:,1);
    ind2=find(norm==0);
    a(ind2)=1.0;
    b(ind2)=0.0;
    
    % update rotation matrix
    q1=q(:,1); q2=q(:,2); q3=q(:,3); q4=q(:,4);
    q(:,1)=a.*q1 - conj(b).*q3;
    q(:,2)=a.*q2 - conj(b).*q4;
    q(:,3)=b.*q1 + conj(a).*q3;
    q(:,4)=b.*q2 + conj(a).*q4;
    
end

% final magnetization
a=q(:,1); b=q(:,3);
mxy=(conj(a).*conj(a))*mxy0 - (b.*b).*conj(mxy0) + (2.0*conj(a).*b).*mz0;
mz=(-conj(a).*conj(b))*mxy0 - (a.*b).*conj(mxy0) + (a.*conj(a)-b.*conj(b))*mz0;

mxf=real(mxy);
myf=imag(mxy);
mzf=real(mz);











