function traj=ZOptAna(Amax,Vmax,SampRate,V0,V1,dist,Verbose)
%SampRate: samples per ms
Xtarget=dist;
V1=-V1;

T1=-V0/Amax+1/Amax*sqrt((V0^2+V1^2)/2+Amax*Xtarget);
Vhat=V0+Amax*T1;
if abs(Vhat)>Vmax
    T1=(Vmax-V0)/Amax;
end    
Xi=V0*T1+1/2*Amax*(T1)^2;
Vi=V0+Amax*T1;

T3=V1/Amax+1/Amax*sqrt(abs((V1^2+V0^2)/2+Amax*Xtarget));
Vhat=V1-Amax*T3;
if abs(Vhat)>Vmax
    T3=abs((-Vmax-V1)/Amax);    
end  
Xk=Vi*T3-1/2*Amax*(T3)^2;

T2=(Xtarget-Xi-Xk)/Vi;   

if sum([T1,T2,T3]<-1/SampRate)>0
    %if Verbose
        disp('velocity change and distance change are incompatible')
    %end
    traj=NaN;
    return
end

Ts1=linspace(0,T1,SampRate*T1);
Ts2=linspace(0,T2,SampRate*T2);
Ts3=linspace(0,T3,SampRate*T3);

F1= @(t) V0.*t+1/2.*Amax.*t.^2;
F2= @(t) Xi+Vi.*t;
F3= @(t) Xtarget-Xk+Vi.*t-1/2.*Amax.*t.^2;

traj=[F1(Ts1),F2(Ts2),F3(Ts3)];
if Verbose
    if isempty(Ts1)
        Inc1=0;
    else
        Inc1=Ts1(end);
    end
    if isempty(Ts2)
        Inc2=0;
    else
        Inc2=Ts2(end);
    end
    plot([Ts1,Inc1+Ts2,Inc1+Inc2+Ts3],traj)
end