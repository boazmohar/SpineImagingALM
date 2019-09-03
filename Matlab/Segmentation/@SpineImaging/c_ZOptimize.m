function c_ZOptimize(obj)
%centering scan fields according to max Z span
%     Center = obj.OptionsStruct.ZmaxTravel./2;
%     CenterZ = ( max(obj.All.z) + min(obj.All.z))/2;
%     obj.Offset = Center - CenterZ;
FieldCenters(:,1)=obj.All.x;
FieldCenters(:,2)=obj.All.y;
FieldCenters(:,3)=obj.All.z+obj.OptionsStruct.offset;

maxAccel=obj.OptionsStruct.maxAccel;
maxVel=obj.OptionsStruct.maxVel;
LineRate=obj.OptionsStruct.lineRate;
linesPerField=obj.linesPerField;
accuracy=obj.OptionsStruct.targetAccuracy;
%CenterZ: Z-position of fields in um
%maxAccel: max acceleration in um/ms^2
%maxVel: max velocity in um/ms
%LineRate: imagaing line rate in Hz
%linesPerField: lines including fly lines
%accuracy: maximum error from ideal position in um


maxAccelOrig=maxAccel;
maxVelOrig=maxVel;
SampRate=LineRate;
maxAccel=maxAccelOrig./((SampRate/1000)^2);
interpMethod='cubic';

%sort fields according to Z position
CenterZ=FieldCenters(:,3);
[CenterZ2,I]=sort(CenterZ,'ascend');
FieldCenters=FieldCenters(I,:);

%Pad begining and end
CenterZ2=[CenterZ2(1)*ones(1,1),CenterZ2',CenterZ2(end)*ones(1,1)]';

hA=obj.Axes2Handle;
plot(hA,CenterZ2)
hold(hA,'on')


fieldDur=linesPerField/LineRate;
fieldDurSamp=round(fieldDur.*SampRate);
fieldDurHalfSamp=round(fieldDurSamp/2);
fieldIntervals=ones(1,length(CenterZ2)-1).*fieldDur;
fieldIntervalsSamp=fieldIntervals*SampRate;
fieldTimesSampPos=round(cumsum([fieldDurHalfSamp,fieldIntervalsSamp]));
%only care about real fields for error checking
RealFieldsIdx=zeros(1,length(fieldTimesSampPos));
RealFieldsIdx(2:(end-1))=1;
fieldTimesCheck=fieldTimesSampPos(RealFieldsIdx==1);

xx=linspace(fieldDur/2+1/SampRate,sum(fieldIntervals)+fieldDur/2,SampRate*sum(fieldIntervals));
F=griddedInterpolant(cumsum([fieldDur/2,fieldIntervals]),CenterZ2,interpMethod);
PosIdeal=[ones(1,fieldDurHalfSamp)*CenterZ2(1),F(xx),ones(1,fieldDurHalfSamp)*CenterZ2(end)];
ZCheck=PosIdeal(fieldTimesCheck);
VelIdeal=diff(PosIdeal);

PosCur=PosIdeal;
VelCur=VelIdeal;
filterFail=1;
FailPt=[];
failField=0;
skipAdd=0;
possibleGoBack=0;
counter=0;
CenterZ3=CenterZ2;
while filterFail
    stop=0;
    counter=counter+1;
    if counter>100
        out=PosCur(fieldTimesCheck);
        disp('step limit reached')
        title(hA,'did not solve')
        break;
    end    
    if ~isempty(FailPt)
        plot(hA,CenterZ3(2:(end-1)))
        pause(0.1)
        
        V0=(CenterZ3(failField)-CenterZ3(failField-1))/fieldDur/1000;
        V1=(CenterZ3(failField+2)-CenterZ3(failField+1))/fieldDur/1000;
        dist=CenterZ3(failField+1)-CenterZ3(failField);
        traj=obj.c_ZOptAna(maxAccelOrig/1.1,maxVelOrig/1.1,SampRate/1000,V0,V1,dist,0);
        traj=traj+CenterZ3(failField);
        nFields=ceil(length(traj)/fieldDurSamp)-1;
        if nFields<1 && ~possibleGoBack
            if (failField+2)<length(CenterZ3)
                failField=failField+1;
            else
                filterFail=0;
            end
            skipAdd=1;
        else
            if possibleGoBack
                failField=failField-1;
                V0=(CenterZ3(failField)-CenterZ3(failField-1))/fieldDur/1000;
                V1=(CenterZ3(failField+2)-CenterZ3(failField+1))/fieldDur/1000;
                dist=CenterZ3(failField+1)-CenterZ3(failField);
                traj=obj.c_ZOptAna(maxAccelOrig/1.1,maxVelOrig/1.1,SampRate/1000,V0,V1,dist,0);
                traj=traj+CenterZ3(failField);
                nFields=ceil(length(traj)/fieldDurSamp)-1;
                skipAdd=0;
            end
        end
        possibleGoBack=0;
        
        if ~skipAdd
            traj=resample(traj,(nFields+1).*fieldDurSamp,length(traj));
            %remove high-frequency garbage produced by resample
            traj=smooth(traj,round(SampRate/1000))';
            newFieldTimes=cumsum(ones(1,nFields)*fieldDurSamp);
            newCenterZ=traj(newFieldTimes);
            CenterZ3=[CenterZ3(1:failField)',newCenterZ,CenterZ3((failField+1):end)']';
            RealFieldsIdx=[RealFieldsIdx(1:failField),zeros(1,nFields),RealFieldsIdx((failField+1):end)];
            fieldIntervals=[fieldIntervals,ones(1,nFields).*fieldDur];
            fieldIntervalsSamp=fieldIntervals*SampRate;
            fieldTimesSampPos=round(cumsum([fieldDurHalfSamp,fieldIntervalsSamp]));
            fieldTimesCheck=fieldTimesSampPos(RealFieldsIdx==1);
        end
        skipAdd=0;
        
        xx=linspace(fieldDur/2+1/SampRate,sum(fieldIntervals)+fieldDur/2,SampRate*sum(fieldIntervals));
        F=griddedInterpolant(cumsum([fieldDur/2,fieldIntervals]),CenterZ3,interpMethod);
        PosCur=[ones(1,fieldDurHalfSamp)*CenterZ3(1),F(xx),ones(1,fieldDurHalfSamp)*CenterZ3(end)];
        VelCur=diff(PosCur);
    end
   while stop~=1
        Accel=diff(VelCur);
        tooFast=find(abs(Accel(fieldDurSamp:(end-fieldDurSamp)))>maxAccel,1,'first')+fieldDurSamp;
        if isempty(tooFast) || tooFast>(length(Accel)-2*round(fieldDurSamp))
            stop=1;
            filterFail=0;
        else
            filterFail=1;
            for i=(2:round(length(VelCur)/8))  
                Istart=max(1,(tooFast-i));
                Istop=min((tooFast+i),length(VelCur));
                VelCurTemp=VelCur;
                seg=VelCur(Istart:Istop);
                %segFilt=moving_average2(seg,round(i/2));
                segFilt=smooth(seg,round(i/2));
                VelCurTemp(Istart:Istop)=segFilt;
                AccelCurTemp=diff(VelCurTemp);
                nextTooFast=find(abs(AccelCurTemp(fieldDurSamp:(end-fieldDurSamp)))>maxAccel,1,'first')+fieldDurSamp;
                if isempty(nextTooFast)
                    nextTooFast=Inf;
                end
                PosCurTemp=cumsum([CenterZ2(1),VelCurTemp]);
                if nextTooFast>tooFast && max(abs(PosCurTemp(fieldTimesCheck)-ZCheck))<accuracy
                    VelCur(Istart:Istop)=segFilt;
                    filterFail=0;
                    break;
                end
            end
            if filterFail
                failFieldNew=find((fieldTimesSampPos-tooFast)>0,1,'first');
                if (failFieldNew+2)>length(CenterZ3)
                    filterFail=0;
                else
                    
                    if failFieldNew>failField
                        failField=failFieldNew;
                        
                    else
                        if failFieldNew==failField;
                            disp('hi')
                            possibleGoBack=1;    
                        end
                        failField=failField;    
                    end   
                    FailPt=tooFast;
                end
                stop=1;
            end
        end
    end
end
plot(hA,CenterZ3(2:(end-1)))
pause(0.1)

PosCur=cumsum([CenterZ2(1),VelCur((fieldDurSamp+1):(end-fieldDurSamp))]);
CenterZ3=CenterZ3(2:(end-1));
RealFieldsIdx=RealFieldsIdx(2:(end-1));

%Interpolate XY position of "fake" fields to maximize likelihood of extra
%data

NewFieldCenters=zeros(length(CenterZ3),3);
NewFieldCenters(:,3)=CenterZ3;
NewFieldCenters(RealFieldsIdx==1,1:2)=FieldCenters(:,1:2);
if sum(RealFieldsIdx==0)>0
    NewCentersXY(:,1)=interp1(find(RealFieldsIdx==1)./length(RealFieldsIdx),FieldCenters(:,1),find(RealFieldsIdx==0)./length(RealFieldsIdx));
    NewCentersXY(:,2)=interp1(find(RealFieldsIdx==1)./length(RealFieldsIdx),FieldCenters(:,2),find(RealFieldsIdx==0)./length(RealFieldsIdx));
    NewFieldCenters(RealFieldsIdx==0,1:2)=NewCentersXY(:,1:2);
end


%Add Flyback
% V0=-(CenterZ3(end)-CenterZ3(end-1))/fieldDur/1000;
V0= -(PosCur(end)-PosCur(end-1))/(1/LineRate)/1000;
% V1=-(CenterZ3(2)-CenterZ3(1))/fieldDur/1000;
V1= -(PosCur(2)-PosCur(1))/(1/LineRate)/1000;
dist=(PosCur(end)-PosCur(1));
% traj=CenterZ3(end)-obj.c_ZOptAna(maxAccelOrig,maxVelOrig,SampRate/1000,V0,V1,dist,0);
traj=PosCur(end)-obj.c_ZOptAna(maxAccelOrig,maxVelOrig,SampRate/1000,V0,V1,dist,0);

traj2 = [repmat(traj(1),1,100) traj repmat(traj(end),1,100)];
nFields=ceil(length(traj)/fieldDurSamp)-1;
nFields2=ceil(length(traj2)/fieldDurSamp)-1;
traj2=resample(traj2,(nFields2+1).*fieldDurSamp,length(traj2));
traj2=smooth(traj2,round(SampRate/1000));
factor = round(100*(nFields+1).*fieldDurSamp./length(traj));
traj3=  traj2(101:end-100);
newFieldTimes=cumsum(ones(1,nFields)*fieldDurSamp);
newCenterZ=traj(newFieldTimes);
CenterZ3=[CenterZ3',newCenterZ]';
%RealFieldsIdx=[RealFieldsIdx(2:(end-1)),zeros(1,nFields)];
RealFieldsIdx=[RealFieldsIdx,zeros(1,nFields)];
fieldIntervals=[fieldIntervals(2:(end-1)),ones(1,nFields).*fieldDur];
fieldIntervalsSamp=fieldIntervals*SampRate;
fieldTimesSampPos=round(cumsum([fieldDurHalfSamp,fieldIntervalsSamp]));
fieldTimesCheck=fieldTimesSampPos(RealFieldsIdx==1);
PosCur=[PosCur,traj];

FlybackFieldCenters(:,1)=interp1([0,1],FieldCenters([end,1],1),[1:nFields]/nFields);
FlybackFieldCenters(:,2)=interp1([0,1],FieldCenters([end,1],2),[1:nFields]/nFields);
FlybackFieldCenters(:,3)=newCenterZ;
NewFieldCenters((end+1):(end+nFields),:)=FlybackFieldCenters;
hold(hA,'off')
% ploting
if obj.OptionsStruct.ShowZoptoFields
    x = NewFieldCenters(~RealFieldsIdx,1);
    y = NewFieldCenters(~RealFieldsIdx,2);
    z = NewFieldCenters(~RealFieldsIdx,3)-obj.OptionsStruct.offset;
    axes(obj.AxesHandle)
    scatter3(x,y,-z,'filled', 'MarkerEdgeColor','[0.8 .8 .8]',...
        'MarkerFaceColor',[0.8 .8 .8])
end

%Output
obj.Zwave=PosCur;
obj.TrueFieldsMask=RealFieldsIdx;
obj.TrueFieldsCenterSamp=fieldTimesCheck;
obj.All.Id(RealFieldsIdx==1) = obj.All.Id(I);
obj.All.Id(RealFieldsIdx==0) = -1;
obj.All.x = NewFieldCenters(:,1);
obj.All.y = NewFieldCenters(:,2);
obj.All.z = NewFieldCenters(:,3)-obj.OptionsStruct.offset;