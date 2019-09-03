function [Zwave,NewFieldCenters,TrueFieldsMask,TrueFieldsCenterSamp]=ZOptimize(FieldCenters,Options)

maxAccel=Options.maxAccel;
maxVel=Options.maxVel;
LineRate=Options.lineRate;
linesPerField=Options.linesPerField;
accuracy=Options.targetAccuracy;
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
counter=0;
CenterZ3=CenterZ2;
while filterFail
    stop=0;
    counter=counter+1;
    if counter>50
        out=PosCur(fieldTimesCheck);
        disp('step limit reached')
        break;
    end    
    if ~isempty(FailPt)
        failField=find((fieldTimesSampPos-FailPt)>0,1,'first');
        V0=(CenterZ3(failField)-CenterZ3(failField-1))/fieldDur/1000;
        V1=(CenterZ3(failField+2)-CenterZ3(failField+1))/fieldDur/1000;
        dist=CenterZ3(failField+1)-CenterZ3(failField);
        traj=ZOptAna(maxAccelOrig,maxVelOrig,SampRate/1000,V0,V1,dist,0);
        traj=traj+CenterZ3(failField);
        nFields=ceil(length(traj)/fieldDurSamp)-1;
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
                FailPt=tooFast;
                stop=1;
            end
        end
    end
end
PosCur=cumsum([CenterZ2(1),VelCur((fieldDurSamp+1):(end-fieldDurSamp))]);
CenterZ3=CenterZ3(2:(end-1));
RealFieldsIdx=RealFieldsIdx(2:(end-1));

%Interpolate XY position of "fake" fields to maximize likelihood of extra
%data
NewCentersXY=interparc(find(RealFieldsIdx==0)./length(RealFieldsIdx),FieldCenters(:,1),FieldCenters(:,2),'linear');
NewFieldCenters=zeros(length(CenterZ3),3);
NewFieldCenters(:,3)=CenterZ3;
NewFieldCenters(RealFieldsIdx==1,1:2)=FieldCenters(:,1:2);
NewFieldCenters(RealFieldsIdx==0,1:2)=NewCentersXY(:,1:2);

%Add Flyback
V0=-(CenterZ3(end)-CenterZ3(end-1))/fieldDur/1000;
V1=-(CenterZ3(2)-CenterZ3(1))/fieldDur/1000;
dist=(PosCur(end)-PosCur(1));
traj=CenterZ3(end)-ZOptAna(maxAccelOrig,maxVelOrig,SampRate/1000,V0,V1,dist,0);
nFields=ceil(length(traj)/fieldDurSamp)-1;
traj=resample(traj,(nFields+1).*fieldDurSamp,length(traj));
traj=smooth(traj,round(SampRate/1000));
newFieldTimes=cumsum(ones(1,nFields)*fieldDurSamp);
newCenterZ=traj(newFieldTimes);
CenterZ3=[CenterZ3',newCenterZ']';
RealFieldsIdx=[RealFieldsIdx(2:(end-1)),zeros(1,nFields)];
fieldIntervals=[fieldIntervals(2:(end-1)),ones(1,nFields).*fieldDur];
fieldIntervalsSamp=fieldIntervals*SampRate;
fieldTimesSampPos=round(cumsum([fieldDurHalfSamp,fieldIntervalsSamp]));
fieldTimesCheck=fieldTimesSampPos(RealFieldsIdx==1);
PosCur=[PosCur,traj'];

FlybackFieldCenters=interparc([1:nFields]/nFields,FieldCenters([end,1],1),FieldCenters([end,1],2),'linear');
FlybackFieldCenters(:,3)=newCenterZ;
NewFieldCenters((end+1):(end+nFields),:)=FlybackFieldCenters;

%Output
Zwave=PosCur;
TrueFieldsMask=RealFieldsIdx;
TrueFieldsCenterSamp=fieldTimesCheck;
