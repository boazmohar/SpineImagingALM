function c_PiezoOptIterative(obj)

SampRate=obj.OptionsStruct.lineRate;
accuracy=obj.OptionsStruct.targetAccuracy;
volPeriodAdj=obj.OptionsStruct.volPeriodAdj;
iter=obj.OptionsStruct.zOptIter;

cmdOffset=-0.0025;
um2volt=1/1000;
scaling=1;
Target=obj.Zwave*um2volt;
accuracy=accuracy*um2volt;
Target(Target<0)=0;
CmdWave=Target*scaling+cmdOffset;
%flatten period at end of waveform to model early end of waveform
volPeriodAdjSamp=round(volPeriodAdj./(10^6).*SampRate);
%CmdWave((end-volPeriodAdjSamp+1):end)=CmdWave(end-volPeriodAdjSamp);
CmdWave((end-volPeriodAdjSamp+1):end)=CmdWave(1);
CmdWaveHist(:,1)=CmdWave;
CmdWave=[CmdWave,CmdWave,CmdWave];
SampPts=obj.TrueFieldsCenterSamp;
%increase command points to prevent oscillations
linesPerField=obj.OptionsStruct.ImagingLines;
quartField=round(linesPerField/4);
if (isfield(obj.OptionsStruct, 'ZoptDampen') && obj.OptionsStruct.ZoptDampen) %&& 0 % replace with || 1 to use it
    SampPts=[SampPts, SampPts-quartField, SampPts+quartField];
    SampPts=sort(SampPts);
end
CmdPts=SampPts;



%flatten period at end of waveform to model early end of waveform
volPeriodAdjSamp=round(volPeriodAdj./(10^6).*SampRate);

curdir=pwd;
cd('C:\Vidrio\ScanImage2015')
import dabs.ni.daqmx.* 
% preallocation
maxError = zeros(1,iter);
AllampShift = cell(1,iter);
for i=1:iter
    dur=length(CmdWave)*2;
    x=ones(1,dur)*min(CmdWave);
    x2=zeros(size(x));
    x((round(dur/4)+1):(round(dur/4)+length(CmdWave)))=CmdWave;
    %x((round(dur*SampRate/2)+1):(round(dur*SampRate/2)+10))=0;
    x2(1:10)=5;
    outputData=[x' x2'];
    try
        hTaskAO = Task('TaskAO');
        hTaskAO.createAOVoltageChan('Dev1-SI_FastZ',0:1,[],-10,10);
        hTaskAO.writeAnalogData([cmdOffset 0]);
        hTaskAO.cfgSampClkTiming(SampRate,'DAQmx_Val_FiniteSamps',length(x),'OnboardClock');
        hTaskAI = Task('TaskAI');
        hTaskAI.createAIVoltageChan('Dev3-SI_Beams',[4:5],[],-10,10);
        hTaskAI.cfgSampClkTiming(SampRate,'DAQmx_Val_FiniteSamps',length(x),'OnboardClock');
        hTaskAI.cfgDigEdgeStartTrig('PFI3','DAQmx_Val_Rising');
        hTaskAI.start();
        hTaskAO.cfgOutputBuffer(length(x));
        hTaskAO.writeAnalogData(outputData, 5)
        hTaskAO.start();
        hTaskAO.waitUntilTaskDone(25);
        hTaskAI.waitUntilTaskDone;
        Praw = hTaskAI.readAnalogData([],'scaled',0); % read all available data
        
        Praw(:,1)=-(Praw(:,1)/10-0.0025);
        
        hTaskAO.stop();
        delete(hTaskAO);
        clear hTaskAO;
        
        hTaskAI.stop();
        delete(hTaskAI);
        clear hTaskAI;
        
        
    catch err % clean up task if error occurs
        if exist('hTaskAO','var')
            delete(hTaskAO);
            clear hTaskAO;
        end
        if exist('hTaskAI','var')
            delete(hTaskAI);
            clear hTaskAI;
        end
        rethrow(err);
    end
    
    Psens=Praw((round(dur/4)+1):(round(dur/4)+length(CmdWave)),1);
    [~, I]=max(xcorr([Target,Target,Target],Psens));
    phaseShift=I-length(CmdWave);
    
    midSens=Psens((length(Target)+1):(length(Target)*2));
    sensHist(:,i)=midSens;
    ampShift=Target(SampPts)-midSens(SampPts)';
    maxError(i)=max(abs(ampShift));
    AllampShift(i) = {abs(ampShift)};
    ampShift(abs(ampShift)<accuracy)=0;
    ampShift=ampShift./3;
    CmdZ=CmdWave(length(Target)+CmdPts)+ampShift*scaling;
    CmdPts=CmdPts+phaseShift;
    F=griddedInterpolant([CmdPts,CmdPts+length(Target),CmdPts+(2*length(Target))],...
        [CmdZ,CmdZ,CmdZ],'spline');
    CmdWave=F(1:(length(Target)*3));
    CmdWave=CmdWave((length(Target)+1):(length(Target)*2));
    CmdWave((end-volPeriodAdjSamp+1):end)=CmdWave(end-volPeriodAdjSamp);
    CmdWaveHist(:,i+1)=CmdWave;
    CmdWave=[CmdWave,CmdWave,CmdWave];
    CmdWave(CmdWave<0)=0;
    CmdWave(CmdWave>1)=1;
    
end

[maxError, BestIter]=min(maxError);
obj.maxErrorMicron=maxError/um2volt;
CmdWave=CmdWaveHist(:,BestIter);
out(:,1)=(CmdWave-cmdOffset)./scaling;
out(:,2)=Target;
out(:,3)=sensHist(:,BestIter);
CmdWave=CmdWave(1:(end-volPeriodAdjSamp));
global UserWaveform
UserWaveform=CmdWave;
obj.CmdWave = CmdWave;
axes(obj.Axes2Handle);
cla;
X = (1:length(Target))./obj.OptionsStruct.lineRate*1000;
X2 = [X;X;X]';
Error = sensHist(:,BestIter)- Target';
global ZError
ZError=Error/um2volt;
obj.ZError=ZError;
hAx = plotyy(X2,out/um2volt,SampPts/SampRate*1000,...
    Error(SampPts)/um2volt) ;
title('Iter Z Optimization')
xlabel('Time (ms)')
ylim(hAx(1),[0 obj.OptionsStruct.ZmaxTravel]);
% ylabel(hAx(1),'Waveform(V)') % left y-axis
ylabel(hAx(2),'Error(\mum)') % right y-axis
l=legend('Command','Target','Sensor','Diff \mum');
l.Location = 'best';
cd(curdir)