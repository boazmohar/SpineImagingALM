function [CmdWave,maxErrorMicron]=PiezoOptIterative(waveform,RefPts,Options)

SampRate=Options.lineRate;
accuracy=Options.targetAccuracy;
volPeriodAdj=Options.volPeriodAdj;
iter=Options.zOptIter;

cmdOffset=1.4286;
um2volt=10/144;
scaling=10/14;
Target=waveform*um2volt;
accuracy=accuracy*um2volt;
Target(Target<0)=0;
CmdWave=Target*scaling+cmdOffset;
%flatten period at end of waveform to model early end of waveform
volPeriodAdjSamp=round(volPeriodAdj./(10^6).*SampRate);
CmdWave((end-volPeriodAdjSamp+1):end)=CmdWave(end-volPeriodAdjSamp);
CmdWaveHist(:,1)=CmdWave;
CmdWave=[CmdWave,CmdWave,CmdWave];
CmdPts=RefPts;

%flatten period at end of waveform to model early end of waveform
volPeriodAdjSamp=round(volPeriodAdj./(10^6).*SampRate);

curdir=pwd;
cd('C:\Vidrio\ScanImage2015')
import dabs.ni.daqmx.* 
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
        hTaskAO.waitUntilTaskDone(10);
        hTaskAI.waitUntilTaskDone;
        Praw = hTaskAI.readAnalogData([],'scaled',0); % read all available data
        
        Praw(:,1)=Praw(:,1)*1.998;
        
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
    [r I]=max(xcorr([Target,Target,Target],Psens));
    phaseShift=I-length(CmdWave);
    
    midSens=Psens((length(Target)+1):(length(Target)*2));
    sensHist(:,i)=midSens;
    ampShift=Target(RefPts)-midSens(RefPts)';
    maxError(i)=max(abs(ampShift));
    
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
    CmdWave(CmdWave>10)=10;
    
end

[maxError, BestIter]=min(maxError);
maxErrorMicron=maxError/um2volt;
CmdWave=CmdWaveHist(:,BestIter);
out(:,1)=(CmdWave-cmdOffset)./scaling;
out(:,2)=Target;
out(:,3)=sensHist(:,BestIter);
CmdWave=CmdWave(1:(end-volPeriodAdjSamp));
figure;plot(out);
cd(curdir)