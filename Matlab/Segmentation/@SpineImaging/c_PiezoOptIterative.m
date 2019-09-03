function c_PiezoOptIterative(obj)

global fit_um_to_V
if (isempty(fit_um_to_V))
    data = csvread('C:\Vidrio\ScanImage2016\calibZ_20161212.csv');
    x = data(:,2);
    y = data(:,1);
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Lower = [-Inf -Inf -Inf 0];
    opts.Upper = [Inf Inf Inf 0];
    ft = fittype( 'poly3' );
    [fit_um_to_V, ~] = fit( x, y, ft, opts);
end

global fit_V_to_um
if (isempty(fit_V_to_um))
    data = csvread('C:\Vidrio\ScanImage2016\calibZ_20161212.csv');
    y = data(:,2);
    x = data(:,1);
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Lower = [-Inf -Inf -Inf 0];
    opts.Upper = [Inf Inf Inf 0];
    ft = fittype( 'poly3' );
    [fit_V_to_um, ~] = fit( x, y, ft, opts);
end

global fit_Sensor_V_to_um
if (isempty(fit_Sensor_V_to_um))
    data = csvread('C:\Vidrio\ScanImage2016\Calib_Sensor_20161212_v2_rev.csv');
    x = data(:,2);
    y = data(:,1);
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Lower = [-Inf -Inf -Inf 0];
    opts.Upper = [Inf Inf Inf 0];
    ft = fittype( 'poly3' );
    [fit_Sensor_V_to_um, ~] = fit( x, y, ft, opts);

end
            
SampRate=obj.OptionsStruct.lineRate;
accuracy=obj.OptionsStruct.targetAccuracy;
% accuracy=fit_V_to_um(accuracy);
volPeriodAdj=obj.OptionsStruct.volPeriodAdj;
iter=obj.OptionsStruct.zOptIter;
rep = obj.OptionsStruct.zOptRep;
cmdOffset=0;
% um2volt=1/obj.OptionsStruct.ZmaxTravel;
scaling=1;
Target_um = obj.Zwave;
Target=fit_um_to_V(Target_um-10);

% figure(111)
% clf
% x = 1:length(obj.Zwave);
% plot(x,Target_um)
% hold on
% plot(Target)
% Target=fit_V_to_um(obj.Zwave*-1);

% Target=obj.Zwave*um2volt;
% accuracy=accuracy*um2volt;
%Target(Target<0)=0;
CmdWave=Target*scaling+cmdOffset;
%flatten period at end of waveform to model early end of waveform
volPeriodAdjSamp=round(volPeriodAdj./(10^6).*SampRate);
%CmdWave((end-volPeriodAdjSamp+1):end)=CmdWave(end-volPeriodAdjSamp);
CmdWave((end-volPeriodAdjSamp+1):end)=CmdWave(1);
CmdWaveHist(:,1)=CmdWave;
CmdWave=repmat(CmdWave, rep,1);
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

% curdir=pwd;
% cd('C:\Vidrio\ScanImage2015')
import dabs.ni.daqmx.* 
% preallocation
maxError = zeros(1,iter);
AllampShift = cell(1,iter);
for i=1:iter
    figure()
    plot(Target_um)
    hold on
    CmdPts=SampPts;
    buff = round(0.05*SampRate);
    x = [ones(1,buff)*min(CmdWave(1)) CmdWave' ones(1,buff)*min(CmdWave(end))];
    x2=zeros(size(x));
    x2(1:10)=5;
    outputData=[x' x2'];
    try
        hTaskAO = Task('TaskAO');
        hTaskAO.createAOVoltageChan('Dev1-SI_FastZ',0:1,[],-10,10);
        %   terminalConfig: (OPTIONAL) One of {'DAQmx_Val_Cfg_Default', 'DAQmx_Val_RSE', 'DAQmx_Val_NRSE', 'DAQmx_Val_Diff', 'DAQmx_Val_PseudoDiff'}. Specifies the input terminal configuration for the channel.
        % If omitted/blank, 'DAQmx_Val_Cfg_Default' is used, NI-DAQmx to choose the default terminal configuration for the channel.
        hTaskAO.writeAnalogData([cmdOffset 0]);
        hTaskAO.cfgSampClkTiming(SampRate,'DAQmx_Val_FiniteSamps',length(x),'OnboardClock');
        hTaskAI = Task('TaskAI');
        hTaskAI.createAIVoltageChan('Dev3-SI_Beams',[4],[],-1,1, 'DAQmx_Val_Volts','','DAQmx_Val_NRSE');
        hTaskAI.createAIVoltageChan('Dev3-SI_Beams',[5],[],-10,10);
        hTaskAI.cfgSampClkTiming(SampRate,'DAQmx_Val_FiniteSamps',length(x),'OnboardClock');
        hTaskAI.cfgDigEdgeStartTrig('PFI3','DAQmx_Val_Rising');
        hTaskAI.start();
        hTaskAO.cfgOutputBuffer(length(x));
        hTaskAO.writeAnalogData(outputData, 5)
        hTaskAO.start();
        hTaskAO.waitUntilTaskDone(25);
        hTaskAI.waitUntilTaskDone;
        Praw = hTaskAI.readAnalogData([],'scaled',0); % read all available data
        
        Praw(:,1)=Praw(:,1)*-1;
        
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
    
    Psens_V=Praw(buff:end-buff-1,1);
    %sensor V --> um
    Psens_um =fit_Sensor_V_to_um(Psens_V)*-1;
%     Psens3 = fit_V_to_um(Psens2)*-1;
    Pmat = reshape(Psens_um, length(Target_um), rep);
    avgTarget = mean(Pmat,2);
    plot(avgTarget)
    [~, I]=max(xcorr(Target_um,avgTarget));
    phaseShift=I-length(Target_um);
    
%     midSens=Psens((length(Target)+1):(length(Target)*2));
    midSens = avgTarget;
    sensHist(:,i)=midSens;
    ampShift=Target_um(SampPts)'-midSens(SampPts);
    maxError(i)=max(abs(ampShift));
    AllampShift(i) = {abs(ampShift)};
    ampShift(abs(ampShift)<accuracy)=0;
    ampShift=ampShift./2;
    if (size(CmdWave,1)==1)
        CmdWave=CmdWave';
    end
    CmdPts=CmdPts+phaseShift;
    CmdZ=fit_V_to_um(CmdWave(CmdPts))+ampShift.*scaling;
   
    points = [];
    for j=0:rep-1
        points = [points CmdPts+j*length(Target)];
    end
    CmdZAll = repmat(CmdZ, rep,1);
    F=griddedInterpolant(points,CmdZAll,'spline');
    CmdWave=F(1:(length(Target)*rep));
    CmdWave=CmdWave((length(Target)+1):(length(Target)*2));
%     CmdWave((end-volPeriodAdjSamp+1):end)=CmdWave(end-volPeriodAdjSamp);
    CmdWaveHist(:,i+1)=CmdWave;
    plot(CmdWave);
    plot(CmdPts, CmdWave(CmdPts),'.');
    legend('Target','sensor','newcommand','points')
    title(maxError(i))
    pause()
    CmdWave=repmat(CmdWave,1,rep);
    CmdWave = fit_um_to_V(CmdWave);
    CmdWave(CmdWave<-10)=-10;
    CmdWave(CmdWave>10)=10;
    
end
[maxError2, BestIter]=min(maxError);
obj.maxErrorMicron=maxError2;
CmdWave=CmdWaveHist(:,BestIter);
out(:,1)=(CmdWave-cmdOffset)./scaling;
out(:,2)=Target;
out(:,3)=sensHist(:,BestIter);
CmdWave=CmdWave(1:(end-volPeriodAdjSamp));
% CmdWave=CmdWave*-1;
global UserWaveform
UserWaveform=fit_um_to_V(CmdWave);
obj.CmdWave = CmdWave;
figure(111)
plot(CmdWave)
a = sensHist(:,BestIter);
plot(a)
plot(SampPts, a(SampPts),'*')
% figure; plot(hSI.scannerAO.ao_volts.G)
legend('ZWave','Target', 'CmdWave', 'sensHist', 'fields')
axes(obj.Axes2Handle);
cla;
X = (1:length(Target))./obj.OptionsStruct.lineRate*1000;
X2 = [X;X;X]';
Error = sensHist(:,BestIter)- Target;
global ZError
ZError=fit_um_to_V(Error);
obj.ZError=ZError;
hAx = plotyy(X2,out,SampPts/SampRate*1000,...
   fit_um_to_V( Error(SampPts))) ;
title('Iter Z Optimization')
xlabel('Time (ms)')
ylim(hAx(1),[0 obj.OptionsStruct.ZmaxTravel]);
% ylabel(hAx(1),'Waveform(V)') % left y-axis
ylabel(hAx(2),'Error(\mum)') % right y-axis
l=legend('Command','Target','Sensor','Diff \mum');
l.Location = 'best';
% cd(curdir)