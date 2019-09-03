function c_SetPower(obj)
% Sets power for flyback and pading frames to 0
% has the option of increasing power with increasnig z
% this function assumes that the firs scan field is the top one!
Power = obj.OptionsStruct.Power;
Scatter = obj.OptionsStruct.ZPowerAdjLength;
if  ~isempty(obj.OptionsStruct.ZPowerAdjLength)
    % z depth of all field (excepth the first) in relation to the first one
    Zdepth = obj.All.z(2:end) - min(obj.All.z);
    %             first     second and onwards
    obj.Powers = [Power;    Power*exp(Zdepth./Scatter)];
else
    obj.Powers = ones(length(obj.All.x))*Power;
end
% zero for false fields
obj.Powers(~obj.TrueFieldsMask) = 0;