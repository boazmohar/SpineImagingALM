%% what to run when adding a new experiment
%% 0: create a new folder in database and copy over:
% 1: expanded_new.tif
% 2: Sp.mat
%% 1: up sample to create top_vol_up.tif
upSampleAll();
%% 2: in imageJ saturate 0.1% of px and save as 8 bit: top_vol_up_8bit.tif
%% 3: in Neuromantic trace dendrites save as dendrite.swc
%% 4: prepare masks for FOV aligment
prepareMasksAll();
prepareMasksAllSWC();
%% 5: fix the Sp table to enable python to read it
fixTableAll();
%% 6: load the correct top_vol_up.mat and trace using:
dend = DendriteSeg(top_vol_up);
%% 7: load the dendrites from file-->load dendrite and trace spines
%% 8: save the mask using file-->save mask
%% 9: copy over the *_python.mat file over to the session folder and rename
% to: Run1mask.mat 
dend2 = DendriteSeg(filt_vol_up);