%% fn_movie

%% Syntax
%  fn_movie(M.Y[,M.opt])
%  fn_movie(M.Y,option1,value1,...)

%% Description
%  available options are
%  - temporal manipulations:
%        tnorm       normalize by average M.fr (overlays with normal M.fr)
%        tspec       normalize by a local average (using frames t-tspec:t+tspec)
%        tbin        binning
%  - spatial manipulations:
%        xhigh       sigma for high-pass
%        xlow        sigma for low-pass
%        xbin        binning
%  - colors (M.fr denotes the M.fr at the instant an option is cM.hanged)
%        brightness  bias in the center of clipping range ((max(M.fr)-min(M.fr))*brightness)
%        contrast    amplitude of clipping range ((max(M.fr)-min(M.fr))/contrast)
%        cmap        choose colormap
%        overlay     overlay computed frames on top of the average M.fr in
%                    grayscale
%  - movie playing
%        start       M.fr from which to start playing
%        end         M.fr until which to play
%        speed       number of frames per second

%% Source
% Thomas Deneux
%
% Copyright 2008-2012
%
