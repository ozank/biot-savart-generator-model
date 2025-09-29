function [BSmag] = BSmag_init()
%---------------------------------------------------
%  NAME:      BSmag_init.m
%  WHAT:      Initializes a Biot-Savart magnetostatic analysis.
%  REQUIRED:  BSmag Toolbox 20150407
%  AUTHOR:    20150407, L. Queval (loic.queval@gmail.com)
%  COPYRIGHT: 2015, Loic Quéval, BSD License (http://opensource.org/licenses/BSD-3-Clause).
%
%  USE:
%    [BSmag] = BSmag_init()
%
%  INPUTS:
%
%  OUTPUTS:
%    BSmag   = Initialized BSmag data structure
%      BSmag.Nfilament      = Number of filaments
%---------------------------------------------------

BSmag.Nfilament = 0; %Number of source filament 

% Open default figure to plot source points and field points

%Only plot figures if plot_figures = 1 
%if the variable is not defined it the code will generate figures by
%default plot_figures can be set to 0 to reduce computation time

global plot_figures

if   exist('plot_figures','var') && ~plot_figures
    %Do not plot figures
  
else

    %Plot Figures
    figure(1), hold on, grid on, box on, axis equal
    xlabel('x [m]'), ylabel('y [m]'), zlabel('z [m]')
    view(3), axis tight
end


