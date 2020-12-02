function [bo] = compactFACETIIrunner(P1,P2,Beam0)
% Scripts runs beam through short lattice, accounting for needed manual
% changes

% Defines global variables
global BEAMLINE WF KLYSTRON PS

load FACET2_shortLattice

% Changes phases for parameter scan
BEAMLINE{1,1}.Phase = P1;
BEAMLINE{11,1}.Phase = P2;

% Tracks beam, correcting for change in average position
[~,bo]=TrackThru(1,10,Beam0,1,1);
bo.Bunch.x(5,:) = bo.Bunch.x(5,:) - mean(bo.Bunch.x(5,:));
bo.Bunch.x(isinf(bo.Bunch.x))= 0;
[~,bo]=TrackThru(11,20,bo,1,1);
bo.Bunch.x(5,:) = bo.Bunch.x(5,:) - mean(bo.Bunch.x(5,:));
bo.Bunch.x(isinf(bo.Bunch.x))= 0;


[~,bo]=TrackThru(21,40,bo,1,1);     % CHANGE TO length(BEAMLINE) when debugged
bo.Bunch.x(5,:) = bo.Bunch.x(5,:) - mean(bo.Bunch.x(5,:));
bo.Bunch.x(isinf(bo.Bunch.x))= 0;

% Gives data on beam
% set(0,'DefaultFigureVisible','off');
% beamdata = beamImage(bo);
% rmslength = beamdata.rmsz;
% peakcurrent = beamdata.pkI;
end
