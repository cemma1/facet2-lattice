function match_L3(Initial)
%MATCH_L3 Matching routines for L3
%match_L3(Initial)
%  Initial: Initial structure for start of FACET2e beamline
%  FACET2e BEAMLINE assumed to be pre-loaded
global BEAMLINE PS

% Assumed PS not assigned in base lattice
if ~isempty(PS)
  error('PS not empty')
end

% - Match into Sector 20
betax=3.17323882507; alphax=0.757843430229 ;
betay=4.001700149031; alphay=0.771049709497 ;
imatch = findcells(BEAMLINE,'Name','BEGBC20') ;
% qm={'Q19501', 'Q19601', 'Q19701', 'Q19801',  'Q19851', 'Q19871'} ;
qm={'Q19701', 'Q19801',  'Q19851', 'Q19871'} ;
M=Match;
M.beam=MakeBeam6DGauss(Initial,1e3,3,1);
M.iInitial=1;
M.initStruc=Initial;
M.verbose=false; % see optimizer output or not
M.optim='fminsearch';
M.optimDisplay='iter';
M.useParallel=false;
M.addMatch(imatch,'beta_x',betax,0.001);
M.addMatch(imatch,'beta_y',betay,0.001);
M.addMatch(imatch,'alpha_x',alphax,0.001);
M.addMatch(imatch,'alpha_y',alphay,0.001);
for iq=1:length(qm)
  iele=findcells(BEAMLINE,'Name',qm{iq});
  AssignToPS( iele, iq ) ;
  M.addVariable('PS',iq,'Ampl',0.1,10);
end
M.doMatch();
display(M);
for iq=1:length(qm)
  for iele=findcells(BEAMLINE,'Name',qm{iq})
    RenormalizePS(iq);
    BEAMLINE{iele}.PS=0;
  end
end
PS=[];
TwissPlot(1,length(BEAMLINE),Initial,[1 1 0]);
clight=2.99792458e8; % speed of light (m/sec)
Cb=1e9/clight;       % rigidity conversion (T-m/GeV)
for iq=1:length(qm)
  iele=findcells(BEAMLINE,'Name',qm{iq});
  fprintf('K%s := %g\n',BEAMLINE{iele(1)}.Name,BEAMLINE{iele(1)}.B/(BEAMLINE{iele(1)}.L*Cb*BEAMLINE{iele(1)}.P));
end
