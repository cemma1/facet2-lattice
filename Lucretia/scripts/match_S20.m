function match_S20(Initial)
%MATCH_S20 Matching routines for Sector 20
%match_S20(Initial)
%  Initial: Initial structure for start of FACET2e beamline
%  FACET2e BEAMLINE assumed to be pre-loaded
global BEAMLINE PS

% Assumed PS not assigned in base lattice
if ~isempty(PS)
  error('PS not empty')
end

% Main IP beta function (round)
ipbeta=0.05;

% Required beamline indices
pent=findcells(BEAMLINE,'Name','PENT');
iscr=findcells(BEAMLINE,'Name','PDUMP');

% Match IP
if ipbeta>0.5
  lval=[-44 -44 -44 0 0]; uval=[0 44 0 44 20.3];
  optim='fminsearch';
else
  lval=[-44 -44 -44 -44 -20.3]; uval=[44 44 44 44 20.3];
  optim='lsqnonlin';
end
M=Match;
qm={'Q0FF' 'Q1FF' 'Q2FF' 'Q3FF' 'Q4FF' 'Q5FF'};
iele=[findcells(BEAMLINE,'Name','Q0FF') findcells(BEAMLINE,'Name','Q2FF')];
AssignToPS( iele, 1 ) ;
M.addVariable('PS',1,'Ampl',lval(1),uval(1)); 
for ips=[2 4 5 6]
  iele=findcells(BEAMLINE,'Name',qm{ips});
  AssignToPS( iele, length(PS)+1 ) ;
  M.addVariable('PS',length(PS),'Ampl',lval(length(PS)),uval(length(PS))); 
end
MovePhysicsVarsToPS(1:length(PS));
M.beam=MakeBeam6DGauss(Initial,1e3,3,1);
M.iInitial=1;
M.initStruc=Initial;
M.verbose=false; % see optimizer output or not
M.optim=optim;
M.optimDisplay='iter';
M.addMatch(pent,'alpha_x',0,0.0001);
M.addMatch(pent,'alpha_y',0,0.0001);
M.addMatch(pent,'beta_x',ipbeta,0.0001);
M.addMatch(pent,'beta_y',ipbeta,0.0001);
M.doMatch();
disp(M);

% Match dump optics from plasma exit or IP
qd=[findcells(BEAMLINE,'Name','Q0D') findcells(BEAMLINE,'Name','Q2D')];
qf=findcells(BEAMLINE,'Name','Q1D');
for iele=qf;BEAMLINE{iele}.B=1; end
for iele=qd;BEAMLINE{iele}.B=-1; end
AssignToPS(qf,length(PS)+1); psqf=length(PS);
AssignToPS(qd,length(PS)+1); psqd=length(PS);
MovePhysicsVarsToPS([psqf psqd]);
qmax=100;
M=Match;
M.beam=MakeBeam6DGauss(Initial,1e3,3,1);
M.iInitial=1;
M.initStruc=Initial;
M.verbose=false; % see optimizer output or not
M.optim='lsqnonlin';
M.optimDisplay='iter';
M.addMatch(iscr,'R',0,1e-9,'12');
M.addMatch(iscr,'R',0,1e-9,'34');
M.addVariable('PS', psqf,'Ampl',0,qmax*2);
M.addVariable('PS', psqd,'Ampl',-qmax*2,0);
M.doMatch();
disp(M);

% Display matched quads values
qm=[qm {'Q0D' 'Q1D' 'Q2D'}];
for ips=1:length(PS)
  for iele=PS(ips).Element
    RenormalizePS(ips);
    BEAMLINE{iele}.PS=0;
  end
end
PS=[];
TwissPlot(1,length(BEAMLINE),Initial,[1 1 0]);
clight=2.99792458e8; % speed of light (m/sec)
Cb=1e9/clight;       % rigidity conversion (T-m/GeV)
for ips=1:length(qm)
  iele=findcells(BEAMLINE,'Name',qm{ips});
  fprintf('K%s := %g\n',BEAMLINE{iele(1)}.Name,BEAMLINE{iele(1)}.B/(BEAMLINE{iele(1)}.L*Cb*BEAMLINE{iele(1)}.P));
end

