function compactFACETIIgenerator
% Make a compact form of FACET-II lattice
global BEAMLINE WF KLYSTRON PS

% Simulation switches
dosr=false; % Switch on ISR & CSR effects in bends
mb_dec=1; % macro-particle gen decimation factor - set this >1 to reduce number of tracked particles (original #=1e6) (e.g. mb_dec=1000 would track 1e3 particles)
% Load design model
load FACET2e_baseline BEAMLINE bstore WF Twiss Initial KLYSTRON PS
% WF(1) = s-band WF(2) = x-band
WFx=WF; WFx.ZSR(1)=[]; WFx.TSR(1)=[];
WF.ZSR(2)=[]; WF.TSR(2)=[];
WFs=WF;

% Store some required response matrices from original design
i1=findcells(BEAMLINE,'Name','BEGL1F');  i2=findcells(BEAMLINE,'Name','BCX11314A')-1; % region from start of L1 to first bend in BC11 compression chicane
[~,R_L1]=RmatAtoB(i1,i2); % response matrix through L1 accelerating section
i12=findcells(BEAMLINE,'Name','BCX11355B')+1;  i22=findcells(BEAMLINE,'Name','BCX14720A')-1; % region from start of L2 to first bend in BC14 compression chicane
[~,R_L2]=RmatAtoB(i12,i22); % response matrix through L2 accelerating section
i13=i22+1;  i23=findcells(BEAMLINE,'Name','ENDL3F'); % region from start of L2 to first bend in BC14 compression chicane
[~,R_L3]=RmatAtoB(i13,i23); % response matrix through L2 accelerating section
%[~,Beam_L1]=TrackThru(i1,i2,bstore.L1,1,1); % Beam at start of L1
%[~,Beam_BC11]=TrackThru(i1,findcells(BEAMLINE,'Name','BCX11355B'),bstore.L1,1,1); % Beam at exit of last bend in BC11 chicane
%[~,Beam_BC14]=TrackThru(findcells(BEAMLINE,'Name','BCX11355B')+1,findcells(BEAMLINE,'Name','BCX14883B'),Beam_BC11,1,1); % Beam at exit of last bend in BC14 chicane
I = TwissToInitial(Twiss,i1,Initial); % Lucretia Initial structure (including Twiss parameters) at L1 entrance
clear global PS KLYSTRON


% FOR TESTING
%[~,COMPTEST] = TrackThru(i1,findcells(BEAMLINE,'Name','BCX14720A')-1,bstore.L1,1,1);
%[~,COMPTEST2] = TrackThru(i1,findcells(BEAMLINE,'Name','ENDL3F'),bstore.L1,1,1);
%[~,COMPTEST3] = TrackThru(i1,1469,bstore.L1,1,1);
% END TESTING


% Copy original BEAMLINE design and start with new one
BL1=BEAMLINE;
BEAMLINE={};

% Replace L1 with single structure
i1=findcells(BL1,'Name','BEGL1F'); i2=findcells(BL1,'Name','BCX11314A')-1; % beamline from start of L1 to first bend of BC11 chicane
L=0;
V=0;
Egain=0;
for iele=findcells(BL1,'Class','LCAV',i1,i2)
  if BL1{iele}.Freq==2856
    L=L+BL1{iele}.L;
    V=V+BL1{iele}.Volt;
    aper=BL1{iele}.aper;
    Egain=Egain+BL1{iele}.Egain;
    phi=BL1{iele}.Phase;                              % Phase Parameter to be optimized
    Kloss=BL1{iele}.Kloss;
  end
end
% Structure with 1m length
BEAMLINE{1}=RFStruc(1,V,phi,2856,0,0,0,aper,'L1');
BEAMLINE{1}.Egain=Egain;
BEAMLINE{1}.Kloss=Kloss*L;
BEAMLINE{1}.Phase=phi;
BEAMLINE{1}.P=BL1{i1}.P;
BEAMLINE{1}.S=0;

% Scale wakefield to total linac length
WF.ZSR(1).K=WF.ZSR(1).K*L;
WF.TSR(1).K=WF.TSR(1).K*L;
BEAMLINE{1}.Wakes=[1 1]; 
BEAMLINE{1}.TrackFlag.SRWF_Z=1; % Switch on just longitudinal wake component
BEAMLINE{1}.TrackFlag.Aper=1; % Switch on checking of aperture excursions for tracked particles

% Add a TMAP element to transversely match the beam to the entrance of BC11
BEAMLINE{2,1} = TMapStruc('L1Match',0);
BEAMLINE{2}.S=BEAMLINE{1}.S+BEAMLINE{1}.L;
BEAMLINE{2}.P=BL1{i2}.P;
[~,R]=RmatAtoB(1,1); % response matrix through our lattice
R2 = R_L1(1:4,1:4) / R(1:4,1:4) ; % Required TMAT R matrix coefficients to "Match" transverse optics
BEAMLINE{2}.R(1:4,1:4) = R2(1:4,1:4) ; % Output of RmatAtoB(1,2) should now match R_L1

% Make the starting beam @ beginning of L1
Beam0=bstore.L1;
% - Reduce number of macro-particles
origlen=length(Beam0.Bunch.Q);
Beam0.Bunch.x=Beam0.Bunch.x(:,1:mb_dec:end);
Beam0.Bunch.stop=Beam0.Bunch.stop(1:mb_dec:end);
Beam0.Bunch.Q=Beam0.Bunch.Q(1:mb_dec:end).*(origlen/ceil(origlen/mb_dec));

% Track a beam through to end of L1 section and compare to original model
% Expect small differences in transverse plane (e.g. from space charge not being modeled), should be good to a few % though
%[~,bo]=TrackThru(1,2,Beam0,1,1);
%beamImage(bo);
%beamImage(Beam_L1);

% Now add BC11 bunch compressor bends
ibc1=findcells(BL1,'Name','BCX11314A'); ibc2=findcells(BL1,'Name','BCX11355B')+1;
ibend=findcells(BL1,'Class','SBEN',ibc1,ibc2); % there are only 4 bends, but 8 indices as each bend split in 2 pieces in model, we don't need that here so we will put them back together
Lb=BL1{ibend(1)}.L*2;
B=BL1{ibend(1)}.B*2;
ANG=BL1{ibend(1)}.Angle*2;
FINT=BL1{ibend(1)}.FINT + BL1{ibend(2)}.FINT ;
% We also need to add drift elements between bends to have the optics work out correctly - get those required distances now from the original model
L0 = BL1{ibend(1)}.S - BL1{ibc1}.S ;
L1 = BL1{ibend(3)}.S - BL1{ibend(2)+1}.S;
L2 = BL1{ibend(5)}.S - BL1{ibend(4)+1}.S;
L3 = BL1{ibend(7)}.S - BL1{ibend(6)+1}.S;
% Layout the BC11 elements now
BEAMLINE{3}=DrifStruc(L0,'D0BC11'); BEAMLINE{3}.P=BL1{ibc1}.P;
BEAMLINE{4}=SBendStruc(Lb, B, ANG, [0 ANG], 0, BL1{ibend(1)}.HGAP, FINT, 0, 'BCX11314') ; BEAMLINE{4}.P=BL1{ibend(1)}.P;
BEAMLINE{5}=DrifStruc(L1,'D1BC11'); BEAMLINE{5}.P=BL1{ibend(2)+1}.P;
BEAMLINE{6}=SBendStruc(Lb, -B, -ANG, [-ANG 0], 0, BL1{ibend(1)}.HGAP, FINT, 0, 'BCX11331') ; BEAMLINE{6}.P=BL1{ibend(3)}.P;
BEAMLINE{7}=DrifStruc(L2,'D2BC11'); BEAMLINE{7}.P=BL1{ibend(4)+1}.P;
BEAMLINE{8}=SBendStruc(Lb, -B, -ANG, [-ANG 0], 0, BL1{ibend(1)}.HGAP, FINT, 0, 'BCX11338') ; BEAMLINE{8}.P=BL1{ibend(5)}.P;
BEAMLINE{9}=DrifStruc(L3,'D3BC11'); BEAMLINE{9}.P=BL1{ibend(6)+1}.P;
BEAMLINE{10}=SBendStruc(Lb, B, ANG, [0 ANG], 0, BL1{ibend(1)}.HGAP, FINT, 0, 'BCX11355') ; BEAMLINE{10}.P=BL1{ibend(7)}.P;
SetSPositions(1,length(BEAMLINE),0);
% Switch on ISR and CSR effects in bends?
for iele=findcells(BEAMLINE,'Class','SBEN')
  if dosr
    BEAMLINE{iele}.TrackFlag.SynRad=2;
    BEAMLINE{iele}.TrackFlag.CSR=-1;
    BEAMLINE{iele}.TrackFlag.CSR_SmoothFactor=1;
    BEAMLINE{iele}.TrackFlag.CSR_DriftSplit=25;
    BEAMLINE{iele}.TrackFlag.Split=25;
  else
    BEAMLINE{iele}.TrackFlag.SynRad=0;
    BEAMLINE{iele}.TrackFlag.CSR=0;
    BEAMLINE{iele}.TrackFlag.Split=0;
  end
end

% We can check the linear Twiss parameters propogate correctly through our lattice so far...
[~,T]=GetTwiss(1,10,I.x.Twiss,I.y.Twiss);
fprintf('@ ENDBC11: betax = %g (%g) alphax = %g (%g) betay = %g (%g) alphay = %g (%g)\n',T.betax(end),Twiss.betax(ibc2),T.alphax(end),Twiss.alphax(ibc2),T.betay(end),Twiss.betay(ibc2),T.alphay(end),Twiss.alphay(ibc2))
% - there will be small differences, but should match at the few % level for beta-functions and <0.1 differences for alphas

% Now check particle tracking and compare to exit of BC11 in the original model
% [~,bo]=TrackThru(1,length(BEAMLINE),Beam0,1,1);
% beamImage(bo);
%beamImage(Beam_BC11);

% Here, add L2 linac section
i1=findcells(BL1,'Name','BEGL2F'); i2=findcells(BL1,'Name','BCX14720A')-1; % beamline from start of L2 to first bend of BC14 chicane
L=0;
V=0;
Egain=0;
for iele=findcells(BL1,'Class','LCAV',i1,i2)
  if BL1{iele}.Freq==2856
    L=L+BL1{iele}.L;
    V=V+BL1{iele}.Volt;
    aper=BL1{iele}.aper;
    Egain=Egain+BL1{iele}.Egain;
    phi=BL1{iele}.Phase;                              % Phase Parameter to be optimized
    Kloss=BL1{iele}.Kloss;
  end
end
% L2 Structure with 1m length
BEAMLINE{end+1}=RFStruc(1,V-178,phi,2856,0,0,0,aper,'L2');
BEAMLINE{end}.Egain=Egain;
BEAMLINE{end}.Kloss=Kloss*L;
BEAMLINE{end}.Phase=phi;
BEAMLINE{end}.P=BL1{i1}.P;
SetSPositions(1,length(BEAMLINE),0);

% Scale wakefield to total linac length
WF.ZSR(2)=WF.ZSR(1); WF.TSR(2)=WF.TSR(1);
WF.ZSR(2).K=WFs.ZSR.K*L;
WF.TSR(2).K=WFs.TSR.K*L;
BEAMLINE{end}.Wakes=[2 2]; 
BEAMLINE{end}.TrackFlag.SRWF_Z=1; % Switch on just longitudinal wake component
BEAMLINE{end}.TrackFlag.Aper=1; % Switch on checking of aperture excursions for tracked particles

% FOR TESTING
%[~,bo]=TrackThru(1,length(BEAMLINE),Beam0,1,1);
%beamImage(bo);
%beamImage(COMPTEST);
% END TESTING

% RAFI SECTION
% Add a TMAP element to transversely match the beam to the entrance of BC14
BEAMLINE{12,1} = TMapStruc('L2Match',0);
BEAMLINE{12}.S=BEAMLINE{11}.S+BEAMLINE{11}.L;
BEAMLINE{12}.P=BL1{i2}.P;
[~,R]=RmatAtoB(11,11); % response matrix through our lattice
R3 = R_L2(1:4,1:4) / R(1:4,1:4) ; % Required TMAT R matrix coefficients to "Match" transverse optics
BEAMLINE{12}.R(1:4,1:4) = R3 ; % Output of RmatAtoB(11,12) should now match R_L2

% Now add BC14 bunch compressor bends
ibc1=findcells(BL1,'Name','BCX14720A'); ibc2=findcells(BL1,'Name','BCX14883B');
ibend=findcells(BL1,'Class','SBEN',ibc1,ibc2); % there are only 4 bends, but 8 indices as each bend split in 2 pieces in model, we don't need that here so we will put them back together
Lb=BL1{ibend(1)}.L*2;
B=BL1{ibend(1)}.B*2;
ANG=BL1{ibend(1)}.Angle*2;
FINT=BL1{ibend(1)}.FINT + BL1{ibend(2)}.FINT ;
% We also need to add drift elements between bends to have the optics work out correctly - get those required distances now from the original model
L0 = BL1{ibend(1)}.S - BL1{ibc1}.S ;
L1 = BL1{ibend(3)}.S - BL1{ibend(2)+1}.S;
L2 = BL1{ibend(5)}.S - BL1{ibend(4)+1}.S;
L3 = BL1{ibend(7)}.S - BL1{ibend(6)+1}.S;
% Layout the BC14 elements now
BEAMLINE{13}=DrifStruc(L0,'D0BC14'); BEAMLINE{13}.P=BL1{ibc1}.P;
BEAMLINE{14}=SBendStruc(Lb, B, ANG, [0 ANG], 0, BL1{ibend(1)}.HGAP, FINT, 0, 'BCX14720') ; BEAMLINE{14}.P=BL1{ibend(1)}.P;
BEAMLINE{15}=DrifStruc(L1,'D1BC14'); BEAMLINE{15}.P=BL1{ibend(2)+1}.P;
BEAMLINE{16}=SBendStruc(Lb, -B, -ANG, [-ANG 0], 0, BL1{ibend(1)}.HGAP, FINT, 0, 'BCX14796') ; BEAMLINE{16}.P=BL1{ibend(3)}.P;
BEAMLINE{17}=DrifStruc(L2,'D2BC14'); BEAMLINE{17}.P=BL1{ibend(4)+1}.P;
BEAMLINE{18}=SBendStruc(Lb, -B, -ANG, [-ANG 0], 0, BL1{ibend(1)}.HGAP, FINT, 0, 'BCX14808') ; BEAMLINE{18}.P=BL1{ibend(5)}.P;
BEAMLINE{19}=DrifStruc(L3,'D3BC14'); BEAMLINE{19}.P=BL1{ibend(6)+1}.P;
BEAMLINE{20}=SBendStruc(Lb, B, ANG, [0 ANG], 0, BL1{ibend(1)}.HGAP, FINT, 0, 'BCX14883') ; BEAMLINE{20}.P=BL1{ibend(7)}.P;
SetSPositions(1,length(BEAMLINE),0);
% Switch on ISR and CSR effects in bends?
for iele=findcells(BEAMLINE,'Class','SBEN')
  if dosr
    BEAMLINE{iele}.TrackFlag.SynRad=2;
    BEAMLINE{iele}.TrackFlag.CSR=-1;
    BEAMLINE{iele}.TrackFlag.CSR_SmoothFactor=1;
    BEAMLINE{iele}.TrackFlag.CSR_DriftSplit=25;
    BEAMLINE{iele}.TrackFlag.Split=25;
  else
    BEAMLINE{iele}.TrackFlag.SynRad=0;
    BEAMLINE{iele}.TrackFlag.CSR=0;
    BEAMLINE{iele}.TrackFlag.Split=0;
  end
end

% Now check particle tracking and compare to exit of BC14 in the original model
% [~,bo]=TrackThru(1,length(BEAMLINE),Beam0,1,1);
% beamImage(bo);
%beamImage(Beam_BC14);

% Here, add L3 linac section
i1=findcells(BL1,'Name','BEGL3F'); i2=findcells(BL1,'Name','BEGBC20*'); % beamline from start of L2 to first bend of BC20 chicane
L=0;
V=0;
Egain=0;
for iele=findcells(BL1,'Class','LCAV',i1,i2)
  if BL1{iele}.Freq==2856
    L=L+BL1{iele}.L;
    V=V+BL1{iele}.Volt;
    aper=BL1{iele}.aper;
    Egain=Egain+BL1{iele}.Egain;
    phi=BL1{iele}.Phase;
    Kloss=BL1{iele}.Kloss;
  end
end
% L3 Structure with 1m length
BEAMLINE{end+1}=RFStruc(1,V-423,phi,2856,0,0,0,aper,'L3');
BEAMLINE{end}.Egain=Egain;
BEAMLINE{end}.Kloss=Kloss*L;
BEAMLINE{end}.Phase=0;
BEAMLINE{end}.P=BL1{i1}.P;
SetSPositions(1,length(BEAMLINE),0);

% Scale wakefield to total linac length
WF.ZSR(3)=WF.ZSR(1); WF.TSR(3)=WF.TSR(1);
WF.ZSR(3).K=WFs.ZSR.K*L;
WF.TSR(3).K=WFs.TSR.K*L;
BEAMLINE{end}.Wakes=[3 3]; 
BEAMLINE{end}.TrackFlag.SRWF_Z=1; % Switch on just longitudinal wake component
BEAMLINE{end}.TrackFlag.Aper=1; % Switch on checking of aperture excursions for tracked particles

% TEST SECTION
% [~,bo]=TrackThru(1,length(BEAMLINE),Beam0,1,1);
% beamImage(bo);
% beamImage(COMPTEST2);
% END TEST SECTION

% Add a TMAP element to transversely match the beam to the entrance of BC20
BEAMLINE{22,1} = TMapStruc('L3Match',0);
BEAMLINE{22}.S=BEAMLINE{21}.S+BEAMLINE{21}.L;
BEAMLINE{22}.P=BL1{i2}.P;
[~,R]=RmatAtoB(21,21); % response matrix through our lattice
R4 = R_L3(1:4,1:4) / R(1:4,1:4) ; % Required TMAT R matrix coefficients to "Match" transverse optics
BEAMLINE{22}.R(1:4,1:4) = R4 ; % Output of RmatAtoB(21,22) should now match R_L3

% Now add BC20, this one is more complicated than BC11 and BC14, just add all magnetic elements (bends, quadrupoles & sextupoles) and drifts from BEGBC20 to ENDBC20
ldBC20E=load('BC20.mat','BEAMLINE','Initial');
for ibl=1:length(ldBC20E.BEAMLINE)
  if isfield(ldBC20E.BEAMLINE{ibl},'PS')
    ldBC20E.BEAMLINE{ibl}.PS=0;
  end
end
BEAMLINE=[BEAMLINE(:);ldBC20E.BEAMLINE];
wigtrim=findcells(BEAMLINE,'Name','YCWIGE'); for iele=wigtrim; BEAMLINE{iele}.B=0; end
SetSPositions(1,length(BEAMLINE),0);

% - Match into optics
iele=findcells(BEAMLINE,'Name','L3Match');
fminsearch(@(x) S20match(x,'x',iele,I,ldBC20E.Initial.x.Twiss.beta,ldBC20E.Initial.x.Twiss.alpha),BEAMLINE{iele}.R(1:2,1:2),optimset('Display','iter'));
fminsearch(@(x) S20match(x,'y',iele,I,ldBC20E.Initial.y.Twiss.beta,ldBC20E.Initial.y.Twiss.alpha),BEAMLINE{iele}.R(3:4,3:4));

% Switch on ISR and CSR effects in bends, quadrupoles, and sextupoles?
for iele=findcells(BEAMLINE,'Class','SBEN')
  if dosr
    BEAMLINE{iele}.TrackFlag.SynRad=2;
    BEAMLINE{iele}.TrackFlag.CSR=-1;
    BEAMLINE{iele}.TrackFlag.CSR_SmoothFactor=1;
    BEAMLINE{iele}.TrackFlag.CSR_DriftSplit=25;
    BEAMLINE{iele}.TrackFlag.Split=25;
  else
    BEAMLINE{iele}.TrackFlag.SynRad=0;
    BEAMLINE{iele}.TrackFlag.CSR=0;
    BEAMLINE{iele}.TrackFlag.Split=0;
  end
end

% FINAL TEST
TwissPlot(1,findcells(BEAMLINE,'Name','ENDBC20'),I,[1 1 0]);
[~,bo]=TrackThru(1,findcells(BEAMLINE,'Name','ENDBC20'),Beam0,1,1);
beamImage(bo);
%beamImage(COMPTEST3);
% END FINAL TEST

% END RAFI SECTION
save FACET2_shortLattice BEAMLINE WF % Save your new short-version of the lattice

function o=S20match(x,dim,ind,I,bmatch,amatch)
global BEAMLINE

if dim=='x'
  di=0;
else
  di=2;
end
BEAMLINE{ind}.R(1+di,1+di)=x(1);
BEAMLINE{ind}.R(2+di,1+di)=x(2);
BEAMLINE{ind}.R(1+di,2+di)=x(3);
BEAMLINE{ind}.R(2+di,2+di)=x(4);
[stat,T]=GetTwiss(1,ind,I.x.Twiss,I.y.Twiss);
if stat{1}~=1
  o=1e6;
elseif dim=='x'
  o = (T.betax(end)-bmatch)^2/bmatch^2 + (T.alphax(end)-amatch)^2/amatch^2 ;
else
  o = (T.betay(end)-bmatch)^2/bmatch^2 + (T.alphay(end)-amatch)^2/amatch^2 ;
end