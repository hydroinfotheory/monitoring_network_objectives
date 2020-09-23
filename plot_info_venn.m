function [H,h,S]=plot_info_venn(Hs,Hc,Hf,JE,JEtot,h);
%%
if nargin<6
    h=figure();
else
    figure(h);
end

% clf
 xlim([-1.5 2.5])
 axis equal, axis off
 
if nargin<1
      % inputs: JE, JEtot,  Hc, TCF,
      JEtot=4.10;
      JE=3.5;
      Hc=1.9;
      Hs=2;
      Hf=3;
      %checks Hc<Hf Hc<JE Hs<JE Hc+Hs>JE Hs+Hf >JEtot
end

Tsf=Hs+Hf-JEtot;
Tsc=Hs+Hc-JE;
Tscf=Tsc;
A=[Hs,Hf,Hc];
I=[Tsf-Tscf,Tscf,Hc,Tscf];
I=[Tsf,Tscf,Hc,Tscf];
%Z=[Hs,Hf,Hc,Tsf,Tsc,Tcf,Tscf];
fminopts=optimset('fminbnd');
%fminopts.MaxFunEvals=100000;
%fminopts.MaxIter=100000;
ccolors={'r','y','c'};
cAlpha={0.3,0.3,0.3};
[H,S]=venn_robust(A,I,fminopts,'ErrMinMode','ChowRodgers','FaceColor',ccolors,'FaceAlpha',cAlpha,'EdgeColor','black');

fmtstr={'color','w','FontWeight','bold'};
  
%hoih=text(S.Position(1,1) ,S.Position(1,2)-0.9*S.Radius(1),  '<H_S>',fmtstr{:});
%text(S.ZoneCentroid(2,1), S.ZoneCentroid(2,2), 'H_F_c');
%text(S.Position(2,1) ,S.Position(2,2)-0.9*S.Radius(2),  '<H_F>',fmtstr{:});
%text(S.Position(3,1) ,S.Position(3,2)-0.9*S.Radius(3),  '<H_F_c>',fmtstr{:});
%text(S.ZoneCentroid(3,1), S.ZoneCentroid(3,2), 'H_F');
%text(S.ZoneCentroid(6,1), S.ZoneCentroid(6,2), 'H_F_c|S',fmtstr{:});

%% make legend values

ccolors={'r','y','c'};
cAlpha={0.3,0.3,0.3};

mAlpha={cAlpha{[1,2,3,1,2,1,3,2,3,1,2,3]}};
mColors={ccolors{[1,2,3,1,2,1,3,2,3,1,2,3]}};
xloc=[1,1,1,1,1,1,1]-2;
yloc=[1.4:-0.2:0]-1;
mxloc=xloc([1,2,3,4,4,5,5,6,6,7,7,7]);
myloc=yloc([1,2,3,4,4,5,5,6,6,7,7,7]);
legend_text={'H_S','H_F','H_F_c','T(S,F)','','hoi13','','hoi23','','hoi123','',''}
radii=(0*myloc)+0.1;
%axis equal

%drawCirclesvenn(xloc, yloc, radii,{'r','y','c','c','c'},{0.3,0.3,0.3,0.3,0.3},{'hoi1','hoi2','hoi3','hoi2','hoi3'});
%drawCirclesvenn(mxloc, myloc, radii,mColors,mAlpha,legend_text);
