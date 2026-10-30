%% R1_C6_R3_C1_Fig9_VGPP_Workflow_R2024b.m
% Publication-oriented VGPP Fig. 9. Standalone MATLAB R2024b (base MATLAB).
% Improves legibility and alignment while preserving seed construction,
% relaxation, growth/backtracking, reduced-count fallback, and strict auditing.
% NOTE: A reduced-count fallback does NOT meet the original target.
clc; clear; close all;

projectDir = '.';
picsDir = fullfile(projectDir,'pics');
if ~isfolder(picsDir), mkdir(picsDir); end
pngFile = fullfile(picsDir,'Fig9_VGPP_workflow_R1_C6_R3_C1_final_adjusted.png');

F = figure('Color','w','Units','pixels','Position',[80 40 2800 1750], ...
           'Toolbar','none','MenuBar','none');
set(F,'PaperPositionMode','auto');
font = 'Times New Roman';
navy=[0.10 0.21 0.36]; blue=[0.12 0.37 0.68];
green=[0.13 0.42 0.26]; red=[0.69 0.21 0.22]; amber=[0.65 0.36 0.10];
lightBlue=[0.96 0.98 1.00]; lightGreen=[0.965 0.990 0.969];
lightRed=[1.00 0.965 0.958]; lightAmber=[1.00 0.983 0.955];

%% Stage panels: equal height, simple headings, no redundant oversized title
annotation(F,'rectangle',[.018 .040 .375 .927],'Color',blue,'FaceColor',[.984 .991 1], ...
    'LineStyle','--','LineWidth',1.4);
annotation(F,'rectangle',[.408 .040 .574 .927],'Color',green,'FaceColor',[.986 .995 .987], ...
    'LineStyle','--','LineWidth',1.4);
node([.058 .899 .296 .050],'STAGE 1 | SEED CONSTRUCTION',[1 1 1],blue,17.5,blue);
node([.558 .885 .275 .050],'STAGE 2 | FIXED-RADIUS DENSIFICATION',[1 1 1],green,17.0,green);

%% Seed construction: regular spacing, centered within left panel
s1=[.065 .782 .281 .067];
s2=[.065 .663 .281 .067];
s3=[.065 .544 .281 .067];
s4=[.065 .425 .281 .067];
s5=[.065 .299 .281 .074];
node(s1,{'Parameter setup','Domain, target area fraction, radius, PBC'},lightBlue,blue,15.0,navy);
node(s2,{'Define target and seed','Set target N and reduced seed radius'},lightBlue,blue,15.0,navy);
node(s3,{'Void-guided candidate insertion','Map voids; prioritize admissible centers'},lightBlue,blue,14.7,navy);
node(s4,{'Stochastic rearrangement','Update void map when insertion stalls'},lightBlue,blue,14.7,navy);
node(s5,{'Terminate seed-insertion phase','Target count reached or stopping budget exhausted'},lightBlue,blue,14.4,navy);
arr([cx(s1) low(s1)],[cx(s2) high(s2)],navy);
arr([cx(s2) low(s2)],[cx(s3) high(s3)],navy);
arr([cx(s3) low(s3)],[cx(s4) high(s4)],navy);
arr([cx(s4) low(s4)],[cx(s5) high(s5)],navy);
route([left(s4) cy(s4);.044 cy(s4);.044 cy(s3);left(s3) cy(s3)],blue,'--');
tag([.019 .493 .053 .020],'REMAP',blue,11);

%% Main growth/relaxation column
A=[.462 .787 .232 .067];
B=[.462 .687 .232 .067];
D1=[.480 .570 .196 .093];
C=[.462 .480 .232 .067];
D2=[.480 .365 .196 .092];
E=[.462 .274 .232 .067];
G=[.462 .185 .232 .064];
D3=[.755 .175 .192 .089];
OK=[.462 .066 .232 .063];
BAD=[.755 .041 .192 .068];
R=[.757 .570 .190 .093];
BACK=[.757 .782 .190 .067];
FALL=[.757 .374 .190 .083];
node(A,{'Propose incremental radius growth','Limit trial radius to prescribed final value'},lightGreen,green,14.2,navy);
node(B,{'Collective periodic relaxation','Resolve overlaps at the trial radius'},lightGreen,green,14.4,navy);
decision(D1,{'Relaxation','successful?'},amber,14.4,navy);
node(C,{'Accept growth increment','Update centers and current radius'},lightGreen,green,14.4,navy);
decision(D2,{'Final radius','attained?'},amber,14.4,navy);
node(E,{'Exact-final-radius relaxation','Resolve overlaps at prescribed final radius'},lightGreen,green,14.0,navy);
node(G,{'Final periodic geometric audit','Check count, radius and overlap tolerance'},lightGreen,green,14.1,navy);
decision(D3,{'All acceptance','criteria satisfied?'},amber,14.0,navy);
node(OK,{'ACCEPT','Save validated configuration'},lightGreen,green,14.5,navy);
node(BAD,{'NOT ACCEPTED','Original count unmet or audit failed'},lightRed,red,13.7,navy);
decision(R,{'Retry budget','remaining?'},amber,14.1,navy);
node(BACK,{'Backtrack','Reduce increment; repeat relaxation'},lightRed,red,13.7,navy);
node(FALL,{'Reduced-count fallback','Remove particle; audit reduced geometry'},lightRed,red,13.0,navy);

% Seed -> growth: dedicated external lane, away from all boxes
route([right(s5) cy(s5);.402 cy(s5);.402 cy(A);left(A) cy(A)],navy,'-');
tag([.402 .835 .053 .018],'PROCEED',navy,11);

% Main flow, with right-angled connectors only
arr([cx(A) low(A)],[cx(B) high(B)],navy);
arr([cx(B) low(B)],[cx(D1) high(D1)],navy);
arr([cx(D1) low(D1)],[cx(C) high(C)],green);
tag([.615 .542 .040 .021],'YES',green,11.5);
arr([cx(C) low(C)],[cx(D2) high(D2)],navy);
arr([cx(D2) low(D2)],[cx(E) high(E)],green);
tag([.616 .343 .040 .021],'YES',green,11.5);
arr([cx(E) low(E)],[cx(G) high(G)],navy);
% G and D3 have slightly different center heights: use a stepped connector.
route([right(G) cy(G);.730 cy(G);.730 cy(D3);left(D3) cy(D3)],navy,'-');

% NO: final radius not yet attained -> return to radius growth
route([left(D2) cy(D2);.438 cy(D2);.438 cy(A);left(A) cy(A)],blue,'-');
tag([.417 .458 .035 .020],'NO',red,11.3);

% Relaxation failure -> retry decision
arr([right(D1) cy(D1)],[left(R) cy(R)],red);
tag([.697 .626 .035 .020],'NO',red,11.3);
% Retry YES -> backtrack, then return to relaxation from right edge
arr([cx(R) high(R)],[cx(BACK) low(BACK)],green);
tag([.876 .684 .039 .020],'YES',green,11.3);
route([left(BACK) cy(BACK);.722 cy(BACK);.722 cy(B);right(B) cy(B)],red,'-');
% Retry NO -> fallback via right-side external lane
route([right(R) cy(R);.970 cy(R);.970 cy(FALL);right(FALL) cy(FALL)],red,'-');
tag([.929 .542 .034 .020],'NO',red,11.3);

% Final acceptance exits: BOTH arrows start at the decision diamond.
% YES exits from its LOWER vertex and bends left just above the ACCEPT box.
% NO exits from its RIGHT vertex, using a separate lower route to the lowered NOT ACCEPTED box.
% These paths neither branch from the incoming audit line nor cross the audit box.
route([cx(D3) low(D3);cx(D3) .136;cx(OK) .136;cx(OK) high(OK)],green,'-');
tag([.705 .149 .036 .020],'YES',green,11.2);
route([right(D3) cy(D3);.955 cy(D3);.955 .118;.859 .118;.859 high(BAD)],red,'-');
tag([.946 .189 .031 .020],'NO',red,11.2);

% Fallback includes its own geometric audit; even if that audit passes,
% the original particle-count target is not attained. The red arrow
% therefore leads to the NOT ACCEPTED (original target) terminal state.
% Routed on far-right boundary and enters the rejection box from the right.
route([right(FALL) cy(FALL);.978 cy(FALL);.978 cy(BAD);right(BAD) cy(BAD)],red,'-');

% Compact footnote, outside stage panels (do not imply success after fallback)
annotation(F,'textbox',[.175 .003 .650 .030], ...
    'String','Outputs: coordinates, attained area fraction, periodic audit and target-attainment status', ...
    'LineStyle','none','FontName',font,'FontSize',14.0,'FontWeight','bold', ...
    'Color',navy,'HorizontalAlignment','center','VerticalAlignment','middle');

exportgraphics(F,pngFile,'Resolution',600,'BackgroundColor','white');
fprintf('Figure saved: %s\n',pngFile);

%% Local functions (standalone; no toolboxes required)
function node(r,t,bg,ec,fs,tc)
annotation('textbox',r,'String',t,'FitBoxToText','off','FontName','Times New Roman', ...
    'FontSize',fs,'FontWeight','bold','Color',tc,'BackgroundColor',bg, ...
    'EdgeColor',ec,'LineWidth',1.45,'HorizontalAlignment','center', ...
    'VerticalAlignment','middle','Interpreter','tex');
end
function decision(r,t,ec,fs,tc)
x=r(1);y=r(2);w=r(3);h=r(4);
P=[x+w/2 y+h;x+w y+h/2;x+w/2 y;x y+h/2;x+w/2 y+h];
for k=1:4
    annotation('line',P(k:k+1,1),P(k:k+1,2),'Color',ec,'LineWidth',1.7);
end
annotation('textbox',[x+.18*w y+.16*h .64*w .68*h],'String',t, ...
    'LineStyle','none','FontName','Times New Roman','FontSize',fs, ...
    'FontWeight','bold','Color',tc,'HorizontalAlignment','center', ...
    'VerticalAlignment','middle');
end
function arr(p,q,col)
assert(abs(p(1)-q(1))<1e-8 || abs(p(2)-q(2))<1e-8,'Nonorthogonal segment');
annotation('arrow',[p(1),q(1)],[p(2),q(2)],'Color',col, ...
    'LineWidth',1.65,'HeadLength',8,'HeadWidth',8);
end
function route(P,col,style)
for k=1:size(P,1)-2
 p=P(k,:);q=P(k+1,:);
 assert(abs(p(1)-q(1))<1e-8 || abs(p(2)-q(2))<1e-8,'Nonorthogonal route');
 if norm(p-q)>1e-8
  annotation('line',[p(1),q(1)],[p(2),q(2)],'Color',col, ...
   'LineWidth',1.65,'LineStyle',style);
 end
end
arr(P(end-1,:),P(end,:),col);
end
function tag(r,t,c,fs)
annotation('textbox',r,'String',t,'LineStyle','none', ...
    'BackgroundColor','white','FontName','Times New Roman', ...
    'FontSize',fs,'FontWeight','bold','Color',c, ...
    'HorizontalAlignment','center','VerticalAlignment','middle');
end
function x=left(r),x=r(1);end
function x=right(r),x=r(1)+r(3);end
function y=high(r),y=r(2)+r(4);end
function y=low(r),y=r(2);end
function x=cx(r),x=r(1)+r(3)/2;end
function y=cy(r),y=r(2)+r(4)/2;end
