%% test greedy ok -> not always, but works for LSM data, mostly close
close all
clear all

% add paths to relevant toolboxes
addbox('infotheory')  ;          % adding hydro-info-theory-lab toolboxes from COMMON folder 
%for export
addbox('plot2svg')    ;          % adding hydro-info-theory-lab toolboxes from COMMON folder
addbox('latex_export');         % add the path for various latex export tools
addbox('venn');
outdir=('./output');

%settings
pareto_plots=1;  % choose whether to plot pareto front plots for all combinations
venn_plots=1;

%seed random number generator for predictable outputs for certain examples
%rng(1);  % 3 replacements at 3,4,5 with gencor
%rng(3);  % 3 replacements more spread/intereting with gencor incl place 2
%rng(4);  % 6 replacements more spread/intereting with gencor from place 3
%rng(6);  % replacement in 2 and 3
%rng(1);  % replacement in 2 and 3
rng(7);   %very interesting case with lots of replacements including at 2 and a double one at 4
% alternative data loading options (overwrite ddata)
gencorr=0;                      % setting for replacing Brazos data with correlated gaussian
genrand=0;                      % setting for replacing Brazos data with random uniform uncorrelated

load ./input/origddata.mat      % load data


%optional: generate normally distributed correlated data
if gencorr==1                   % generate random series with correlation matrix equal to original
    load input/corrQbrazos; %load corr matrix for Q (corrQ)
    
    series=generate_random_correlated(corrQ, 860, 1);
    a=0.4;
    ddata=1+floor((2*series+a)/(2*a)); %series of bin numbers (for joint)
end

%optional: generate random uniform series of same dimensions as original
if genrand==1
    a=0.2;
    series=rand(240,12);
    ddata=1+floor((2*series+a)/(2*a)); %series of bin numbers (for joint)
end

%% this finds the true optimum by exhaustive optimization
prevset=[];                 % initialize
for k=1:12                  % loop over different sub-network sizes  
    tic                     % start measuring time
    combis=nchoosek(1:12,k);% generate all combinations of k out of 12
    nr_c=size(combis,1);    % number of combinations to test
    JE=zeros(nr_c,1);       % initialize
    TC=JE;                  % initialize
    for ii=1:nr_c          % for all combinations of k stations... 
        [JE(ii),TC(ii)]=merging_stats_fast(ddata,combis(ii,:));  %calculate information statistics
    end
    % identify best set of k sensors for both minTC and maxJE objectives
    [JEmax(k),idx_maxJE]=max(JE); %does not identify multiple optima
    [TCmin(k),idx_minTC]=min(TC); 
    
    %identify single best combination of sensors (not robust against ties)
    JEcombi{k}=combis(idx_maxJE,:);
    TCcombi{k}=combis(idx_minTC,:);
    
    %identify sets of sensors dropped and added for each step in k
    JE_optim_added{k}=setdiff(JEcombi{k},prevset);  %note: this approach is not robust against multiple optima
    JE_optim_dropped{k}=setdiff(prevset,JEcombi{k});
    
    try
        JEorder(k)=setdiff(JEcombi{k},prevset);  % try identifying a single added station and add it to ranking vector
    catch
        JEorder(k)=-1;                          % if it fails, replacement occurred, -1 marker added
    end
    prevset=JEcombi{k};                         % update previous set for next iteration
    timetake(k)=toc                             % measure time taken
    
    % pareto plots
    if pareto_plots                             % switch for plotting
        figure();
        %title(num2str(k));
        plot(JE,TC,'r.');                       %plot info measures for all combinations 
        xlabel('Joint Entropy (bits)');
        ylabel('Total Correlation (bits)');
        title(['Information measures for all combinations of ' num2str(k) ' out of 12 locations']);
    end 
    
    
    % all combination's results exhaustive
    allsets(k).combis=combis;                   %Record sets for all combinations
    allsets(k).JE=JE;                           %Record JE for all combinations
    allsets(k).TC=TC;                           %Record TC for all combinations
end
order_JE_optim=JEorder;                         %record for tables
JE_results_optim=JEmax;                         %
total_time=sum(timetake);                       %

%% now test a greedy algorithm maximizing joint entropy by adding one at a
% time
%initialize
idx_left=1:12;              %stations left to select from, not in selected pool
idx_selected=[];            %stations selected in sub network (initially empty)
for k=1:12                  %loop over network sizes
    fi=0;                   %initialize counter
    JE=zeros(1,12);         %initialize JE vector
    for ff=[idx_left]       %loop over candidates for addition
        fi=fi+1;            %counter for candidate sets
        [JE(fi),TC(fi)]=merging_stats_fast(ddata,[idx_selected ff]); % evaluate info metrics for candidate expanded set

    end
    [JEmax(k),idx_maxJE]=max(JE);       % among candiate networks, find expanded network that retains max JE, 
    idx_new_select=idx_left(idx_maxJE);   %get the station nr for added station
    idx_selected=[idx_selected idx_new_select]   % update the selected set to be the previous network plus the added station
    idx_left=setdiff(idx_left,idx_selected);     % update the pool of non-selected candidate stations for next steps.
end
JE_results_add=JEmax;                       % Vector of JE for tables. 
order_greedyJE_ADD=idx_selected;            % Vector for 
JEorder
        

%% Greedy Drop: now test a greedy algorithm maximizing joint entropy by dropping one at a time
%backwards optim
%initialize
idx_left=[]; 
idx_selected=1:12;  %alls stations are selected initially
nr_mon=length(idx_selected);
[JEmax(1),TC_JEmax(1)]=merging_stats_fast(ddata,idx_selected);
for k=2:nr_mon                %loop over network sizes  
    fi=0;                     %initialize counter
    JE=zeros(1,nr_mon);       %initialize JE vector
    for ff=[idx_selected]     % candidates for dropping : loop over all monitors in current selected se
        fi=fi+1;              % counter
        idx_keep_potential=setdiff(idx_selected,ff);    % generate potential reduced set
        [JE(fi),TC(fi)]=merging_stats_fast(ddata,idx_keep_potential);   % evaluate Joint entropy and total correlation for a potential reduced networks

    end
    [JEmax(k),idx_maxJE]=max(JE);    % among candiate networks, find reduced network that retains max JE, 
    idx_drop=idx_selected(idx_maxJE);   %get the station nr for dropped station
    idx_selected=setdiff(idx_selected,idx_drop);  % update the selected set to be the previous network minus the dropped station.
    idx_left=[idx_left idx_drop];                 %keep track of order of dropped stations
end
%add last remaining station
idx_left=[idx_left idx_selected]; 
JE_results_drop=fliplr(JEmax);          %Flip vector of JE, to align with greedy add results for comparison
order_greedyJE_DROP=fliplr(idx_left);   %Flip vector of JE, to align with greedy add results for comparison


%% Ridolfi HESS algorithm : min MI(selected set; new station)
% time
%initialize first selection of max entropy station manually
idx_left=1:11;              
idx_selected=[12];

for k=2:12      % for other stations
    fi=0;       % init counter
    MIlast=0*idx_left; %initialize
    for ff=[idx_left]  %loop over unselected stations
        fi=fi+1;        %counter
        %calculate the Mutual information of the candidate station with the selected set
        [MIlast(fi)]=MI_last_added(ddata,[idx_selected ff]);  

    end
    [minMIlast(k),idx_minMIlast]=min(MIlast);  
    idx_new_select=idx_left(idx_minMIlast);
    idx_selected=[idx_selected idx_new_select]
    idx_left=setdiff(idx_left,idx_selected);
end
MI_results_Ridolfi=minMIlast; 
gauge_order_Ridolfi=idx_selected;

%add gauge orders taken from respective papers
gauge_order_MIMR=[12,6,1,8,2,3,4,7,5,9,10,11];
gauge_order_WMP=[12,9,7,6,5,4,3,2,1,8,11,10];

% Gather goude orders in table
%Gorder=order_JE_optim;
%Gorder=gauge_order;%
GO(1,:)=gauge_order_MIMR;
GO(2,:)=gauge_order_WMP;
GO(3,:)=gauge_order_Ridolfi;
GO(4,:)=order_JE_optim;


%% Plot matrix of Venn Diagrams: 4 methods, multiple steps of k
if venn_plots
    figure('position',[100,100,300,350]);
    hold on;
end
    
%choice of which values of k to plot Venn diagrams for
its=2:6;

% specific row and column for corner of multi-venn plot
c=4;
r=6;

cc=0; %current colummn counter
titles={'MIMR','WMP','minT(S,Fc)','maxJE'};
for ordr=1:4 % loop over different methods
    
    %load gauge addition order from variable
    Gorder=GO(ordr,:);
    % check non greedy xxx SW
    if min(GO(4,:))<1
        GO(4,:)= order_greedyJE_ADD;
    end
    
    cr=0; %current row
    cc=cc+1; %current_col
    % loop over values of k)
    for it=its
        plotnr=it+12*(ordr-1); 
        %subplot(r,c,plotnr,'align');
        cr=cr+1; % next row
        if venn_plots
            subplot('position',[(cc-1)/c,0.95-(cr/r),1/c,0.9/r]);   %positioning of subplot
            text(1.6,0,num2str(Gorder(it)));                        %print added station number
            % add title to top row
            if it==2
                title(titles(ordr));                                %add method name to top row
            end
        end %if venn_plots
        [JEtot,TCtot,Hmargopt]=merging_stats_fast(ddata,Gorder);
        [JEopt,TCopt]=merging_stats_fast(ddata,Gorder(1:it));
        [JEoptF,TCoptF]=merging_stats_fast(ddata,Gorder(it:end));
        [Hs]=merging_stats_fast(ddata,Gorder(1:it-1));

        JEresults=[JE_results_optim;JE_results_add;JE_results_drop]
        JEorder=[order_JE_optim;order_greedyJE_ADD;order_greedyJE_DROP]


        JEtot=JEtot;
        JE=JEopt;
        Hc=Hmargopt(it);
        Hf=JEoptF;
        %Hs=JEresults(1,it-1);
        if venn_plots
            try
                plot_info_venn(Hs,Hc,Hf,JE,JEtot,gcf);
            catch
                text(-1,0,'draw error');
            end
            box on;
            grid on;
        end %if venn_plots
    end
end
%subplot('position',[0,0,1,2/r]);
%plot_info_venn(Hs,Hc,Hf,JE,JEtot,gcf);
%plot2svg('SVG6x4a.svg',gcf)
%svg2pdf([pwd '\SVG6x4a.svg'])
%%
%infile='D:\Users\weijs\Documents\matlab\li_singh_mishra_comment\SVG6x4.svg'
%svg2pdf(infile);

%% single venn plot
ordr=4;it=3; %select JE method, 3 sensors added
Gorder=GO(ordr,:); %select JE method
% collect stats for 3way Venn diagram
[JEtot,TCtot,Hmargopt]=merging_stats_fast(ddata,Gorder);
[JEopt,TCopt]=merging_stats_fast(ddata,Gorder(1:it));
[JEoptF,TCoptF]=merging_stats_fast(ddata,Gorder(it:end));
[Hs]=merging_stats_fast(ddata,Gorder(1:it-1));
JEtot=JEtot;
JE=JEopt;
Hc=Hmargopt(it);
Hf=JEoptF;

if venn_plots
% start plotting
figure1 = figure('InvertHardcopy','off','Color',[1 1 1],'position',[100,100,520,400]);

% Create axes
axes('Visible','off','Parent',figure1,...
    'PlotBoxAspectRatio',[1.2678936605317 1 1],...
    'DataAspectRatio',[1 1 1],'position',[0 0 1 1]);


%draw Venn Circles
[~,~,Vc]=plot_info_venn(Hs,Hc,Hf,JE,JEtot,gcf);


%plot JE circle line to highlight 
draw_arc(Vc,1,3,0  ,1);
%note xxx : if red lines on Venn plot look off try switching to the
%commented out line of the 2 below
draw_arc(Vc,3,1,0,0); %option 1
%draw_arc(Vc,3,1,0,1);  %option 2

%set boundaries of figure
xlim([-1.2,4]);
ylim([-2 2]);
%Annotation
% Create textarrow
annotation(figure1,'textarrow',[0.257692307692308 0.23346153846154],...
    [0.9325 0.65],'TextEdgeColor','none','HorizontalAlignment','left',...
    'String',{'redundant info added by candidate sensor'});

% Create textarrow
annotation(figure1,'textarrow',[0.471153846153846 0.411538461538464],...
    [0.885 0.655],'TextEdgeColor','none','HorizontalAlignment','left',...
    'String',{'new info added by candidate sensor'});

% Create textarrow
annotation(figure1,'textarrow',[0.492307692307692 0.42967032967033],...
    [0.16 0.355833333333333],'TextEdgeColor','none',...
    'HorizontalAlignment','left',...
    'String',{'Non-captured information','(uncertainty left,','room for improvement)'});

% Create textarrow
annotation(figure1,'textarrow',[0.148076923076923 0.178846153846154],...
    [0.14 0.28],'TextEdgeColor','none','HorizontalAlignment','left',...
    'String',{'net info collected','by current set of sensors'},...
    'HeadStyle','diamond');



xshift=1.5;
yshift=-1;
radius=0.13;
ccolors={'r','y','c'};
cAlpha={0.3,0.3,0.3};


mAlpha={cAlpha{[1,2,3,1,2,1,3,2,3,1,2,3]}};
mColors={ccolors{[1,2,3,1,2,1,3,2,3,1,2,3]}};
xloc=[1,1,1,1,1,1,1]+xshift;
yloc=[1.4:-0.2:0]+yshift;
mxloc=xloc([1,2,3,4,4,5,5,6,6,7,7,7]);
myloc=yloc([1,2,3,4,4,5,5,6,6,7,7,7]);
legend_text={'full:H(S), free:H(S|F)','full:H(F), free:H(F|S,F_c)','full:H(F_c), free: 0','T(S;F)-T(S;F_c)','','0, because F_c \subset{} F','','H(F_c|S)','','T(S;F_c)','',''}
radii=(0*myloc)+0.1;


%% draw legend
drawCirclesvenn([1-0.3*radius,1,1+0.3*radius]+xshift, [-1.7*radius,-1*radius,-1.7*radius]+yshift+2, [1,1,1]*radius,ccolors,cAlpha,{'','color legend','',''});
drawCirclesvenn(mxloc, myloc, radii,mColors,mAlpha,legend_text);
drawCirclesvenn(mxloc-0.2, myloc, 0*radii,mColors,mAlpha,{'1','2','3','','4','','5','','6','','','7'});

%save as svg and then as pdf
plot2svg('single_venn_fixed1.svg',gcf)
svg2pdf([pwd '\single_venn_fixed1.svg']);
end % if venn_plots

%% develop subset of JE optimal solutions
prevset=[];
for k=1:12
    %Find all max JE's
    
    %[JEmax(k),idx_maxJE]=max(allsets(k).JE); %xxx not robust against multiple optima
    %find all indices of JE values equal to max.
    idx_max_JE=find(allsets(k).JE==max(allsets(k).JE));
    allsets_opt_JE(k).JE=allsets(k).JE(idx_max_JE);
    allsets_opt_JE(k).TC=allsets(k).TC(idx_max_JE);
    allsets_opt_JE(k).combis=allsets(k).combis(idx_max_JE,:);
    
    idx_min_TC=find(allsets(k).TC==min(allsets(k).TC));
    allsets_opt_TC(k).JE=allsets(k).JE(idx_min_TC);
    allsets_opt_TC(k).TC=allsets(k).TC(idx_min_TC);
    allsets_opt_TC(k).combis=allsets(k).combis(idx_min_TC,:);  
    
    
    %JE_optim_added{k}=setdiff(allsets_opt_JE(k).combis
    
%     [TCmin(k),idx_minTC]=min(TC);
%     JEcombi{k}=combis(idx_maxJE,:);
%     TCcombi{k}=combis(idx_minTC,:);
    JEcombi{k}=allsets_opt_JE(k).combis;
    JE_optim_added{k}=setdiff(JEcombi{k},prevset);
    JE_optim_dropped{k}=setdiff(prevset,JEcombi{k});
    try
        JEorder(k)=setdiff(JEcombi{k},prevset);
        %make cell
        %order_JE_exh{k}=(setdiff(JEcombi{k},prevset));
        order_JE_exh{k}=[JE_optim_added{k} -JE_optim_dropped{k}];
        
    catch
        JEorder(k)=-1;
        order_JE_exh{k}=[JE_optim_added{k} -JE_optim_dropped{k}];
        
    end
    prevset=JEcombi{k};

end

%% test greedy ok : test of exhaustive and greedy algorithms are robustly the same, taking into account ties.
% for each optimal set for k, see if there at least one optimal
% set for k+1 that includes it. If this is true for all optimal
% sets, then greedy add is OK.
has_no_dead_ends(1:11)=0;
for k=1:11
    nr_optk=length(allsets_opt_JE(k).JE);
    nr_optk1=length(allsets_opt_JE(k+1).JE);
    child_found(1:nr_optk)=0;
    for ii=1:nr_optk % for all current sets
        for jj=1:nr_optk1 % check all next sets
            setk=allsets_opt_JE(k).combis(ii,:);
            setk1=allsets_opt_JE(k+1).combis(jj,:);
            % check if there is no sensor in k that is not in k+1
            % if there is no such sensor, then set k,ii has a
            % follow up set ("child")
            if isempty(setdiff(setk,setk1))
                child_found(ii)=1;
                break
            end
        end
    end
    % check if parent found for all ii (optimal combinations of k sensors)
    if min(child_found)==1
        %now we know there are follow up optimal sets of size k+1...
        % for each optimal set in k
        has_no_dead_ends(k)=1;
    else
        max_idx_deadend(k)=ii;
    end
end
% if this is true of all sizes of sets, then the greedy add stategy
% is optimal for this example
if min(has_no_dead_ends)==1
    disp('greedy add is optimal');
else
    disp('greedy add is non-optimal')
end

%print all optimal sets from exhaustive optimization.

allsets_opt_JE.combis        

%do check if back-forth sweep algorithm would be optimal
% first: simple test on JE values, could compare sets for each K to be more
% robust
JE_loss_greedy_add=JEresults(1,:)-JEresults(2,:)
JE_loss_greedy_drop=JEresults(1,:)-JEresults(3,:)
JE_loss_greedy_back_forth=min([JE_loss_greedy_add;JE_loss_greedy_drop])

         
%% Export tables 
maindir=cd();
cd(outdir);
% convert the relevant matrices to latex tables

% 1)Gauge orders

matrix2latex(GO,'testGO1.txt','rowlabels',{'MIMR','WMP1/2','minT','maxJE'},'columnlabels', [1:12])

% 2) Joint entropy for 4 methods
matrix2latex(JE,'testJE1.txt','rowlabels',{'M1','M2','minT','maxJE'},'columnlabels', [1:12])

% 3) Exhaustive, greedy add, greedy drop
% convert everything to cells for case where cells are not single numbers
rw1=order_JE_exh;
%rw2=mat2cell(order_greedyJE_ADD,1,ones(1,12));
%rw3=mat2cell(order_greedyJE_DROP,1,ones(1,12));
rw2=num2cell(order_greedyJE_ADD);
rw3=num2cell(order_greedyJE_DROP);
JEorderTable1=[rw1;rw2;rw3];
JEorderTable=cellfun(@num2str,JEorderTable1,'UniformOutput',false);


%gauge order tables
%matrix2latex(JEresults,'JE_greedy_compare.txt','rowlabels',{'Exhausitve','Greedy Add','Greedy Drop'},'columnlabels', [1:12])
%matrix2latex(JEorderTable,'order_greedy_compare.txt','rowlabels',{'Exhausitve','Greedy Add','Greedy Drop'},'columnlabels', [1:12],'format','%-1.0f')
matrix2latex(JEorderTable,'order_greedy_compare.txt','rowlabels',{'Exhausitve','Greedy Add','Greedy Drop'},'columnlabels', [1:12])

%JE for orders
matrix2latex(JEresults,'JE_greedy_compare.txt','rowlabels',{'Exhausitve','Greedy Add','Greedy Drop'},'columnlabels', [1:12])

% format files names

% write tables to latex files

 %% write all optimal combis in big table
 rows_total=[0 [1:12]]
 for k=1:12
     nr_comb=size(allsets_opt_JE(k).combis,1);
     rows_sel=[666*ones(nr_comb,13)];
     for ii=1:nr_comb
         rows_sel(ii,1)=k;
         rows_sel(ii,1+allsets_opt_JE(k).combis(ii,:))=777;
     end
     rows_total=[rows_total; rows_sel]
 end
 
 %note 666 and 777 can be replaced by dots and empty spaces in text editor
 %in latex out file
 matrix2latex(rows_total,'optimal_combis.txt','columnlabels', [1 1:12])
     
cd(maindir);
            
            
