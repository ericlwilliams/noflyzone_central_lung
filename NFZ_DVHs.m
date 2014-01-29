function NFZ_DVHs
tic; close all;
screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 screen_size(4)/2];

do_print = true;
fig_loc = 'Z:/elw/MATLAB/nfz_analy/slides/figures/latest/';


mLymanN = -1:0.1:1;

%structures = {'ILUNG' 'ESOPHAGUS' 'HEART' 'NFZ' 'PBT' 'LUNGS'};
%structures = {'GTV'}
structures = {'ESOPHAGUS'}
%toxicities = {'rp','pultox','esotox};

%toxicities = {'lclfail'};
toxicities = {'esotox'};
%toxicities = {'pultox'};

fp = 'Z:\elw\MATLAB\nfz_analy\meta_data\';

%a2b = {'Inf' '3' '10'};
a2b = {'10'};

for i=1:length(toxicities)

    for j=1:length(structures)
    cur_fig_ctr = (10*i)+j-1;
    
    fprintf('\n');
    disp(['Tox: ',toxicities{i}]);
    disp(['Struct: ',structures{j}]);    
    disp(['Counter: ',num2str(cur_fig_ctr)]);
    fprintf('\n');
    
     fig_basename = [fig_loc,'nfz_',...
                        structures{j},'_',...
                        toxicities{i},'_a2b',...
                        a2b{1}];
    
    %% load data
    fn = ['NFZ_',structures{j},'_',toxicities{i},'_a2b',a2b{1},'_data.mat'];
    load(strcat(fp,fn),'CGobj_org');
    CGobj = CGobj_org;
    LymanN = log10(CGobj.mLymanN);
    CGobj.mLymanN = LymanN;        
    f = [CGobj.mGrp.mFlgCensor];


    %% DVHS
    cur_fig=figure(cur_fig_ctr); hold on;
    set(gcf,'Position',ss_four2three)
    CGobj.fDVHCurvesSummary_DVH(); hold off;    
    
    %title(['DVHs: ',structures{j}],'FontSize',15);
    xlabel('Dose [Gy]','FontSize',15);
    ylabel('Volume [cc]','FontSize',15);
   
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,...
            [fig_basename,'_avg_dvhs'],'-pdf');
        disp(['Saving ',fig_basename,'_avg_dvhs.pdf']);

    end
    

    cur_fig=figure(cur_fig_ctr+100); hold on;
    set(gcf,'Position',ss_four2three);
    % DVHs of censored patients
    hold on;
    for k=1:length(f),
        if f(k), % cens
            dvh_color = 'b';
        else
            dvh_color = 'r';
        end
        
        dosebins = CGobj.mGrp(k).mDoseBins_LQ;
        volcum = CGobj.mGrp(k).mVolCum;
        plot(dosebins(1:end-1), volcum(1:end-1),dvh_color);
        
    end
    %title(['DVHs: ',structures{j}],'FontSize',15);
    set(gca,'xminortick','on','yminortick','on');
    xlabel('Dose [Gy]','FontSize',18);
    ylabel('Volume [cc]','FontSize',18);
    %ylim([0 40]);
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,...
            [fig_basename,'_dvhs'],'-pdf');
        disp(['Saving ',fig_basename,'_dvhs.pdf']);

    end
    
     
    avg_eud = [CGobj.mGrp.mEUD]';
    
    avg_cens_eud = mean(avg_eud(f,:));
    std_cens_eud = std(avg_eud(f,:));
    
    avg_comp_eud = mean(avg_eud(~f,:));
    std_comp_eud = std(avg_eud(~f,:));
    
    
    fig_avg_euds=figure(cur_fig_ctr+200); hold on; grid on;
    set(gcf,'Position',ss_four2three);
    %if k==1, set(gcf,'Position',ss_four2three);end;
    %subplot(2,2,k);
    hold on;
    
    h(1)=plot(avg_comp_eud,mLymanN,'-r*','LineWidth',2);
    h(2)=plot(avg_comp_eud+std_comp_eud,mLymanN,'--r','LineWidth',0.5);
    plot(avg_comp_eud-std_comp_eud,mLymanN,'--r','LineWidth',0.5);
    
    h(3)=plot(avg_cens_eud,mLymanN,'-b*','LineWidth',2);
    h(4)=plot(avg_cens_eud+std_cens_eud,mLymanN,'--b','LineWidth',0.5);
    plot(avg_cens_eud-std_cens_eud,mLymanN,'--b','LineWidth',0.5);
    
    ylim([-1 1]);
    xlim([0 max([max(avg_cens_eud+std_cens_eud) max(avg_comp_eud+std_comp_eud)])+0.5]);
    
   %title(['Avg. gEUDs: ',structures{j}],'FontSize',15);
    set(gca,'YTickLabel',1:-0.2:-1);
    legend(h,'Avg. with comp','Central 68%','Avg. no comp','Central 68%',...
        'Location','NorthEast');
    ylabel('log_1_0(a)','FontSize',18);
    xlabel('gEUD [Gy]','FontSize',18);
    grid on;
  
      
    if do_print,
        set(fig_avg_euds,'Color','w');

        export_fig(fig_avg_euds,...
            [fig_basename,'_avg_euds'],'-pdf');
        disp(['Saving ',fig_basename,'_avg_euds.pdf']);

    end;
    
    % spaghetti euds
    
    
    fig_avg_euds=figure(cur_fig_ctr+300); hold on; grid on;
    set(gcf,'Position',ss_four2three);
    % DVHs of censored patients
    first_comp=true;
    first_cens=true;
    max_eud=-1;
    for l=1:length(f),
        if f(l), % cens
            dvh_color = 'b';
        else
            dvh_color = 'r';
        end
        
        euds = CGobj.mGrp(l).mEUD';
        if max(euds)>max_eud
            max_eud=max(euds);
        end
        if first_comp && ~f(l)
            g(1)=plot(euds,mLymanN,dvh_color);
            first_comp=false;
        elseif first_cens && f(l)
            g(2)=plot(euds,mLymanN,dvh_color);
            first_cens=false;
        else
            plot(euds,mLymanN,dvh_color);
        end
                      
    end
    ylim([-1 1]);
        %xlim([0 max([max(avg_cens_eud+std_cens_eud) max(avg_comp_eud+std_comp_eud)])+0.5]);
    xlim([0 max_eud+0.5]);
    %title(['gEUDs: ',structures{j}],'FontSize',15);
    set(gca,'YTickLabel',1:-0.2:-1);
    set(gca,'xminortick','on','yminortick','on');
    legend(g,'Complication','No Complication',...
        'Location','NorthEast');
    
    ylabel('log_1_0(a)','FontSize',15);
    xlabel('gEUD [Gy]','FontSize',15);

    
end
end
end