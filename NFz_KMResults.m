function NFz_KMResults
tic; close all;
screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check flags!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do_print = true;
do_gd3_exclude = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('**** START Flags ****');
disp([num2str(do_print),' - do_print']);
disp('**** END Flags ****');
disp(sprintf('\n'));

fig_loc = 'Z:/elw/MATLAB/nfz_analy/slides/figures/latest/';

%a2b = {'Inf' '3'};
a2b = {'10'};

%structures = {'PBT' 'ILUNG' 'ESOPHAGUS' 'HEART' 'NFZ' 'LUNGS'};
structures = {'ESOPHAGUS'};

%toxicities = {'rp','pultox','esotox'};
toxicities = {'esotox'};

fp = 'Z:\elw\MATLAB\nfz_analy\meta_data\';

for i=1:length(toxicities)

    for j=1:length(structures)
        cur_fig_ctr = (10*i)+j-1;
        
        fprintf('\n');
        disp(['Tox: ',toxicities{i}]);
        disp(['Struct: ',structures{j}]);
        disp(['Counter: ',num2str(cur_fig_ctr)]);
        fprintf('\n');
        
        %% load data

    if do_gd3_exclude
        fn = ['NFZ_',structures{j},'_',toxicities{i},'_a2b',a2b{1},'_nogd3_data.mat'];
    else
        fn = ['NFZ_',structures{j},'_',toxicities{i},'_a2b',a2b{1},'_data.mat'];
    end
        disp(['']);
        disp(['Loading ',fn]);
        disp(['']);
        load(strcat(fp,fn),'CGobj_org');
        CGobj = CGobj_org;
        
        if j==1  %one complication curve for all structures
            
            %% Complication incidence curve
            cur_fig=figure(i);clf reset;
            set(cur_fig,'Position',ss_four2three);
            % complication incidence curve
            
            sa = CGobj.mKaplanMeierCompOverall;
            sa_times = sa.mSurvivalTimeSorted{1};
            sa_times(sa_times<0)=0;
            sa_curve = 1-sa.mSurvivalCurve{1};
            stairs(sa_times,sa_curve,'b','LineWidth',2);hold on;
            plot(sa_times(sa.mCensorStatistics{1}(:,1)),...
                sa_curve(sa.mCensorStatistics{1}(:,1)),...
                'b+','MarkerSize',20);hold off;
            xlim([0 max(sa_times)+5]);
            ylim([0 max(sa_curve)+0.02]);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'FontSize',18);
            set(gca,'box','on');
            grid on;
            set(gca,'GridLineStyle','--')
            
            xlabel('Months','fontsize',20);
            ylabel('Probability of Complication','fontsize',20);
            
            if do_print,
                set(cur_fig,'Color','w');
                export_fig(cur_fig,...
                    [fig_loc,'nfz_',...
                    toxicities{i},'_km_overall'],'-pdf');
                %system('pdfcrop fname fname');% to remove borders
                disp(['Saving ',...
                    fig_loc,'nfz_',...
                    toxicities{i},'_km_overall.pdf']);
            end
        end
    end
end