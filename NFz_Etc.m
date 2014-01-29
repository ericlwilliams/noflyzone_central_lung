function NFz_Etc
tic; close all;
screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check flags!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do_print = true;
do_debug = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('**** START Flags ****');
disp([num2str(do_print),' - do_print']);
disp([num2str(do_debug),' - do_debug']);
disp('**** END Flags ****');
fprintf('\n');

fig_loc = 'Z:/elw/MATLAB/nfz_analy/slides/figures/latest/';

%a2b = {'Inf' '3'};
%a2b = {'10'};
a2b = {'3'};

structures = {'ILUNG' 'ESOPHAGUS' 'HEART' 'LUNGS' 'NFZ' 'PBT'};
%structures = {'LUNGS','ILUNG'};

toxicities = {'rp','pultox','esotox'};
%toxicities = {'pultox'};

fp = 'Z:\elw\MATLAB\nfz_analy\meta_data\';

for i=1:length(toxicities)
    
    spearman_corrs = inf(length(structures),1);
    logreg_pvals = inf(length(structures),1);
    logreg_pos_corrs = inf(length(structures),1);
    rnksum_pvals = inf(length(structures),1);
    rnksum_pos_corrs = inf(length(structures),1);
    
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
        disp(['']);
        disp(['Loading ',fn]);
        disp(['']);
        load(strcat(fp,fn),'CGobj_org');
        CGobj = CGobj_org;
        clear CGObj_org

        pttotal = ones(CGobj.mNumInGrp,1);
        ptcomp = ones(CGobj.mNumInGrp,1); ptcomp([CGobj.mGrp.mFlgCensor])=0;
        
        max_vol=max([CGobj.mGrp.mVolCum])';
        
        %% Max volume histograms
        
         cur_fig=figure(cur_fig_ctr);clf reset;
        set(cur_fig,'Position',ss_four2three);
        cur_fig_ctr=cur_fig_ctr+1;
        
      [n_vol,x_vol]=hist(max_vol(logical(~ptcomp)),0:max(max_vol)/10:max(max_vol));
      [n_vol_comp,x_vol_comp]=hist(max_vol(logical(ptcomp)),0:max(max_vol)/10:max(max_vol));
    
        b_vol_cens=bar(x_vol,n_vol,'b');hold on;
        b_vol_comp=bar(x_vol_comp,n_vol_comp,'r');
        xlim([0 max(max_vol)+max(max_vol)/10]);
        lgnd_vol=legend([b_vol_cens b_vol_comp],'No Complication',' Complication',...
            'Location','Best');
        set(lgnd_vol,'FontSize',16);    
        set(gca,'FontSize',15);
    
        xlabel(['V_{',structures{j},'} [cc]'],'FontSize',18);
        ylabel('Frequency','FontSize',18);
        if do_print,
           set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'h_nfz_',structures{j},'_',toxicities{i},'_volume'],'-pdf');
        end
    
        %% look at large volumes
        if isequal(structures{j},'LUNGS')
            is_large = max_vol > 6000;
        elseif isequal(structures{j},'ILUNG')
            is_large = max_vol > 3000;
        else
            is_large = max_vol > 6000;
        end
        large_vols = max_vol(is_large);
        euds = [CGobj.mGrp.mEUD]';
        lrg_vol_mean_doses = euds(:,11);
        lrg_vol_mean_doses = lrg_vol_mean_doses(is_large);
        
        lrg_vol_ptcomp = ptcomp(is_large);
        
        [sort_lrg_vols,sort_idx] = sort(large_vols);
        sort_lrg_vol_mean_doses = lrg_vol_mean_doses(sort_idx);
        sort_lrg_vol_ptcomp = lrg_vol_ptcomp(sort_idx);
        
        
        cur_fig=figure(cur_fig_ctr);clf reset;
        set(cur_fig,'Position',ss_four2three);
        cur_fig_ctr=cur_fig_ctr+1;
        
        sorted_xaxis = 0:length(sort_lrg_vol_mean_doses)+1;% 0-4
        h_srt_comp=plot(sorted_xaxis(logical(sort_lrg_vol_ptcomp))+1,...
            sort_lrg_vol_mean_doses(logical(sort_lrg_vol_ptcomp)),...
            'o','MarkerSize',20,'MarkerFaceColor','r');hold on;
        h_srt_cens=plot(sorted_xaxis(logical(~sort_lrg_vol_ptcomp))+1,...
            sort_lrg_vol_mean_doses(logical(~sort_lrg_vol_ptcomp)),...
            'o','MarkerSize',20,'MarkerFaceColor','b');
        cur_lgnd=legend([h_srt_comp h_srt_cens],'Complication','No Complication');
        set(cur_lgnd,'FontSize',18);
        set(gca,'XTick',[1 2 3]);
        set(gca,'XTickLabel',num2str(sort_lrg_vols,5));

%       set(gca,'YTickLabel',num2str(sort_lrg_vol_mean_doses,2))
        set(gca,'FontSize',18);
        if isequal(structures{j},'LUNGS'),
            xlabel(['Ranked ',structures{j},' volume [cc]',10,'For V>6000cc'],'FontSize',20)
        elseif isequal(structures{j},'ILUNG')
             xlabel(['Ranked ',structures{j},' volume [cc]',10,'For V>3000cc'],'FontSize',20)
        else
             xlabel(['Ranked ',structures{j},' volume [cc]',10,'For V>6000cc'],'FontSize',20)
        end
        ylabel(['Mean Dose to ',structures{j}],'FontSize',20);
        grid on;
        
        if do_print,
           set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'nfz_',structures{j},'_',toxicities{i},'_vol_ranked_dmean'],'-pdf');
        end
    
        %% Volume/toxicity correlations
        disp(['Correlations: ',structures{j},' volume and ',toxicities{i}]); 
            
        
              
        %% spearman
        spearman_corrs(j) = corr(max_vol,ptcomp,'type','Spearman');
        disp(['']);
        disp(['Spearman \rho: ',...
            num2str(spearman_corrs(j),4)]);
        
        %% rank sum
        [rnksum_pvals(j),~,rnksum_stats] = ranksum(max_vol(logical(ptcomp)),max_vol(logical(~ptcomp)));
        rnksum_pos_corrs(j) = logical(rnksum_stats.zval > 0);
        disp(['']);
        disp(['Rank Sum p-value: ',...
            num2str(rnksum_pvals(j),4)]);
        
        %% logistic
        [b,~,s]=glmfit(max_vol,[ptcomp pttotal],'binomial','link','logit');
        logreg_pvals(j) = s.p(2);
        logreg_pos_corrs(j) = logical(b(2) > 0);
        disp(['']);
        disp(['Logistic regression beta, p-value: ',...
            num2str(b(2),4),', ',num2str(logreg_pvals(j),4)]);
        
    end
end