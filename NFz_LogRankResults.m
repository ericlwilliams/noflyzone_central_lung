function NFz_LogRankResults
%% Also includes logistic regression results for range of dv and vd
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
a2b = {'Inf'};

%structures = {'ILUNG' 'ESOPHAGUS' 'HEART' 'LUNGS' 'NFZ' 'PBT'};
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
          
        pttotal = ones(CGobj.mNumInGrp,1);
        ptcomp = ones(CGobj.mNumInGrp,1); ptcomp([CGobj.mGrp.mFlgCensor])=0;
        
        volume_bins=CGobj.mBinsVol;
        dose_bins=CGobj.mBinsDose;
    
        %% Get Logrank results
        f_lr = cellfun(@(x) strcmpi('VDx',x),CGobj.mLogRank(:,1));
        pvx = CGobj.mLogRank{f_lr,2};


          %% Get Vx results
        med_vx_p = nan(length(dose_bins),1);
        med_vx_hr = nan(length(dose_bins),1);
            
        lr_vx_p = nan(length(dose_bins),1);
        is_corr = logical(ones(length(dose_bins),1));
        
        for d=1:length(dose_bins)
            
            cur_vx_area = squeeze(pvx(d,:,6));
            cur_vx_p = squeeze(pvx(d,:,5));
            cur_vx_hr = squeeze(pvx(d,:,7));
            
            cur_pos_corrs = (cur_vx_area==0);
            
            if sum(cur_pos_corrs)==0,
                continue;
            end
            
            %Find median Dx at this volume
            Vx=zeros(CGobj.mNumInGrp,1);
            for k=1:CGobj.mNumInGrp
                Vx(k) = CGobj.mGrp(k).fVolAtDose( dose_bins(d));
            end
            
            [cur_b,~,cur_lr_s]=glmfit(Vx,[ptcomp pttotal],'binomial','link','logit');
            is_corr(d) = cur_b(2)>0;
            lr_vx_p(d) = cur_lr_s.p(2);
             
             
            cur_med_vx = median(Vx);
            
            [~,fvol] = min(abs(volume_bins - cur_med_vx));
            
            if cur_pos_corrs(fvol),
                med_vx_p(d) = cur_vx_p(fvol);
                med_vx_hr(d) = cur_vx_hr(fvol);
            end
            
        end       
        
        

        
        %% Logistic Regression Vx results
        cur_fig=figure(cur_fig_ctr+1000);clf reset;
        set(cur_fig,'Position',ss_four2three);
        h_lr_corr=semilogy(dose_bins(is_corr),lr_vx_p(is_corr),'.','MarkerSize',20);hold on;
        h_lr_acorr=semilogy(dose_bins(~is_corr),lr_vx_p(~is_corr),'r.','MarkerSize',20);
        h_lr_sig=semilogy([0 max(dose_bins)],[0.05 0.05],'g--','LineWidth',2);
        ylim([0 0.38]);
        xlabel('(D_{V}) Volume [cc]','FontSize',22);
        ylabel('Logistic Regression p-value','FontSize',22);
        set(gca,'FontSize',18)
        h_lr_lgnd = legend([h_lr_corr h_lr_acorr h_lr_sig],...
            'Positive Corr.','Negative Corr.','p = 0.05','Location','Best');
        set(h_lr_lgnd,'FontSize',18);
            
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,...
                [fig_basename,'_logreg_vd_pvals'],'-pdf');
            disp(['Saving ',...
                fig_basename,'_logreg_vd_pvals.pdf']);
        end

        
         %% VDx logrank p-value and HR plots
        
        %% VDx results
        
        % only look for minimum p-values in 'reasonable range' (first 70% of range to
        %ingore outliers)
        [min_vx_pval, min_vx_pval_idx] = min(med_vx_p(1:round(length(dose_bins)*0.7)));
        min_vx_pval_hr = med_vx_hr(min_vx_pval_idx);
        
        
        
        %% VDx logrank p-value and HR plots
        cur_fig=figure(cur_fig_ctr);clf reset;
        set(cur_fig,'Position',ss_four2three);
        
        
        [ax,h1,h2]=plotyy(dose_bins,med_vx_p,...
            dose_bins,med_vx_hr,@semilogy);hold on;
        
        semilogy([0 max(dose_bins)],[0.05 0.05],'r--','LineWidth',2);
        h_best_vx_pval_hr=semilogy([dose_bins(min_vx_pval_idx) dose_bins(min_vx_pval_idx)],...
            ylim, 'b--','LineWidth',2);
        
        hold off; % grid on;
        
        cur_lgnd=legend(h_best_vx_pval_hr,['Best p-value: ',num2str(min_vx_pval,'%3.1e'),10,...
            'with HR: ',num2str(min_vx_pval_hr,'%3.1f')],...
            'Location','Best');
        set(cur_lgnd,'FontSize',18);
        
        set(h1,'LineWidth',2);
        set(h2,'LineWidth',2);
        cur_ylim = [ylim];
        set(ax(1),'YLim',[cur_ylim(1) 1]);
        set(get(ax(1),'Ylabel'),'String','Logrank p-value','FontSize',20);
        
        set(get(ax(2),'Ylabel'),'String','Hazard Ratio (HR)','FontSize',20);
        set(ax(1),'XLim',[0 max(dose_bins)]);
        set(ax(2),'XLim',[0 max(dose_bins)]);
        set(ax,'FontSize',18);
        %set(gca,'box','on');
        xlabel('(V_{D}) Dose [Gy]','fontsize',20);
        
        if do_print,
            
            set(cur_fig,'Color','w');
            export_fig(cur_fig,...
                [fig_basename,'_lgrk_vd_pvals_hrs'],'-pdf');
            disp(['Saving ',fig_basename,'_lgrk_vd_pvals_hrs.pdf']);
            %system('pdfcrop fname fname');% to remove borders
        end
        
        
        
        
        %% Best median split KM curve
        dose=dose_bins(min_vx_pval_idx);
        split=-1;
        
        [cur_fig, ~, ~]=fPlotKaplanMeierCurve_VDx({CGobj},dose,split);
        
        set(cur_fig,'Position',ss_four2three);
        grid on;
        set(gca,'GridLineStyle','--')
        
        if do_print,
            
            set(cur_fig,'Color','w');
            export_fig(cur_fig,...
                [fig_basename,'_km_vd_best_med_split'],'-pdf');
            disp(['Saving ',...
                fig_basename,'_km_vd_best_med_split.pdf']);
            %system('pdfcrop fname fname');% to remove borders
        end
        
                
      %% Best COX median split KM curve
      %dose=87.5;%a2b=3
      %dose=44.8;%a2b=10
      dose=25.6;%a2b=Inf
      split=-1;
        
        [cur_fig, ~, ~]=fPlotKaplanMeierCurve_VDx({CGobj},dose,split);
        
        set(cur_fig,'Position',ss_four2three);
        grid on;
        set(gca,'GridLineStyle','--')
        
        if do_print,
            
            set(cur_fig,'Color','w');
            export_fig(cur_fig,...
                [fig_basename,'_km_vd_best_cox_med_split'],'-pdf');
            disp(['Saving ',...
                fig_basename,'_km_vd_best_cox_med_split.pdf']);
            %system('pdfcrop fname fname');% to remove borders
        end
                
        %% Get DVx results
        med_dx_p = nan(length(volume_bins),1);
        med_dx_hr = nan(length(volume_bins),1);
                
        lr_dx_p = nan(length(volume_bins),1);
        is_corr = logical(ones(length(volume_bins),1));
        
        for v=1:length(volume_bins)
            
            cur_dx_area = squeeze(pvx(:,v,6));
            cur_dx_p = squeeze(pvx(:,v,5));
            cur_dx_hr = squeeze(pvx(:,v,7));
            
            cur_pos_corrs = (cur_dx_area==0);
            
            if sum(cur_pos_corrs)==0,
                continue;
            end
            
            %Find median Dx at this volume
            Dx=zeros(CGobj.mNumInGrp,1);
            for k=1:CGobj.mNumInGrp
                Dx(k) = CGobj.mGrp(k).fDoseAtVol( volume_bins(v));
            end
            
            [cur_b,~,cur_lr_s]=glmfit(Dx,[ptcomp pttotal],'binomial','link','logit');
            is_corr(v) = cur_b(2)>0;
            lr_dx_p(v) = cur_lr_s.p(2);
            
            cur_med_dx = median(Dx);
            
            [~,fdose] = min(abs(dose_bins - cur_med_dx));
            
            if cur_pos_corrs(fdose),
                med_dx_p(v) = cur_dx_p(fdose);
                med_dx_hr(v) = cur_dx_hr(fdose);
            end
            
        end
        
        
         %% Logistic Regression Vx results
        cur_fig=figure(cur_fig_ctr+10000);clf reset;
        set(cur_fig,'Position',ss_four2three);
        h_lr_corr=semilogy(volume_bins(is_corr),lr_dx_p(is_corr),'.','MarkerSize',20);hold on;
        h_lr_acorr=semilogy(volume_bins(~is_corr),lr_dx_p(~is_corr),'r.','MarkerSize',20);
        h_lr_sig=semilogy([0 max(volume_bins)],[0.05 0.05],'g--','LineWidth',2);
        ylim([0 0.38]);
        xlabel('(V_{D}) Dose [Gy]','FontSize',22);
        ylabel('Logistic Regression p-value','FontSize',22);
        set(gca,'FontSize',18)
        h_lr_lgnd = legend([h_lr_corr h_lr_acorr h_lr_sig],...
            'Positive Corr.','Negative Corr.','p = 0.05','Location','Best');
        set(h_lr_lgnd,'FontSize',18);
        
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,...
                [fig_basename,'_logreg_dv_pvals'],'-pdf');
            disp(['Saving ',...
                fig_basename,'_logreg_dv_pvals.pdf']);
        end
        
        [min_dx_pval, min_dx_pval_idx] = ...%min(med_dx_p)
            min(med_dx_p(1:round(length(volume_bins)*0.7)));
        min_dx_pval_hr = med_dx_hr(min_dx_pval_idx);
       
        
         %% DVx logrank p-value and HR plots
        cur_fig=figure(cur_fig_ctr+100);clf reset;
        set(cur_fig,'Position',ss_four2three);

        [ax,h1,h2]=plotyy(volume_bins,med_dx_p,...
            volume_bins,med_dx_hr,@semilogy);hold on;
        
        semilogy([0 max(volume_bins)],[0.05 0.05],'r--','LineWidth',2);
        h_best_dx_pval_hr=semilogy([volume_bins(min_dx_pval_idx) volume_bins(min_dx_pval_idx)],...
            ylim, 'b--','LineWidth',2);
        hold off; % grid on;
        
        cur_lgnd=legend(h_best_dx_pval_hr,['Best p-value: ',num2str(min_dx_pval,'%3.1e'),10,...
            'with HR: ',num2str(min_dx_pval_hr,'%3.1f')],...
            'Location','Best');
        set(cur_lgnd,'FontSize',18);
        
        set(h1,'LineWidth',2);
        set(h2,'LineWidth',2);
        cur_ylim = [ylim];
        set(ax(1),'YLim',[cur_ylim(1) 1]);
        set(get(ax(1),'Ylabel'),'String','Logrank p-value','FontSize',20);
        set(get(ax(2),'Ylabel'),'String','Hazard Ratio (HR)','FontSize',20);
        
        set(ax(1),'XLim',[0 max(volume_bins)]);
        set(ax(2),'XLim',[0 max(volume_bins)]);
        set(ax,'FontSize',18);
        xlabel('(D_{V}) Volume [cc]','fontsize',20);
        
        if do_print,
            
            set(cur_fig,'Color','w');
            export_fig(cur_fig,...
                [fig_basename,'_lgrk_dv_pvals_hrs'],'-pdf');
            disp(['Saving ',...
                fig_basename,'_lgrk_dv_pvals_hrs.pdf']);
            %system('pdfcrop fname fname');% to remove borders
        end

        
    
             %% Best median split KM curve
        volume=volume_bins(min_dx_pval_idx);
        split=-1;
        
        [cur_fig, ~, ~]=fPlotKaplanMeierCurve_DVx({CGobj},volume,split);
        
        set(cur_fig,'Position',ss_four2three);
        grid on;
        set(gca,'GridLineStyle','--')
        
        if do_print,
            
            set(cur_fig,'Color','w');
            export_fig(cur_fig,...
                [fig_basename,'_km_dv_best_med_split'],'-pdf');
            disp(['Saving ',...
                fig_basename,'_km_dv_best_med_split.pdf']);
                
            %system('pdfcrop fname fname');% to remove borders
        end
        
        
        %% D3.5
        volume=3.5;
        split=-1;
        
        [cur_fig, ~, ~]=fPlotKaplanMeierCurve_DVx({CGobj},volume,split);
        
        ylim([0 0.3]);
        set(cur_fig,'Position',ss_four2three);
        grid on;
        set(gca,'GridLineStyle','--')
        
        if do_print,
            
            set(cur_fig,'Color','w');
            export_fig(cur_fig,...
                [fig_basename,'_km_d35_med_split'],'-pdf');
            disp(['Saving ',...
                fig_basename,'_km_d35_med_split.pdf']);
                
            %system('pdfcrop fname fname');% to remove borders
        end
        
           %% D05
        volume=5.0;
        split=-1;
        
        [cur_fig, ~, ~]=fPlotKaplanMeierCurve_DVx({CGobj},volume,split);
        
        ylim([0 0.3]);
        set(cur_fig,'Position',ss_four2three);
        grid on;
        set(gca,'GridLineStyle','--')
        
        if do_print,
            
            set(cur_fig,'Color','w');
            export_fig(cur_fig,...
                [fig_basename,'_km_d05_med_split'],'-pdf');
            disp(['Saving ',...
                fig_basename,'_km_d05_med_split.pdf']);
                
            %system('pdfcrop fname fname');% to remove borders
        end
    
              %% Dmax
        volume=-1;% max dose
        split=-1;
        
        [cur_fig, ~, ~]=fPlotKaplanMeierCurve_DVx({CGobj},volume,split);
        
        ylim([0 0.3]);
        set(cur_fig,'Position',ss_four2three);
        grid on;
        set(gca,'GridLineStyle','--')
        
        if do_print,
            
            set(cur_fig,'Color','w');
            export_fig(cur_fig,...
                [fig_basename,'_km_dmax_med_split'],'-pdf');
            disp(['Saving ',...
                fig_basename,'_km_dmax_med_split.pdf']);
                
            %system('pdfcrop fname fname');% to remove borders
        end
    
        
    end

    
    fig_nfz_basename = [fig_loc,'nfz_',...
                        toxicities{i}];
                    

    %% NFZ split
    
    sa=classKaplanMeierCurve(); % initialize a survivalanalysis obj
        
                % survival/complication time
    f2 = ~cellfun('isempty',{CGobj.mGrp.mDateComp}); % patients with no complication date
    f3 = ~cellfun('isempty',{CGobj.mGrp.mDateLastFollowup}); % patients with no last follow up date
    compdate = inf(CGobj.mNumInGrp,1);
    lastfollowup = inf(CGobj.mNumInGrp,1);
    compdate(f2) = ([CGobj.mGrp(f2).mDateComp] - [CGobj.mGrp(f2).mDateBaseline])' / 30;
    lastfollowup(f3) = ([CGobj.mGrp(f3).mDateLastFollowup] - [CGobj.mGrp(f3).mDateBaseline])' / 30;
    compdate = min( lastfollowup, compdate );
    flgcensor = [CGobj.mGrp.mFlgCensor]';

    flg_nfz = [CGobj.mGrp.mFlgNFz]'; % 1 - in NFz, 0 - not in NFz
    flg_nfz = logical(flg_nfz);
     cox_beta=coxphfit(flg_nfz,compdate,'baseline',0,'censoring',flgcensor);
     cox_hr = exp(cox_beta);
     disp(['Hazard Ratio: ',num2str(cox_hr)]);
%                 % assign properties of object sa
     survivedate={compdate(~flg_nfz); compdate(flg_nfz)}; % survive time of each group
     fcensor={flgcensor(~flg_nfz); flgcensor(flg_nfz)}; % censor flag for each group
     sa.mSurvivalTime=survivedate;
     sa.mFlgCensor=fcensor;
     sa.mHR = cox_hr;

%                 % compute survival curves and compare them
     sa=sa.fCalculateSurvivalCurve();
     sa=sa.fCombineSurvivalTime();
     sa=sa.fCompareSurvivalByLogrank();
     pdx=sa.mpValue;
     %
     cur_fig=figure(cur_fig_ctr+200);clf reset;hold on;
     set(cur_fig,'Position',ss_four2three);
     
     sa_times_below=sa.mSurvivalTimeSorted{1};
     sa_curve_below=1-sa.mSurvivalCurve{1};
     sa_times_below(sa_times_below<0)=0;%set negative onset to t=0
     
     sa_times_above=sa.mSurvivalTimeSorted{2};
     sa_curve_above=1-sa.mSurvivalCurve{2};
     sa_times_above(sa_times_above<0)=0;%set negative onset to t=0
     
     h_km(1)=stairs(sa_times_below,sa_curve_below,'LineWidth',2);
     plot(sa_times_below(sa.mCensorStatistics{1}(:,1)),...
         sa_curve_below(sa.mCensorStatistics{1}(:,1)),'+','MarkerSize',20);
     
     h_km(2)=stairs(sa_times_above,sa_curve_above,'r','LineWidth',2);
     plot(sa_times_above(sa.mCensorStatistics{2}(:,1)),...
         sa_curve_above(sa.mCensorStatistics{2}(:,1)),'r+','MarkerSize',20);
     
     str_pval2 = ['Logrank p-value = ',num2str(pdx,'%3.2f\n'),10,...
         'HR = ',num2str(cox_hr,'%3.1f')];
     %text(38,0.25,str_pval2,'FontSize',18,'Location',);
     
     lgnd=legend(h_km,...
         'Not in NFz',...
         'In NFz',...
         'Location','Best');
     set(lgnd,'FontSize',18);
     h_lgnd=legend;
     set(h_lgnd,'interpreter','latex');
     h_lgnd_title = get(h_lgnd,'title');
     set(h_lgnd_title,'string','Tumor Location','FontSize',18);
   
     textbp(str_pval2,'FontSize',22);
     
     set(gca,'xminortick','on','yminortick','on');
     set(gca,'FontSize',18);
     xlabel(['Months'],'fontsize',20);
     ylabel(['Probability of Complication'],'fontsize',20);
    
     if do_print,
            
            set(cur_fig,'Color','w');
            export_fig(cur_fig,...
                [fig_nfz_basename,'_nfz_split'],'-pdf');
            disp(['Saving ',...
                fig_nfz_basename,'_nfz_split.pdf']);
                
            %system('pdfcrop fname fname');% to remove borders
        end
end
end