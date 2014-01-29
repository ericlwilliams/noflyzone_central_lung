function NFz_Alpha2BetaAnaly
tic;close all;
% prepare

% var_str controls what to look at (e.g. max, mean, min)
% these have been saved in files in meta_data/a2b_data/

%var_str = 'max'; 
var_str = 'd35'; 

screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

do_plot_km = false;
do_print = true;
fig_loc = 'Z:/elw/MATLAB/nfz_analy/slides/figures/latest/';

%structures = {'PBT' 'ILUNG' 'ESOPHAGUS' 'HEART' 'NFZ' 'LUNGS'};
structures = {'ESOPHAGUS'};

%toxicities = {'rp','pultox','esotox'};
toxicities = {'esotox'};

fp = 'Z:\elw\MATLAB\nfz_analy\meta_data\';

for i=1:length(toxicities)
    
    for j=1:length(structures)
        tic;
        cur_fig_ctr = (10*i)+j-1;
        
        fprintf('\n');
        disp(['Tox: ',toxicities{i}]);
        disp(['Struct: ',structures{j}]);
        disp(['Counter: ',num2str(cur_fig_ctr)]);
        fprintf('\n');
        
        fig_basename = [fig_loc,'nfz_',...
            structures{j},'_',...
            toxicities{i}];
        %% load data
        if exist('CGobj','var')==1,
            clear CGobj;
        end
        fn = ['NFZ_',structures{j},'_',toxicities{i},'_a2bInf_data.mat'];
        disp(['']);
        disp(['Loading ',fn]);
        disp(['']);
        load(strcat(fp,fn),'CGobj_org');
        CGobj = CGobj_org;
        clear CGobj_org;
        
        %% load a2b data
        fn2 = ['a2b_data\NFZ_',structures{j},'_',toxicities{i},'_a2b_data.mat'];
        load(strcat(fp,fn2)); % a2b_dmax, a2b_range, a2b_dmean, a2b_d05, a2b_d35
        
        
        % survival/complication time
        f2 = ~cellfun('isempty',{CGobj.mGrp.mDateComp}); % patients with no complication date
        f3 = ~cellfun('isempty',{CGobj.mGrp.mDateLastFollowup}); % patients with no last follow up date
        compdate = inf(CGobj.mNumInGrp,1);
        lastfollowup = inf(CGobj.mNumInGrp,1);
        compdate(f2) = ([CGobj.mGrp(f2).mDateComp] - [CGobj.mGrp(f2).mDateBaseline])' / 30;
        lastfollowup(f3) = ([CGobj.mGrp(f3).mDateLastFollowup] - [CGobj.mGrp(f3).mDateBaseline])' / 30;
        compdate = min( lastfollowup, compdate );
        flgcensor = [CGobj.mGrp.mFlgCensor]';
        
        
        
        
        pttotal = ones(CGobj.mNumInGrp,1);
        ptcomp = ones(CGobj.mNumInGrp,1);
        ptcomp([CGobj.mGrp.mFlgCensor])=0;
        
        %a2b_max_doses = [cellfun(@(x) max(x),a2b_doses)];
        
        sa=classKaplanMeierCurve(); % initialize a survivalanalysis obj
        
        if isequal(var_str,'dmax')
            a2b_var = a2b_dmax;
        elseif isequal(var_str,'dmean')
            a2b_var = a2b_dmean;
        elseif isequal(var_str,'d05')
            a2b_var = a2b_d05;
        elseif isequal(var_str,'d35')
            a2b_var = a2b_d35;
        else
            disp(['Unknown variable: ' var_str]);
        end
            
        
        lr_pvals = Inf(size(a2b_var,2),1);
        km_hrs = zeros(size(a2b_var,2),1);
        
        pvals = Inf(size(a2b_var,2),1);
        llhds = Inf(size(a2b_var,2),1);
        
        fig_ctr=1;
        %loop over all a2b valus
        for k=1:size(a2b_var,2)
            cur_a2b = a2b_range(k);
            cur_a2b_var = [a2b_var(:,k)];
            
            [b,dev,s]=glmfit(cur_a2b_var,[ptcomp pttotal],'binomial','link','logit');
            pvals(k) = s.p(2);
            
            B0 = b(1);
            B1 = b(2);
            pr = exp(B0+B1*cur_a2b_var);
            
%             %tmp
%             dpf = dev; % deviations
%             df = [s.dfe]; % degree of freedom
%             dpf = dpf./df; % deviations per degree of freedom
%             llhds(k) = -0.5*dpf;
            
            pr = pr./(1+pr); % logistic probability
            pr(flgcensor) = 1-pr(flgcensor); % non-complication patients
            pr = log(pr); % log likelihood of each patients
            llhds(k) = sum(pr); % loglikelihood of all
            
            
%             %% KM split
%             
%             med_a2b_var = median(cur_a2b_var);
%             flg_below=cur_a2b_var<=med_a2b_var;
%             
%             survivedate={compdate(flg_below); compdate(~flg_below)}; % survive time of each group
%             fcensor={flgcensor(flg_below); flgcensor(~flg_below)}; % censor flag for each group
%             sa.mSurvivalTime=survivedate;
%             sa.mFlgCensor=fcensor;
%             % compute survival curves and compare them
%             sa=sa.fCalculateSurvivalCurve();
%             sa=sa.fCombineSurvivalTime();
%             sa=sa.fCompareSurvivalByLogrank();
%             cur_area = sa.mCurveArea;
%             cur_pval=sa.mpValue;
%             
%             %     if sum(~flgcensor(flg_below))
%             %         cox_beta=coxphfit(~flg_below,compdate,'baseline',0,'censoring',flgcensor);
%             %         km_hrs(i) = exp(cox_beta);
%             %     end
%             cox_beta=coxphfit(~flg_below,compdate,'baseline',0,'censoring',flgcensor);
%             cur_hr = exp(cox_beta);
%             if cur_hr<=1,
%                 cur_hr=Inf;
%             end
%             
%             if ~isequal(cur_area,0),
%                 cur_hr=Inf;
%                 cur_pval=Inf;
%             end
%             
%             km_hrs(k) = cur_hr;
%             lr_pvals(k) = cur_pval;
%             %      % plot
%             if (k==1 || k==size(a2b_var,2) || mod(a2b_range(k),1)==0) && do_plot_km,
%                 
%                 cur_fig=figure(100+fig_ctr);clf reset;hold on;
%                 set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
%                 
%                 h_km(1)=stairs(sa.mSurvivalTimeSorted{1},1-sa.mSurvivalCurve{1});
%                 plot(sa.mSurvivalTimeSorted{1}(sa.mCensorStatistics{1}(:,1)),...
%                     1-sa.mSurvivalCurve{1}(sa.mCensorStatistics{1}(:,1)),'+');
%                 h_km(2)=stairs(sa.mSurvivalTimeSorted{2},1-sa.mSurvivalCurve{2},'r');
%                 plot(sa.mSurvivalTimeSorted{2}(sa.mCensorStatistics{2}(:,1)),...
%                     1-sa.mSurvivalCurve{2}(sa.mCensorStatistics{2}(:,1)),'r+');
%                 %         xticks = get(gca,'Xlim'); set(gca,'XTick',0:6:max(xticks));
%                 
%                 str_pval1 = ['Log-Rank p-value = ',num2str(sa.mpValue,2),10,...
%                     '\alpha/\beta = ',num2str(a2b_range(k))];
%                 %text(.6,0.7,str_pval1,...
%                 %    'FontSize',16,'Units','normalized','BackgroundColor','w');
%                 textbp(str_pval1,'FontSize',18);
%                 lgnd=legend(h_km,...
%                     ['D$_{\rm{max}} \leq',num2str(med_a2b_var,4),...
%                     '~\rm{Gy}_{',num2str(a2b_range(k)),'}$'],...
%                     ['D$_{\rm{max}} > ',num2str(med_a2b_var,4),...
%                     '~\rm{Gy}_{',num2str(a2b_range(k)),'}$'],...
%                     'Location','SouthEast');
%                 set(lgnd,'FontSize',18);
%                 h=legend;
%                 set(h,'interpreter','latex');
%                 
%                 set(gca,'xminortick','on','yminortick','on');
%                 xlabel(['Months'],'fontsize',20);
%                 ylabel(['Probability of Complication'],'fontsize',20);
%                 %title('Brachial Plexopathy Incidence','fontsize',14);
%                 
%                 
%                 if do_print,
%                     set(cur_fig,'Color','w');
%                     export_fig(cur_fig,[fig_basename,'_a2b_km-',num2str(fig_ctr)],'-pdf');
%                     disp(['Saving ',fig_basename,'_a2b_km-',num2str(fig_ctr),'.pdf...']);
%                 end;
%                 
%                 
%                 fig_ctr=fig_ctr+1;
%                 
%             end
        end

        x_axis = a2b_range(1:end-1);
        phys_pval = pvals(end);

        cur_fig=figure(1);clf reset;hold on;
        set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
        plot(x_axis,pvals(1:end-1),'LineWidth',2);% excluding Inf
        h_phys_pval=plot(x_axis,repmat(phys_pval,length(x_axis)),'r--','LineWidth',2);
        [min_pval,min_idx] = min(pvals(1:end-1));
        str_pval2 = ['Min p-value of ',num2str(min_pval,'%3.2e'),10,'at \alpha/\beta = ',num2str(a2b_range(min_idx))];
        %text(.1,0.8,...
        %    ['Min p-value of ',num2str(min_pval),10,'at \alpha/\beta = ',num2str(a2b_range(min_idx))],...
        %    'FontSize',16,'Units','normalized','BackgroundColor','w');
        lgnd_pval=legend(h_phys_pval(1),['Physical dose p-value: ',num2str(phys_pval,'%3.2e')]);
        set(lgnd_pval,'FontSize',18);
        text(0.05,0.9,str_pval2,'FontSize',18,'units','normalized','BackgroundColor','w');

        set(gca,'FontSize',16);
        ylabel('p-value','FontSize',20);
        xlabel('\alpha/\beta value [Gy]','FontSize',20);
        %title(['Logistic Regression',10,'BPx vs. Max BED for given \alpha/\beta'],'FontSize',16);
        grid on;
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_basename,'_a2b_',var_str,'_pvals'],'-pdf');
            disp(['Saving ',fig_basename,'_a2b_',var_str,'_pvals.pdf...']);
        end;
        
        
        cur_fig=figure(2);clf reset;hold on;
        set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
        phys_llhd = llhds(end);
        plot(x_axis,llhds(1:end-1),'LineWidth',2);% excluding Inf
        hold on;
        [max_llhd,max_idx] = max(llhds(1:end-1));
        lowCI68 = max_llhd - 0.5; % 68% confidence

        plot([0 max(x_axis)],[lowCI68 lowCI68],'g--','LineWidth',2);
        plot([x_axis(max_idx) x_axis(max_idx)],ylim,'b--','LineWidth',2);
        h_phys_llhd=plot(x_axis,repmat(phys_llhd,length(x_axis)),'r--','LineWidth',2);
        str_ll = ['Max LogL of ',num2str(max_llhd,'%3.2e'),10,'at \alpha/\beta = ',num2str(a2b_range(max_idx))];
        %text(.1,0.4,...
        %    ['Max LogL of ',num2str(max_llhd),10,'at \alpha/\beta = ',num2str(a2b_range(max_idx))],...
        %    'FontSize',16,'Units','normalized','BackgroundColor','w');
        text(0.05,0.9,str_ll,'FontSize',18,'units','normalized','BackgroundColor','w');
        
        lgnd_llhd=legend(h_phys_llhd(1),['Physical dose',10,'Log-likelihood: ',num2str(phys_llhd,'%3.3g')]);
        set(lgnd_llhd,'FontSize',18);
        set(gca,'FontSize',16);
        ylabel('Log-Likelihood','FontSize',20);
        xlabel('\alpha/\beta value [Gy]','FontSize',20);
        %title(['Logistic regression',10,'BPx vs. Max BED for given \alpha/\beta'],'FontSize',16);
        grid on;
        
        
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_basename,'_a2b_',var_str,'_llhds'],'-pdf');
            disp(['Saving ',fig_basename,'_a2b_',var_str,'_llhds.pdf...']);
        end;
        
        
%        
%         
%         cur_fig=figure(3);clf reset;hold on;
%         set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
%         
%         [ax,h1,h2]=plotyy(x_axis,lr_pvals(1:end-1),x_axis,log(km_hrs(1:end-1)),@semilogy);hold on;
%         set(h1,'LineWidth',2);
%         set(h2,'LineWidth',2);
%         %set(ax(1),'YLim',[0.001 1]);
%         set(get(ax(1),'Ylabel'),'String','p-value','FontSize',16);
%         [min_pval,min_idx] = min(lr_pvals(1:end-1));
%         str_pval3 = ['Min p-value of ',num2str(min_pval,3),10,'at \alpha/\beta = ',num2str(a2b_range(min_idx))];
%         %text(.5,0.6,...
%         %    ['Min p-value of ',num2str(min_pval,3),10,'at \alpha/\beta = ',num2str(a2b_range(min_idx))],...
%         %    'FontSize',16,'Units','normalized','BackgroundColor','w');
%         textbp(str_pval3,'FontSize',18);
%         %set(ax(2),'YLim',[ylim(1) max(log(km_hrs(1:end-1)))*1.2]);
%         set(get(ax(2),'Ylabel'),'String','ln(Hazard Ratio)','FontSize',20);
%         %semilogy(x_axis,repmat(0.05,length(x_axis)),'r--','LineWidth',1);
%         set(ax(1),'XLim',[min(x_axis) max(x_axis)]);
%         set(ax(2),'XLim',[min(x_axis) max(x_axis)]);
%         set(ax,'FontSize',16);
%         
%         hold off; % grid on;
%         xlabel('\alpha/\beta value [Gy]','fontsize',20);
%         
%         
%         if do_print,
%             set(cur_fig,'Color','w');
%             export_fig(cur_fig,[fig_basename,'_a2b_lr_pvals_km_hr'],'-pdf');
%             disp(['Saving ',fig_basename,'_a2b_lr_pvals_km_hr.pdf...']);
%         end;
        
    end
end
end
