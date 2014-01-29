function NFZ_gEUD
tic;close all;

%% TMP
do_special_tox=false;


screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];
do_print=true;

fig_loc = 'Z:/elw/MATLAB/nfz_analy/slides/figures/latest/';

%structures = {'ILUNG' 'ESOPHAGUS' 'HEART' 'NFZ' 'PBT' 'LUNGS'};
%str_colors = {'k' 'm' 'b' 'r' 'g' 'c'};

structures = {'ESOPHAGUS'};
str_colors = {'k'};

%toxicities = {'pultox'};
%toxicities = {'rp','pultox','esotox'};
toxicities = {'esotox'};

comp_rate = 0.2;

fp = 'Z:\elw\MATLAB\nfz_analy\meta_data\';

%a2b = {'3'};
%cur_a2b = 3;
a2b = {'10'};
cur_a2b = 10;


for i=1:length(toxicities)
    
    %% pvalues
    figure(i); clf reset;
    set(gcf,'Position',ss_four2three);
    %title(['Logistic Regression P-values',10,tox_titles{i}],'FontSize',18);
    
    %% loglikelihoods
    figure(i+100); clf reset;
    set(gcf,'Position',ss_four2three);
    %title(['Logistic Regression Log-likelihoods',10,tox_titles{i}],'FontSize',18);
    fig_tox_basename = [fig_loc,'nfz_',toxicities{i},'_a2b',a2b{1}];
    
      %% tvalues
    figure(i+1000); clf reset;
    set(gcf,'Position',ss_four2three);
    h_pvals = inf(length(structures),1);
    h_tvals = inf(length(structures),1);
    h_llhds = inf(length(structures),1);
    
    for j=1:length(structures)
        
        cur_fig_ctr = (10*i)+j-1;
        
        fprintf('\n');
        disp(['Tox: ',toxicities{i}]);
        disp(['Struct: ',structures{j}]);
        
        fprintf('\n');
        
        fig_basename = [fig_loc,'nfz_',...
            structures{j},'_',...
            toxicities{i},'_a2b',...
            a2b{1}];
        
        
        %% load data
        if ~do_special_tox,
            fn = ['NFZ_',structures{j},'_',toxicities{i},'_a2b',a2b{1},'_data.mat'];
        else
            fn = ['NFZ_',structures{j},'_',toxicities{i},'_a2b',a2b{1},'_geq3_data.mat'];
        end
        disp(['Loading: ',fn,'...']);
        
        load(strcat(fp,fn),'CGobj_org');
        CGobj = CGobj_org;
        LymanN = log10(CGobj.mLymanN);
        CGobj.mLymanN = LymanN;
        
        
        figure(i);hold on;
        [h_pvals(j) min_lr_pvals(j)]=CGobj.fLogisticRegressionPvalueExactFig_a_EUD(strcat(str_colors{j},'s--'),2);
        %%loglikelihoods
         figure(i+1000);hold on;
        h_tvals(j) =CGobj.fLogisticRegressionTvalueExactFig_a_EUD(strcat(str_colors{j},'s--'),2);
        
        figure(i+100);
        [h_llhds(j) min_lr_llhds(j) max_lr_llhds(j)]=CGobj.fLogisticRegressionLikelyhoodExactFig_a_EUD('loga',strcat(str_colors{j},'s--'),2);
        

            
        %% Mean Dose
        cur_fig=figure(cur_fig_ctr+200);
        set(gcf,'Position',ss_four2three);
        %title(['Logistic Regression Mean-Dose response',10,tox_titles{i}],'FontSize',18);
        [loga,pval,~] = CGobj.fLogisticRegressionGridRespondingCurveExactFig_a_EUD(0,'k',1);
        
        set(gca,'xminortick','on','yminortick','on');
        set(gca,'box','on');
        ylabel('Complication probability');
        
        %%tmp
        %[loga,pval] = CGobj.fLogisticRegressionRespondingCurveExactFig_a_EUD(0,'k',1);
        %[loga,pval] = CGobj.fLogisticRegressionRespondingCurveExactFig_a_EUD(0,'r',1,'MM');
        %%tmp
        
        CGobj.fComplicationObservedFig_EUD(loga,4,'k*',2);
        ylim([0 0.5]);
        
        xlabel('Mean Dose [Gy]','FontSize',20);
        loga_str = ['log_1_0(a) = %s',10,'p-val: %s',10,'%s'];
        %str = sprintf(loga_str,num2str(loga),num2str(pval,2),structures{j});
        str = sprintf(loga_str,num2str(loga),num2str(pval,'%3.2e'));
        text(0.65,0.85,str,'FontSize',24,'Units','normalized');
       
      
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,...
                [fig_basename,'_dmean_lr'],'-pdf');
            disp(['Saving ',fig_basename,'_dmean_lr.pdf']);
            
        end
        
           %tmp

        
        
        %% Best loga
        cur_fig=figure(cur_fig_ctr+300);
        set(gcf,'Position',ss_four2three);
        %title(['Logistic Regression Best gEUD response',10,tox_titles{i}],'FontSize',18);
        
        [loga,pval,~] = CGobj.fLogisticRegressionGridRespondingCurveExactFig_a_EUD('loga','k',1);
        CGobj.fComplicationObservedFig_EUD(loga,4,'k*',2);
        
        ylim([0 0.5]);
        set(gca,'xminortick','on','yminortick','on');
        set(gca,'box','on');
        ylabel('Complication probability');
        
        loga_str = ['log_1_0(a) = %s',10,'p-val: %s',10,'%s'];
        %str = sprintf(loga_str,num2str(loga),num2str(pval,2),structures{j});
        str = sprintf(loga_str,num2str(loga),num2str(pval,'%3.2e'));
        text(0.65,0.85,str,'FontSize',24,'Units','normalized');
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,...
                [fig_basename,'_geud_lr'],'-pdf');
            disp(['Saving ',fig_basename,'_geud_lr.pdf']);
        end
        
        
        %% D_Max
        
        cur_fig=figure(cur_fig_ctr+400);
        set(gcf,'Position',ss_four2three);
        
        pttotal = ones(CGobj.mNumInGrp,1);
        ptcomp = ones(CGobj.mNumInGrp,1); ptcomp([CGobj.mGrp.mFlgCensor])=0;
    
        d_bins = [CGobj.mGrp.mDoseBins_LQ];
        v_bins = [CGobj.mGrp.mVolCum];
    
        %dmax = [max([CGobj.mGrp.mDoseBins_LQ])]';
        
        d05s = zeros(size(v_bins,2),1);
        d35s = zeros(size(v_bins,2),1);
        dmax = zeros(size(v_bins,2),1);
    
        for k=1:size(v_bins,2) % for each patient
            vol = v_bins(:,k);
            
            %vol = vol./max(vol);
            %v05 = vol<0.05;
            
            
            %v05 = vol<5; %% < 5 cc
            %v05 = vol<3.5; %% < 5 cc
    
            v05 = vol<5; %% < 5 cc
            d05_inds = find(v05);
            dose = d_bins(:,k);
            d05 = min(dose(d05_inds));
            d05s(k) = d05;
            
            v35 = vol<3.5; %% < 5 cc
            d35_inds = find(v35);
            %dose = d_bins(:,j);
            d35 = min(dose(d35_inds));
            d35s(k) = d35;
            
            vol(~vol)=nan;
            nan_inds = find(isnan(vol));
            if ~isempty(nan_inds)
                min_ind = nan_inds(1)-1;
            else
                min_ind = length(vol);
            end
            %[~,min_ind] = min(vol);
            dmax(k) = dose(min_ind);
        end
        
        
        
        % regression using exact EUD
        [b,dev,s]=glmfit(dmax,[ptcomp pttotal],'binomial','link','logit');
        %CGobj.mLogisticRegressionMat(k).b=b;
        %            CGobj.mLogisticRegressionMat(k).dev=dev;
        %            CGobj.mLogisticRegressionMat(k).stats=s;
        
        pvalue = [s.p];
        pval = pvalue(2); % the p-value corresponding to gEUD
        
        %tmp
        doses = (0:max(dmax))'; % doses (gEUD) whose RP probability will be computed
        %doses = (0:90)';
        [rpb,rplo,rphi] = glmval(b, doses,'logit',s); % the responding function values at doses
                
        % plot
        hold on;
        plot(doses,rpb,'k','LineWidth',2); % responding function
        plot(doses,rpb-rplo,'k','LineWidth',1); % low CI curve
        plot(doses,rpb+rphi,'k','LineWidth',1); % high CI curve
        
        
        
        dmax_lim = interp1(rpb,doses,comp_rate);
        dmax_limhi = interp1(rpb-rplo,doses,comp_rate);
        dmax_limlo = interp1(rpb+rphi,doses,comp_rate);
        
        if (~isnan(dmax_lim) && ~isnan(dmax_limhi) && ~isnan(dmax_limlo)) &&...
                (~isinf(dmax_lim) && ~isinf(dmax_limhi) && ~isinf(dmax_limlo)),
            
            dmax_phys = [max(roots([1 cur_a2b*3 -cur_a2b*3*dmax_lim])),...
                max(roots([1 cur_a2b*4 -cur_a2b*4*dmax_lim])),...
                max(roots([1 cur_a2b*5 -cur_a2b*5*dmax_lim]))];
            dmax_phys_hi = [max(roots([1 cur_a2b*3 -cur_a2b*3*dmax_limhi])),...
                max(roots([1 cur_a2b*4 -cur_a2b*4*dmax_limhi])),...
                max(roots([1 cur_a2b*5 -cur_a2b*5*dmax_limhi]))];
            dmax_phys_lo = [max(roots([1 cur_a2b*3 -cur_a2b*3*dmax_limlo])),...
                max(roots([1 cur_a2b*4 -cur_a2b*4*dmax_limlo])),...
                max(roots([1 cur_a2b*5 -cur_a2b*5*dmax_limlo]))];
        else
            dmax_phys =[0 0 0];
            dmax_phys_hi =[0 0 0];
            dmax_phys_lo =[0 0 0];
        end
        
        
        disp('******');
        disp(['D_{max} corresponding to ',num2str(comp_rate*100),'% rate: ',10,...
            num2str(dmax_lim),...
            ' [',num2str(dmax_limlo),', ',num2str(dmax_limhi),']']);
        disp(['= Physical doses =']);
        disp(['Nfx = 3',10,...
            num2str(dmax_phys(1)),...
            ' [',num2str(dmax_phys_lo(1)),', ',num2str(dmax_phys_hi(1)),']']);
        disp(['Nfx = 4',10,...
            num2str(dmax_phys(2)),...
            ' [',num2str(dmax_phys_lo(2)),', ',num2str(dmax_phys_hi(2)),']']);
        disp(['Nfx = 5',10,...
            num2str(dmax_phys(3)),...
            ' [',num2str(dmax_phys_lo(3)),', ',num2str(dmax_phys_hi(3)),']']);   
        disp([]);
        disp('******');
        disp([]);

        flg=[CGobj.mGrp.mFlgCensor]; % censor flags of patients
        
        [medianeud,~,~,binlow,binhigh,numcomp,numtotal,betainv84,betainv16] = EventObserved(flg,dmax,4);
        prob = numcomp./numtotal;
        % plot
        errorbar(medianeud,prob,max(0,prob-betainv16),max(0,betainv84-prob),'k*','LineWidth',1);
        errorbar_x(medianeud,prob,(medianeud-binlow),(binhigh-medianeud),'k*');
        
        %tmp
        ylim([0 0.5]);
        xlim([0 max(binhigh)*1.1]);
        
        %xlims = xlim;
        %xmax = xlims(2);
        %xlim([0 xmax]);
        %ylim([0 0.5]);
        
        loga_str = ['p-val: %s'];
        %str = sprintf(loga_str,num2str(loga),num2str(pval,2),structures{j});
        str = sprintf(loga_str,num2str(pval,'%3.2e'));
        %text(0.65,0.85,str,'FontSize',24,'Units','normalized');
        textbp(str,'FontSize',24);
        
        set(gca,'xminortick','on','yminortick','on');
        set(gca,'box','on');
        set(gca,'FontSize',16);
        xlabel('D_{max} [Gy_{10}]','FontSize',25);
        ylabel('Complication rate observed','FontSize',22);
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,...
                [fig_basename,'_lr_dmax'],'-pdf');
            disp(['Saving ',fig_basename,'_lr_dmax.pdf']);
        end
        
        %% rate vs dmax
        rates = 0:0.01:.3;
        
          rate_vs_dmax = inf(length(rates),1);
        rate_vs_dmax_hi = inf(length(rates),1);
        rate_vs_dmax_lo = inf(length(rates),1);
        
        rate_vs_phys_dmax = inf(length(rates),3);
        rate_vs_phys_dmax_hi = inf(length(rates),3);
        rate_vs_phys_dmax_lo = inf(length(rates),3);
        
        for l=1:length(rates)
            rate_vs_dmax(l) = interp1(rpb,doses,rates(l));
            rate_vs_dmax_hi(l) = interp1(rpb-rplo,doses,rates(l));
            rate_vs_dmax_lo(l) = interp1(rpb+rphi,doses,rates(l));
            
            if ~isnan(rate_vs_dmax(l)),
            rate_vs_phys_dmax(l,1) = max(roots([1 cur_a2b*3 -cur_a2b*3*rate_vs_dmax(l)]));
            rate_vs_phys_dmax(l,2) =  max(roots([1 cur_a2b*4 -cur_a2b*4*rate_vs_dmax(l)]));
            rate_vs_phys_dmax(l,3) =  max(roots([1 cur_a2b*5 -cur_a2b*5*rate_vs_dmax(l)]));
            end
            
            if ~isnan(rate_vs_dmax_hi(l)),
            rate_vs_phys_dmax_hi(l,1) = max(roots([1 cur_a2b*3 -cur_a2b*3*rate_vs_dmax_hi(l)]));
            rate_vs_phys_dmax_hi(l,2) =  max(roots([1 cur_a2b*4 -cur_a2b*4*rate_vs_dmax_hi(l)]));
            rate_vs_phys_dmax_hi(l,3) =  max(roots([1 cur_a2b*5 -cur_a2b*5*rate_vs_dmax_hi(l)]));
            end
            
            if ~isnan(rate_vs_dmax_lo(l)),
            rate_vs_phys_dmax_lo(l,1) = max(roots([1 cur_a2b*3 -cur_a2b*3*rate_vs_dmax_lo(l)]));
            rate_vs_phys_dmax_lo(l,2) =  max(roots([1 cur_a2b*4 -cur_a2b*4*rate_vs_dmax_lo(l)]));
            rate_vs_phys_dmax_lo(l,3) =  max(roots([1 cur_a2b*5 -cur_a2b*5*rate_vs_dmax_lo(l)]));
            end
            
        end
        
        cur_fig=figure(cur_fig_ctr+450);
        set(gcf,'Position',ss_four2three);
        hold on;
        h_dmax(1)=plot(rates,rate_vs_dmax,'LineWidth',2);
        h_dmax(2)=plot(rates,rate_vs_dmax_hi,'--');
        plot(rates,rate_vs_dmax_lo,'--');
        grid on;
        set(gca,'GridLineStyle','--')
        set(gca,'xminortick','on','yminortick','on');
        set(gca,'box','on');
        set(gca,'FontSize',18);
        ylabel(['D_{\rm{max}} Gy_{',a2b{1},'}'],'FontSize',20);
        xlabel('Complication Rate (%)','FontSize',20);
        [~,ten_pct_ind] = min(abs(rates-0.1));
        [~,fftn_pct_ind] = min(abs(rates-0.15));
        [~,twnty_pct_ind] = min(abs(rates-0.2));
        
        text(0.05,0.9,['10% Comp. rate -> ',...
            num2str(rate_vs_dmax(ten_pct_ind),'%3.1f'),' Gy_{',num2str(a2b{1}),'} [',...
            num2str(rate_vs_dmax_lo(ten_pct_ind),'%3.1f'),' - ',...
            num2str(rate_vs_dmax_hi(ten_pct_ind),'%3.1f'),']',10,...
            '15% Comp. rate -> ',...
            num2str(rate_vs_dmax(fftn_pct_ind),'%3.1f'),' Gy_{',num2str(a2b{1}),'} [',...
            num2str(rate_vs_dmax_lo(fftn_pct_ind),'%3.1f'),' - ',...
            num2str(rate_vs_dmax_hi(fftn_pct_ind),'%3.1f'),']',10,...
            '20% Comp. rate -> ',...
            num2str(rate_vs_dmax(twnty_pct_ind),'%3.1f'),' Gy_{',num2str(a2b{1}),'} [',...
            num2str(rate_vs_dmax_lo(twnty_pct_ind),'%3.1f'),' - ',...
            num2str(rate_vs_dmax_hi(twnty_pct_ind),'%3.1f'),']'],...
            'units','normalized','FontSize',18,'BackgroundColor','w');
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,...
                [fig_basename,'_rate_vs_bed_dmax'],'-pdf');
            disp(['Saving ',fig_basename,'_rate_vs_bed_dmax.pdf']);
        end
        
        
        %% rate vs dmax physical doses
         cur_fig=figure(cur_fig_ctr+4500);
        set(gcf,'Position',ss_four2three);
        hold on;
        h_dmax_nfx3=plot(rates,rate_vs_phys_dmax(:,1),'LineWidth',2);
        plot(rates,rate_vs_phys_dmax_hi(:,1),'--');
        plot(rates,rate_vs_phys_dmax_lo(:,1),'--');
        h_dmax_nfx4=plot(rates,rate_vs_phys_dmax(:,2),'g','LineWidth',2);
        plot(rates,rate_vs_phys_dmax_hi(:,2),'g--');
        plot(rates,rate_vs_phys_dmax_lo(:,2),'g--');
        h_dmax_nfx5=plot(rates,rate_vs_phys_dmax(:,3),'r','LineWidth',2);
        plot(rates,rate_vs_phys_dmax_hi(:,3),'r--');
        plot(rates,rate_vs_phys_dmax_lo(:,3),'r--');
        grid on;
        legend([h_dmax_nfx3 h_dmax_nfx4 h_dmax_nfx5],...
            {'3 Fx','4 Fx','5 Fx'},'FontSize',18);
            
        set(gca,'GridLineStyle','--')
        set(gca,'xminortick','on','yminortick','on');
        set(gca,'box','on');
        set(gca,'FontSize',18);
        ylabel(['D_{\rm{max}} Gy_{PHYS}'],'FontSize',20);
        xlabel('Complication Rate (%)','FontSize',20);

        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,...
                [fig_basename,'_rate_vs_phys_dmax'],'-pdf');
            disp(['Saving ',fig_basename,'_rate_vs_phys_dmax.pdf']);
        end
        
        %% D_05
        
        cur_fig=figure(cur_fig_ctr+500);
        set(gcf,'Position',ss_four2three);
        
        %pttotal = ones(CGobj.mNumInGrp,1);
        %ptcomp = ones(CGobj.mNumInGrp,1); ptcomp([CGobj.mGrp.mFlgCensor])=0;
        
        
       
        
        % regression using exact EUD
        [b,dev,s]=glmfit(d05s,[ptcomp pttotal],'binomial','link','logit');
        %CGobj.mLogisticRegressionMat(k).b=b;
        %            CGobj.mLogisticRegressionMat(k).dev=dev;
        %            CGobj.mLogisticRegressionMat(k).stats=s;
        
        pvalue = [s.p];
        pval = pvalue(2); % the p-value corresponding to gEUD
        
        %tmp
        doses = (0:max(d05s))'; % doses (gEUD) whose RP probability will be computed
        %doses = (0:90)';
        [rpb,rplo,rphi] = glmval(b, doses,'logit',s); % the responding function values at doses
        
        % plot
        hold on;
        plot(doses,rpb,'k','LineWidth',2); % responding function
        plot(doses,rpb-rplo,'k','LineWidth',1); % low CI curve
        plot(doses,rpb+rphi,'k','LineWidth',1); % high CI curve
        
       
        
        d05_lim = interp1(rpb,doses,comp_rate);
        d05_limhi = interp1(rpb-rplo,doses,comp_rate);
        d05_limlo = interp1(rpb+rphi,doses,comp_rate);
         
        if (~isnan(d05_lim) && ~isnan(d05_limhi) && ~isnan(d05_limlo)) &&...
                (~isinf(d05_lim) && ~isinf(d05_limhi) && ~isinf(d05_limlo)),
            
        d05_phys = [max(roots([1 cur_a2b*3 -cur_a2b*3*d05_lim])),...
                    max(roots([1 cur_a2b*4 -cur_a2b*4*d05_lim])),...
                    max(roots([1 cur_a2b*5 -cur_a2b*5*d05_lim]))];
        d05_phys_hi = [max(roots([1 cur_a2b*3 -cur_a2b*3*d05_limhi])),...
                        max(roots([1 cur_a2b*4 -cur_a2b*4*d05_limhi])),...
                        max(roots([1 cur_a2b*5 -cur_a2b*5*d05_limhi]))];
        d05_phys_lo = [max(roots([1 cur_a2b*3 -cur_a2b*3*d05_limlo])),...
                        max(roots([1 cur_a2b*4 -cur_a2b*4*d05_limlo])),...
                        max(roots([1 cur_a2b*5 -cur_a2b*5*d05_limlo]))];
                
           else
            d05_phys =[0 0 0];
            d05_phys_hi =[0 0 0];
            d05_phys_lo =[0 0 0];
        end
        
        
        disp([]);
        disp('******');
        disp(['D_{05 cc} corresponding to ',num2str(comp_rate*100),'% rate: ',10,...
            num2str(d05_lim),...
            ' Gy_{',num2str(a2b{1}),'} [',num2str(d05_limlo),', ',num2str(d05_limhi),']']);
        disp(['= Physical doses =']);
        disp(['Nfx = 3',10,...
            num2str(d05_phys(1)),...
            ' Gy_{',num2str(a2b{1}),'} [',num2str(d05_phys_lo(1)),', ',num2str(d05_phys_hi(1)),']']);
        disp(['Nfx = 4',10,...
            num2str(d05_phys(2)),...
            ' Gy_{',num2str(a2b{1}),'} [',num2str(d05_phys_lo(2)),', ',num2str(d05_phys_hi(2)),']']);
        disp(['Nfx = 5',10,...
            num2str(d05_phys(3)),...
            ' Gy_{',num2str(a2b{1}),'} [',num2str(d05_phys_lo(3)),', ',num2str(d05_phys_hi(3)),']']);
                
            
        disp('******');
        disp([]);
        flg=[CGobj.mGrp.mFlgCensor]; % censor flags of patients
        
        [medianeud,~,~,binlow,binhigh,numcomp,numtotal,betainv84,betainv16] = EventObserved(flg,d05s,4);
        prob = numcomp./numtotal;
        % plot
        errorbar(medianeud,prob,max(0,prob-betainv16),max(0,betainv84-prob),'k*','LineWidth',1);
        errorbar_x(medianeud,prob,(medianeud-binlow),(binhigh-medianeud),'k*');
        
        %tmp
        ylim([0 0.5]);
        xlim([0 max(binhigh)*1.1]);
        
        %xlims = xlim;
        %xmax = xlims(2);
        %xlim([0 xmax]);
        %ylim([0 0.5]);
        
        loga_str = ['p-val: %s'];
        %str = sprintf(loga_str,num2str(loga),num2str(pval,2),structures{j});
        str = sprintf(loga_str,num2str(pval,'%3.2e'));
        %text(0.65,0.85,str,'FontSize',24,'Units','normalized');
        textbp(str,'FontSize',24);
        
        set(gca,'xminortick','on','yminortick','on');
        set(gca,'box','on');
        set(gca,'FontSize',16);
        %xlabel('D_{05 cc} [Gy]','FontSize',25);
        xlabel('D_{5 cc} [Gy_{10}]','FontSize',25);
        ylabel('Complication rate observed','FontSize',22);
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,...
                [fig_basename,'_lr_d05'],'-pdf');
            disp(['Saving ',fig_basename,'_lr_d05.pdf']);
        end
        
        %% rate vs d05
          rate_vs_d05 = inf(length(rates),1);
        rate_vs_d05_hi = inf(length(rates),1);
        rate_vs_d05_lo = inf(length(rates),1);
        
        rate_vs_phys_d05 = inf(length(rates),3);
        rate_vs_phys_d05_hi = inf(length(rates),3);
        rate_vs_phys_d05_lo = inf(length(rates),3);
        
        for l=1:length(rates)
            rate_vs_d05(l) = interp1(rpb,doses,rates(l));
            rate_vs_d05_hi(l) = interp1(rpb-rplo,doses,rates(l));
            rate_vs_d05_lo(l) = interp1(rpb+rphi,doses,rates(l));
            
        if ~isnan(rate_vs_d05(l)),
            rate_vs_phys_d05(l,1) = max(roots([1 cur_a2b*3 -cur_a2b*3*rate_vs_d05(l)]));
            rate_vs_phys_d05(l,2) =  max(roots([1 cur_a2b*4 -cur_a2b*4*rate_vs_d05(l)]));
            rate_vs_phys_d05(l,3) =  max(roots([1 cur_a2b*5 -cur_a2b*5*rate_vs_d05(l)]));
            end
            
            if ~isnan(rate_vs_d05_hi(l)),
            rate_vs_phys_d05_hi(l,1) = max(roots([1 cur_a2b*3 -cur_a2b*3*rate_vs_d05_hi(l)]));
            rate_vs_phys_d05_hi(l,2) =  max(roots([1 cur_a2b*4 -cur_a2b*4*rate_vs_d05_hi(l)]));
            rate_vs_phys_d05_hi(l,3) =  max(roots([1 cur_a2b*5 -cur_a2b*5*rate_vs_d05_hi(l)]));
            end
            
            if ~isnan(rate_vs_d05_lo(l)),
            rate_vs_phys_d05_lo(l,1) = max(roots([1 cur_a2b*3 -cur_a2b*3*rate_vs_d05_lo(l)]));
            rate_vs_phys_d05_lo(l,2) =  max(roots([1 cur_a2b*4 -cur_a2b*4*rate_vs_d05_lo(l)]));
            rate_vs_phys_d05_lo(l,3) =  max(roots([1 cur_a2b*5 -cur_a2b*5*rate_vs_d05_lo(l)]));
            end
            
        end
        
         cur_fig=figure(cur_fig_ctr+550);
        set(gcf,'Position',ss_four2three);
        hold on;
        h_d05(1)=plot(rates,rate_vs_d05,'LineWidth',2);
        h_d05(2)=plot(rates,rate_vs_d05_hi,'--');
        plot(rates,rate_vs_d05_lo,'--');
        grid on;
        set(gca,'GridLineStyle','--')
        set(gca,'xminortick','on','yminortick','on');
        set(gca,'box','on');
        set(gca,'FontSize',18);
        ylabel(['D_{5.0\rm{cc}} Gy_{',a2b{1},'}'],'FontSize',20);
        xlabel('Complication Rate (%)','FontSize',20);
        [~,ten_pct_ind] = min(abs(rates-0.1));
        [~,fftn_pct_ind] = min(abs(rates-0.15));
        [~,twnty_pct_ind] = min(abs(rates-0.2));
        
        text(0.05,0.9,['10% Comp. rate -> ',...
            num2str(rate_vs_d05(ten_pct_ind),'%3.1f'),' Gy_{',num2str(a2b{1}),'} [',...
            num2str(rate_vs_d05_lo(ten_pct_ind),'%3.1f'),' - ',...
            num2str(rate_vs_d05_hi(ten_pct_ind),'%3.1f'),']',10,...
            '15% Comp. rate -> ',...
            num2str(rate_vs_d05(fftn_pct_ind),'%3.1f'),' Gy_{',num2str(a2b{1}),'} [',...
            num2str(rate_vs_d05_lo(fftn_pct_ind),'%3.1f'),' - ',...
            num2str(rate_vs_d05_hi(fftn_pct_ind),'%3.1f'),']',10,...
            '20% Comp. rate -> ',...
            num2str(rate_vs_d05(twnty_pct_ind),'%3.1f'),' Gy_{',num2str(a2b{1}),'} [',...
            num2str(rate_vs_d05_lo(twnty_pct_ind),'%3.1f'),' - ',...
            num2str(rate_vs_d05_hi(twnty_pct_ind),'%3.1f'),']'],...
            'units','normalized','FontSize',18,'BackgroundColor','w');
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,...
                [fig_basename,'_rate_vs_bed_d05'],'-pdf');
            disp(['Saving ',fig_basename,'_rate_vs_bed_d05.pdf']);
        end
        
        %% rate vs d05 physical doses
         cur_fig=figure(cur_fig_ctr+5500);
        set(gcf,'Position',ss_four2three);
        hold on;
        h_d05_nfx3=plot(rates,rate_vs_phys_d05(:,1),'LineWidth',2);
        plot(rates,rate_vs_phys_d05_hi(:,1),'--');
        plot(rates,rate_vs_phys_d05_lo(:,1),'--');
        h_d05_nfx4=plot(rates,rate_vs_phys_d05(:,2),'g','LineWidth',2);
        plot(rates,rate_vs_phys_d05_hi(:,2),'g--');
        plot(rates,rate_vs_phys_d05_lo(:,2),'g--');
        h_d05_nfx5=plot(rates,rate_vs_phys_d05(:,3),'r','LineWidth',2);
        plot(rates,rate_vs_phys_d05_hi(:,3),'r--');
        plot(rates,rate_vs_phys_d05_lo(:,3),'r--');
        grid on;
        legend([h_d05_nfx3 h_d05_nfx4 h_d05_nfx5],...
            {'3 Fx','4 Fx','5 Fx'},'FontSize',18);
            
        set(gca,'GridLineStyle','--')
        set(gca,'xminortick','on','yminortick','on');
        set(gca,'box','on');
        set(gca,'FontSize',18);
        ylabel(['D_{5.0\rm{cc}} Gy_{PHYS}'],'FontSize',20);
        xlabel('Complication Rate (%)','FontSize',20);

        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,...
                [fig_basename,'_rate_vs_phys_d05'],'-pdf');
            disp(['Saving ',fig_basename,'_rate_vs_phys_d05.pdf']);
        end
        
        %% D_35
        
        cur_fig=figure(cur_fig_ctr+600);
        set(gcf,'Position',ss_four2three);
       
        
        % regression using exact EUD
        [b,~,s]=glmfit(d35s,[ptcomp pttotal],'binomial','link','logit');
        
        pvalue = [s.p];
        pval = pvalue(2); % the p-value corresponding to gEUD
        
        %tmp
        doses = (0:max(d35s))'; % doses (gEUD) whose RP probability will be computed
        [rpb,rplo,rphi] = glmval(b, doses,'logit',s); % the responding function values at doses
        
        % plot
        hold on;
        plot(doses,rpb,'k','LineWidth',2); % responding function
        plot(doses,rpb-rplo,'k','LineWidth',1); % low CI curve
        plot(doses,rpb+rphi,'k','LineWidth',1); % high CI curve
        
        
        
        
        
        
        
        d35_lim = interp1(rpb,doses,comp_rate);
        d35_limhi = interp1(rpb-rplo,doses,comp_rate);
        d35_limlo = interp1(rpb+rphi,doses,comp_rate);
        
     if (~isnan(d35_lim) && ~isnan(d35_limhi) && ~isnan(d35_limlo)) &&...
                (~isinf(d35_lim) && ~isinf(d35_limhi) && ~isinf(d35_limlo)),
            
        d35_phys = [max(roots([1 cur_a2b*3 -cur_a2b*3*d35_lim])),...
                    max(roots([1 cur_a2b*4 -cur_a2b*4*d35_lim])),...
                    max(roots([1 cur_a2b*5 -cur_a2b*5*d35_lim]))];
        d35_phys_hi = [max(roots([1 cur_a2b*3 -cur_a2b*3*d35_limhi])),...
                        max(roots([1 cur_a2b*4 -cur_a2b*4*d35_limhi])),...
                        max(roots([1 cur_a2b*5 -cur_a2b*5*d35_limhi]))];
        d35_phys_lo = [max(roots([1 cur_a2b*3 -cur_a2b*3*d35_limlo])),...
                        max(roots([1 cur_a2b*4 -cur_a2b*4*d35_limlo])),...
                        max(roots([1 cur_a2b*5 -cur_a2b*5*d35_limlo]))];
                
           else
            d35_phys =[0 0 0];
            d35_phys_hi =[0 0 0];
            d35_phys_lo =[0 0 0];
        end
        
                
        
        disp([]);
        disp('******');
        disp(['D_{3.5 cc} corresponding to ',num2str(comp_rate*100),'% rate: ',10,...
            num2str(d35_lim),...
            ' Gy_{',num2str(a2b{1}),'} [',num2str(d35_limlo),', ',num2str(d35_limhi),']']);
        disp(['= Physical doses =']);
        disp(['Nfx = 3',10,...
            num2str(d35_phys(1)),...
            ' Gy_{',num2str(a2b{1}),'} [',num2str(d35_phys_lo(1)),', ',num2str(d35_phys_hi(1)),']']);
        disp(['Nfx = 4',10,...
            num2str(d35_phys(2)),...
            ' Gy_{',num2str(a2b{1}),'} [',num2str(d35_phys_lo(2)),', ',num2str(d35_phys_hi(2)),']']);
        disp(['Nfx = 5',10,...
            num2str(d35_phys(3)),...
            ' Gy_{',num2str(a2b{1}),'} [',num2str(d35_phys_lo(3)),', ',num2str(d35_phys_hi(3)),']']);
                
            
        disp('******');
        disp([]);
        flg=[CGobj.mGrp.mFlgCensor]; % censor flags of patients
        
        [medianeud,~,~,binlow,binhigh,numcomp,numtotal,betainv84,betainv16] = EventObserved(flg,d35s,4);
        prob = numcomp./numtotal;
        % plot
        errorbar(medianeud,prob,max(0,prob-betainv16),max(0,betainv84-prob),'k*','LineWidth',1);
        errorbar_x(medianeud,prob,(medianeud-binlow),(binhigh-medianeud),'k*');
        
        %tmp
        ylim([0 0.5]);
        xlim([0 max(binhigh)*1.1]);
        
        %xlims = xlim;
        %xmax = xlims(2);
        %xlim([0 xmax]);
        %ylim([0 0.5]);
        
        loga_str = ['p-val: %s'];
        %str = sprintf(loga_str,num2str(loga),num2str(pval,2),structures{j});
        str = sprintf(loga_str,num2str(pval,'%3.2e'));
        %text(0.65,0.85,str,'FontSize',24,'Units','normalized');
        textbp(str,'FontSize',24);
        
        set(gca,'xminortick','on','yminortick','on');
        set(gca,'box','on');
        set(gca,'FontSize',16);
        %xlabel('D_{05 cc} [Gy]','FontSize',25);
        xlabel('D_{3.5 cc} [Gy_{10}]','FontSize',25);
        ylabel('Complication rate observed','FontSize',22);
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,...
                [fig_basename,'_lr_d35'],'-pdf');
            disp(['Saving ',fig_basename,'_lr_d35.pdf']);
        end
        
        
        %% rate vs d35
         rate_vs_d35 = inf(length(rates),1);
        rate_vs_d35_hi = inf(length(rates),1);
        rate_vs_d35_lo = inf(length(rates),1);
        
        rate_vs_phys_d35 = inf(length(rates),3);
        rate_vs_phys_d35_hi = inf(length(rates),3);
        rate_vs_phys_d35_lo = inf(length(rates),3);
        
        for l=1:length(rates)
            rate_vs_d35(l) = interp1(rpb,doses,rates(l));
            rate_vs_d35_hi(l) = interp1(rpb-rplo,doses,rates(l));
            rate_vs_d35_lo(l) = interp1(rpb+rphi,doses,rates(l));
            
           if ~isnan(rate_vs_d35(l)),
            rate_vs_phys_d35(l,1) = max(roots([1 cur_a2b*3 -cur_a2b*3*rate_vs_d35(l)]));
            rate_vs_phys_d35(l,2) =  max(roots([1 cur_a2b*4 -cur_a2b*4*rate_vs_d35(l)]));
            rate_vs_phys_d35(l,3) =  max(roots([1 cur_a2b*5 -cur_a2b*5*rate_vs_d35(l)]));
            end
            
            if ~isnan(rate_vs_d35_hi(l)),
            rate_vs_phys_d35_hi(l,1) = max(roots([1 cur_a2b*3 -cur_a2b*3*rate_vs_d35_hi(l)]));
            rate_vs_phys_d35_hi(l,2) =  max(roots([1 cur_a2b*4 -cur_a2b*4*rate_vs_d35_hi(l)]));
            rate_vs_phys_d35_hi(l,3) =  max(roots([1 cur_a2b*5 -cur_a2b*5*rate_vs_d35_hi(l)]));
            end
            
            if ~isnan(rate_vs_d35_lo(l)),
            rate_vs_phys_d35_lo(l,1) = max(roots([1 cur_a2b*3 -cur_a2b*3*rate_vs_d35_lo(l)]));
            rate_vs_phys_d35_lo(l,2) =  max(roots([1 cur_a2b*4 -cur_a2b*4*rate_vs_d35_lo(l)]));
            rate_vs_phys_d35_lo(l,3) =  max(roots([1 cur_a2b*5 -cur_a2b*5*rate_vs_d35_lo(l)]));
            end
                        
        end
        
        cur_fig=figure(cur_fig_ctr+650);
        set(gcf,'Position',ss_four2three);
        hold on;
        h_d35(1)=plot(rates,rate_vs_d35,'LineWidth',2);
        h_d35(2)=plot(rates,rate_vs_d35_hi,'--');
        plot(rates,rate_vs_d35_lo,'--');
        grid on;
        set(gca,'GridLineStyle','--')
        set(gca,'xminortick','on','yminortick','on');
        set(gca,'box','on');
        set(gca,'FontSize',18);
        ylabel(['D_{3.5\rm{cc}} Gy_{',a2b{1},'}'],'FontSize',20);
        xlabel('Complication Rate (%)','FontSize',20);
        [~,ten_pct_ind] = min(abs(rates-0.1));
        [~,fftn_pct_ind] = min(abs(rates-0.15));
        [~,twnty_pct_ind] = min(abs(rates-0.2));
        
        text(0.05,0.9,['10% Comp. rate -> ',...
            num2str(rate_vs_d35(ten_pct_ind),'%3.1f'),' Gy_{',num2str(a2b{1}),'} [',...
            num2str(rate_vs_d35_lo(ten_pct_ind),'%3.1f'),' - ',...
            num2str(rate_vs_d35_hi(ten_pct_ind),'%3.1f'),']',10,...
            '15% Comp. rate -> ',...
            num2str(rate_vs_d35(fftn_pct_ind),'%3.1f'),' Gy_{',num2str(a2b{1}),'} [',...
            num2str(rate_vs_d35_lo(fftn_pct_ind),'%3.1f'),' - ',...
            num2str(rate_vs_d35_hi(fftn_pct_ind),'%3.1f'),']',10,...
            '20% Comp. rate -> ',...
            num2str(rate_vs_d35(twnty_pct_ind),'%3.1f'),' Gy_{',num2str(a2b{1}),'} [',...
            num2str(rate_vs_d35_lo(twnty_pct_ind),'%3.1f'),' - ',...
            num2str(rate_vs_d35_hi(twnty_pct_ind),'%3.1f'),']'],...
            'units','normalized','FontSize',18,'BackgroundColor','w');
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,...
                [fig_basename,'_rate_vs_bed_d35'],'-pdf');
            disp(['Saving ',fig_basename,'_rate_vs_bed_d35.pdf']);
        end
        
        %% rate vs d35 physical doses
         cur_fig=figure(cur_fig_ctr+6500);
        set(gcf,'Position',ss_four2three);
        hold on;
        h_d35_nfx3=plot(rates,rate_vs_phys_d35(:,1),'LineWidth',2);
        plot(rates,rate_vs_phys_d35_hi(:,1),'--');
        plot(rates,rate_vs_phys_d35_lo(:,1),'--');
        h_d35_nfx4=plot(rates,rate_vs_phys_d35(:,2),'g','LineWidth',2);
        plot(rates,rate_vs_phys_d35_hi(:,2),'g--');
        plot(rates,rate_vs_phys_d35_lo(:,2),'g--');
        h_d35_nfx5=plot(rates,rate_vs_phys_d35(:,3),'r','LineWidth',2);
        plot(rates,rate_vs_phys_d35_hi(:,3),'r--');
        plot(rates,rate_vs_phys_d35_lo(:,3),'r--');
        grid on;
        legend([h_d35_nfx3 h_d35_nfx4 h_d35_nfx5],...
            {'3 Fx','4 Fx','5 Fx'},'FontSize',18);

        ylim([0 30]);
        set(gca,'GridLineStyle','--')
        set(gca,'xminortick','on','yminortick','on');
        set(gca,'box','on');
        set(gca,'FontSize',18);
        ylabel(['D_{3.5\rm{cc}} Gy_{PHYS}'],'FontSize',20);
        xlabel('Complication Rate (%)','FontSize',20);

        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,...
                [fig_basename,'_rate_vs_phys_d35'],'-pdf');
            disp(['Saving ',fig_basename,'_rate_vs_phys_d35.pdf']);
        end
        
        
    end
    
    %% Print and save
    % Logistic regression pvalues
    cur_fig=figure(i);
    
    ylim([min([0.04 min(min_lr_pvals)*0.75]) 1]);
    set(gca,'YScale','log');
    lgnd_pvals=legend(h_pvals,structures,'Location','Best');
    set(lgnd_pvals,'FontSize',16);
    
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,...
            [fig_tox_basename,'_lr_eud_pvals'],'-pdf');
        disp(['Saving ',fig_tox_basename,'_lr_eud_pvals.pdf']);
        
    end
    
    % Logistic regression loglikelihoods
    cur_fig=figure(i+100);
    if isequal(toxicities{i},'pultox')
        ylim([-0.59 -0.574]);
    end
    
    lgnd_llhds=legend(h_llhds,structures,'Location','Best');
    set(lgnd_llhds,'FontSize',18);
    
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,...
            [fig_tox_basename,'_lr_eud_llhds'],'-pdf');
        disp(['Saving ',fig_tox_basename,'_lr_eud_llhds.pdf']);
    end
    
     % t-values
    cur_fig=figure(i+1000);hold on;
    lgnd_tvals=legend(h_tvals,structures,'Location','Best');
    set(lgnd_tvals,'FontSize',18);
    grid off;
    plot([-1 1],[0 0],'--k','LineWidth',0.4);
    if isequal(toxicities{i},'pultox')
        ylim([-2.2 2.2]);
    end
    
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,...
            [fig_tox_basename,'_lr_eud_tvals'],'-pdf');
        disp(['Saving ',fig_tox_basename,'_lr_eud_tvals.pdf']);
    end
    
end
toc;