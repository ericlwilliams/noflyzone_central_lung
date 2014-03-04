function NFz_CoxResults
tic; close all;
warning('off'); % ignore negative data warning
screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check flags!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do_print = true;
do_tcp_exclude = false; %exclude 10 patients with low Rx (<= 3000 cGy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
disp('**** START Flags ****');
disp([num2str(do_print),' - do_print']);
disp([num2str(do_tcp_exclude),' - do_tcp_exclude']);
disp('**** END Flags ****');
disp(sprintf('\n'));

fig_loc = 'Z:/elw/MATLAB/nfz_analy/slides/figures/latest/';

%structures = {'ILUNG' 'ESOPHAGUS' 'HEART' 'LUNGS' 'NFZ' 'PBT'};
structures = {'ESOPHAGUS'};
%structures = {'PTV' 'GTV'};
%toxicities = {'rp','pultox','esotox','lclfail};
%toxicities = {'lclfail'};
toxicities = {'esotox'};

fp = 'Z:\elw\MATLAB\nfz_analy\meta_data\';

%a2b = {'Inf' '3' '10'};
a2b = {'Inf'};

for i=1:length(toxicities)
    
    for j=1:length(structures)
        cur_fig_ctr = (10*i)+j-1;
        
        fprintf('\n');
        disp(['Tox: ',toxicities{i}]);
        disp(['Struct: ',structures{j}]);    
        disp(['Counter: ',num2str(cur_fig_ctr)]);
        fprintf('\n');
        
        if do_tcp_exclude
         fig_basename = [fig_loc,'nfz_',...
                        structures{j},'_',...
                        toxicities{i},'_a2b',...
                        a2b{1},'_lowrx'];
        else
            fig_basename = [fig_loc,'nfz_',...
                        structures{j},'_',...
                        toxicities{i},'_a2b',...
                        a2b{1}];
        end
     %% load data
    if do_tcp_exclude,
        fn = ['NFZ_',structures{j},'_',toxicities{i},'_a2b',a2b{1},'_lowrx_data.mat'];
    else
        fn = ['NFZ_',structures{j},'_',toxicities{i},'_a2b',a2b{1},'_data.mat'];
    end
        disp(['']);
        disp(['Loading ',fn]);
        disp(['']);
        load(strcat(fp,fn),'CGobj_org');
        CGobj = CGobj_org;
        LymanN = log10(CGobj.mLymanN);
        CGobj.mLymanN = LymanN;
        
        %get compdate and censoring info for cox analy
        f2 = ~cellfun('isempty',{CGobj.mGrp.mDateComp}); % patients with no complication date
        f3 = ~cellfun('isempty',{CGobj.mGrp.mDateLastFollowup}); % patients with no last follow up date
        compdate = inf(CGobj.mNumInGrp,1);
        lastfollowup = inf(CGobj.mNumInGrp,1);
        compdate(f2) = ([CGobj.mGrp(f2).mDateComp] - [CGobj.mGrp(f2).mDateBaseline])' / 30;
        lastfollowup(f3) = ([CGobj.mGrp(f3).mDateLastFollowup] - [CGobj.mGrp(f3).mDateBaseline])' / 30;
        compdate = min( lastfollowup, compdate );
        flgcensor = [CGobj.mGrp.mFlgCensor]';
        
        %% DVx Cox PH Model results
        [DVxCox,flgCox,flganti] = CGobj.fCoxParameter_DVH('DVx'); % find availabe Cox models
        flgCox(flganti)=false; % anti-correlations were not be considered
        
        infFlg = isinf([DVxCox.beta])';

        dx_no_corr = (sum(flgCox)<2); % no correlations found, shouldn't happen?
    
        % remove infinites
        flgCox = flgCox(~infFlg);
        DVxCox = DVxCox(~infFlg);
        x_dvx=CGobj.mBinsVol(~infFlg);
        
       
        
        
        %% DVx llhds
        logl = [DVxCox.logl]'; %logl(~flgCox) = -inf; % log likelihood of Cox model, anti-correlation points not counted
       
        
        cur_fig=figure(cur_fig_ctr); clf reset;
        set(gcf,'Position',ss_four2three);
        
        if ~dx_no_corr,
           
           
            corr_col = 'b';
           if isequal(toxicities{i},'lclfail'), % find max for *anti-corr*
               [mx,doseloc]=max(logl(~flgCox)); % the best fitting of Cox model
                x_pos_dvx = x_dvx(~flgCox);
                corr_col = 'r';
           else
               [mx,doseloc]=max(logl(flgCox)); % the best fitting of Cox model
                x_pos_dvx = x_dvx(flgCox);
           end
           lowCI68 = mx - 0.5; % 68% confidence
           %lowCI95 = mx - 1.96; % 95% confidence
             best_cox_dv = DVxCox(doseloc).data_exposure;
        
            loglog(x_dvx,logl','k-','LineWidth',1);hold on;
            h_corr=loglog(x_dvx(flgCox), logl(flgCox)','.','MarkerSize',25);
            h_acorr=loglog(x_dvx(~flgCox), logl(~flgCox)','r.','MarkerSize',25);
            h_low_68cl=loglog([x_dvx(2) max(x_dvx)],[lowCI68 lowCI68],'g--','LineWidth',2);
            h_mx_logl=loglog([x_pos_dvx(doseloc) x_pos_dvx(doseloc)],ylim,strcat(corr_col,'--'),'LineWidth',2);
            hold off;
            xlim([0 max(x_dvx)]);
            set(gca,'fontsize',18);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('(D_{V}) Volume [cc]','fontsize',20);
            ylabel('Cox model log-likelihood','fontsize',20);
            logl_str = {'Positive Corr.',...
                        'Negative Corr.',...
                        ['Max LL = ', num2str(mx,3),10,...
                'at D_{',num2str(x_pos_dvx(doseloc),3),'cc}'],...
                '68% CL'};
            legend([h_corr h_acorr h_mx_logl h_low_68cl],logl_str,...
                'Location','Best');
        else  %still printing empty figure for slides
            
            plot(x_dvx,logl','k-','LineWidth',1);hold on;
            %h_corr=plot(x_dvx(flgCox), logl(flgCox)','.','MarkerSize',25);
            h_acorr=plot(x_dvx(~flgCox), logl(~flgCox)','r.','MarkerSize',25);
            hold off;
            xlim([0 max(x_dvx)]);
            set(gca,'fontsize',18);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('(D_{V}) Volume [cc]','fontsize',20);
            ylabel('Cox model log-likelihood','fontsize',20);
            logl_str = {'Negative Corr.'};
            legend([h_acorr],logl_str,...
                'Location','Best');
            
        end
        
        
        if do_print,
            set(cur_fig,'Color','w');
         export_fig(cur_fig,...
            [fig_basename,'_cox_dv_llhds'],'-pdf');
        disp(['Saving ',fig_basename,'_cox_dv_llhds.pdf']);
        end
        
        
        %% DVx Correlatinos
        if ~dx_no_corr,
            dvx_corrs = corr([DVxCox.data_exposure],[DVxCox.data_exposure]).^2;
            figure(cur_fig_ctr+100);
            imagesc(x_dvx,x_dvx,dvx_corrs);
            set(gca,'YDir','normal');
            colorbar;
            title('D_{V} R^2 Correlations','FontSize',14);
            xlabel('(D_{V}) Volume [cc]','FontSize',14)
            ylabel('(D_{V}) Volume [cc]','FontSize',14)
        end
        
        %% DVx p-vals
        p = ones(size(DVxCox));
        p = [DVxCox.p]';
        
        
       if isequal(toxicities{i},'lclfail'), % find max for *anti-corr*
           [min_p,ploc] = min(p(~flgCox));
       else
            [min_p,ploc] = min(p(flgCox));
       end


        
        cur_fig=figure(cur_fig_ctr+200); clf reset;
        set(gcf,'Position',ss_four2three);
        
        if ~dx_no_corr, %still print empty figure for slides
            loglog(x_dvx,p','k-','LineWidth',1);hold on;
            h_pos_pval=loglog(x_dvx(flgCox),p(flgCox)','.','MarkerSize',25);hold on;
            h_neg_pval=loglog(x_dvx(~flgCox),p(~flgCox)','r.','MarkerSize',25);
            
            h_min_pval = loglog([x_pos_dvx(ploc) x_pos_dvx(ploc)],ylim,strcat(corr_col,'--'),'LineWidth',2);
            
            disp(['Cox model significant for D_V, with V <',num2str(interp1(p(1:20),x_dvx(1:20),0.05))]);
            xlim([0 max(x_dvx)]);
            ylim([min(min(p),0.05)-0.001 1]);
            h_sig=loglog([min(x_dvx(x_dvx>0)) max(x_dvx)],[0.05 0.05],'g--','LineWidth',2);hold off;
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            pval_str = {'Positive Corr.',...
                        'Negative Corr.',...
                        ['Min $p$-val = ', num2str(min_p,2),10,...
                'at $D_{',num2str(x_pos_dvx(ploc),3),' cc}$'],...
                '$p = 0.05$'};

           cur_lgnd = legend([h_pos_pval h_neg_pval h_min_pval h_sig],pval_str,...
                'Location','SouthEast');
            set(cur_lgnd,'FontSize',21);
            set(cur_lgnd,'interpreter','latex');
            xlabel('(D_{V}) Volume [cc]','fontsize',24);
            ylabel('Cox model p-value','fontsize',24);
            set(gca,'FontSize',22);
        else  %still printing empty figure for slides
            semilogy(x_dvx,p','k-','LineWidth',1);hold on;
            h_neg_pval=semilogy(x_dvx(~flgCox),p(~flgCox)','r.','MarkerSize',25);
            h_sig=semilogy([0 max(x_dvx)],[0.05 0.05],'g--','LineWidth',2);
            hold off;
            xlim([0 max(x_dvx)]);
            ylim([min(min(p),0.05)-0.001 1]);
            set(gca,'fontsize',18);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            pval_str = {'Negative Corr.' 'p = 0.05'};
           legend([h_neg_pval h_sig],pval_str,...
                    'Location','Best');
            %set(cur_lgnd,'FontSize',20);
            xlabel('(D_{V}) Volume [cc]','fontsize',20);
            ylabel('Cox model p-value','fontsize',20);
        end
        
        if do_print,
            set(cur_fig,'Color','w');
         export_fig(cur_fig,...
            [fig_basename,'_cox_dv_pvals'],'-pdf');
        disp(['Saving ',fig_basename,'_cox_dv_pvals.pdf']);
        
        end
        
        
        
        
        
          %% VDx Cox PH Model results
        [VDxCox,flgCox,flganti] = CGobj.fCoxParameter_DVH('VDx'); % find availabe Cox models
        
        %[VDxCox,flgCox,flganti] = CGobj.fCoxParameter_DVH('Vx'); % Vx -> zeros excluded!
        
        flgCox(flganti)=false; % anti-correlations were not be considered
        
        infFlg = isinf([VDxCox.beta])';

        %no_corr = (sum(~infFlg)<2); % no correlations found, shouldn't happen?
        vx_no_corr = (sum(flgCox)<2); % no correlations found, shouldn't happen?
        
        % remove infinites
        flgCox = flgCox(~infFlg);
        VDxCox = VDxCox(~infFlg);
        x_vdx=CGobj.mBinsDose(~infFlg);
        

        
        %% VDx llhds
        
        
        %loglikelihood = [VDxCox.logl]'; %logl(~flgCox) = -inf; % log likelihood of Cox model, anti-correlation points not counted
     
        %df = [cellfun(@(x) length(x), {VDxCox(:).data_exposure})]';
         
        %logl = loglikelihood./df;
        logl = [VDxCox.logl]'; %logl(~flgCox) = -inf; % log likelihood of Cox model, anti-correlation points not counted   
            
            
            
        corr_col = 'b';
        if isequal(toxicities{i},'lclfail'), % find max for *anti-corr*
           [mx,doseloc]=max(logl(~flgCox)); % the best fitting of Cox model
           x_pos_vdx = x_vdx(~flgCox);
           corr_col = 'r';
        else
            [mx,doseloc]=max(logl(flgCox)); % the best fitting of Cox model
            x_pos_vdx = x_vdx(flgCox);
        end
        
        best_cox_vd = VDxCox(doseloc).data_exposure;
        
        lowCI68 = mx - 0.5; % 68% confidence
        lowCI95 = mx - 1.96; % 95% confidence
        
        cur_fig=figure(cur_fig_ctr+300); clf reset;
        set(gcf,'Position',ss_four2three);
        
        if ~vx_no_corr,
            plot(x_vdx,logl','k-','LineWidth',1);hold on;
            h_corr=plot(x_vdx(flgCox), logl(flgCox)','.','MarkerSize',25);
            h_acorr=plot(x_vdx(~flgCox), logl(~flgCox)','r.','MarkerSize',25);
            h_low_68cl=plot([0 max(x_vdx)],[lowCI68 lowCI68],'g--','LineWidth',2);
            h_mx_logl=plot([x_pos_vdx(doseloc) x_pos_vdx(doseloc)],ylim,strcat(corr_col,'--'),'LineWidth',2);
            hold off;
            xlim([0 max(x_vdx)]);
            set(gca,'fontsize',18);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('(V_{D}) Dose [Gy]','fontsize',20);
            ylabel('Cox model log-likelihood / degree of freedom','fontsize',20);
            logl_str = {'Positive Corr.',...
                        'Negative Corr.',...
                        ['Max LL = ', num2str(mx,3),10,...
                'at V_{',num2str(x_pos_vdx(doseloc),3),'Gy.}'],...
                '68% CL'};
            legend([h_corr h_acorr h_mx_logl h_low_68cl],logl_str,...
                'Location','Best');
            
            

            
        else  %still printing empty figure for slides
           plot(x_vdx,logl','k-','LineWidth',1);hold on;
           h_acorr=plot(x_vdx(~flgCox), logl(~flgCox)','r.','MarkerSize',25);
            hold off;
            xlim([0 max(x_vdx)]);
            set(gca,'fontsize',18);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('(V_{D}) Dose [Gy]','fontsize',20);
            ylabel('Cox model log-likelihood','fontsize',20);
            logl_str = {'Negative Corr.'};
            legend([h_acorr],logl_str,...
                'Location','Best');
            
        end


 
              
            if do_print,
                set(cur_fig,'Color','w');
                
                export_fig(cur_fig,...
                    [fig_basename,'_cox_vd_llhds'],'-pdf');
                disp(['Saving ',fig_basename,'_cox_vd_llhds.pdf']);
                
            end
            
            
            
            
        
        %% VDx p-vals
        p = ones(size(VDxCox));
        p = [VDxCox.p]';
        
        
        if isequal(toxicities{i},'lclfail'), % find max for *anti-corr*
            [min_p,ploc] = min(p(~flgCox));
        else
            [min_p,ploc] = min(p(flgCox));
        end

        
        cur_fig=figure(cur_fig_ctr+500); clf reset;
        set(gcf,'Position',ss_four2three);
        
        if ~vx_no_corr, %still print empty figure for slides
            semilogy(x_vdx,p','k-','LineWidth',1);hold on;
            h_pos_pval=semilogy(x_vdx(flgCox),p(flgCox)','.','MarkerSize',25);hold on;
            h_neg_pval=semilogy(x_vdx(~flgCox),p(~flgCox)','r.','MarkerSize',25);
            h_sig=semilogy([0 max(x_vdx)],[0.05 0.05],'g--','LineWidth',2);
            h_min_pval = semilogy([x_pos_vdx(ploc) x_pos_vdx(ploc)],ylim,strcat(corr_col,'--'),'LineWidth',2);
            hold off;
            disp(['Cox model significant for V_D, for D >',num2str(interp1(p(flgCox),x_vdx(flgCox),0.05))]);
            xlim([0 max(x_vdx)]);
            ylim([min(min(p),0.05)-0.001 1]);
            set(gca,'fontsize',22);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
             pval_str = {'Positive Corr.',...
                        'Negative Corr.',...
                        ['Min $p$-val = ', num2str(min_p,2),10,...
                'at $V_{',num2str(x_pos_vdx(ploc),3),'~Gy_{10}}$'],...
                '$p = 0.05$'};


           cur_lgnd=legend([h_pos_pval h_neg_pval h_min_pval h_sig],pval_str,...
                'Location','SouthWest');
            
            set(cur_lgnd,'FontSize',21);
            set(cur_lgnd,'interpreter','latex');
%             legend_best_fit(gca);
            %set(cur_lgnd,'FontSize',20);
            xlabel('(V_{D}) Dose [Gy_{10}]','fontsize',24);
            ylabel('Cox model p-value','fontsize',24);
        else  %still printing empty figure for slides
            semilogy(x_vdx,p','k-','LineWidth',1);hold on;
            h_neg_pval=semilogy(x_vdx(~flgCox),p(~flgCox)','r.','MarkerSize',25);
            h_sig=semilogy([0 max(x_vdx)],[0.05 0.05],'g--','LineWidth',2);
            hold off;
            xlim([0 max(x_vdx)]);
            ylim([min(min(p),0.05)-0.001 1]);
            set(gca,'fontsize',22);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            pval_str = {'Negative Corr.',...
                        'p = 0.05'};
           legend([h_neg_pval],pval_str,...
                'Location','SouthWest');
            legend_best_fit(gca);
            %set(cur_lgnd,'FontSize',20);
            xlabel('(V_{D}) Dose [Gy]','fontsize',22);
            ylabel('Cox model p-value','fontsize',22);
        end
        
        if do_print,
            set(cur_fig,'Color','w');

         export_fig(cur_fig,...
            [fig_basename,'_cox_vd_pvals'],'-pdf');
        disp(['Saving ',fig_basename,'_cox_vd_pvals.pdf']);

        end
      
        %% VDx Correlatinos
        if ~vx_no_corr,
            vdx_corrs = corr([VDxCox.data_exposure],[VDxCox.data_exposure]).^2;
            figure(cur_fig_ctr+400);
            imagesc(x_vdx,x_vdx,vdx_corrs);
            set(gca,'YDir','normal');
            colorbar;
            title('V_{D} R^2 Correlations','FontSize',14);
            xlabel('(V_{D}) Dose [Gy]','FontSize',14)
            ylabel('(V_{D}) Dos [Gy]','FontSize',14)
        end
          
        %% Best DVx VDx Cox model
    if ~vx_no_corr && ~dx_no_corr
        [~,~,~,cur_stats]=...
                 coxphfit([best_cox_dv best_cox_vd],compdate,'baseline',0,'censoring',flgcensor);
        comb_dxvx_pvals = [cur_stats.p]
         disp(['Best D_V + V_D Cox model pvalues: ',num2str([cur_stats.p]')]);
    end
% %         
       
if isequal(a2b{i},'Inf'),


 %% DVx +  Fx Cox PH Model results
        [DVxCox,orgFlgCox,flganti] = CGobj.fCoxParameter_DVH('DVxFx'); % find availabe Cox models

        orgFlgCox(flganti)=false; % anti-correlations were not be considered
        flgCox = orgFlgCox;
        
        infFlg = isinf([DVxCox.beta])';

        infFlg = [isinf(infFlg(:,1))].*[isinf(infFlg(:,2))];
        
        dx_no_corr = (sum(flgCox)<2); % no correlations found, shouldn't happen?
    
        % remove infinites
         flgCox = flgCox(~infFlg);
         DVxCox = DVxCox(~infFlg);
         x_dvx=CGobj.mBinsVol(~infFlg);
        
       
        
        
        %% DVx llhds
        logl = [DVxCox.logl]'; %logl(~flgCox) = -inf; % log likelihood of Cox model, anti-correlation points not counted
       
        
        cur_fig=figure(cur_fig_ctr); clf reset;
        set(gcf,'Position',ss_four2three);
        
        if ~dx_no_corr,
           
           
            corr_col = 'b';
           if isequal(toxicities{i},'lclfail'), % find max for *anti-corr*
               [mx,doseloc]=max(logl(~flgCox)); % the best fitting of Cox model
                x_pos_dvx = x_dvx(~flgCox);
                corr_col = 'r';
           else
               [mx,doseloc]=max(logl(flgCox)); % the best fitting of Cox model
                x_pos_dvx = x_dvx(flgCox);
           end
           lowCI68 = mx - 0.5; % 68% confidence
           %lowCI95 = mx - 1.96; % 95% confidence
             best_cox_dv = DVxCox(doseloc).data_exposure;
        
            loglog(x_dvx,logl','k-','LineWidth',1);hold on;
            h_corr=loglog(x_dvx(flgCox), logl(flgCox)','.','MarkerSize',25);
            h_acorr=loglog(x_dvx(~flgCox), logl(~flgCox)','r.','MarkerSize',25);
            h_low_68cl=loglog([x_dvx(2) max(x_dvx)],[lowCI68 lowCI68],'g--','LineWidth',2);
            h_mx_logl=loglog([x_pos_dvx(doseloc) x_pos_dvx(doseloc)],ylim,strcat(corr_col,'--'),'LineWidth',2);
            hold off;
            xlim([0 max(x_dvx)]);
            set(gca,'fontsize',18);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('(D_{V}) Volume [cc]','fontsize',20);
            ylabel('Cox model log-likelihood','fontsize',20);
            logl_str = {'Positive Corr.',...
                        'Negative Corr.',...
                        ['Max LL = ', num2str(mx,3),10,...
                'at D_{',num2str(x_pos_dvx(doseloc),3),'cc}'],...
                '68% CL'};
            legend([h_corr h_acorr h_mx_logl h_low_68cl],logl_str,...
                'Location','Best');
        else  %still printing empty figure for slides
            
            plot(x_dvx,logl','k-','LineWidth',1);hold on;
            %h_corr=plot(x_dvx(flgCox), logl(flgCox)','.','MarkerSize',25);
            h_acorr=plot(x_dvx(~flgCox), logl(~flgCox)','r.','MarkerSize',25);
            hold off;
            xlim([0 max(x_dvx)]);
            set(gca,'fontsize',18);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            xlabel('(D_{V}) Volume [cc]','fontsize',20);
            ylabel('Cox model log-likelihood','fontsize',20);
            logl_str = {'Negative Corr.'};
            legend([h_acorr],logl_str,...
                'Location','Best');
            
        end
        
        
        if do_print,
            set(cur_fig,'Color','w');
         export_fig(cur_fig,...
            [fig_basename,'_cox_dvfx_llhds'],'-pdf');
        disp(['Saving ',fig_basename,'_cox_dvfx_llhds.pdf']);
        end
        
        
        %% DVx Correlatinos
        if ~dx_no_corr,
            dvx_corrs = corr([DVxCox.data_exposure],[DVxCox.data_exposure]).^2;
            figure(cur_fig_ctr+100);
            imagesc(x_dvx,x_dvx,dvx_corrs);
            set(gca,'YDir','normal');
            colorbar;
            title('D_{V} R^2 Correlations','FontSize',14);
            xlabel('(D_{V}) Volume [cc]','FontSize',14)
            ylabel('(D_{V}) Volume [cc]','FontSize',14)
        end
        
        %% DVx p-vals
        p = ones(size(DVxCox));
        p = [DVxCox.p]';
        
        
       if isequal(toxicities{i},'lclfail'), % find max for *anti-corr*
           [min_p,ploc] = min(p(~flgCox));
       else
            [min_p,ploc] = min(p(flgCox));
       end

    p_dv = p(:,1);
    p_fx = p(:,2);
        
        cur_fig=figure(cur_fig_ctr+200); clf reset;
        set(gcf,'Position',ss_four2three);
        
        if ~dx_no_corr, %still print empty figure for slides
            loglog(x_dvx,p','k-','LineWidth',1);hold on;
            h_pos_pval=loglog(x_dvx(flgCox),p_dv(flgCox)','.','MarkerSize',25);hold on;
            loglog(x_dvx(flgCox),p_fx(flgCox)','.','MarkerSize',25);hold on;
            h_neg_pval=loglog(x_dvx(~flgCox),p_dv(~flgCox)','r.','MarkerSize',25);
            loglog(x_dvx(~flgCox),p_fx(~flgCox)','r.','MarkerSize',25);
            
            h_min_pval = loglog([x_pos_dvx(ploc) x_pos_dvx(ploc)],ylim,strcat(corr_col,'--'),'LineWidth',2);
            
            disp(['Cox model significant for D_V, with V <',num2str(interp1(p(1:20),x_dvx(1:20),0.05))]);
            xlim([0 max(x_dvx)]);
            ylim([min(min(min(p)),0.05)-0.001 1]);
            h_sig=loglog([min(x_dvx(x_dvx>0)) max(x_dvx)],[0.05 0.05],'g--','LineWidth',2);hold off;
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            pval_str = {'Positive Corr.',...
                        'Negative Corr.',...
                        ['Min $p$-val = ', num2str(min_p,2),10,...
                'at $D_{',num2str(x_pos_dvx(ploc),3),' cc}$'],...
                '$p = 0.05$'};

           cur_lgnd = legend([h_pos_pval h_neg_pval h_min_pval h_sig],pval_str,...
                'Location','SouthEast');
            set(cur_lgnd,'FontSize',21);
            set(cur_lgnd,'interpreter','latex');
            xlabel('(D_{V}) Volume [cc]','fontsize',24);
            ylabel('Cox model p-value','fontsize',24);
            set(gca,'FontSize',22);
        else  %still printing empty figure for slides
            semilogy(x_dvx,p','k-','LineWidth',1);hold on;
            h_neg_pval=semilogy(x_dvx(~flgCox),p(~flgCox)','r.','MarkerSize',25);
            h_sig=semilogy([0 max(x_dvx)],[0.05 0.05],'g--','LineWidth',2);
            hold off;
            xlim([0 max(x_dvx)]);
            ylim([min(min(p),0.05)-0.001 1]);
            set(gca,'fontsize',18);
            set(gca,'xminortick','on','yminortick','on');
            set(gca,'box','on');
            pval_str = {'Negative Corr.' 'p = 0.05'};
           legend([h_neg_pval h_sig],pval_str,...
                    'Location','Best');
            %set(cur_lgnd,'FontSize',20);
            xlabel('(D_{V}) Volume [cc]','fontsize',20);
            ylabel('Cox model p-value','fontsize',20);
        end
        
        if do_print,
            set(cur_fig,'Color','w');
         export_fig(cur_fig,...
            [fig_basename,'_cox_dvfx_pvals'],'-pdf');
        disp(['Saving ',fig_basename,'_cox_dvfx_pvals.pdf']);
        
        end
end
        

    end
end
end
