function NFz_LowBED
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
structures = {'LUNGS'};
%structures = {'PTV' 'GTV'};
%toxicities = {'rp','pultox','esotox','lclfail};
%toxicities = {'lclfail'};
toxicities = {'pultox'};

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
         
        ptcomp = ones(CGobj.mNumInGrp,1); 
        ptcomp([CGobj.mGrp.mFlgCensor])=0;
         
        rxs = [CGobj.mGrp.mDoseTx]./100;
        fxs = [CGobj.mGrp.mFxNum];
        rx_beds = rxs.*(1+rxs./(10*fxs));rx_beds = rxs.*(1+rxs./(10*fxs));
        % prescription BED histogram

        cur_fig=figure(cur_fig_ctr);
        set(gcf,'Position',ss_four2three);
        x_bins = 0:15:200;
        
        h2=hist(rx_beds(fxs==2),x_bins);
        h3=hist(rx_beds(fxs==3),x_bins);
        h4=hist(rx_beds(fxs==4),x_bins);
        h5=hist(rx_beds(fxs==5),x_bins);
        b_plot=bar(x_bins,[h2;h3;h4;h5]','stacked');
        
        lgnd=legend(b_plot,'Nfx = 2','Nfx = 3','Nfx = 4','Nfx = 5');
        set(lgnd,'FontSize',21);
        xlim([0 200]);
        xlabel('Prescription BED ($\alpha/\beta = 10$Gy)','FontSize',22,'interpreter','latex');
        ylabel('Number of patients','FontSize',22,'interpreter','latex');
        set(gca,'FontSize',20);

        cur_fig_ctr = cur_fig_ctr+1;
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,...
                [fig_basename,'_bed_hist'],'-pdf');
            disp(['Saving ',fig_basename,'_bed_hist.pdf']);
        end
        
        
        %plot ranked prescription BED with complications
        [sorted_beds,sorted_inds] = sort(rx_beds);
        sorted_ptcomp = ptcomp(sorted_inds);

        cur_fig=figure(cur_fig_ctr);
        set(gcf,'Position',ss_four2three);
        
        x_range = 1:length(sorted_beds);
         h_cens=plot(x_range(~logical(sorted_ptcomp)),sorted_beds(~logical(sorted_ptcomp)),'b+','LineWidth',2,'MarkerSize',12);%'MarkerSize',8,'MarkerFaceColor','b');
        hold on;
        h_comp=plot(x_range(logical(sorted_ptcomp)),sorted_beds(logical(sorted_ptcomp)),'r+','LineWidth',2,'MarkerSize',12);%,'MarkerFaceColor','r');
        set(gca,'FontSize',20);
        xlabel('BED Rank (\#)','FontSize',22,'interpreter','latex');
        ylabel('BED ($\alpha/\beta = 10$Gy)','FontSize',22,'interpreter','latex');
        ylim([30 190]);
        xlim([0 126]);
        lgnd = legend([h_comp h_cens],...
            [toxicities{i},'~$\geq 2$, Med. BED: ',num2str(median(rx_beds(logical(ptcomp))),3)],...
            [toxicities{i},'~$< 2$, Med. BED: ',num2str(median(rx_beds(~logical(ptcomp))),3)],'Location','NorthWest')
        set(lgnd,'FontSize',22);
        set(lgnd,'interpreter','latex');
       
                cur_fig_ctr = cur_fig_ctr+1;
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,...
                [fig_basename,'_rank_bed'],'-pdf');
            disp(['Saving ',fig_basename,'_rank_bed.pdf']);
        end
    cur_fig_ctr = cur_fig_ctr+1;        
    end
end
