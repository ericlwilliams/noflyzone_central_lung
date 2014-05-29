function NFz_BoxPlots

tic; %close all;

screen_size=get(0,'ScreenSize');
ss_two2two = [screen_size(3)/2 0 screen_size(4) screen_size(4)];
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

do_print = false;
%fig_loc = 'Z:/elw/MATLAB/nfz_analy/slides/figures/latest/';
fig_loc = 'V:/cwdvhs/elw/MATLAB/nfz_analy/slides/figures/latest/';
%fig_loc = 'pensmph6/MpcsResearch1/cwdvhs/elw/MATLAB/nfz_analy/slides/figures/latest/';

%structures = {'ILUNG' 'ESOPHAGUS' 'HEART' 'NFZ' 'PBT' 'LUNGS'};
structures = {'ESOPHAGUS'};
%toxicities = {'rp','pultox','esotox'};
toxicities = {'esotox'};

%fp = 'Z:\elw\MATLAB\nfz_analy\meta_data\';
%fp = 'pensmph6\MpcsResearch1\cwdvhs\elw\MATLAB\nfz_analy\meta_data\';
fp = 'V:\cwdvhs\elw\MATLAB\nfz_analy\meta_data\';
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
        
        grp = CGobj.mGrp;
        comps = ~f;
        eud = [grp(:).mEUD];
        
        %% Mean Boxplots
        
        eud_mean_loc=11;
        
        %msk
        eud_cens = eud(:,~comps);
        eud_comps = eud(:,comps);
        mean_cens = eud_cens(eud_mean_loc,:);
        mean_comps = eud_comps(eud_mean_loc,:);
        group = [ones(size(mean_cens'));...
            ones(size(mean_comps'))+1];
        means = cat(1,mean_cens',mean_comps');
        
        
        cur_fig=figure(cur_fig_ctr);clf reset;
        set(gcf,'Position',ss_four2three);
        
        boxplot(means,group);
        set(gca,'XTick',1:2,'XTickLabel',...
            {'W/out comp.','With comp.'},...
            'FontSize',20);
        ylabel('Mean Dose [Gy]','FontSize',20);
        %title([structures{j},': ',tox_titles{i}],'FontSize',18);
        
        if do_print,
            
         set(cur_fig,'Color','w');
         export_fig(cur_fig,...
            [fig_basename,'_box_euds_mld'],'-png');
        disp(['Saving ',fig_basename,'_box_euds_mld.pdf']);
        end
        
        
        d_bins = [CGobj.mGrp.mDoseBins_LQ];
        v_bins = [CGobj.mGrp.mVolCum];

        dmax = zeros(size(v_bins,2),1);
    
        for k=1:size(v_bins,2) % for each patient
            vol = v_bins(:,k);
            dose = d_bins(:,k);
            
            vol(~vol)=nan;
            nan_inds = find(isnan(vol));
            if ~isempty(nan_inds)
                min_ind = nan_inds(1)-1;
            else
                min_ind = length(vol);
            end
            dmax(k) = dose(min_ind);
        end

        
        fx_grp = [CGobj.mGrp.mFxNum];
        
        cur_fig=figure(cur_fig_ctr+10);clf reset;
        set(gcf,'Position',ss_four2three);
        boxplot(dmax,fx_grp);
        set(findobj(gca,'Type','text'),'FontSize',14);
        set(gca,'FontSize',14);
        xlabel('Num Fx','FontSize',18);
        ylabel('Physical D_{max} (Gy)','FontSize',18);
        
         if do_print,
            
         set(cur_fig,'Color','w');
         export_fig(cur_fig,...
            [fig_basename,'_box_dmax_vs_nfx'],'-pdf');
        disp(['Saving ',fig_basename,'_box_dmax_vs_nfx.pdf']);
        end
       
        
        
    end
end
end
