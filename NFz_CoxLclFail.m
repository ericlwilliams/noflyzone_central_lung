function NFz_CoxLclFail
tic; close all;
warning('off'); % ignore negative data warning
screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check flags!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do_print = true;
do_tcp_exclude = true; %exclude 10 patients with low Rx (<= 3000 cGy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
excluded_mrns = {'35178082' '35093359' '35170866-1' '35190407'...
                '35176922' '35273833' '35235354' '35237016'...
                '35229199' '35233998'};

disp('**** START Flags ****');
disp([num2str(do_print),' - do_print']);
disp([num2str(do_tcp_exclude),' - do_tcp_exclude']);
disp('**** END Flags ****');
disp(sprintf('\n'));

fig_loc = 'Z:/elw/MATLAB/nfz_analy/slides/figures/latest/';

%structures = {'PBT' 'ILUNG' 'ESOPHAGUS' 'HEART' 'NFZ' 'LUNGS' 'GTV' 'PTV'};
%structures = {'PTV' 'GTV'};
structures = {'PTV' 'GTV'};
%toxicities = {'rp','pultox','esotox','lclfail};
toxicities = {'lclfail'};

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
        
        mrns = {CGobj.mGrp.mID};
        
        excluded_pts = ismember(mrns,excluded_mrns);
        
        %get compdate and censoring info for cox analy
        f2 = ~cellfun('isempty',{CGobj.mGrp.mDateComp}); % patients with no complication date
        f3 = ~cellfun('isempty',{CGobj.mGrp.mDateLastFollowup}); % patients with no last follow up date
        compdate = inf(CGobj.mNumInGrp,1);
        lastfollowup = inf(CGobj.mNumInGrp,1);
        compdate(f2) = ([CGobj.mGrp(f2).mDateComp] - [CGobj.mGrp(f2).mDateBaseline])' / 30;
        lastfollowup(f3) = ([CGobj.mGrp(f3).mDateLastFollowup] - [CGobj.mGrp(f3).mDateBaseline])' / 30;
        compdate = min( lastfollowup, compdate );
        flgcensor = [CGobj.mGrp.mFlgCensor]';
        flgcomp = ~flgcensor;
        
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
        
        logl = [VDxCox.logl]';
        
        [~,doseloc]=max(logl(~flgCox)); % the best fitting of Cox model
        VDxCox = VDxCox(~flgCox);
        x_vdx = x_vdx(~flgCox);
        
        best_vx = VDxCox(doseloc).data_exposure;
        
        disp(['Best V_{D} at V_{',num2str(x_vdx(doseloc)),'}']);
        
        [sort_best_vx,sort_best_vx_inds] = sort(best_vx);
        sort_excluded_pts = excluded_pts(sort_best_vx_inds)';
        
        flg_excluded_comp = logical(flgcomp.*sort_excluded_pts);
        flg_excluded_cens = logical(flgcensor.*sort_excluded_pts);
        
        flg_not_excluded_comp = logical(flgcomp.*~sort_excluded_pts);
        flg_not_excluded_cens = logical(flgcensor.*~sort_excluded_pts);
        
        x_axis = [1:length(sort_best_vx)];
        cur_fig=figure(cur_fig_ctr); clf reset;
        set(gcf,'Position',ss_four2three); hold on;
        %not excluded censored (blue empty)
        h_not_excl_cens = plot(x_axis(flg_not_excluded_cens),sort_best_vx(flg_not_excluded_cens),'o','MarkerSize',10);
        % not excluded complication (red empty)
        h_not_excl_comp = plot(x_axis(flg_not_excluded_comp),sort_best_vx(flg_not_excluded_comp),'ro','MarkerSize',10);
        % excluded cens
        h_excl_cens = plot(x_axis(flg_excluded_cens),sort_best_vx(flg_excluded_cens),'o','MarkerSize',10,'MarkerFaceColor','b');
        % excluded ccomp
        h_excl_comp = plot(x_axis(flg_excluded_comp),sort_best_vx(flg_excluded_comp),'ro','MarkerSize',10,'MarkerFaceColor','r');
        
        xlabel(['V_{',num2str(x_vdx(doseloc),3),'} Rank'],'FontSize',18);
        ylabel(['V_{',num2str(x_vdx(doseloc),3),'} (%)'],'FontSize',18);
        set(gca,'FontSize',16);
        title([structures{j},' V_{',num2str(x_vdx(doseloc),3),'}'],'FontSize',20);
    
        if do_tcp_exclude
            
        h_lgnd=legend([h_not_excl_cens h_not_excl_comp],'Not excluded','Not excluded w/ lcl fail',...
            'Location','Best');
        else
        h_lgnd=legend([h_excl_cens h_excl_comp h_not_excl_cens h_not_excl_comp],'Excluded','Excluded w/ lcl fail','Not excluded','Not excluded w/ lcl fail',...
            'Location','Best');
        end
        
        set(h_lgnd,'FontSize',16);
        
        if do_print,
            set(cur_fig,'Color','w');
         export_fig(cur_fig,...
            [fig_basename,'_best_vds'],'-pdf');
        disp(['Saving ',fig_basename,'_best_vds.pdf']);
        end
        
    end
end
end