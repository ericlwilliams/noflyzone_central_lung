function NFz_LoadAlpha2BetaDoses
close all;
% prepare

%fn = {'BPx_DiVj_DVHs_fx-1_a2bInf.mat'};
a2b_corr='BED';
    


screen_size=get(0,'ScreenSize');

%structures = {'PBT' 'ILUNG' 'ESOPHAGUS' 'HEART' 'NFZ' 'LUNGS'};
structures = {'ESOPHAGUS'};

%toxicities = {'rp','pultox','esotox','lclfail'};
toxicities = {'esotox'};

fp = 'Z:\elw\MATLAB\nfz_analy\meta_data\';
fp_save = 'Z:\elw\MATLAB\nfz_analy\meta_data\a2b_data\';

for i=1:length(toxicities)
    
    for j=1:length(structures)
        tic;
        cur_fig_ctr = (10*i)+j-1;
        
        fprintf('\n');
        disp(['Tox: ',toxicities{i}]);
        disp(['Struct: ',structures{j}]);
        disp(['Counter: ',num2str(cur_fig_ctr)]);
        fprintf('\n');
        
        save_basename = [fp_save,'NFZ_',...
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
        
    for n=1:length([CGobj.mGrp]);
        CGobj.mGrp(n).mBeta2AlphaCorrection = a2b_corr;
    end
    
    a2b_range = [.1:.1:30];
    % #pts x #a2b values
    a2b_doses = cell(length(CGobj.mGrp),length(a2b_range)+1);
    a2b_vols = cell(length(CGobj.mGrp),length(a2b_range)+1);
    
    %last entry is a2b=Inf
    a2b_doses(:,end) = {CGobj.mGrp.mDoseBins_LQ}';
    a2b_vols(:,end) = {CGobj.mGrp.mVolCum}';
    
    n_pt = CGobj.mNumInGrp;
    
    for k=1:length(a2b_range),
        
        tmpCGobj = CGobj;       
        tmpCGobj.mBeta2Alpha = 1/a2b_range(k);
        a2b_doses(:,k) = {tmpCGobj.mGrp.mDoseBins_LQ}';
        a2b_vols(:,k) ={tmpCGobj.mGrp.mVolCum}'; 
        clear tmpCGobj;
    end
   
     a2b_range = [a2b_range Inf];
  
     
      a2b_d05 = inf(n_pt,length(a2b_range));
      a2b_d35 = inf(n_pt,length(a2b_range));
      a2b_dmax = inf(n_pt,length(a2b_range));
      a2b_dmean = inf(n_pt,length(a2b_range));

    
        for l=1:length(a2b_range)%loop over all a/b dose bins
            cur_doses = [a2b_doses{:,l}];
            cur_vols = [a2b_vols{:,l}]; %501 x 125
        
            for m=1:n_pt
                vol = cur_vols(:,m);
                dose = cur_doses(:,m);
                a2b_dmean(m,l) = mean(dose);
                
                v05 = vol<5; %% < 5 cc
                d05_inds = find(v05);
                d05 = min(dose(d05_inds));
                a2b_d05(m,l) = d05;
            
                v35 = vol<3.5; %% < 5 cc
                d35_inds = find(v35);
                d35 = min(dose(d35_inds));
                a2b_d35(m,l) = d35;
            
                vol(~vol)=nan;
                nan_inds = find(isnan(vol));
                if ~isempty(nan_inds)
                    min_ind = nan_inds(1)-1;
                else
                    min_ind = length(vol);
                end
            %[~,min_ind] = min(vol);
               a2b_dmax(m,l) = dose(min_ind);
            end
    end
        
     %a2b_max_doses = [cellfun(@(x) max(x),a2b_doses)];
 %a2b_max_doses = dmax;
     %a2b_mean_doses = [cellfun(@(x) mean(x),a2b_doses)];
     
     
     clear a2b_doses;
     
     
     
     disp(['Saving ' save_basename '_a2b_dosebins.mat...']);
     save([save_basename '_a2b_data.mat'],...
         'a2b_dmax',...
         'a2b_dmean',...
          'a2b_d05',...
          'a2b_d35',...
         'a2b_range');
     clear a2b_max_doses;
     clear a2b_mean_doses;
     toc;
    end
end
end