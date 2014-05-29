function NFz_Anonymize

% structures/toxicity doesn't matter, anonymization on cohort
%structures = {'ILUNG' 'ESOPHAGUS' 'HEART' 'NFZ' 'PBT' 'LUNGS'};
structures = {'ESOPHAGUS'};
%toxicities = {'rp','pultox','esotox'};
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
      
     
     %% load data
        fn = ['NFZ_',structures{j},'_',toxicities{i},'_a2b',a2b{1},'_data.mat'];
        load(strcat(fp,fn),'CGobj_org');

          CGgrp = CGobj_org.mGrp;

    % encrypt mID 
   anon_mrns = cell(length(CGgrp),2);
   for j=1:length(CGgrp)
       
       
       tmp_id = CGgrp(j).mID;
        tmp_md5 = DataHash(tmp_id);
        anon_mrns{j,1} = tmp_id;
        anon_mrns{j,2} = tmp_md5;
       
        CGgrp(j).mID = tmp_md5;
    
   end
   
   [~,md5_inds] = sort({anon_mrns{:,2}});
   new_ids = mat2cell(md5_inds,1,ones(1,size(md5_inds,2)));
   anon_mrns(:,2) = new_ids';
   
   for k=1:length(CGgrp)
       CGgrp(k).mID = num2str(anon_mrns{k,2});
   end
   
   CGobj_org.mGrp = CGgrp;
   xlswrite([fp, 'anon_data\CWp_mrn_encrypt.xlsx'],anon_mrns)
   anon_fn=[fp,'anon_data\',strrep(fn,'_data','_anon_data')];
   
   save(anon_fn,'CGobj_org');
end
end