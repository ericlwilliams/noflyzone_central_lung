function NFz_LoadDVHs
% DVH is saved in 3 cell columns, the first stores the patient name,
% the second the original dose, original differential volume, cumulative volume from oriignal diffirential volume
% the third the resampled dose, differential volume, and cumulative volume
tic;
dose_step=0.50;
pt_added=0;

do_tcp_exclude=true;
%tmp
if do_tcp_exclude
    filepathname='Z:/elw/MATLAB/original_data/NFZ/Central_tumor_dataset_tcp_4-22-13'
    pathname='Z:/elw/MATLAB/original_data/NFZ/nfz_tcp_dvhs/';
else
    filepathname='Z:/elw/MATLAB/original_data/NFZ/Central_tumor_dataset_4-22-13'
    pathname='Z:/elw/MATLAB/original_data/NFZ/nfz_dvhs/';
end


[~,~,xlsraw]=xlsread(filepathname);
save(filepathname,'xlsraw');

% read data from each .TXT file

%structures = {'CLUNG' 'ESOPHAGUS' 'GTV' 'HEART'...
%            'ILUNG' 'NFZ' 'PBT' 'PTV' 'TRACHEA' 'LUNGS'};

structures = {'PTV' 'GTV'};

% to_exclude = {'Byrne_35140428_CLUNG', 'Caruso_35359984_CLUNG',...
%     'DEGRASSE_35061243_CLUNG',...
%     'Rota_35367615_CLUNG',...
%     'Rota_35367615_ESOPHAGUS',...
%     'Rota_35367615_GTV',...
%     'Rota_35367615_HEART',...
%     'Rota_35367615_ILUNG',...
%     'Rota_35367615_LUNGS',...
%     'Rota_35367615_NFZ',...
%     'Rota_35367615_PBT',...
%     'Rota_35367615_PTV',...
%     'Rota_35367615_TRACHEA',...
%     'COPPINS_00997860_GTV',...
%     'MACK_35304193_ILUNG',...
%     'Mutsch_35256267_ILUNG',...
%     'Fabiani_35273833_NFZ',...
%     'Magan_35171817_NFZ',...
%     'Solomon_35186653_NFZ',...
%     'Caruso_35359984_PBT',...
%     'Kokoris_00195114_PBT',...
%     'ARENA_00509113_TRACHEA',...
%     'Harmonay_00191050_TRACHEA',...
%     'Matofsky_00518051_TRACHEA'};
to_exclude = {''};

 % load DVH names in d
d=dir(pathname);
if isempty(d)
    disp(['!! Bad DVH location name', pathname]);
    return
end
   
f=false(size(d,1),1);
for k=1:size(d,1)
    if isempty(findstr(d(k).name,'.TXT')) || isdir(d(k).name)
        f(k)=true;
    end
end
d(f)=[];
dvh_file_names = {d.name}';
% patient last name
f=cellfun(@(x) strcmp(x,'Patient Last Name'),xlsraw(1,:)); f=find(f); f=f(1); % found the column of "Patient Last Name"
   
ptnames=xlsraw(2:end,f); % patient last name


%% get MRNs
g=cellfun(@(x) strcmp(x,'MRN'),xlsraw(1,:)); g=find(g); g=g(1); % found the column of "MRN"
mrns=xlsraw(2:end,g);

% convert all mrn data to strings
containsNumbers = cellfun(@isnumeric,mrns);
mrns(containsNumbers) = cellfun(@num2str,mrns(containsNumbers),'UniformOutput',false);
orig_mrns=mrns;
mrns=regexprep(mrns,'-.','');



for i=1:length(structures)

    
    DVH=cell(length(ptnames),4); 
    flg_nofile=false(size(ptnames,1),1);  

structure = structures{i};
disp(['Loading ',structure]);


    for k=1:size(ptnames,1)
        f=cellfun(@(x) strcmp(x,mrns{k}), mrns); f=find(f);
    
        ptname=ptnames{k};
        mrn=mrns(k);
        % remove any '-' from mrn
        dash_ind = strfind(mrn,'-');
        if ~isempty(dash_ind{1}),
            mrn = mrn{1}(1:dash_ind{1}-1);
        else
        mrn=mrn{1};
        end
    
        dvh_file_name = strcat(ptname,'_',mrn,'_',structure);
    
        excluding=false;
        for j=1:length(to_exclude)
            if isequal(to_exclude{j},dvh_file_name)
                
                %% TMP
                if ~isempty(findstr('CLUNG',dvh_file_name))
                    dvh_file_name='Asseoff_35154734_CLUNG';
                elseif ~isempty(findstr('ESOPHAGUS',dvh_file_name))
                    dvh_file_name='Asseoff_35154734_ESOPHAGUS';
                elseif ~isempty(findstr('GTV',dvh_file_name))
                    dvh_file_name='Asseoff_35154734_GTV';
                elseif ~isempty(findstr('PTV',dvh_file_name))
                    dvh_file_name='Asseoff_35154734_PTV';
                elseif ~isempty(findstr('HEART',dvh_file_name))
                    dvh_file_name='Asseoff_35154734_HEART';
                elseif ~isempty(findstr('NFZ',dvh_file_name))
                    dvh_file_name='Asseoff_35154734_NFZ';
                elseif ~isempty(findstr('PBT',dvh_file_name))
                    dvh_file_name='Asseoff_35154734_PBT';
                elseif ~isempty(findstr('ILUNG',dvh_file_name))
                    dvh_file_name='Asseoff_35154734_ILUNG';
                elseif ~isempty(findstr('TRACHEA',dvh_file_name))
                    dvh_file_name='Asseoff_35154734_TRACHEA';
                elseif ~isempty(findstr('LUNGS',dvh_file_name))
                    dvh_file_name='Asseoff_35154734_LUNGS';
                else
                    excluding=true;
                end
                break;
            end
        end
        if excluding
            disp(['Excluding ',dvh_file_name]);
            continue;
        end
        
        
        
        if length(f)>1
        %disp(['repeated last names: ',ptnames{k}]);
        if length(f)>2
            disp([ptnames{k},' repeated ',num2str(length(f)),' times, quitting...']);
        end
        
        loc_str = CWLocationString(ptname,txsites{k});
        
        dvh_file_name = strcat(dvh_file_name,'_',loc_str);
        end
    % No duplicates, look up DVH with last name and MRN, strip any extra '-' first
    
    
    cur_dvh = strfind(dvh_file_names,dvh_file_name);

    disp(['Loading ',dvh_file_name,'...']);
    dvh_name=dvh_file_names{~cellfun(@isempty,cur_dvh)};
    
    
    dvh_loc = strcat(pathname,dvh_name);
    %disp(['Loading ',dvh_loc,'...']);
    
    
    % original DVH
    try
        %df=textread(strcat(pathname,ptnames{k},fnSuffix),'%s'); % read the file into cells df (data from file)
        df=textread(dvh_loc,'%s'); % read the file into cells df (data from file)
    catch ME
        disp(strcat(pathname,ptnames{k},fnSuffix)); disp(ME.message);
        %             for errcount=1:length(ME.stack)
        %                 disp(ME.stack(errcount));
        %             end
        flg_nofile(i,k)=true; continue;
    end
    
    DVH{k,1}=ptnames{k};
    DVH{k,4}=orig_mrns{k};
    %DVH{k,i,5}=structure;    
    
    df(1)=[]; % the first cell is the description text of the data, remove it
    data=cellfun(@(x) str2num(x), df);
    
    % original DVH as step function
    dose_org=data(1:4:end)/100; vol_org=data(2:4:end); % /100 to change cGy to Gy
    % cumulative dvh
    vol_cumulative=vol_org;
    %         vol_cumulative(1:end-1)=vol_cumulative(1:end-1).*diff(dose_org);
    vol_cumulative=flipud(vol_cumulative); vol_cumulative=cumsum(vol_cumulative); vol_cumulative=flipud(vol_cumulative);
    DVH{k,2}=[dose_org,vol_org,vol_cumulative];
    
    % regenerate the differentail DVH by adding new bins (which is the expected bins)
    dose_interpolation=(0:dose_step:max(DVH{k,2}(:,1))+dose_step)'; % locations of new bins
    vol_interpolation=zeros(size(dose_interpolation,1),1);
    f=ismember(dose_interpolation, dose_org); % the doses alreadyin the original dose vector should be removed in new bins
    dose_interpolation(f)=[]; vol_interpolation(f)=[];
    dose_new=[dose_org; dose_interpolation]; % add the new bins and sort it so new bins and old ones are arranged properly
    vol_new=[vol_org; vol_interpolation];
    [dose_new,g]=sort(dose_new); vol_new=vol_new(g);
    for n=1:size(dose_interpolation,1)-1 % the end point can not be interpolated, and should be zero since the last dose "bin" should correspond no volume
        f=find(dose_interpolation(n)>dose_new); % the location of interpolated dose in the DVH with combined bins, which is applied to determine the left end point (d(i-1)), d(i-1)=f(end)
        g=find(dose_interpolation(n)<dose_org); % the location of interpolated dose in the DVH with original bins, which is applied to determine the right end point (d(i)), d(i)=g(1)
        vol_new(f(end)+1) = (dose_org(g(1))-dose_interpolation(n)) / (dose_org(g(1))-dose_new(f(end))) * vol_new(f(end)); % v(j)=(d(i)-d(j))/(d(i)-d(i-1))*v(i-1)
        vol_new(f(end)) = (dose_interpolation(n)-dose_new(f(end))) / (dose_org(g(1))-dose_new(f(end))) * vol_new(f(end)); % v(i-1)=(d(j)-d(i-1))/(d(i)-d(i-1))*v(i-1)
    end
    % regenerate the differential DVH using the new bin
    dose_interpolation=(0:dose_step:max(DVH{k,2}(:,1))+dose_step)'; % locations of new bins
    vol_interpolation=zeros(size(dose_interpolation,1),1);
    for n=1:size(dose_interpolation,1)-1
        f = ( dose_interpolation(n)<=dose_new & dose_interpolation(n+1)>dose_new); % the range of new bins in the interpolated bin
        vol_interpolation(n) = sum(vol_new(f));
    end
    % cumulative dvh
    vol_cumulative=vol_interpolation;% vol_cumulative(1:end-1)=vol_cumulative(1:end-1).*diff(dose_interpolation);
    vol_cumulative=flipud(vol_cumulative); vol_cumulative=cumsum(vol_cumulative); vol_cumulative=flipud(vol_cumulative);
    DVH{k,3}=[dose_interpolation,vol_interpolation,vol_cumulative];
    
    
    
    % plot the result
    figure(1); a=data(1:2:end)/100; b=data(2:2:end); plot(a,b); hold on;
    stairs(DVH{k,2}(:,1),DVH{k,2}(:,2),'r');
    stairs(DVH{k,3}(:,1),DVH{k,3}(:,2),'k');
  
    grid on; title([ptnames{k}, ' ',structure,' ',num2str(k)]); hold off;
    
    figure(2); stairs(DVH{k,2}(:,1),DVH{k,2}(:,3)); hold on;
    stairs(DVH{k,3}(:,1),DVH{k,3}(:,3),'r'); hold off;
    grid on; pause(0);
    if i==1,
        pt_added=pt_added+1;
    end
end


 if do_tcp_exclude
    fn=['Z:\elw\MATLAB\nfz_analy\meta_data\NFZ_TCP_',structure,'_DVHs.mat'];
 else
    fn=['Z:\elw\MATLAB\nfz_analy\meta_data\NFZ_',structure,'_DVHs.mat'];
 end
 
if 0
    disp(['Saving ',fn]);
    save(fn,'DVH');
end

end

disp(' ');

disp(['total number of patients in .xls file: ',num2str(pt_added)]);
disp(['total number of patients with the DVH files: ', num2str(sum(~flg_nofile))]);
disp(['total number of patients without the DVH files: ', num2str(sum(flg_nofile))]);
toc;

end