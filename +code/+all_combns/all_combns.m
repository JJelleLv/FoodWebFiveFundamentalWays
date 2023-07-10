function foldername = all_combns(NRspec,resultFolder)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Function that stores all combinations of species %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% make new folder or go back if already exists
foldername = ...
    sprintf('%s%sALL_comb_S%d%s',resultFolder,filesep,NRspec,filesep);
if (~exist(foldername, 'dir'))
    mkdir(foldername);
% else
%     return
end

for SpecNR=1:NRspec
            
    ALL_comb = nchoosek(1:NRspec,SpecNR);
    ALL_comb = round(ALL_comb);
    NRcomb=length(ALL_comb(:,1));   %#ok<NASGU>
    
    save(sprintf('%sALL_comb_%d', foldername, SpecNR), ...
        'ALL_comb', 'NRcomb','-v7.3')
    
    %% print progress        
    fprintf('Spec COMBs: %d\n', SpecNR);

end
