% Set MOSEK as the solver for LP problems
changeCobraSolver('mosek', 'LP');

% Set MOSEK as the solver for QP problems
changeCobraSolver('mosek', 'QP');

% Set MOSEK as the solver for MILP problems
changeCobraSolver('mosek', 'MILP');

% List of model names and corresponding biomass reactions
ModelNames = {'iAF987.mat', 'iEK1008.mat','iIT341.mat','iJN1463.mat','iML1515.mat', 'iMM904.mat','iNJ661.mat','iSBO_1134.mat', 'iYO844.mat', 'STM_v1_0.mat', };
BiomassRxns = {'BIOMASS_Gm_GS15_core_79p20M', 'BIOMASS__2','BIOMASS_HP_published','BIOMASS_KT2440_WT3','BIOMASS_Ec_iML1515_core_75p37M','BIOMASS_SC5_notrace', 'BIOMASS_Mtb_9_60atp','BIOMASS_Ec_iJO1366_core_53p95M', 'BIOMASS_BS_10', 'BIOMASS_iRR1083'};

dietConstraintsIndividual = {
    'EX_pi[u]', '-1000', '1000';  
    'EX_co2[u]', '-1000', '1000';  
    'EX_fe3[u]', '-1000', '1000';  
    'EX_h[u]', '-1000', '1000';    
    'EX_mn2[u]', '-1000', '1000';  
    'EX_fe2[u]', '-1000', '1000';  
    'EX_glc__D[u]', '-10', '1000';
    'EX_zn2[u]', '-1000', '1000';  
    'EX_mg2[u]', '-1000', '1000';  
    'EX_ca2[u]', '-1000', '1000';  
    'EX_ni2[u]', '-1000', '1000';  
    'EX_cu2[u]', '-1000', '1000';  
    'EX_sel[u]', '-1000', '1000';  
    'EX_cobalt2[u]', '-1000', '1000'; 
    'EX_h2o[u]', '-1000', '1000';  
    'EX_mobd[u]', '-1000', '1000'; 
    'EX_so4[u]', '-1000', '1000';  
    'EX_nh4[u]', '-1000', '1000';  
    'EX_k[u]', '-1000', '1000';    
    'EX_na1[u]', '-1000', '1000';  
    'EX_cl[u]', '-1000', '1000';  
    'EX_o2[u]', '-1000', '1000';   
    'EX_tungs[u]', '-1000', '1000';  
    'EX_slnt[u]', '-1000', '1000'; 
    'EX_h2co3[u]', '-1000', '1000';
    'EX_cit[u]', '-1000', '1000';
    'EX_pime[u]', '-1000', '1000';
    'EX_pheme[u]', '-1000', '1000';
    'EX_cbl1[u]', '-1000', '1000';
    'EX_nac[u]', '-1000', '1000';
    'EX_thm[u]', '-1000', '1000';
    'EX_glyc[u]', '-1000', '1000';

};

% Folder containing the pairwise models
pairwiseModelFolder = 'ExampleFolder';

% Step 1: Join models pairwise
joinModelsPairwiseFromList(ModelNames, pairwiseModelFolder, 'biomasses', BiomassRxns);

% Step 2: Simulate pairwise interactions with diet constraints
[pairwiseInteractions, pairwiseSolutions] = simulatePairwiseInteractions(pairwiseModelFolder, 'inputDiet', dietConstraintsIndividual);

% Display the pairwise interactions
disp(pairwiseInteractions);

% Step 3: Extract only the pairings that include iML1515
iML1515_index = find(contains(ModelNames, 'iML1515'));
filteredInteractions = pairwiseInteractions(contains(pairwiseInteractions(:, 1), 'iML1515') | contains(pairwiseInteractions(:, 2), 'iML1515'), :);

% Assuming filteredInteractions(1, :) are not proper variable names, define them manually
variableNames = {'pairedModelID', 'ModelID1', 'ModelID2', 'pairedGrowth_Model1', 'pairedGrowth_Model2', 'singleGrowth_Model1', 'singleGrowth_Model2', 'Outcome_Model1', 'Outcome_Model2', 'Total_Outcome'};

% If filteredInteractions(1, :) are already correct names, you can use:
% variableNames = filteredInteractions(1, :);

% Convert the cell array to a table using the valid variable names
filteredTable = cell2table(filteredInteractions(2:end, :), 'VariableNames', variableNames);

% Save the table to a CSV file
writetable(filteredTable, 'iML1515_pairwise.csv');

% Step 5: Generate a pie chart of interaction types for pairings involving iML1515
figure('rend', 'painters', 'pos', [10 10 900 600])

typesIA = unique(filteredInteractions(2:end, 10)); % Unique interaction types
dat = zeros(length(typesIA), 1);

for j = 1:length(typesIA)
    dat(j) = sum(strcmp(filteredInteractions(:, 10), typesIA{j}));
end

pie(dat)
set(gca, 'FontSize', 10)
title('Pairwise interactions involving iML1515')

legend1 = legend(typesIA);
set(legend1, 'Position', [0.42 0.45 0.2 0.2], 'FontSize', 12)

