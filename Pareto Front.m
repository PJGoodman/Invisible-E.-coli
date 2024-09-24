
% Set MOSEK as the solver for LP, QP, and MILP problems
changeCobraSolver('mosek', 'LP');
changeCobraSolver('mosek', 'QP');
changeCobraSolver('mosek', 'MILP');

% List of model names and corresponding biomass reactions, can add more
% models as desired
ModelNames = {'new_reduced_iML1515.mat', 'iML1515.mat', 'iYO844.mat'};
BiomassRxns = {'BIOMASS_Ec_iML1515_core_75p37M','BIOMASS_Ec_iML1515_core_75p37M','BIOMASS_BS_10'};

% Diet constraints
dietConstraints = {
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

% Index of iML1515 in ModelNames
iML1515_index = find(strcmp(ModelNames, 'iYO844.mat'));

% Iterate over the pairs of models involving iML1515 and perform the Pareto optimality analysis
dinc = 0.1; % Increment for Pareto frontier steps
for i = 1:length(ModelNames)
    if i == iML1515_index
        continue; % Skip the loop when i equals iML1515_index (to avoid iML1515 with itself)
    end
    models = {};
    bioID = {};
    nameTagsModels = {};

    % Load iML1515 as the first model
    model1 = readCbModel(ModelNames{iML1515_index});
    models{1,1} = model1;
    bioID{1,1} = BiomassRxns{iML1515_index};
    nameTagsModels{1,1} = strcat(ModelNames{iML1515_index}, '_');

    % Load the second model
    model2 = readCbModel(ModelNames{i});
    models{2,1} = model2;
    bioID{2,1} = BiomassRxns{i};
    nameTagsModels{2,1} = strcat(ModelNames{i}, '_');

    % Create the paired model
    [pairedModel] = createMultipleSpeciesModel(models, 'nameTagsModels', nameTagsModels);

    % Couple the biomass reactions to ensure simultaneous optimization
    pairedModel = coupleRxnList2Rxn(pairedModel, pairedModel.rxns(strmatch(nameTagsModels{1,1}, pairedModel.rxns)), strcat(ModelNames{iML1515_index}, '_', bioID{1,1}));
    pairedModel = coupleRxnList2Rxn(pairedModel, pairedModel.rxns(strmatch(nameTagsModels{2,1}, pairedModel.rxns)), strcat(ModelNames{i}, '_', bioID{2,1}));

    % Apply diet constraints
    pairedModel = useDiet(pairedModel, dietConstraints);

    % Compute Pareto optimality frontier
    ParetoFrontier = computeParetoOptimality(pairedModel, strcat(ModelNames{iML1515_index}, '_', bioID{1,1}), strcat(ModelNames{i}, '_', bioID{2,1}), 'dinc', dinc, 'FVAflag', false);
end
