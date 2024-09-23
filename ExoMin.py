##
import numpy as np
def knockout_and_check_viability(model, genes):
    with model:
        cobra.manipulation.knock_out_model_genes(model, genes)
        solution = model.slim_optimize()
        return solution

def build_pathways(model, metabolite, max_layers):
    if not model.metabolites.has_id(metabolite.id):
        return {}

    q = Queue()
    q.put((metabolite, [metabolite.id], []))

    pathways = {}

    while not q.empty():
        current_metabolite, current_path, current_reactions = q.get()

        if len(current_path) >= max_layers or (
                current_metabolite.id in pathways and current_reactions == pathways[current_metabolite.id][0]):
            continue

        visited_metabolites = set(current_path)

        for reaction in current_metabolite.reactions:
            if reaction.id in biomass_reaction_names:
                continue

            if current_metabolite.id in common_metabolites_to_exclude:
                next_reactions = current_reactions + [reaction.id]
                pathways[current_metabolite.id] = (current_reactions, current_path)
                continue

            next_reactions = current_reactions + [reaction.id]

            if current_metabolite in reaction.products:
                for metabolite_in_reaction, coefficient in reaction.metabolites.items():
                    if coefficient < 0 and metabolite_in_reaction.id not in visited_metabolites:
                        next_metabolite = metabolite_in_reaction
                        next_path = current_path + [metabolite_in_reaction.id]
                        q.put((metabolite_in_reaction, next_path, next_reactions))
                        visited_metabolites.add(metabolite_in_reaction.id)

            if reaction.reversibility:
                for metabolite_in_reaction, coefficient in reaction.metabolites.items():
                    if coefficient > 0 and metabolite_in_reaction.id not in visited_metabolites:
                        next_metabolite = metabolite_in_reaction
                        next_path = current_path + [metabolite_in_reaction.id]
                        q.put((metabolite_in_reaction, next_path, next_reactions))
                        visited_metabolites.add(metabolite_in_reaction.id)

        pathways[current_metabolite.id] = (current_reactions, current_path)

    return pathways


def generate_pathways(model, metabolite_ids, max_layers):
    pathways_dict = {}

    for metabolite_id in metabolite_ids:
        try:
            metabolite = model.metabolites.get_by_id(metabolite_id)
            pathways_dict[metabolite_id] = build_pathways(model, metabolite, max_layers)
        except KeyError:
            continue

    return pathways_dict


def load_model_and_get_extracellular_metabolites(model_path, common_metabolites_to_exclude):
    """
    Load the model and get the list of extracellular metabolites converted to their cytosolic counterparts.

    Parameters:
    model_path (str): Path to the SBML model file.
    common_metabolites_to_exclude (set): Set of common metabolites to exclude.

    Returns:
    list: List of cytosolic counterparts of extracellular metabolites, excluding common metabolites.
    """
    # Load the model
    model = read_sbml_model(model_path)

    # Define the list of extracellular metabolites to consider
    extracellular_metabolites = [metabolite.id for metabolite in model.metabolites if metabolite.compartment == 'e']

    # Convert to cytosolic counterparts and exclude common metabolites
    extracellular_metabolites_c = []
    for met in extracellular_metabolites:
        new_met = met[:-1] + 'c'
        if new_met not in common_metabolites_to_exclude:
            extracellular_metabolites_c.append(new_met)

    return extracellular_metabolites_c

from cobra.exceptions import Infeasible
from warnings import warn

def test_carbon_sources(model, carbon_sources, min_growth_rate, pathways, anaerobic):
    KO_Dicts = []
    # Set a desired growth rate or objective value

    med = model.medium

    # Get the value of the original carbon source
    original_carbon_value = med.get("EX_glc__D_e", 0)
    med["EX_glc__D_e"] = 0

    # Remove the original carbon source from the medium
    if anaerobic:
        med['EX_o2_e'] = 0

    for carbon_source in carbon_sources:
        with model:
            try:
                medium = med.copy()
                medium[carbon_source] = original_carbon_value
                model.medium = medium

                # Optimize the model to get the growth rate
                growth_rate = model.slim_optimize()

                if not np.isnan(growth_rate) and growth_rate > min_growth_rate:
                    KO_Dict = create_gene_KO_list(model, pathways)
                    KO_Dicts.append(KO_Dict)
                    print(f'Sufficient growth rate on {carbon_source}: {growth_rate}')
                else:
                    print(f'Low growth rate on {carbon_source}: {growth_rate}')
            except Infeasible:
                print(f'Model is infeasible for carbon source {carbon_source}.')
            except Exception as e:
                print(f'An error occurred for carbon source {carbon_source}: {e}')

    return KO_Dicts
def check_gene_relationship(reaction):
    """
    Determine the gene relationship for a given reaction.

    Parameters:
    reaction (cobra.core.Reaction): The reaction to analyze.

    Returns:
    str: The type of gene relationship ('AND relationship', 'OR relationship', 'Unknown or single gene relationship').
    """
    rule = reaction.gene_reaction_rule
    if ' and ' in rule:
        return 'AND relationship'
    elif ' or ' in rule:
        return 'OR relationship'
    else:
        return 'Unknown or single gene relationship'


def perform_direct_knockout_analysis(extracellular_metabolites, common_metabolites_to_exclude,
                                   min_growth_rate, model):
    """
    Perform gene knockout analysis for metabolites in the given model.

    Parameters:
    extracellular_metabolites (list): List of extracellular metabolites to consider.
    common_metabolites_to_exclude (set): Set of common metabolites to exclude.

    Returns:
    Direct_KO_Genes_Dict
    """
    Direct_KO_Genes_Dict = {}

    for met_id in extracellular_metabolites:
        with model:
            metabolite_genes = []

            try:
                met = model.metabolites.get_by_id(met_id)
                reactions = [reaction for reaction in met.reactions if met in reaction.products]
                if len(reactions) == 0:
                    continue
                Direct_KO_Genes_Dict.setdefault(met.id, {})

                for reaction in reactions:
                    reaction_id = reaction.id
                    if reaction_id.endswith('pp'):
                        continue
                    AND = (check_gene_relationship(reaction) == 'AND relationship')
                    genes = reaction.genes
                    if len(genes) == 0:
                        continue
                    Direct_KO_Genes_Dict[met.id].setdefault(reaction.id, [])
                    if AND:
                        gene_impact_list = []
                        low_impact_gene = None
                        for gene in genes:
                            reaction_length = len(gene.reactions)
                            gene_impact_list.append(reaction_length)
                            solution = knockout_and_check_viability(model, [gene.id])
                            if reaction_length == min(gene_impact_list) and solution > min_growth_rate and not np.isnan(solution):
                                low_impact_gene = gene
                        if low_impact_gene:
                            metabolite_genes.append(low_impact_gene)
                            Direct_KO_Genes_Dict[met.id][reaction.id].append(low_impact_gene)
                    else:
                        solution = knockout_and_check_viability(model,[gene.id for gene in genes])
                        if solution > min_growth_rate and not np.isnan(solution):
                            for gene in genes:
                                Direct_KO_Genes_Dict[met.id][reaction.id].append(gene)
            except:
                continue

    return Direct_KO_Genes_Dict
def evaluate_knockout_impact(model, reaction1_id, reaction2_id, biomass_threshold=0.1):
    with model:
        # Knock out genes associated with the first reaction
        reaction1 = model.reactions.get_by_id(reaction1_id)
        genes_to_knockout = [gene.id for gene in reaction1.genes]
        cobra.manipulation.knock_out_model_genes(model, genes_to_knockout)

        biomass_flux = model.slim_optimize()
        is_lethal = biomass_flux < biomass_threshold or np.isnan(biomass_flux)

        # Evaluate if knockout disrupted the second reaction
        model.objective = reaction2_id
        model.reactions.get_by_id(reaction2_id).upper_bound = 1000.
        second_reaction_flux = model.optimize().objective_value
        is_disrupted = second_reaction_flux == 0 or np.isnan(second_reaction_flux)

        return is_lethal, is_disrupted

non_effective_upstream_KOs = {}
def perform_upstream_knockout_analysis(pathways, model, min_growth_rate, Direct_KO_Genes_Dict):
    Upstream_KO_Dict = {}

    for metabolite_id, paths in pathways.items():
        Upstream_KO_Dict.setdefault(metabolite_id,{})
        if metabolite_id not in Direct_KO_Genes_Dict:
            continue

        for path_name, path_data in paths.items():
            if path_name == metabolite_id:
                continue

            reactions = path_data[0]
            terminal_reaction = model.reactions.get_by_id(reactions[-1])
            terminal_reaction_id = terminal_reaction.id
            if terminal_reaction_id.endswith('pp'):
                continue

            non_effective_upstream_KOs.setdefault(reactions[0], [])
            if terminal_reaction in non_effective_upstream_KOs[reactions[0]]:
                continue

            KO_Genes = terminal_reaction.genes

            if len(KO_Genes) == 0:
                continue

            AND = (check_gene_relationship(terminal_reaction) == 'AND relationship')

            if AND:
                gene_impact_list = []
                low_impact_gene = None

                for gene in KO_Genes:
                    reaction_length = len(gene.reactions)
                    gene_impact_list.append(reaction_length)
                    solution = knockout_and_check_viability(model, [gene.id])
                    if reaction_length == min(gene_impact_list) and solution > min_growth_rate and not np.isnan(solution):
                        low_impact_gene = gene

                if low_impact_gene:
                    KO_Genes = [low_impact_gene]

            if reactions[0] in Direct_KO_Genes_Dict[metabolite_id]:
                if len(Direct_KO_Genes_Dict[metabolite_id][reactions[0]]) != 0:
                    if len(KO_Genes) >= len(Direct_KO_Genes_Dict[metabolite_id][reactions[0]]):
                        non_effective_upstream_KOs[reactions[0]].append(reactions[-1])
                        continue

            is_lethal, is_disrupted = evaluate_knockout_impact(model, reactions[0], reactions[-1])

            if not is_lethal and is_disrupted:
                Upstream_KO_Dict[metabolite_id].setdefault(reactions[0],[])
                for gene in KO_Genes:
                    Upstream_KO_Dict[metabolite_id][reactions[0]].append(gene)
            elif not is_disrupted:
                non_effective_upstream_KOs[reactions[0]].append(reactions[-1])

    return Upstream_KO_Dict


def update_KO_Genes(Upstream_KO_Dict, Direct_KO_Genes_Dict):
    """
    Update KO_Genes_Dict with genes from Upstream_KO_Dict.

    Parameters:
    Upstream_KO_Dict (dict): Dictionary containing metabolite IDs and their pathways with knockout genes.
    Direct_KO_Genes_Dict (dict): Dictionary containing metabolite IDs, reactions, and their associated genes.

    Returns:
    dict: Updated KO_Genes_Dict.
    """
    for metabolite_id, knockout_data in Upstream_KO_Dict.items():
        if metabolite_id not in Direct_KO_Genes_Dict:
            Direct_KO_Genes_Dict[metabolite_id] = {}

        for reaction_id, gene_list in knockout_data.items():
            if reaction_id not in Direct_KO_Genes_Dict[metabolite_id]:
                Direct_KO_Genes_Dict[metabolite_id][reaction_id] = []

            # Append genes from upstream KO data to direct KO data
            for gene in gene_list:
                if gene not in Direct_KO_Genes_Dict[metabolite_id][reaction_id]:
                    Direct_KO_Genes_Dict[metabolite_id][reaction_id].append(gene)

    return Direct_KO_Genes_Dict


def intersect_KO_Genes_Dicts(dicts):
    """
    Intersects a list of dictionaries with the structure of updated_KO_Genes_Dict.

    Parameters:
    dicts (list of dict): The list of dictionaries to intersect.

    Returns:
    dict: A new dictionary containing only the entries found in all input dictionaries.
    """
    if not dicts:
        return {}

    # Start with the intersection dictionary as the first dictionary
    intersection_dict = dicts[0]

    # Iterate over the rest of the dictionaries to find the intersection
    for other_dict in dicts[1:]:
        new_intersection_dict = {}

        for metabolite_id in intersection_dict:
            if metabolite_id in other_dict:
                intersection_reactions = {}
                for reaction_id in intersection_dict[metabolite_id]:
                    if reaction_id in other_dict[metabolite_id]:
                        # Add the reaction and its genes to the intersection dictionary
                        intersection_reactions[reaction_id] = intersection_dict[metabolite_id][reaction_id]

                if intersection_reactions:
                    new_intersection_dict[metabolite_id] = intersection_reactions

        # Update the intersection dictionary
        intersection_dict = new_intersection_dict

    return intersection_dict


from queue import Queue
from cobra.io import read_sbml_model
import cobra
import itertools

from collections import defaultdict



def find_all_combinations_intersections(labeled_KO_Dicts):
    intersections = {}
    # Generate combinations of dictionaries of all sizes
    for r in range(2, len(labeled_KO_Dicts) + 1):
        for combination in itertools.combinations(labeled_KO_Dicts, r):
            labels = [label for label, _ in combination]
            dicts = [d for _, d in combination]
            intersection_label = ' - '.join(labels)
            intersection_dict = intersect_KO_Genes_Dicts(dicts)
            intersections[intersection_label] = intersection_dict
    return intersections

labels = [
    'Aerobic Lactose', 'Aerobic Glucose', 'Aerobic Fructose', 'Aerobic Maltose',
    'Aerobic Glycerol', 'Aerobic Glutamate', 'Aerobic Succinate',
    'Anaerobic Lactose', 'Anaerobic Glucose', 'Anaerobic Fructose', 'Anaerobic Maltose',
    'Anaerobic Glycerol', 'Anaerobic Glutamate', 'Anaerobic Succinate'
]

# Function to label each dictionary in KO_Dicts
def label_KO_Dicts(KO_Dicts, labels):
    labeled_KO_Dicts = [(label, ko_dict) for label, ko_dict in zip(labels, KO_Dicts)]
    return labeled_KO_Dicts

def create_value_set_lengths(intersections):
    value_set_lengths = {}
    for subdict_name, subdict in intersections.items():
        value_set = set()
        for values in subdict.values():
            for genes in values.values():
                value_set.update(gene.id for gene in genes)
        value_set_lengths[subdict_name] = len(value_set)
    return value_set_lengths


biomass_reaction_names = {'BIOMASS_Ec_iML1515_core_75p37M', 'BIOMASS_Ec_iML1515_WT_75p37M'}
common_metabolites_to_exclude = {'gtp_c', 'nadh_c', 'h2o_e', 'h2o_c', 'h2o_p', 'h_p', 'h_c', 'h_e', 'o2', 'hco3_c',
                                 'co2_c', 'co2_e', 'co2_p', 'nad_c', 'nadp_c', 'nadph_c', 'atp_c', 'adp_c', 'amp_c',
                                 'o2_e', 'o2_p', 'o2_c', 'pi_e', 'pi_c', 'pi_p'}


# Set the minimum growth rate
model = read_sbml_model(r'input_path')
growth_rate_initial = model.slim_optimize()
min_growth_rate = 0.45 * growth_rate_initial

# List of carbon source exchange reactions to set to zero, add or remove carbon sources if needed
carbon_sources = [
    'EX_lcts_e',  # Lactose
    'EX_glc__D_e',  # Glucose
    'EX_fru_e',  # Fructose
    'EX_malt_e',  # Maltose
    'EX_glyc_e'  # Glycerol
]

# Generate Extracellular Metabolites
extracellular_metabolites_c = load_model_and_get_extracellular_metabolites(r'input_path', common_metabolites_to_exclude)
print('Extracellular Metabolites Collected')
# Generate Extracellular Pathways
pathways = generate_pathways(model, extracellular_metabolites_c, 5)
print('Pathways Generated')


def create_gene_KO_list(model, pathways):

    # Check Direct KO Results
    Direct_KO_Dict= perform_direct_knockout_analysis(extracellular_metabolites_c,common_metabolites_to_exclude, min_growth_rate,model)
    print('Direct KO Candidates Gathered')
    # Check Upstream KO Results
    Upstream_KO_Dict = perform_upstream_knockout_analysis(pathways, model, min_growth_rate, Direct_KO_Dict)
    print('Upstream KO Candidates Gathered')
    # Combine Upstream / Direct Dictionaries
    updated_KO_Genes_Dict = update_KO_Genes(Direct_KO_Dict, Upstream_KO_Dict)

    return updated_KO_Genes_Dict

import logging.config
logging.config.dictConfig({
    'version': 1,
    'disable_existing_loggers': True,
})


#Create aerobic KO-List
print('Testing Aerobic Conditions')
KO_Dicts_ae = test_carbon_sources(model,carbon_sources, min_growth_rate, pathways, anaerobic=False)

#Create anerboic KO List
print('Testing Anaerobic Conditions')
KO_Dicts_an = test_carbon_sources(model,carbon_sources, min_growth_rate, pathways, anaerobic=True)
#Merge Dictionaries
KO_Dicts = KO_Dicts_ae + KO_Dicts_an

# Apply the labels to KO_Dicts
labeled_KO_Dicts = label_KO_Dicts(KO_Dicts, labels)

# Intersect KO genes across conditions
intersected_dict = find_all_combinations_intersections(labeled_KO_Dicts)

# Write the output to a file
with open('Output_Dictionary.txt', 'wt') as filehandler:
    filehandler.write(str(intersected_dict))

dictionary_lengths = create_value_set_lengths(intersected_dict)
print(dictionary_lengths)

print('Intersected Dict')
print(intersected_dict)

def find_intersection_of_most_conditions(intersected_dict):
    """
    Find the list of genes that are common in the intersection of the most conditions.

    Parameters:
    intersected_dict (dict): Dictionary containing intersected KO genes.

    Returns:
    list: List of genes that are common in the intersection of the most conditions.
    """
    max_conditions = 0
    common_genes_set = set()

    for subdict_name, subdict in intersected_dict.items():
        # Flatten the genes in the subdict
        current_genes_set = set(gene.id for reactions in subdict.values() for genes in reactions.values() for gene in genes)
        print(f"Processing {subdict_name} with genes: {current_genes_set}")  # Debugging statement
        if len(subdict_name.split(' - ')) > max_conditions:
            max_conditions = len(subdict_name.split(' - '))
            common_genes_set = current_genes_set
        elif len(subdict_name.split(' - ')) == max_conditions:
            common_genes_set &= current_genes_set

    return list(common_genes_set)

common_genes = find_intersection_of_most_conditions(intersected_dict)

print("List of genes in the intersection of the most conditions:", common_genes)

##

from itertools import combinations, product, chain

def create_expanded_deletion_segments(genes):
    n = len(genes)
    segments = {
        "100%": genes,
        "90%A": genes[:int(0.9 * n)],
        "90%B": genes[int(0.1 * n):],
        "80%A": genes[:int(0.8 * n)],
        "80%B": genes[int(0.2 * n):],
        "70%A": genes[:int(0.7 * n)],
        "70%B": genes[int(0.3 * n):],
        "60%A": genes[:int(0.6 * n)],
        "60%B": genes[int(0.4 * n):],
        "50%A": genes[:int(0.5 * n)],
        "50%B": genes[int(0.5 * n):],
        "33%A": genes[:int(0.33 * n)],
        "33%B": genes[int(0.33 * n):int(0.66 * n)],
        "33%C": genes[int(0.66 * n):],
        "25%A": genes[:int(0.25 * n)],
        "25%B": genes[int(0.25 * n):int(0.5 * n)],
        "25%C": genes[int(0.5 * n):int(0.75 * n)],
        "25%D": genes[int(0.75 * n):],
        "12.5%A": genes[:int(0.125 * n)],
        "12.5%B": genes[int(0.125 * n):int(0.25 * n)],
        "12.5%C": genes[int(0.25 * n):int(0.375 * n)],
        "12.5%D": genes[int(0.375 * n):int(0.5 * n)],
        "12.5%E": genes[int(0.5 * n):int(0.625 * n)],
        "12.5%F": genes[int(0.625 * n):int(0.75 * n)],
        "12.5%G": genes[int(0.75 * n):int(0.875 * n)],
        "12.5%H": genes[int(0.875 * n):]
    }
    return segments

def evaluate_segments(model, segments):
    viable_segments = {}
    for segment_name, genes in segments.items():
        viability = knockout_and_check_viability(model, genes)
        print(f"Evaluating segment {segment_name}: {viability}")  # Debug statement
        if viability > min_growth_rate and not np.isnan(viability):
            viable_segments[segment_name] = genes
    return viable_segments

def powerset(iterable):
    """Generates a powerset, all possible unique combinations of a set."""
    s = list(set(iterable))  # Ensure no duplicate elements
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s) + 1))

def outputToLists(components, boardnames, genecodes, color):
    """Uses powerset() to generate combinations and records the names of segments in the combinations,
    matching with the combined deleted genes."""
    successfulcomponents = [component for component in components if component in boardnames]
    powersetcombos_key = []
    powersetcombos = []

    for combo in powerset(successfulcomponents):
        combo_key = ' '.join(combo)
        powersetcombos_key.append(combo_key)

        genesdeleted = []
        for item in combo:
            if item.strip():
                x = boardnames.index(item)
                genesdeleted.extend(genecodes[x])

        powersetcombos.append(genesdeleted)

    return list(zip(powersetcombos_key, powersetcombos))

def third_stage(model, viable_segments):
    # Identify the three largest viable segments
    segment_names = list(viable_segments.keys())
    segment_names.sort(key=lambda x: len(viable_segments[x]), reverse=True)
    top_segments = segment_names[:3]

    combined_results = {}

    # For each of the top three segments, combine with structured combinations of smaller segments
    boardnames = list(viable_segments.keys())
    genecodes = list(viable_segments.values())

    matching_logic = {
        '80%A': ['80%A', '12.5%G', '12.5%H'],
        '80%B': ['80%B', '12.5%A', '12.5%B'],
        '70%A': ['70%A', '12.5%F', '12.5%G', '12.5%H'],
        '70%B': ['70%B', '12.5%A', '12.5%B', '12.5%C'],
        '60%A': ['60%A', '33%C', '12.5%E', '12.5%F', '12.5%G', '12.5%H'],
        '60%B': ['60%B', '33%A', '12.5%A', '12.5%B', '12.5%C', '12.5%D'],
        '50%A': ['50%A', '33%C', '12.5%E', '12.5%F', '12.5%G', '12.5%H'],
        '50%B': ['50%B', '33%A', '12.5%A', '12.5%B', '12.5%C', '12.5%D'],
        '33%A': ['33%A', '33%B', '33%C', '12.5%D', '12.5%E', '12.5%F', '12.5%G', '12.5%H'],
        '33%B': ['33%B', '33%A', '33%C', '12.5%A', '12.5%B', '12.5%F', '12.5%G', '12.5%H'],
        '33%C': ['33%C', '33%A', '33%B', '12.5%A', '12.5%B', '12.5%C', '12.5%D', '12.5%E'],
        '25%A': ['25%A', '12.5%C', '12.5%D', '12.5%E', '12.5%F', '12.5%G', '12.5%H'],
        '25%B': ['25%B', '12.5%A', '12.5%B', '12.5%E', '12.5%F', '12.5%G', '12.5%H'],
        '25%C': ['25%C', '12.5%A', '12.5%B', '12.5%C', '12.5%D', '12.5%G', '12.5%H'],
        '25%D': ['25%D', '12.5%A', '12.5%B', '12.5%C', '12.5%D', '12.5%E', '12.5%F'],
        '12.5%A': ['12.5%A', '12.5%B', '12.5%C', '12.5%D', '12.5%E', '12.5%F', '12.5%G', '12.5%H'],
        '12.5%B': ['12.5%A', '12.5%B', '12.5%C', '12.5%D', '12.5%E', '12.5%F', '12.5%G', '12.5%H'],
        '12.5%C': ['12.5%A', '12.5%B', '12.5%C', '12.5%D', '12.5%E', '12.5%F', '12.5%G', '12.5%H'],
        '12.5%D': ['12.5%A', '12.5%B', '12.5%C', '12.5%D', '12.5%E', '12.5%F', '12.5%G', '12.5%H'],
        '12.5%E': ['12.5%A', '12.5%B', '12.5%C', '12.5%D', '12.5%E', '12.5%F', '12.5%G', '12.5%H'],
        '12.5%F': ['12.5%A', '12.5%B', '12.5%C', '12.5%D', '12.5%E', '12.5%F', '12.5%G', '12.5%H'],
        '12.5%G': ['12.5%A', '12.5%B', '12.5%C', '12.5%D', '12.5%E', '12.5%F', '12.5%G', '12.5%H'],
        '12.5%H': ['12.5%A', '12.5%B', '12.5%C', '12.5%D', '12.5%E', '12.5%F', '12.5%G', '12.5%H']
    }

    for color in top_segments:
        if color in matching_logic:
            components = matching_logic[color]
            combos_and_genes = outputToLists(components, boardnames, genecodes, color)
            for combo_key, combined_genes in combos_and_genes:
                print(f"Evaluating combination {color}_{combo_key} with genes: {combined_genes}")  # Debug statement
                solution =  knockout_and_check_viability(model, combined_genes)
                if solution > min_growth_rate and not np.isnan(solution):
                    combined_results[f"{color}_{combo_key}"] = combined_genes

    return combined_results

def fourth_stage(model, initial_largest_combination, remaining_genes):
    largest_combination = initial_largest_combination

    remaining_genes = list(remaining_genes)  # Convert set to list

    while len(remaining_genes) > 8:
        groups = [remaining_genes[i::8] for i in range(8)]
        found_viable_combination = False
        combined_results = {}

        # Generate powerset for the groups
        for combination in chain.from_iterable(combinations(groups, r) for r in range(1, len(groups) + 1)):
            subset = [gene for group in combination for gene in group]
            combined_genes = list(largest_combination) + subset
            solution = knockout_and_check_viability(model, combined_genes)
            if solution > min_growth_rate and not np.isnan(solution):
                combined_results[f"combination_{combination}"] = combined_genes
                largest_combination = combined_genes
                found_viable_combination = True

        if not found_viable_combination:
            # Try single gene knockouts
            single_gene_results = {}
            for gene in remaining_genes:
                combined_genes = list(largest_combination) + [gene]
                solution = knockout_and_check_viability(model, combined_genes)
                if solution > min_growth_rate and not np.isnan(solution):
                    single_gene_results[f"single_{gene}"] = combined_genes
                    largest_combination = combined_genes
                    found_viable_combination = True

            if not found_viable_combination:
                break  # Exit if no viable single gene knockout is found

            remaining_genes = list(set(remaining_genes) - set(chain(*single_gene_results.values())))

        remaining_genes = list(set(remaining_genes) - set(largest_combination))

    if len(remaining_genes) <= 8:
        combined_results = {}
        groups = [remaining_genes[i::8] for i in range(8)]
        for combination in product([0, 1], repeat=8):
            subset = [gene for i, include in enumerate(combination) for gene in groups[i] if include]
            combined_genes = list(largest_combination) + subset
            solution = knockout_and_check_viability(model, combined_genes)
            if solution > min_growth_rate and not np.isnan(solution):
                combined_results[f"final_combination_{combination}"] = combined_genes
                largest_combination = combined_genes

    return combined_results, largest_combination

# Stage 2: Create and evaluate deletion segments
segments = create_expanded_deletion_segments(common_genes)
print(f"Segments: {segments}")  # Debug statement
viable_segments = evaluate_segments(model, segments)
print(f"Viable segments: {viable_segments}")  # Debug statement

# Stage 3: Process the three largest viable segments
third_stage_results = third_stage(model, viable_segments)
print(f"Third stage results: {third_stage_results}")  # Debug statement

# Find the largest combination from the third stage results
if third_stage_results:
    largest_combination = max(third_stage_results.values(), key=len)
    remaining_genes = set(common_genes) - set(chain(*third_stage_results.values()))

    # Stage 4: Split the remaining genes and evaluate combinations
    fourth_stage_results, final_largest_combination = fourth_stage(model, largest_combination, remaining_genes)
else:
    fourth_stage_results, final_largest_combination = {}, []

print(set(final_largest_combination))
print(len(set(final_largest_combination)))

gene_list = final_largest_combination

##
def get_genes_for_reaction(model, reaction_id):
    """Returns the set of genes associated with a given reaction."""
    try:
        reaction = model.reactions.get_by_id(reaction_id)
        return set(gene.id for gene in reaction.genes)
    except KeyError:
        return set()

def find_metabolites_with_genes_in_pathways(model, pathways_dict, gene_list):
    """Returns a list of metabolites whose pathways include reactions encoded by genes in the gene_list."""
    GL_Met = []

    for metabolite, pathways in pathways_dict.items():
        for path_reactions, path_metabolites in pathways.values():
            for reaction_id in path_reactions:
                reaction_genes = get_genes_for_reaction(model, reaction_id)
                if not reaction_genes.isdisjoint(gene_list):
                    GL_Met.append(metabolite)
                    break  # No need to check further reactions in this pathway if one reaction already matched

    return list(set(GL_Met))  # Remove duplicates

def extract_metabolites_with_genes(model, pathways_dict, gene_list):
    """Main function to extract metabolites whose pathways involve genes in the gene_list."""
    GL_Met = find_metabolites_with_genes_in_pathways(model, pathways_dict, gene_list)
    return GL_Met

GL_Met = extract_metabolites_with_genes(model, pathways, gene_list)
print(GL_Met)

import cobra
from cobra.io import write_sbml_model

from cobra.manipulation.delete import remove_genes

# Remove the specified genes and their associated reactions
remove_genes(model, gene_list, remove_reactions=True)

# Save the modified model
write_sbml_model(model, r'output_path')

print("Modified model saved.")
