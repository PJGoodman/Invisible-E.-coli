import cometspy as c
import cobra
import pandas as pd
import matplotlib.pyplot as plt
import shutup; shutup.please()
from matplotlib import rc

# Load the models
ecoli = c.model(cobra.io.read_sbml_model("ecoli_path"))
ecoli.id = 'ecoli'

bsubtilis = c.model(cobra.io.read_sbml_model(r"bsubtilis_path"))
bsubtilis.id = 'bsubtilis'

# Set initial biomass, 5e-6 gr at coordinate [0,0]
ecoli.initial_pop = [0, 0, 5e-6]
bsubtilis.initial_pop = [0, 0, 5e-6]

# Remove bounds from all exchange reactions for both species
def remove_bounds_from_exchanges(comets_model):
    for reaction_id in comets_model.reactions:
        if reaction_id.startswith('EX_'):  # Identifying exchange reactions by prefix 'EX_'
            comets_model.change_bounds(reaction_id, -1000, 1000)  # Set bounds to be unrestricted

remove_bounds_from_exchanges(ecoli)
remove_bounds_from_exchanges(bsubtilis)

# Create an empty layout
test_tube = c.layout()

# Add the models to the test tube
test_tube.add_model(ecoli)
test_tube.add_model(bsubtilis)

# Remove bounds from all exchange reactions for non-static metabolites for both species
def change_bounds_for_non_static_exchanges(comets_model, static_metabolites):
    for reaction_id in comets_model.reactions:
        if reaction_id.startswith('EX_'):  # Identifying exchange reactions by prefix 'EX_'
            metabolite = reaction_id.split('_')[1]  # Extract the metabolite part
            if metabolite + '_e' not in static_metabolites:  # Check if it's non-static
                comets_model.change_bounds(reaction_id, -1000, 1000)  # Set bounds to be unrestricted


# Add typical trace metabolites and oxygen coli as static
trace_metabolites = ['cbl1_e', 'cobalt2_e', 'cu2_e', 'fe2_e', 'fe3_e', 'h_e', 'h2o_e', 'k_e',
                     'mn2_e', 'mobd_e', 'ni2_e', 'o2_e', 'pi_e', 'sel_e', 'slnt_e', 'tungs_e', 'zn2_e']

change_bounds_for_non_static_exchanges(ecoli, trace_metabolites)
change_bounds_for_non_static_exchanges(bsubtilis, trace_metabolites)

for i in trace_metabolites:
    test_tube.set_specific_metabolite(i, 1000)
    test_tube.set_specific_static(i, 1000)

# Add Media
test_tube.set_specific_metabolite('glc__D_e', 5.6)   # Glucose concentration to track
test_tube.set_specific_metabolite('na1_e', 9.3)
test_tube.set_specific_metabolite('k_e', 2.21)
test_tube.set_specific_metabolite('nh4_e', 0.0009347)  # Ammonium concentration to track
test_tube.set_specific_metabolite('ca2_e', 0.1)
test_tube.set_specific_metabolite('mg2_e', 0.01)
test_tube.set_specific_metabolite('cl_e', 2.92)
test_tube.set_specific_metabolite('so4_e', 0.01)
test_tube.set_specific_metabolite('pi_e', 4.22)
test_tube.set_specific_metabolite('urea_e', 0.0)  # Urea concentration to track
test_tube.set_specific_metabolite('arg__L_e', 0.009347)

# Define parameters and run the simulation
comp_params = c.params()
comp_params.set_param('maxCycles', 70)
comp_params.set_param('writeMediaLog', True)  # Enable media logging to track metabolites

comp_assay = c.comets(test_tube, comp_params)
comp_assay.run()

# Extract biomass data
biomass = comp_assay.total_biomass
biomass['t'] = biomass['cycle'] * comp_assay.parameters.all_params['timeStep']

# Extract glucose concentration data from the media log
media = comp_assay.media.copy()
glucose_concentration = media[media['metabolite'] == 'glc__D_e'][['cycle', 'conc_mmol']]
glucose_concentration['t'] = glucose_concentration['cycle'] * comp_assay.parameters.all_params['timeStep']

# Extract urea concentration data from the media log
urea_concentration = media[media['metabolite'] == 'urea_e'][['cycle', 'conc_mmol']]
urea_concentration['t'] = urea_concentration['cycle'] * comp_assay.parameters.all_params['timeStep']

# Extract NH4 concentration data from the media log
nh4_concentration = media[media['metabolite'] == 'nh4_e'][['cycle', 'conc_mmol']]
nh4_concentration['t'] = nh4_concentration['cycle'] * comp_assay.parameters.all_params['timeStep']

# Extract Arginine concentration data from the media log
arg_concentration = media[media['metabolite'] == 'arg__L_e'][['cycle', 'conc_mmol']]
arg_concentration['t'] = arg_concentration['cycle'] * comp_assay.parameters.all_params['timeStep']

# Plot biomass data
plt.figure()
plt.plot(biomass['t'], biomass['ecoli'], label='E. coli Biomass', color='tab:orange')
plt.plot(biomass['t'], biomass['bsubtilis'], label='B. subtilis Biomass', color='tab:blue')
plt.xlabel('Time (hours)')
plt.ylabel('Biomass (gr.)')
plt.title('Biomass Over Time')
plt.legend()
plt.show()

# Plot glucose concentration data
plt.figure()
plt.plot(glucose_concentration['t'], glucose_concentration['conc_mmol'], label='Glucose Concentration', color='tab:green')
plt.xlabel('Time (hours)')
plt.ylabel('Glucose Concentration (mmol)')
plt.title('Glucose Concentration Over Time')
plt.legend()
plt.show()

# Plot urea concentration data
plt.figure()
plt.plot(urea_concentration['t'], urea_concentration['conc_mmol'], label='Urea Concentration', color='tab:red')
plt.xlabel('Time (hours)')
plt.ylabel('Urea Concentration (mmol)')
plt.title('Urea Concentration Over Time')
plt.legend()
plt.show()

# Plot NH4 concentration data
plt.figure()
plt.plot(nh4_concentration['t'], nh4_concentration['conc_mmol'], label='Ammonium Concentration', color='tab:purple')
plt.xlabel('Time (hours)')
plt.ylabel('Ammonium Concentration (mmol)')
plt.title('Ammonium (NH4) Concentration Over Time')
plt.legend()
plt.show()

# Plot Arginine concentration data
plt.figure()
plt.plot(arg_concentration['t'], arg_concentration['conc_mmol'], label='Arginine Concentration', color='tab:purple')
plt.xlabel('Time (hours)')
plt.ylabel('Arginine Concentration (mmol)')
plt.title('Arginine Concentration Over Time')
plt.legend()
plt.show()

# Calculate cfit correctly
final_bsubtilis = biomass.loc[biomass['t'] == 50, 'bsubtilis'].iloc[0]
print(f'Final B. subtilis: {final_bsubtilis}')
final_ecoli = biomass.loc[biomass['t'] == 50, 'ecoli'].iloc[0]
print(f'Final E. coli: {final_ecoli}')

print(f'Ratio: {final_bsubtilis / (final_ecoli + final_bsubtilis)}')
