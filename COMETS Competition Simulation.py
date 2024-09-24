# Start by loading required packages, including the COMETS toolbox
import cometspy as c
import cobra
import pandas as pd
import matplotlib.pyplot as plt
import scipy.io

ecoli = c.model(cobra.io.read_sbml_model('path_to_file'))
ecoli.id = 'ecoli'

bsubtilis = c.model(cobra.io.read_sbml_model('path_to_file')
bsubtilis.id = 'bsubtilis'


# set its initial biomass, 5e-6 gr at coordinate [0,0]
ecoli.initial_pop = [0, 0, 5e-6]
bsubtilis.initial_pop = [0, 0, 5e-6]

# create an empty layout
test_tube = c.layout()

# add the models to the test tube
test_tube.add_model(ecoli)
test_tube.add_model(bsubtilis)

# Add glucose to the media 
test_tube.set_specific_metabolite('glc__D_e', 0.01)

# Add typical trace metabolites and oxygen coli as static
trace_metabolites = ['ca2_e', 'cbl1_e', 'cl_e', 'co2_e', 'cobalt2_e', 'cu2_e', 'fe2_e', 'fe3_e', 'h_e', 'h2o_e', 'k_e',
                     'mg2_e', 'mn2_e', 'mobd_e', 'na1_e', 'nh4_e', 'ni2_e', 'o2_e', 'pi_e', 'sel_e', 'slnt_e', 'so4_e', 'tungs_e', 'zn2_e']

for i in trace_metabolites:
    test_tube.set_specific_metabolite(i, 1000)
    test_tube.set_specific_static(i, 1000)


comp_params = c.params()
comp_params.set_param('maxCycles', 500)

comp_assay = c.comets(test_tube, comp_params)
comp_assay.run()

biomass = comp_assay.total_biomass
biomass['t'] = biomass['cycle'] * comp_assay.parameters.all_params['timeStep']

#Alternatively use dashed lines, useful if overlaying different plots.

#plt.plot(biomass['t'], biomass['ecoli'], linestyle='--', color='orange')  # Set E. coli to orange
#plt.plot(biomass['t'], biomass['bsubtilis'], linestyle='--', color='blue')  # Set B. subtilis to blue

plt.plot(biomass['t'], biomass['ecoli'], color='orange')  # Set E. coli to orange
plt.plot(biomass['t'], biomass['bsubtilis'], color='blue')  # Set B. subtilis to blue



plt.ylim(0, 0.0006)


# Add labels to the axes
plt.ylabel("Biomass (gr.)")
plt.xlabel("Time (t)")

# Save the plot with a transparent background
plt.savefig('path_to_output', transparent=True, dpi = 700)

# Show the plot without the legend or labels for individual lines
plt.show()


# Calculate cfit correctly
initial_bsubtilis = biomass.loc[biomass['t'] == 0, 'bsubtilis'].iloc[0]
final_bsubtilis = biomass.loc[biomass['t'] == 50, 'bsubtilis'].iloc[0]
print(f'Final BSub {final_bsubtilis}')
initial_ecoli = biomass.loc[biomass['t'] == 0, 'ecoli'].iloc[0]
final_ecoli = biomass.loc[biomass['t'] == 50, 'ecoli'].iloc[0]
print(f'Final Ecoli {final_ecoli}')

print(f'Ratio: {final_bsubtilis / (final_ecoli + final_bsubtilis)}')

# Avoiding division by zero and ensuring proper calculation
if initial_bsubtilis != 0 and initial_ecoli != 0:
    bsubtilis_ratio = final_bsubtilis / initial_bsubtilis
    ecoli_ratio = final_ecoli / initial_ecoli
    cfit = final_bsubtilis / initial_ecoli
    print(f'cfit: {cfit}')
else:
    print("Error: Initial biomass values should not be zero.")



