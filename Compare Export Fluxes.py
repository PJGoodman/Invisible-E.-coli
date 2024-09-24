import cobra

# Load your models (assuming you have the models saved as .mat files)
model1 = cobra.io.read_sbml_model('model1')

model2 = cobra.io.read_sbml_model('model2')
# Perform FBA on both models
solution1 = model1.optimize()
solution2 = model2.optimize()

# Extract exchange reactions
exchange_reactions_model1 = [rxn for rxn in model1.reactions if rxn.id.startswith('EX_')]
exchange_reactions_model2 = [rxn for rxn in model2.reactions if rxn.id.startswith('EX_')]

# Create a dictionary to store fluxes for comparison
flux_comparison = {}

for rxn in exchange_reactions_model1:
    if rxn.id in [r.id for r in exchange_reactions_model2]:
        flux_model1 = solution1.fluxes[rxn.id]
        flux_model2 = solution2.fluxes[rxn.id]
        # Only consider reactions where the cell is exporting (flux > 0)
        if flux_model1 > 0 or flux_model2 > 0:
            # Only add to the comparison if the fluxes are different
            if round(flux_model1, 3) != round(flux_model2, 3):
                flux_comparison[rxn.id] = {
                    'Model 1': round(flux_model1, 3), 
                    'Model 2': round(flux_model2, 3)
                }

# Print the comparison results
print("Comparison of Export Exchange Reaction Fluxes Between Model 1 and Model 2 (where fluxes differ):")
print(f"{'Reaction ID':<20}{'Model 1 Flux':<15}{'Model 2 Flux'}")
for rxn_id, fluxes in flux_comparison.items():
    print(f"{rxn_id:<20}{fluxes['Model 1']:<15}{fluxes['Model 2']}")

# Optionally, save the comparison to a file
import pandas as pd
df = pd.DataFrame.from_dict(flux_comparison, orient='index')
df.to_csv('export_flux_comparison_truncated.csv')
