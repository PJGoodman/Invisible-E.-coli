from selenium import webdriver
from selenium.webdriver.firefox.service import Service
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.common.by import By
import cobra

# Load your model
model = cobra.io.load_matlab_model("path_to_file")

# Extract ModelSEED compound IDs from your model
modelseed_ids = [str(met.id[0:-3]) for met in model.metabolites]  # Adjusting to exclude compartment information

# Path to your GeckoDriver
gecko_driver_path = "path_to_GeckoDriver"

# Optional: Specify the path to the Firefox binary if it's not in the default location
firefox_binary_path = "path_to_binary"

# Set up Firefox options (e.g., headless mode)
options = Options()
options.binary_location = firefox_binary_path
options.headless = True  # Run in headless mode (no GUI)

# Create a Service object
service = Service(executable_path=gecko_driver_path)

# Initialize the WebDriver with the Service object and options
driver = webdriver.Firefox(service=service, options=options)

# Base URL for ModelSEED compound pages
base_url = "https://modelseed.org/biochem/compounds/"

# Initialize a list to store the mapping
mapping = []

for compound_id in modelseed_ids:
    url = base_url + compound_id
    driver.get(url)

    try:
        # Wait for the dynamic content to load
        driver.implicitly_wait(3)

        # Look for the BiGG alias link
        bigg_links = driver.find_elements(By.XPATH, "//a[contains(@href, 'bigg.ucsd.edu/universal/metabolites/')]")
        
        bigg_alias = 'N/A'
        if bigg_links:
            bigg_alias = bigg_links[0].text.strip()

        # Store the mapping
        mapping.append({
            'ModelSEED ID': compound_id,
            'BiGG ID': bigg_alias
        })
    
    except Exception as e:
        print(f"Failed to retrieve data for {compound_id}: {str(e)}")

# Close the browser
driver.quit()

# Convert the mapping to a DataFrame
import pandas as pd
df = pd.DataFrame(mapping)

# Save to an Excel file
df.to_excel("modelseed_to_bigg_mapping.xlsx", index=False)

print("Mapping completed and saved to modelseed_to_bigg_mapping.xlsx")
