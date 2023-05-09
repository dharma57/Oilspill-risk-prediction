# Simulation Code based on EPA-540-R-070-002
# Risk Assessment Guidance for Superfund
# Volume 1 : Human Health Evaluation Manual
# Part F, Supplemental Guidance for Inhalation Risk Assessment
# In this project we looking to build a quantitative risk assessment model for east palestine train derialment accident risk assessment 

import numpy as np
import csv
import math

# Gaussian risk assessment
def gaussian_plume_model(Q, u, sigma_y, sigma_z, H, x, y, z):
    #Q = mass_of_chemical_spilled / release_duration_minutes  # g/min
    factor = Q / (2 * math.pi * u * sigma_y * sigma_z)
    exponential_y = math.exp(-0.5 * (y / sigma_y)**2)
    exponential_z1 = math.exp(-0.5 * ((z - H) / sigma_z)**2)
    exponential_z2 = math.exp(-0.5 * ((z + H) / sigma_z)**2)

    concentration = factor * exponential_y * (exponential_z1 + exponential_z2)
    return concentration

# Linear Evaporation 
def linear_evaporation_model(initial_mass, burn_rate, density, alpha):
     # Convert burn rate from mm/min to g/min
    burn_rate_g_min = (burn_rate / 1000) * density *alpha*60

    # Calculate the mass reduction due to burning
    mass_reduction = burn_rate_g_min 

    # Ensure the remaining mass doesn't go below zero
    remaining_mass = initial_mass - mass_reduction
    if remaining_mass < 0:
        remaining_mass = 0

    return remaining_mass

# Chemical Concentration
def calculate_ca(Q, x, y, z,u=1, H=0):
    # Constants and conversions


    # Pasquill-Gifford method for stability class D (estimates)
    sigma_y = 0.08 * x**(0.894)  # m
    sigma_z = 0.056 * x**(0.894)  # m
    #Q = mass_of_chemical_spilled / release_duration_minutes  # g/m
        
    # Calculate CA using the Gaussian plume model
    # Constants and conversions
    CA = gaussian_plume_model(Q, u, sigma_y, sigma_z, H, x, y, z)
    return CA

# Exposure Concentration
def exposure_concentration(CA, ET, EF, ED, AT):
    # Estimating Exposure Concentrations for Assessing Cancer Risks
    if ED <= 1 / 365:
        # Acute Exposures
        EC = CA
        print("exposure_type == acute")
    # Chronic or Subchronic Exposures    
    elif 30 / 365 < ED <= 0.1 * AT / (365 * 24):
        EC = (CA * ET * EF * ED) / AT
        print("exposure_type == sub-chronic")
    else:
        EC = (CA * ET * EF * ED) / AT
        print("exposure_type == chronic")
    return EC

# Calculate BMI for Population
def get_body_weight(age):
    # The relationship between age and body weight here.
    # For simplicity,a linear model with any age-weight relationship that you prefer.
    Bw = 40 + 0.5 * age
    return Bw

# Inhalation Exposure
def inhalation_exposure_dose(EC, IR, Bw, ED, ET, EF, AT):
    # RAGS part A Equation Describing the Estimation of Inhalation Exposure
    D = EC * (IR/Bw) * ET * EF * ED / AT
    return D

# Inhalation Unit Risk
def inhalation_unit_risk(LEC10_HEC):
# Deviavtion of the Inhalation Unit Risk
    IUR = 0.1 / LEC10_HEC
    return IUR

# Cancer Risk
def cancer_risk(IUR, EC):
    # Cancer Risks Characterized an Inhalation Unit Risk
    cancer_risk = IUR * EC
    return cancer_risk

# Hazard Quotient 
def hazard_quotient(EC, TV):
# HQ for inhalation pathway
    HQ = EC / (TV * 1000)
    return HQ

# Constants
'''
LEC10_HEC -  For vinyl chloride, the IUR is 1.0 x 10^(-4) (µg/m³)^(-1). 
'''
LEC10_HEC = 10e3  # Example value
TV = 1  # Toxicity Value (mg/m3) = Inhalation toxicity value
mass_of_chemical_spilled = 1109400  # grams
burn_rate = 4.3  # mm/min
density = 2150  # g/m^3
gallons_to_cubic_meters = 0.00378541 * 5 * 25800
x, y, z = 10, 0.5, 0.5  # coordinates of the point where we want to calculate CA
release_duration_days = 5 # Release time is assumed from Feb 3 - Feb 8
release_duration_minutes = release_duration_days * 24 * 60


# Create the time array
time_array = [i for i in range(0, release_duration_minutes + 1)]  # 0 to 7 days, incrementing over 1 minute

# East palestine population data index 
# Data collected from https://namecensus.com/demographics/ohio/columbiana-county/east-palestine/
population_data = [
    {"age_group": "0 to 9 Years", "avg_age": 4.5, "male": 474, "female": 397},
    {"age_group": "10 to 19 Years", "avg_age": 14.5, "male": 174, "female": 288},
    {"age_group": "20 to 29 Years", "avg_age": 24.5, "male": 234, "female": 268},
    {"age_group": "30 to 39 Years", "avg_age": 34.5, "male": 264, "female": 281},
    {"age_group": "40 to 49 Years", "avg_age": 44.5, "male": 259, "female": 337},
    {"age_group": "50 to 59 Years", "avg_age": 54.5, "male": 231, "female": 192},
    {"age_group": "60 to 69 Years", "avg_age": 64.5, "male": 347, "female": 236},
    {"age_group": "70 Years and over", "avg_age": 75, "male": 299, "female": 458},
]


# Define constants
# Define age groups, inhalation rates (m^3/day) and exposure durations (years)
age_groups = [
    {"name": "0-9 years", "IR": 10* 24 * 60, "ED": 9* 5},
    {"name": "10-19 years", "IR": 15* 24 * 60, "ED": 10* 5},
    {"name": "20-29 years", "IR": 20* 24 * 60, "ED": 10* 5},
    {"name": "30-39 years", "IR": 20* 24 * 60, "ED": 10* 5},
    {"name": "40-49 years", "IR": 20* 24 * 60, "ED": 10* 5},
    {"name": "50-59 years", "IR": 20* 24 * 60, "ED": 10* 5},
    {"name": "60-69 years", "IR": 18* 24 * 60, "ED": 10* 5},
    {"name": "70 years and over", "IR": 15* 24 * 60, "ED": 10* 5},
]
ET = 1 * 60  # Exposure time in hours/day, converted to minutes/day
EF = 1 * 365  # Exposure frequency in days/year, converted to days/minute
'''
Controlled burn rate was implemented by  first responder  
Assumed some rate alpha to control burn
'''
alpha = 0.25
mass_remaining = mass_of_chemical_spilled 
results = []

for elapsed_time in time_array:
    # Update the mass of the chemical
    mass_remaining = linear_evaporation_model(mass_remaining , burn_rate, density, alpha)

    # Update the Q value for the new mass_remaining
    Q = mass_remaining / release_duration_minutes

    CA = calculate_ca(Q , x, y, z, u=1, H=0)
                
    for age_group_data, age_group in zip(population_data, age_groups):

            
        for gender in ["male", "female"]:
            Bw = get_body_weight(age_group_data["avg_age"])
            population = age_group_data[gender]
            IR = age_group["IR"]
            ED = age_group["ED"]
            AT = ED  # Averaging time in days (since ED is already in days)

         
            EC = exposure_concentration(CA, ET, EF, ED, AT)
            D = inhalation_exposure_dose(EC, IR, Bw, ED, ET, EF, AT)
            IUR = inhalation_unit_risk(LEC10_HEC)
            CR = cancer_risk(IUR, EC)
            HQ = hazard_quotient(EC, TV)

            results.append([elapsed_time, age_group["name"], gender, population, Bw,mass_remaining,CA, EC, D, IUR, CR, HQ])



# Save the results to a text file
with open('results.txt', 'w') as f:
    f.write('Elapsed Time (min), Age Group, Gender, Population, Body Weight (Bw),mass_remaining (g), Contaminant Concentration CA(g/m^3), Exposure Concentration (EC) (µg/m3), Intake (g/kg-d), Inhalation Unit Risk (IUR) (µg/m3)^-1), Cancer Risk CR, Hazard Quotient (HQ)\n')
    for row in results:
        f.write(', '.join([str(value) for value in row]) + '\n')

# Save the results to a CSV file
with open('results2.csv', 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(['Elapsed Time (min)', 'Age Group', 'Gender', 'Population', 'Body Weight (Bw)','mass_remaining (g)', 'Contaminant Concentration CA(mg/m^3)', 'Exposure Concentration (EC) (µg/m3)', ' Intake (g/kg-d)', 'Inhalation Unit Risk (IUR) (µg/m3)^-1)', 'Cancer CR', 'Hazard Quotient (HQ)'])
    for row in results:
        csvwriter.writerow(row)