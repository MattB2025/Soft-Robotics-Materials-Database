# -*- coding: utf-8 -*-
"""
Hyperelastic Material Constitutive Law Curve Fitting
Replicates core functions from the Soft Robotics Materials Database app
without web interface. Fits hyperelastic models to uniaxial tensile data from CSV.
"""

import pandas as pd
import numpy as np
from scipy.optimize import least_squares, minimize, NonlinearConstraint, LinearConstraint
from Hyperelastic import Hyperelastic
from HyperelasticStats import HyperelasticStats
import sys
import os

def read_csv_exp_data_files(file_path):
    """Read experimental data from CSV file."""
    try:
        # Read header information
        header = pd.read_csv(file_path, delimiter=';', usecols=["PARAMETER", "INFO", "URL"]).head(15)
        # Read experimental data (skipping header rows)
        data = pd.read_csv(file_path, delimiter=';', skiprows=18,
                          names=['Time (s)', 'True Strain', 'True Stress (MPa)',
                                'Engineering Strain', 'Engineering Stress (MPa)'])
        return data, header
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        sys.exit(1)

def objectiveFun_Callback(parameters, exp_strain, exp_stress, hyperelastic):
    """Cost function for optimization."""
    theo_stress = hyperelastic.ConsitutiveModel(parameters, exp_strain)
    if hyperelastic.fitting_method == 'lm':
        residuals = theo_stress - exp_stress
    elif hyperelastic.fitting_method == 'trust-constr':
        residuals = np.sqrt(sum((theo_stress - exp_stress)**2.0))
    else:
        raise ValueError("Error: please choose either 'lm' or 'trust-constr' as fitting method")
    return residuals

def optimization(model, order, dataframe, data_type, strain_range=None):
    """Perform optimization to fit constitutive model to experimental data."""
    hyperelastic = Hyperelastic(model, np.array([0]), order, data_type)

    # Select data within strain range if specified
    if strain_range:
        min_strain, max_strain = strain_range
        mask = (dataframe[data_type + ' Strain'] >= min_strain) & (dataframe[data_type + ' Strain'] <= max_strain)
        dataframe = dataframe[mask]

    exp_strain = dataframe[data_type + ' Strain'].values
    exp_stress = dataframe[data_type + ' Stress (MPa)'].values

    if hyperelastic.fitting_method == 'trust-constr':
        if hyperelastic.model == 'Ogden':
            const = NonlinearConstraint(hyperelastic.NonlinearConstraintFunction, 0.0, np.inf,
                                      jac=hyperelastic.NonlinearConstraintJacobian, hess='2-point')
        elif hyperelastic.model == 'Mooney Rivlin':
            const = LinearConstraint([[1.0, 1.0, 0.0][:hyperelastic.order],
                                    [0.0, 0.0, 0.0][:hyperelastic.order]], 0.0, np.inf)
        else:
            const = ()

        optim_result = minimize(objectiveFun_Callback, hyperelastic.initialGuessParam,
                              args=(exp_strain, exp_stress, hyperelastic),
                              method='trust-constr', constraints=const, tol=1e-12)
    elif hyperelastic.fitting_method == 'lm':
        optim_result = least_squares(objectiveFun_Callback, hyperelastic.initialGuessParam,
                                   method='lm', gtol=1e-12,
                                   args=(exp_strain, exp_stress, hyperelastic))
    else:
        raise ValueError("Error in fitting method")

    optim_parameters = optim_result.x

    df_model_param = pd.DataFrame(optim_parameters, index=hyperelastic.param_names,
                                columns=[model]).transpose()

    theo_stress = hyperelastic.ConsitutiveModel(optim_parameters, exp_strain)
    data_model = pd.DataFrame({data_type + ' Strain': exp_strain,
                             data_type + ' Stress (MPa)': theo_stress})

    stats = HyperelasticStats(exp_stress, theo_stress, hyperelastic.nbparam)
    aic = stats.aic()
    r_squared = stats.r_squared()

    return df_model_param, data_model, aic, r_squared

def main():
    """Main function to run curve fitting."""
    if len(sys.argv) < 2:
        print("Usage: python hyperelastic_fitting.py <csv_file> [model] [order] [data_type] [min_strain] [max_strain]")
        print("  csv_file: Path to CSV file with experimental data")
        print("  model: Constitutive model (default: Mooney Rivlin)")
        print("         Options: Mooney Rivlin, Ogden, Neo Hookean, Yeoh, Gent, Veronda Westmann, Humphrey")
        print("  order: Model order (default: 3)")
        print("  data_type: True or Engineering (default: True)")
        print("  min_strain, max_strain: Optional strain range for fitting (default: full range)")
        sys.exit(1)

    csv_file = sys.argv[1]
    model = sys.argv[2] if len(sys.argv) > 2 else 'Mooney Rivlin'
    order = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    data_type = sys.argv[4] if len(sys.argv) > 4 else 'True'

    strain_range = None
    if len(sys.argv) > 6:
        min_strain = float(sys.argv[5])
        max_strain = float(sys.argv[6])
        strain_range = (min_strain, max_strain)

    # Validate inputs
    valid_models = ['Mooney Rivlin', 'Ogden', 'Neo Hookean', 'Yeoh', 'Gent', 'Veronda Westmann', 'Humphrey']
    if model not in valid_models:
        print(f"Error: Invalid model '{model}'. Valid options: {', '.join(valid_models)}")
        sys.exit(1)

    if data_type not in ['True', 'Engineering']:
        print("Error: data_type must be 'True' or 'Engineering'")
        sys.exit(1)

    if not os.path.exists(csv_file):
        print(f"Error: CSV file '{csv_file}' not found")
        sys.exit(1)

    # Read data
    print(f"Reading data from {csv_file}...")
    data, header = read_csv_exp_data_files(csv_file)

    # Perform optimization
    print(f"Fitting {model} model (order {order}) to {data_type} data...")
    if strain_range:
        print(f"Using strain range: {strain_range[0]} to {strain_range[1]}")

    df_param, model_data, aic, r_squared = optimization(model, order, data, data_type, strain_range)

    # Output results
    print("\nFitted Parameters:")
    print(df_param.to_string())
    print(".4f")
    print(".4f")

    # Save fitted data
    output_file = os.path.splitext(csv_file)[0] + f"_{model.replace(' ', '_')}_fitted.csv"
    model_data.to_csv(output_file, index=False)
    print(f"\nFitted data saved to: {output_file}")

    # Optional: Save comparison plot data
    comparison_data = pd.DataFrame({
        'Strain': data[data_type + ' Strain'],
        'Experimental_Stress': data[data_type + ' Stress (MPa)'],
        'Model_Stress': np.interp(data[data_type + ' Strain'], model_data[data_type + ' Strain'], model_data[data_type + ' Stress (MPa)'])
    })
    comparison_file = os.path.splitext(csv_file)[0] + f"_{model.replace(' ', '_')}_comparison.csv"
    comparison_data.to_csv(comparison_file, index=False)
    print(f"Comparison data saved to: {comparison_file}")

if __name__ == '__main__':
    main()