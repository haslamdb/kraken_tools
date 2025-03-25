"""
This script demonstrates how to perform a zero-inflated GLMM in Python.

In R, one might use glmmTMB or NBZIMM for zero-inflated negative binomial or Poisson.
In Python, we can use statsmodels' ZeroInflatedPoisson or ZeroInflatedNegativeBinomialP.

We'll run it for each species or do a simple example with one species.

Note: This is quite advanced and might require careful model fitting.

We assume the input has row=samples, columns=metadata + species counts.
We'll do something like:
   ZeroInflatedPoisson.from_formula(
       "Count ~ group + other_metadata",
       data=..., p="logit"
   )

Reference: https://www.statsmodels.org/stable/generated/statsmodels.discrete.count_model.ZeroInflatedPoisson.html

"""
import argparse
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from statsmodels.discrete.count_model import ZeroInflatedPoisson, ZeroInflatedNegativeBinomialP
import os
import logging
import traceback
from kraken_tools.utils.file_utils import check_file_exists_with_logger

def run_glmm_analysis(
    abundance_file,
    sample_key,
    output_dir,
    formula,
    model="negbin",
    group_col="Group",
    min_abundance=0.01,
    min_prevalence=0.1,
    logger=None
):
    """
    Run Generalized Linear Mixed Model analysis on taxonomic data.
    
    Args:
        abundance_file: Path to abundance file
        sample_key: Path to sample metadata file
        output_dir: Directory to save output files
        formula: R-style formula for GLMM (e.g., 'Count ~ Group + (1|Subject)')
        model: Model family for GLMM ("poisson" or "negbin")
        group_col: Column name for grouping
        min_abundance: Minimum relative abundance threshold for inclusion
        min_prevalence: Minimum prevalence threshold for inclusion
        logger: Logger instance
        
    Returns:
        Boolean success flag
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    # Check if files exist
    if not check_file_exists_with_logger(abundance_file, "Abundance file", logger):
        return False
    
    if not check_file_exists_with_logger(sample_key, "Sample key", logger):
        return False
    
    try:
        # Make output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Read data
        abundance_df = pd.read_csv(abundance_file, sep="\t", index_col=0)
        metadata_df = pd.read_csv(sample_key, index_col=None)
        
        # Get sample ID column
        sample_id_col = None
        common_id_names = ["SampleName", "Sample", "SampleID", "Sample_ID", "sample_name", "sample_id"]
        for col in common_id_names:
            if col in metadata_df.columns:
                sample_id_col = col
                break
        
        if sample_id_col is None:
            logger.error("Could not find sample ID column in metadata")
            return False
        
        # Filter low-abundance taxa if needed
        if min_abundance > 0 or min_prevalence > 0:
            before_count = abundance_df.shape[0]
            
            # Calculate mean abundance per taxon
            mean_abundance = abundance_df.mean(axis=1)
            
            # Calculate prevalence (proportion of samples where taxon is present)
            prevalence = (abundance_df > 0).mean(axis=1)
            
            # Apply filters
            abundance_df = abundance_df[(mean_abundance >= min_abundance) & 
                                      (prevalence >= min_prevalence)]
            
            after_count = abundance_df.shape[0]
            logger.info(f"Filtered from {before_count} to {after_count} taxa based on abundance/prevalence thresholds")
        
        # Run GLMM analysis for top taxa
        top_n = 20  # Analyze top 20 taxa by abundance
        top_taxa = abundance_df.mean(axis=1).sort_values(ascending=False).head(top_n).index
        
        # Process each taxon
        success_count = 0
        for taxon in top_taxa:
            try:
                # Create long-format data for modeling
                long_data = []
                for sample in abundance_df.columns:
                    if sample in metadata_df[sample_id_col].values:
                        sample_data = metadata_df[metadata_df[sample_id_col] == sample].iloc[0].to_dict()
                        sample_data['Count'] = abundance_df.loc[taxon, sample]
                        sample_data['Taxon'] = taxon
                        long_data.append(sample_data)
                
                long_df = pd.DataFrame(long_data)
                
                # Run the model using the function from the example script
                logger.info(f"Running GLMM for taxon: {taxon}")
                
                # Call the main function from the glmm_analysis.py example script
                # You can import and use main() from the current script
                from statsmodels.discrete.count_model import ZeroInflatedPoisson, ZeroInflatedNegativeBinomialP
                
                # Create the model
                if model == "poisson":
                    model_class = ZeroInflatedPoisson
                else:
                    model_class = ZeroInflatedNegativeBinomialP
                
                # Fit the model
                glmm_model = model_class.from_formula(formula, data=long_df, exog_infl='1')
                glmm_result = glmm_model.fit(method='bfgs', maxiter=200, disp=False)
                
                # Save results
                output_file = os.path.join(output_dir, f"glmm_{taxon.replace('.', '_')}.txt")
                with open(output_file, 'w') as f:
                    f.write(str(glmm_result.summary()))
                    
                success_count += 1
                logger.info(f"GLMM analysis for {taxon} saved to {output_file}")
                
            except Exception as e:
                logger.error(f"Error in GLMM analysis for {taxon}: {str(e)}")
        
        logger.info(f"GLMM analysis completed for {success_count}/{len(top_taxa)} taxa")
        return success_count > 0
        
    except Exception as e:
        logger.error(f"Error in GLMM analysis: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False


def main():
    parser = argparse.ArgumentParser(description="Zero-inflated GLMM example.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--metadata_cols", type=int, default=5)
    parser.add_argument("--species", default=None,
                        help="If specified, run the zero-inflated model for just this species.")
    parser.add_argument("--model", default="poisson", choices=["poisson", "negbin"],
                        help="Which zero-inflated family to use.")
    parser.add_argument("--formula", default="Count ~ Location + SampleCollectionWeek + PostNatalAntibiotics",
                        help="Right-hand side of model formula.")
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    meta = df.iloc[:, :args.metadata_cols]
    counts = df.iloc[:, args.metadata_cols:].copy()

    # We'll pick one species or do a quick demonstration.

    if args.species is None:
        # pick the first species column for demonstration.
        sp = counts.columns[0]
    else:
        sp = args.species
        if sp not in counts.columns:
            raise ValueError(f"Species {sp} not found in columns.")

    # We'll create a data frame with the columns needed for statsmodels.
    # formula might look like: 'Count ~ group + var1 + var2'
    # We'll do: data['Count'] = counts[sp].values.

    data_for_model = meta.copy()
    data_for_model['Count'] = counts[sp].values

    # We'll do a zero-inflated model.
    # For random effects (GLMM), statsmodels doesn't fully support random effects with zero-inflated.
    # We might do a simpler approach or use another library.

    if args.model == "poisson":
        model_class = ZeroInflatedPoisson
    else:
        model_class = ZeroInflatedNegativeBinomialP

    # We'll build the formula string.
    formula_str = "Count ~ " + args.formula.split("~")[1]

    print(f"Fitting zero-inflated {args.model} for species={sp} with formula={formula_str}")

    # We'll do from_formula, specifying exog_infl='1' means we do logistic model for inflation.
    # Or we can specify something else for zero inflation part.

    # For advanced usage, we can pass offset or exposure.

    model = model_class.from_formula(formula_str,
                                     exog_infl='1',
                                     data=data_for_model)

    result = model.fit(method='bfgs', maxiter=200, disp=False)

    print(result.summary())

    # If the user wanted to do it for multiple species, they'd loop over columns in 'counts'.

if __name__ == "__main__":
    main()
