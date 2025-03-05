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
