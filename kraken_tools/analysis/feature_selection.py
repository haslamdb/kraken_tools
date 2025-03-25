# kraken_tools/analysis/feature_selection.py
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
import scipy.spatial.distance as ssd
import logging
import traceback

from kraken_tools.utils.file_utils import check_file_exists_with_logger
from kraken_tools.analysis.permanova import transform_abundance_data

def encode_categorical_variables(metadata_df, categorical_features, logger=None):
    """
    Encode categorical variables for Random Forest
    
    Args:
        metadata_df: DataFrame with metadata
        categorical_features: List of categorical variables to encode
        logger: Logger instance
        
    Returns:
        DataFrame with encoded variables
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    encoded_df = metadata_df.copy()
    
    for col in categorical_features:
        if col in encoded_df.columns:
            # Fill NAs
            encoded_df[col] = encoded_df[col].fillna('missing')
            
            try:
                # Use label encoding for Random Forest
                encoded_df[col] = LabelEncoder().fit_transform(encoded_df[col])
                logger.info(f"Encoded categorical variable: {col}")
            except Exception as e:
                logger.warning(f"Error encoding {col}, dropping column: {str(e)}")
                encoded_df = encoded_df.drop(columns=[col])
    
    return encoded_df

def scale_numerical_variables(metadata_df, numerical_features, logger=None):
    """
    Scale numerical variables for Random Forest
    
    Args:
        metadata_df: DataFrame with metadata
        numerical_features: List of numerical variables to scale
        logger: Logger instance
        
    Returns:
        DataFrame with scaled variables
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    scaled_df = metadata_df.copy()
    
    # Create scaler
    scaler = StandardScaler()
    
    # Get only numeric columns that exist in the DataFrame
    numeric_cols = [col for col in numerical_features if col in scaled_df.columns 
                   and np.issubdtype(scaled_df[col].dtype, np.number)]
    
    if numeric_cols:
        # Scale numeric columns
        scaled_df[numeric_cols] = scaler.fit_transform(scaled_df[numeric_cols])
        logger.info(f"Scaled numerical variables: {numeric_cols}")
    
    return scaled_df

def create_pairwise_differences(X_scaled, logger=None):
    """
    Create matrix of pairwise differences between samples
    
    Args:
        X_scaled: Scaled/encoded metadata DataFrame
        logger: Logger instance
        
    Returns:
        Tuple of (pairwise differences array, indices)
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    logger.info("Creating pairwise differences matrix")
    
    # Convert to numpy array
    X_array = X_scaled.values
    
    # Create pairwise differences
    X_diff = np.abs(X_array[:, None, :] - X_array[None, :, :])
    
    # Get upper triangle indices
    n_samples = X_array.shape[0]
    indices = np.triu_indices(n_samples, k=1)
    
    # Extract upper triangle (pairwise differences)
    X_diff_upper = X_diff[indices]
    
    logger.info(f"Created {X_diff_upper.shape[0]} pairwise differences with {X_diff_upper.shape[1]} features")
    
    return X_diff_upper, indices

def calculate_microbiome_distances(abundance_df, distance_metric="bray", logger=None):
    """
    Calculate pairwise distances between microbiome samples
    
    Args:
        abundance_df: DataFrame with abundance data
        distance_metric: Distance metric to use
        logger: Logger instance
        
    Returns:
        Array of pairwise distances
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    logger.info(f"Calculating microbiome distances using {distance_metric} metric")
    
    try:
        if distance_metric.lower() == "bray":
            distances = ssd.pdist(abundance_df.values, metric="braycurtis")
        elif distance_metric.lower() == "jaccard":
            distances = ssd.pdist(abundance_df.values, metric="jaccard")
        elif distance_metric.lower() == "euclidean":
            distances = ssd.pdist(abundance_df.values, metric="euclidean")
        else:
            logger.warning(f"Unknown distance metric: {distance_metric}. Using Bray-Curtis.")
            distances = ssd.pdist(abundance_df.values, metric="braycurtis")
        
        logger.info(f"Calculated {len(distances)} pairwise distances")
        return distances
    
    except Exception as e:
        logger.error(f"Error calculating distances: {str(e)}")
        logger.error(traceback.format_exc())
        return None

def train_random_forest(X_diff, y_distances, feature_names, test_size=0.2, 
                       n_estimators=100, max_features="sqrt", random_state=42, logger=None):
    """
    Train Random Forest model to predict microbiome distances from metadata differences
    
    Args:
        X_diff: Pairwise differences in metadata
        y_distances: Pairwise microbiome distances
        feature_names: Names of features
        test_size: Proportion of data for testing
        n_estimators: Number of trees in the forest
        max_features: Max features to consider at each split
        random_state: Random seed for reproducibility
        logger: Logger instance
        
    Returns:
        Dictionary with model, performance metrics, and feature importance
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X_diff, y_distances, test_size=test_size, random_state=random_state
    )
    
    logger.info(f"Training data: {X_train.shape[0]} pairs, Test data: {X_test.shape[0]} pairs")
    
    # Train Random Forest
    logger.info(f"Training Random Forest with {n_estimators} trees")
    rf = RandomForestRegressor(
        n_estimators=n_estimators,
        max_features=max_features,
        random_state=random_state,
        n_jobs=-1
    )
    
    rf.fit(X_train, y_train)
    
    # Predict on test set
    y_pred = rf.predict(X_test)
    
    # Calculate performance metrics
    mse = mean_squared_error(y_test, y_pred)
    rmse = np.sqrt(mse)
    r2 = r2_score(y_test, y_pred)
    mae = mean_absolute_error(y_test, y_pred)
    
    logger.info(f"Model performance: R² = {r2:.4f}, RMSE = {rmse:.4f}, MAE = {mae:.4f}")
    
    # Calculate feature importance
    importance = rf.feature_importances_
    
    # Create feature importance DataFrame
    importance_df = pd.DataFrame({
        'Feature': feature_names,
        'Importance': importance
    })
    
    # Sort by importance
    importance_df = importance_df.sort_values('Importance', ascending=False)
    
    return {
        'model': rf,
        'metrics': {
            'mse': mse,
            'rmse': rmse,
            'r2': r2,
            'mae': mae
        },
        'importance': importance_df,
        'predictions': {
            'y_test': y_test,
            'y_pred': y_pred
        }
    }

def plot_feature_importance(importance_df, output_dir, top_n=15, logger=None):
    """
    Create and save feature importance plots
    
    Args:
        importance_df: DataFrame with feature importance
        output_dir: Directory to save plots
        top_n: Number of top features to include
        logger: Logger instance
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    # Limit to top_n features
    top_features = importance_df.head(top_n)
    
    # Reverse for better visualization
    top_features = top_features.iloc[::-1]
    
    try:
        # Create figure
        plt.figure(figsize=(12, 8))
        
        # Create horizontal line + dot plot
        plt.hlines(
            y=range(len(top_features)), 
            xmin=0, 
            xmax=top_features["Importance"].values,  
            color="skyblue", 
            alpha=0.7, 
            linewidth=2
        )
        
        plt.plot(
            top_features["Importance"].values,  
            range(len(top_features)), 
            "o", 
            markersize=10, 
            color="blue", 
            alpha=0.8
        )
        
        # Add feature names and values
        plt.yticks(range(len(top_features)), top_features["Feature"].values)
        plt.xlabel("Feature Importance Score")
        plt.title("Variables Associated with Microbiome Composition Differences")
        
        # Add values next to dots
        for i, importance in enumerate(top_features["Importance"].values):  
            plt.text(importance + 0.001, i, f"{importance:.4f}", va='center')
        
        plt.tight_layout()
        
        # Save plot
        plot_file = os.path.join(output_dir, "feature_importance_plot.pdf")
        plt.savefig(plot_file, bbox_inches="tight")
        
        # Also save as PNG
        png_file = os.path.join(output_dir, "feature_importance_plot.png")
        plt.savefig(png_file, dpi=300, bbox_inches="tight")
        
        logger.info(f"Feature importance plot saved to {plot_file} and {png_file}")
        plt.close()
        
        # Create a regular bar plot as an alternative visualization
        plt.figure(figsize=(10, 6))
        sns.barplot(x="Importance", y="Feature", data=top_features[::-1], palette="viridis")
        plt.xlabel("Feature Importance Score")
        plt.ylabel("Features")
        plt.title("Top Features Driving Microbiome Differences")
        
        # Save bar plot
        bar_plot_file = os.path.join(output_dir, "feature_importance_barplot.pdf")
        plt.savefig(bar_plot_file, bbox_inches="tight")
        plt.close()
        
    except Exception as e:
        logger.error(f"Error creating feature importance plot: {str(e)}")
        logger.error(traceback.format_exc())

def save_model_results(result_dict, output_dir, logger=None):
    """
    Save model results to files
    
    Args:
        result_dict: Dictionary with model results
        output_dir: Directory to save results
        logger: Logger instance
    """
    if logger is None:
        logger = logging.getLogger('kraken_analysis')
    
    try:
        # Save feature importance
        importance_file = os.path.join(output_dir, "feature_importance.csv")
        result_dict['importance'].to_csv(importance_file, index=False)
        logger.info(f"Feature importance saved to {importance_file}")
        
        # Save model summary
        metrics = result_dict['metrics']
        with open(os.path.join(output_dir, "rf_model_summary.txt"), 'w') as f:
            f.write("Random Forest Model Summary\n")
            f.write("==========================\n\n")
            f.write(f"R² score: {metrics['r2']:.4f}\n")
            f.write(f"RMSE: {metrics['rmse']:.4f}\n")
            f.write(f"MAE: {metrics['mae']:.4f}\n\n")
            f.write("Top 10 features:\n")
            for i, row in result_dict['importance'].head(10).iterrows():
                f.write(f"  {row['Feature']}: {row['Importance']:.4f}\n")
        
        # Save predictions
        pred_df = pd.DataFrame({
            'y_true': result_dict['predictions']['y_test'],
            'y_pred': result_dict['predictions']['y_pred']
        })
        pred_file = os.path.join(output_dir, "rf_predictions.csv")
        pred_df.to_csv(pred_file, index=False)
        
    except Exception as e:
        logger.error(f"Error saving model results: {str(e)}")
        logger.error(traceback.format_exc())

def run_feature_selection(abundance_file, metadata_file, output_dir, predictors=None,
                        n_estimators=100, max_features="sqrt", distance_metric="bray",
                        transform="clr", test_size=0.2, random_state=42, log_file=None):
    """
    Run feature selection analysis to identify important variables
    
    Args:
        abundance_file: Path to abundance file
        metadata_file: Path to metadata file
        output_dir: Directory to save output files
        predictors: List of predictor variables to use
        n_estimators: Number of trees in the Random Forest
        max_features: Max features to consider at each split
        distance_metric: Distance metric to use
        transform: Transformation to apply to abundance data
        test_size: Proportion of data for testing
        random_state: Random seed for reproducibility
        log_file: Path to log file
        
    Returns:
        DataFrame with feature importance
    """
    # Setup logging
    if log_file:
        import logging.handlers
        logger = logging.getLogger('kraken_analysis')
        handler = logging.handlers.RotatingFileHandler(log_file, maxBytes=10_485_760, backupCount=5)
        handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(handler)
    else:
        logger = logging.getLogger('kraken_analysis')
    
    logger.info("Starting feature selection analysis")
    
    # Validate input files
    if not check_file_exists_with_logger(abundance_file, "Abundance file", logger):
        return None
    
    if not check_file_exists_with_logger(metadata_file, "Metadata file", logger):
        return None
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # Read abundance data
        abundance_df = pd.read_csv(abundance_file, sep='\t', index_col=0)
        logger.info(f"Loaded abundance data: {abundance_df.shape[0]} taxa, {abundance_df.shape[1]} samples")
        
        # Read metadata
        metadata_df = pd.read_csv(metadata_file, index_col=0)
        logger.info(f"Loaded metadata: {metadata_df.shape[0]} samples, {metadata_df.shape[1]} variables")
        
        # Ensure indices are strings for consistent joining
        abundance_df.index = abundance_df.index.astype(str)
        metadata_df.index = metadata_df.index.astype(str)
        
        # Identify common samples
        common_samples = set(abundance_df.columns) & set(metadata_df.index)
        if len(common_samples) < 10:
            logger.error(f"Not enough common samples between abundance and metadata: {len(common_samples)}")
            return None
        
        logger.info(f"Found {len(common_samples)} common samples between abundance and metadata")
        
        # Subset to common samples
        abundance_subset = abundance_df[list(common_samples)]
        metadata_subset = metadata_df.loc[list(common_samples)]
        
        # Transform abundance data
        abundance_transformed = transform_abundance_data(abundance_subset, transform, logger)
        
        # Get predictors
        if predictors:
            if isinstance(predictors, str):
                predictors = [p.strip() for p in predictors.split(',')]
            
            # Filter to existing columns
            predictor_cols = [p for p in predictors if p in metadata_subset.columns]
            if not predictor_cols:
                logger.error(f"None of the specified predictors exist in metadata: {predictors}")
                return None
                
            logger.info(f"Using specified predictors: {predictor_cols}")
        else:
            # Use all columns except obviously problematic ones
            exclude_cols = ['sample_id', 'sampleid', 'sample_name', 'samplename']
            predictor_cols = [col for col in metadata_subset.columns if col.lower() not in exclude_cols]
            logger.info(f"Using all available metadata columns as predictors: {predictor_cols}")
        
        # Identify categorical and numerical columns
        categorical_cols = metadata_subset[predictor_cols].select_dtypes(include=['object', 'category']).columns.tolist()
        numerical_cols = [col for col in predictor_cols if col not in categorical_cols]
        
        logger.info(f"Categorical predictors: {categorical_cols}")
        logger.info(f"Numerical predictors: {numerical_cols}")
        
        # Prepare metadata
        encoded_df = encode_categorical_variables(metadata_subset, categorical_cols, logger)
        
        # Scale numerical variables
        if numerical_cols:
            processed_df = scale_numerical_variables(encoded_df, numerical_cols, logger)
        else:
            processed_df = encoded_df
        
        # Get final list of features
        feature_names = [col for col in predictor_cols if col in processed_df.columns]
        if not feature_names:
            logger.error("No valid features available after processing")
            return None
            
        logger.info(f"Final feature set: {feature_names}")
        
        # Create pairwise differences
        X_diff, indices = create_pairwise_differences(processed_df[feature_names], logger)
        
        # Calculate microbiome distances
        y_distances = calculate_microbiome_distances(abundance_transformed.T, distance_metric, logger)
        if y_distances is None:
            logger.error("Failed to calculate microbiome distances")
            return None
        
        # Train Random Forest model
        result_dict = train_random_forest(
            X_diff, y_distances, feature_names,
            test_size=test_size,
            n_estimators=n_estimators,
            max_features=max_features,
            random_state=random_state,
            logger=logger
        )
        
        # Plot feature importance
        plot_feature_importance(result_dict['importance'], output_dir, logger=logger)
        
        # Save results
        save_model_results(result_dict, output_dir, logger)
        
        logger.info("Feature selection analysis completed successfully")
        return result_dict['importance']
        
    except Exception as e:
        logger.error(f"Error in feature selection analysis: {str(e)}")
        logger.error(traceback.format_exc())
        return None
