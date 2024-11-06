import pandas as pd
import numpy as np
from pathlib import Path
import re

class MzTabReader:
    def __init__(self):
        self.metadata = {}
        self.psm_data = None
        
    def read_mztab(self, filepath):
        """
        Read an mzTab file and extract PSM data and metadata.
        
        Parameters:
        -----------
        filepath : str
            Path to the mzTab file
            
        Returns:
        --------
        pandas.DataFrame
            Processed PSM data
        """
        psm_lines = []
        metadata_lines = []
        
        # Read the file and separate metadata from PSM data
        with open(filepath, 'r') as file:
            for line in file:
                if line.startswith('MTD'):
                    metadata_lines.append(line.strip())
                elif line.startswith('PSM'):
                    psm_lines.append(line.strip())
                    
        # Process metadata
        self._process_metadata(metadata_lines)
        
        # Process PSM data
        if psm_lines:
            # First line contains headers
            headers = psm_lines[0].split('\t')
            # Remove 'PSM' prefix from header
            headers = [h[4:] if h.startswith('PSM_') else h for h in headers]
            
            # Process remaining lines
            data = [line.split('\t') for line in psm_lines[1:]]
            self.psm_data = pd.DataFrame(data, columns=headers)
            
        return self.psm_data
    
    def _process_metadata(self, metadata_lines):
        """Process metadata lines into a dictionary."""
        for line in metadata_lines:
            parts = line.split('\t')
            if len(parts) >= 3:
                key = parts[1]
                value = parts[2]
                self.metadata[key] = value

def prepare_precision_recall_data(mztab_file, score_column='search_engine_score[1]', 
                                decoy_prefix='DECOY_', score_threshold=None):
    """
    Prepare data from mzTab file for precision-recall curve analysis.
    
    Parameters:
    -----------
    mztab_file : str
        Path to mzTab file
    score_column : str
        Name of the column containing peptide scores
    decoy_prefix : str
        Prefix used to identify decoy sequences
    score_threshold : float, optional
        Threshold for filtering scores
        
    Returns:
    --------
    tuple
        (y_true, y_scores) arrays for precision-recall curve
    """
    # Read mzTab file
    reader = MzTabReader()
    psm_data = reader.read_mztab(mztab_file)
    
    # Convert score column to numeric, handling any non-numeric values
    psm_data[score_column] = pd.to_numeric(psm_data[score_column], errors='coerce')
    
    # Remove rows with missing scores
    psm_data = psm_data.dropna(subset=[score_column])
    
    # Determine true/decoy status
    psm_data['is_decoy'] = psm_data['sequence'].str.startswith(decoy_prefix)
    
    # Create arrays for precision-recall curve
    y_true = (~psm_data['is_decoy']).astype(int).values
    y_scores = psm_data[score_column].values
    
    # Normalize scores to [0,1] if they aren't already
    if y_scores.max() > 1 or y_scores.min() < 0:
        y_scores = (y_scores - y_scores.min()) / (y_scores.max() - y_scores.min())
    
    # Apply threshold if specified
    if score_threshold is not None:
        mask = y_scores >= score_threshold
        y_true = y_true[mask]
        y_scores = y_scores[mask]
    
    return y_true, y_scores

def analyze_mztab(mztab_file, output_dir=None):
    """
    Analyze mzTab file and generate precision-recall curve with statistics.
    
    Parameters:
    -----------
    mztab_file : str
        Path to mzTab file
    output_dir : str, optional
        Directory to save results
    """
    # Prepare data
    y_true, y_scores = prepare_precision_recall_data(mztab_file)
    
    # Calculate precision-recall curve
    precision, recall, thresholds = precision_recall_curve(y_true, y_scores)
    
    # Create plot
    plt.figure(figsize=(10, 6))
    plt.plot(recall, precision, label='Precision-Recall curve')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall Curve for Peptide Identifications')
    plt.grid(True)
    
    # Add statistics to plot
    total_psms = len(y_true)
    target_psms = sum(y_true)
    decoy_psms = total_psms - target_psms
    
    stats_text = f'Total PSMs: {total_psms}\n'
    stats_text += f'Target PSMs: {target_psms}\n'
    stats_text += f'Decoy PSMs: {decoy_psms}\n'
    plt.text(0.02, 0.02, stats_text, transform=plt.gca().transAxes, 
             bbox=dict(facecolor='white', alpha=0.8))
    
    plt.legend()
    
    # Save results if output directory is specified
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save plot
        plt.savefig(output_dir / 'precision_recall_curve.png')
        
        # Save data points
        results_df = pd.DataFrame({
            'Precision': precision,
            'Recall': recall,
            'Threshold': np.append(thresholds, 1.0)  # Add final threshold
        })
        results_df.to_csv(output_dir / 'precision_recall_data.csv', index=False)
    
    plt.show()
    
    return precision, recall, thresholds

# Example usage
if __name__ == "__main__":
    # Replace with your mzTab file path
    mztab_file = "path/to/your/file.mztab"
    
    try:
        # Analyze single file
        precision, recall, thresholds = analyze_mztab(
            mztab_file,
            output_dir="results"
        )
        
        # Print summary statistics
        print("\nSummary Statistics:")
        print(f"Maximum Precision: {max(precision):.3f}")
        print(f"Maximum Recall: {max(recall):.3f}")
        print(f"Number of threshold points: {len(thresholds)}")
        
    except FileNotFoundError:
        print("Please specify the correct path to your mzTab file.")