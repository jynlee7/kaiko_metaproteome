import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve

# Load your data
tsv_file_path = "/Users/leej741/Desktop/git/comparison.tsv"
data = pd.read_csv(tsv_file_path, sep='\t')

# Print column names to identify the correct columns
print(data.columns)

# Function to plot precision-recall curve
def plot_pr_curve(y_true, y_scores, title):
    precision, recall, _ = precision_recall_curve(y_true, y_scores)
    plt.plot(recall, precision, marker='.')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(title)
    plt.show()

# Filter by charge state and mass threshold and compute precision-recall
charge_states = [3]  # Example charge states
mass_threshold = 1  # Example mass threshold

for charge in charge_states:
    subset = data[(data['charge'] == charge) & (data['pepmass'] > mass_threshold)]
    
    if subset.empty:
        continue
    
    y_true = subset['true_seq'] == subset['casanovo_seq']  # Replace with actual columns for ground truth and predictions
    y_scores = subset['score']
    
    plt.figure()
    plot_pr_curve(y_true, y_scores, f'PR Curve for Charge {charge} and Mass > {mass_threshold}')

# Example to compute precision and recall for a specific cutoff
cutoff = 0.7
for charge in charge_states:
    subset = data[(data['charge'] == charge) & (data['pepmass'] > mass_threshold)]
    
    if subset.empty:
        continue
    
    y_true = subset['true_seq'] == subset['casanovo_seq']  # Replace with actual columns for ground truth and predictions
    y_scores = subset['score']
    
    predictions_above_cutoff = y_scores > cutoff
    M = predictions_above_cutoff.sum()
    M_correct = ((y_true == True) & predictions_above_cutoff).sum()
    N = len(y_true)
    
    precision = M_correct / M if M > 0 else 0
    recall = M_correct / N if N > 0 else 0
    
    print(f'Charge: {charge}, Mass Threshold: {mass_threshold}, Cutoff: {cutoff}')
    print(f'Precision: {precision}, Recall: {recall}')
