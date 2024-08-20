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


exit()
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the data
data = pd.read_csv('/Users/leej741/Desktop/git/comparison.tsv', sep='\t')

# Function to calculate precision and recall at a given cutoff
def calculate_precision_recall(data, cutoff):
    # Apply the cutoff
    data['predicted_class'] = np.where(data['score'] >= cutoff, 1, 0)
    
    # Calculate M (Total predictions above cutoff)
    M = np.sum(data['predicted_class'])
    
    # Calculate M_correct (True Positives)
    M_correct = np.sum((data['predicted_class'] == 1) & (data['actual_class'] == "True"))
    
    # Calculate N (Total True Positives in the dataset)
    N = np.sum(data['actual_class'] == "True")
    
    # Precision = M_correct / M
    precision = M_correct / M if M > 0 else 0
    
    # Recall = M_correct / N
    recall = M_correct / N if N > 0 else 0
    
    return precision, recall

# Define a sequence of cutoff values from 1 to 0
cutoffs = np.arange(1, 0, -0.01)

# Initialize lists to store precision and recall values
precision_values = []
recall_values = []

# Loop through each cutoff to calculate precision and recall
for cutoff in cutoffs:
    precision, recall = calculate_precision_recall(data, cutoff)
    precision_values.append(precision)
    recall_values.append(recall)

# Convert the precision and recall into a dataframe
pr_data = pd.DataFrame({'cutoff': cutoffs, 'precision': precision_values, 'recall': recall_values})

# Plot the Precision-Recall curve
plt.figure(figsize=(10, 6))
plt.plot(pr_data['recall'], pr_data['precision'], color='blue')
plt.title('Precision-Recall Curve')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.grid(True)
plt.show()
