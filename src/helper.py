import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d

def normalize_histogram(data):
    """
    Reads a .hist file, applies min-max normalization to the 'number of distinct k-mers' values,
    and writes the normalized data to a new file.

    Parameters:
        input_path (str): Path to the input .hist file.
        output_dir (str): Directory where the output file will be saved.
        output_file (str): Name of the output file.
    """
    
    # Min-max normalization: (value - min) / (max - min)
    min_val = data['distinct_kmers'].min()
    print(min_val)
    max_val = data['distinct_kmers'].max()
    print(max_val)
    data['distinct_kmers_normalized'] = (data['distinct_kmers'] - min_val) / (max_val - min_val)
    
    return data

def log_normalize_histogram(data):
    """
    Reads a .hist file, applies log normalization to the 'number of distinct k-mers' values,
    and writes the normalized data to a new file.

    Parameters:
        input_path (str): Path to the input .hist file.
        output_dir (str): Directory where the output file will be saved.
        output_file (str): Name of the output file.
    """
    
    # Log normalization: log(value + 1)
    data['distinct_kmers_normalized'] = data['distinct_kmers'].apply(lambda x: np.log(x + 1))
    
    return data

def zscore_normalize_histogram(data):
    """
    Reads a .hist file, applies z-score normalization to the 'number of distinct k-mers' values,
    and writes the normalized data to a new file.

    Parameters:
        input_path (str): Path to the input .hist file.
        output_dir (str): Directory where the output file will be saved.
        output_file (str): Name of the output file.
    """
    
    # Z-score normalization: (value - mean) / std
    mean = data['distinct_kmers'].mean()
    std = data['distinct_kmers'].std()
    data['distinct_kmers_normalized'] = (data['distinct_kmers'] - mean) / std
    
    return data

def plot_kmer_histograms(file_dict, output_dir, output_file, normalize=False, smoothing=False, sigma=2):
    """
    Plots the normalized 27-mer histograms from multiple .hist files.

    Parameters:
        file_dict (dict): Dictionary where keys are labels for the plot and values are paths to the normalized .hist files.
    """
    plt.figure(figsize=(10, 6))
    
    # Iterate over each label and file path in the dictionary
    for label, file_path in file_dict.items():
        
        data = pd.read_csv(file_path, delim_whitespace=True, comment='#', header=None, names=['kmer_frequency', 'distinct_kmers'])
        
        # Check for and drop any rows with missing values
        data.dropna(inplace=True)
        
        # Normalize the data if specified
        if normalize:
            log_normalized_data = log_normalize_histogram(data)
            # for the log normalized data, change the y_data to the normalized values
            log_normalized_data['distinct_kmers'] = log_normalized_data['distinct_kmers_normalized']
            normalized_data = normalize_histogram(log_normalized_data)
            y_data = data['distinct_kmers_normalized']
            ylabel = "Normalized Number of Distinct 27-mers"
        else:
            y_data = data['distinct_kmers']
            ylabel = "Number of Distinct 27-mers"
            
        if smoothing:
            y_data = gaussian_filter1d(y_data, sigma)
            
        # Plot the data
        plt.plot(data['kmer_frequency'], y_data, label=label, alpha=0.5)
    
    # Configure the plot
    #plt.title("27-mer Frequency Spectra")
    plt.xlabel("27-mer Frequency", fontsize=16)
    ylabel = "Number of Distinct 27-mers"
    if normalize:
        ylabel = "Normalized Number of Distinct 27-mers"
    plt.ylabel(ylabel, fontsize=16)
    plt.grid(True)
    
    #y-height max 80
    if not normalize:
        plt.ylim(0, 80)
    else:
        plt.ylim(0, 0.6)
        
    plt.xlim(0, 9999)
        
    #if normalize:
        #use log scale for y-axis
        #plt.yscale('log')
        
    #make the text larger
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    
    #put legend in top left
    plt.legend(loc='upper left', fontsize=14)
    
    
    # Show the plot
    plt.tight_layout()
    plt.show()

def plot_k_mer_histograms_subplots(file_dict, output_dir, output_file, normalize=False):
    """
    Plots the normalized 27-mer histograms from multiple .hist files.

    Parameters:
        file_dict (dict): Dictionary where keys are labels for the plot and values are paths to the normalized .hist files.
    """
    fig, axs = plt.subplots(2, 2, figsize=(15, 10))
    
    # Iterate over each label and file path in the dictionary
    for i, (label, file_path) in enumerate(file_dict.items()):
        
        data = pd.read_csv(file_path, delim_whitespace=True, comment='#', header=None, names=['kmer_frequency', 'distinct_kmers'])
        
        # Check for and drop any rows with missing values
        data.dropna(inplace=True)
        
        if normalize:
            log_normalized_data = log_normalize_histogram(data)
            # for the log normalized data, change the y_data to the normalized values
            log_normalized_data['distinct_kmers'] = log_normalized_data['distinct_kmers_normalized']
            normalized_data = normalize_histogram(log_normalized_data)
            y_data = data['distinct_kmers_normalized']
            ylabel = "Normalized Number of Distinct 27-mers"
        else:
            y_data = data['distinct_kmers']
            ylabel = "Number of Distinct 27-mers"
        #plot histograms
        axs[i//2, i%2].plot(data['kmer_frequency'], y_data, label=label)
        axs[i//2, i%2].fill_between(data['kmer_frequency'], y_data)
        axs[i//2, i%2].set_title(label)
        axs[i//2, i%2].set_ylabel("Number of Distinct 27-mers")
        axs[i//2, i%2].set_xlabel("Frequency")
        axs[i//2, i%2].grid(True)
        if not normalize:
            axs[i//2, i%2].set_ylim(0, 80)
        axs[i//2, i%2].set_xlim(-100, 9999)
        #if normalize:
            #use log scale for y-axis
            #ax.set_yscale('log')
    
    # Show the plot
    plt.tight_layout()
    plt.show()