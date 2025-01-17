# K-mer frequency spectra analysis

To create the k-mer frequency spectra, k-mer analysis toolkit (KAT) [1] was used to count 27-mers, using the following command:

```
kat hist hiv_sample.fastq
```

As datasets (i), (iii) and (iv) had  a coverage of ~200000x and dataset (ii) had a coverage of ~150000x, no normalization over the counts was performed. At low k-mer counts exponentially more distinct k-mers are found compared to higher counts, and a log scaling was performed on the number of distinct k-mers, followed by a min-max normalization. For a more clear visualization, while retaining the original shape of the graph, 1D Gaussian smoothing was applied on the log scaled, normalized number of distinct k-mers, with a standard deviation of 5 [2].

Log-scaling:
$$
y_{i,log\ scaled} = ln(y_i+1)
$$

Min-max normalization:
$$
y_{i,normalized} = \frac{y_{i} - min(y)}{max(y)}
$$

1D Gaussian smoothing:
$$
y_smoothed = scipy.gaussian_filter1d(y, 5)
$$


## Usage
See src/example.ipynb for example usage.

## Data availability

## References
[1] Mapleson, D., Accinelli, G. G., Kettleborough, G., Wright, J., & Clavijo, B. J. (2016). KAT: a K-mer analysis toolkit to quality control NGS datasets and genome assemblies. Bioinformatics, 33(4), 574–576. https://doi.org/10.1093/bioinformatics/btw663

[2] Virtanen, P., Gommers, R., Oliphant, T. E., Haberland, M., Reddy, T., Cournapeau, D., Burovski, E., Peterson, P., Weckesser, W., Bright, J., Van Der Walt, S. J., Brett, M., Wilson, J., Millman, K. J., Mayorov, N., Nelson, A. R. J., Jones, E., Kern, R., Larson, E., . . . Vázquez-Baeza, Y. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. Nature Methods, 17(3), 261–272. https://doi.org/10.1038/s41592-019-0686-2