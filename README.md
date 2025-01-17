# K-mer frequency spectra analysis

To create the k-mer frequency spectra, k-mer analysis toolkit (KAT) [1] was used to count 27-mers, using the following command:

```
kat hist hiv_sample.fastq
```

As datasets (i), (iii) and (iv) had  a coverage of ~200000x and dataset (ii) had a coverage of ~150000x, no normalization over the counts was performed. At low k-mer counts exponentially more distinct k-mers are found compared to higher counts, and a log scaling was performed on the number of distinct k-mers, followed by a min-max normalization. For a more clear visualization, while retaining the original shape of the graph, 1D Gaussian smoothing was applied on the log scaled, normalized number of distinct k-mers, with a standard deviation of 5 [2].

Log-scaling:
```math
y_{i,log\ scaled} = ln(y_i+1)
```

Min-max normalization:
```math
y_{i,normalized} = \frac{y_{i} - min(y)}{max(y)}
```

1D Gaussian smoothing:
```math
y_{smoothed} = scipy.gaussian\_filter1d(y, 5)
```


## Usage
See [Example usage](src/example.ipynb) for example usage.

## Data availability
For the (i) pure HIV strain HIV-1 NL4-3, (ii) a patient sample of HIV-1, (iii) a sample of five pure HIV strains mixed together and (iv) a simulated HIV sample, the data is publicly accessible under accession numbers (i) SRR2976616 [3], (ii) DRR030309 [4], (iii) SRR961514 or dataset (iii) can be found at https://github.com/cbg-ethz/5-virus-mix [5] and dataset (iv) can be found at https://bitbucket.org/jbaaijens/savage-benchmarks [6].

## References
[1] Mapleson, D., Accinelli, G. G., Kettleborough, G., Wright, J., & Clavijo, B. J. (2016). KAT: a K-mer analysis toolkit to quality control NGS datasets and genome assemblies. Bioinformatics, 33(4), 574–576. https://doi.org/10.1093/bioinformatics/btw663

[2] Virtanen, P., Gommers, R., Oliphant, T. E., Haberland, M., Reddy, T., Cournapeau, D., Burovski, E., Peterson, P., Weckesser, W., Bright, J., Van Der Walt, S. J., Brett, M., Wilson, J., Millman, K. J., Mayorov, N., Nelson, A. R. J., Jones, E., Kern, R., Larson, E., . . . Vázquez-Baeza, Y. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. Nature Methods, 17(3), 261–272. https://doi.org/10.1038/s41592-019-0686-2

[3] Bertels, F., Leemann, C., Metzner, K. J., & Regoes, R. R. (2019). Parallel Evolution of HIV-1 in a Long-Term Experiment. Molecular Biology And Evolution, 36(11), 2400–2414. https://doi.org/10.1093/molbev/msz155

[4] Ode, H., Matsuda, M., Matsuoka, K., Hachiya, A., Hattori, J., Kito, Y., Yokomaku, Y., Iwatani, Y., & Sugiura, W. (2015). Quasispecies Analyses of the HIV-1 Near-full-length Genome With Illumina MiSeq. Frontiers in Microbiology, 6. https://doi.org/10.3389/fmicb.2015.01258

[5] Di Giallonardo, F., Töpfer, A., Rey, M., Prabhakaran, S., Duport, Y., Leemann, C., Schmutz, S., Campbell, N. K., Joos, B., Lecca, M. R., Patrignani, A., Däumer, M., Beisel, C., Rusert, P., Trkola, A., Günthard, H. F., Roth, V., Beerenwinkel, N., & Metzner, K. J. (2014). Full-length haplotype reconstruction to infer the structure of heterogeneous virus populations. Nucleic Acids Research, 42(14), e115. https://doi.org/10.1093/nar/gku537

[6] Baaijens, J. A., Aabidine, A. Z. E., Rivals, E., & Schönhuth, A. (2017). De novo assembly of viral quasispecies using overlap graphs. Genome Research, 27(5), 835–848. https://doi.org/10.1101/gr.215038.116
