La Jolla Assembler
==============

### Version: 0.2

La Jolla Assembler(LJA) is a tool for genome assembly from PacBio HiFI reads based on de Bruijn graphs.
LJA uses very high values of k for de Bruijn graph construction thus automatically resolving almost all repeats even in mammalian genomes.
LJA consists of three modules each of which addresses a computational challenge that is unique for HiFi reads:
(i) jumboDBG module for de Bruijn graph construction for arbitrarily large values of k;
(ii) mowerDBG module for almost perfect error correction and
(iii) multiDBG module for using variable values of k in different parts of the genome and using the full length of HiFi reads for repeat resolution.
In addition, since LJA compresses homopolymers in reads to avoid frequent errors in homopolymers, LJA uses LJApolisher tool for uncompressing and polishing final contigs.
For more details please refer to our [paper](https://www.biorxiv.org/content/10.1101/2020.12.10.420448).

Please note that LJA software is still a work in progress.
In current version diploid assembly is an experimental feature and for now we can not combine HiFI reads with other technologies.
We constantly work on improving LJA performance and results.
Please check for the new versions regularly and send any questions and bug reports to [anton.bankevich@gmail.com](mailto:anton.bankevich@gmail.com).  


For LJA installation and running instructions please refer to [LJA manual](docs/lja_manual.md).
We also provide jumboDBG module for de Bruijn graph construction as a separate script.
For jumboDBG running instructions please refer to [jumboDBG manual](docs/jumbodbg_manual.md).

**For running multiplexDBG by itself,** just run the usual LJA installation process (cmake / make).
You should then see a `multiplexDBG` binary in the `bin/` directory:

```
multiplexDBG
Usage: multiplexDBG [options] -o <output-dir> -g <graph> -a <aln> -k <int>

Options:
  -o <file_name> (or --output-dir <file_name>)  Name of output folder. multiplexDBG outputs will be stored here.
  -g <file_name> (or --graph <file_name>)       mowerDBG output GFA graph (i.e. .../01_TopologyBasedCorrection/final_dbg.gfa).
  -a <file_name> (or --aln <file_name>)         mowerDBG output alignment (i.e. .../01_TopologyBasedCorrection/final_dbg.aln).
  -k <int>                                      Big k-mer size that was used in the final mowerDBG step (probably 5001).
  --max-k <int>                                 Default 40000.
  --unique-threshold <int>                      Default 40000.
  --diploid                                     Diploidy flag.
  -h (or --help)                                Print this help message.
  -t <int> (or --threads <int>)                 Number of threads. The default value is 16.
```

License
-------

This tool is distributed under a BSD license. See the [LICENSE file](LICENSE) for details.


Credits
-------
If you use this software in your research please cite our [paper](https://www.biorxiv.org/content/10.1101/2020.12.10.420448) or the upcoming paper in Nature Biotechnology when it is published.
LJA is developed by Anton Bankevich and Andrey Bzikadze in [Pavel Pevzner's lab at UCSD](http://cseweb.ucsd.edu/~ppevzner/)
and Dmitry Antipov in [Center for Algorithmic Biology at SPbSU](https://cab.spbu.ru/).
