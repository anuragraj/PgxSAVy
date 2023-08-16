<div align="center">
<img src="images/logo.jpg" width="350px"/>
</div>

# PgxSAVy
PgxSAVy is a tool for quality control and annotation of variant peptides identified in proteogenomics. It has been rigorously tested and has demonstrated excellent performance across a wide range of datasets- simulated and real-world. PgxSAVy is suitable for rescoring and quality assessment of variant PSMs from proteogenomics search results. PgxSAVy is available as a stand-alone command-line tool, providing freely available, easy-to-use, accessible open source utility.

## Features
PgxSAVy offers quality checks specifically designed for proteogenomics searches, ensuring accurate and reliable results. The tool includes match quality checks to assess the quality and reliability of the identified variant peptides. It also performs isobaric checks, comparing variant peptides with other amino acids to identify potential similarities or discrepancies and using positional variant decoys to assess and validate the positional accuracy of identified variant peptides. The following feature checks are used for calculating the qulaity scores of variant peptides:
-  SAV quality checks from proteogenomics searches
-  Match quality checks
-  Global search parameters (like Search engine counts , if multiple search tools used) and number of sibling PSMs
-  Isobaric checks for other amino acids
-  Isobaric checks for combination of 2 amino acids
-  Isobaric checks for common modifications
-  Checks with wild type peptide quality
-  Checks with positional variant decoys   

## Installation and Download
PgxSAVy is implemented in the cross-platform Perl programming language and is platform independent i.e. can be used on computers running Windows, Linux, or Mac OS X. PgxSAVy is freely available for academic use. You will need to download and setup Perl and some modules to use PgxSAVy. See instructions for [downloading and setting up PgxSAVy](https://github.com/anuragraj/PgxSAVy/wiki/Installation).

## Supported file formats
- ### Tab Separated Text format
PgxSAVy can read inputs and write outputs in tab seperated text formats, making it easy to use. See the [complete documentation](https://github.com/anuragraj/PgxSAVy/wiki), including a list of [Frequently Asked Questions](https://github.com/anuragraj/PgxSAVy/wiki/Frequently-Asked-Questions). Example sample files can be found [here](https://github.com/anuragraj/PgxSAVy/wiki/Examples).

- ### EuGenoSuite Output
PgxSAVy can read [EuGenoSuite](https://github.com/anuragraj/EuGenoSuite) proteogenomic search output directly. 

## Usage
PgxSAVy is a command line tool. It can be used through Command Prompt on Windows and Terminal on Linux and Mac OS X. See the [complete documentation](https://github.com/anuragraj/PgxSAVy/wiki/Overview) and [usage page](https://github.com/anuragraj/PgxSAVy/wiki/Usage) for more details.


## License
PgxSAVy is licensed under the [CC0-1.0 license](https://github.com/anuragraj/PgxSAVy/blob/main/LICENSE).


## Contact
If you have any questions, comments, or suggestions, please contact Anurag Raj at anurag.igib@gmail.com or Dr Amit Kumar Yadav at amit.yadav@thsti.res.in.

## Citation
If you use PgxSAVy in your research, please cite the following publication:  
Raj, Anurag, Suruchi Aggarwal, Amit Kumar Yadav, and Debasis Dash. "[Quality control of variant peptides identified through proteogenomics- catching the (un)usual suspects](https://www.biorxiv.org/content/10.1101/2023.05.31.542998v2.full)". bioRxiv (2023): 2023-05. doi: https://doi.org/10.1101/2023.05.31.542998

