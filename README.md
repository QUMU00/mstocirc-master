# mstocirc-master

## ONE. Usage description
usage: python3 analyze-tool4.py [-h] -p PEPTIDE [-s SEQUENCE] [-j JUNCTION] [-i INFO]
                        [--output OUTPUT] [--annote ANNOTE] [--genome GENOME]
                        [--transcript TRANSCRIPT]
### for example: 
~$ python3 analyze-tool4.py -s seqence.fasta -p peptide.txt -j junction.fasta
<br> ~$ python3 analyze-tool4.py -i circRNA_info3.txt -g GCF_000001735.4_TAIR10.1_genomic.fna -a Arabidopsis_thaliana.TAIR10.47.gtf -p peptide2.txt -j junction2.fasta 


## TWO. input files format description:
 1. peptide.txt is the result, generated by LC-MS/MS software analyse,
 2. junction.fasta is translated by 200nt spanning BSJ sequence by 6-frame,
 3. circRNA_info.txt is downloaded from different online databse, which demenostrates the detailed information of different circRNAs,
 4. sequnce.fasta, circRNA nucleotide sequnce file.
 

## THREE. function presention

project21-10-14/

feature                       | function name         | description                                     | folder
------------------------------|-----------------------|-------------------------------------------------|---------
main function                 | map_corf              | map MS-based on circRNA ORF                     | 1mapp_corf
     2                        | sklearn_codinng       | assess the circRNA-derived peptides             | 1.5skcoding
     3                        | map_junct             | screen MS-based peptide spanning BSJ            | 2mapp_junct
     4                        | peptide_merge         | merge the overlapping MS-based peptide          | 3peptide_merge
     5                        | ires_predict          | predict all the possible IRES elements          | 4ires_predict
     6                        | path_analysis         | GO enrichment and KEGG pathways                 | 5enrich
     7                        | draw_circ             | draw picture of the circRNAs                    | 6draw_circ
optional fucntion             | rem_peptide           | remove MS-based repeat peptide                  | 7option
     2                        | circ_classify         | classify the circRNAs into six groups           | 7option
     3                        | circ_annote           | annote the circRNA by host gene function        | 7option
     4                        | ms_ribo               | add directive evidence Ribo-seq,                | 7option                      
    
<br>generate 'result' directory
 

## FOUR. install requirement:
   Mstocirc is available at https://github.com/QUMU00/mstocirc-master/, and *it can work well without being installed* after you set the parameters properly, saving the time to
install in the oprating system.
   Mstocirc includes one main program and several mini programs，and during statistics analysis, those mini programs are called to perform their specific functions by main 
program. Except enrichment analysis part dependent on package 'clusterProfiler' in R environment，the body of other main and small are programmed by python3 language，meaning 
that *python3 installation for your computer is necessary condition to run mstocirc*. If you want the tool to analyze the GO term and KEGG pathways，R is suggested to be 
installed together. But both python3 and R are the well-used programming languages for science research，therefore it's easy to satisfy the requirements for mstocirc. 
Mstocirc doesn't rely on other tool (e.g. bedtool)，but like other application developed by python3 it needs some dependencies，'matplotlib' and 'scikit-learn' module for 
python3 and package 'clusterProfiler' for R. You can install all the dependencies before running mstocirc，or tempt to run with file in'/ data'. And the tool will tempt firstly to install them automatically, during this period please keep the Internet connected. 
   please keep your important input file format as same as the test input file recognized well by mstocirc. And there are not specific requirements for general input file ，
such as genome sequences file (.fa, .fna format) and genome annotation file (.gtf or .gff format) but we recommend to download them from NCBI website. 
To perform one specific function of small program, here we do not describe them about how to use in detail, and README file in their own directory are available which 
containing the usage information.

   Thank you for your downloading and running this tool. If you have any questions，contact with us by email *(glli@snnu.edu.cn)*.

Edited by Zhou Cao on October 4th, 2021.

Thank you....

