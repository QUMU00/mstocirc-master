# Mass Spectrometry to Translatable CircRNAs (MStoCIRC)
-----
## ONE. Usage description
usage:<br> 
```
~$ python3 analyze-tool4.py [-h] -p PEPTIDE [-s SEQUENCE] [-j JUNCTION] [-i INFO] [--output OUTPUT] [--annote ANNOTE] [--genome GENOME] [--transcript TRANSCRIPT]
<br>
```
<center class="half">
   <img src="https://gimg2.baidu.com/image_search/src=http%3A%2F%2Fnimg.ws.126.net%2F%3Furl%3Dhttp%253A%252F%252Fdingyue.ws.126.net%252F2021%252F0413%252Fa42ac95bj00qrhsip000ec000ht00b5m.jpg%26thumbnail%3D660x2147483647%26quality%3D80%26type%3Djpg&refer=http%3A%2F%2Fnimg.ws.126.net&app=2002&size=f9999,10000&q=a80&n=0&g=0n&fmt=auto?sec=1657941807&t=55a2c10a85257fa1b43145fbdce045e4" width="300">
   <img src="https://img1.baidu.com/it/u=3196272413,620674393&fm=253&fmt=auto&app=138&f=JPG?w=729&h=500" width="300">
</center>
![](https://img1.baidu.com/it/u=1935329263,825017194&fm=253&fmt=auto&app=138&f=PNG?w=1068&h=500)

>>### For example: 
There exits three model for mstocirc to analyse statistics, and below the examples are shown how to set paramaters: 
<br> ~$ python3 analyze-tool4.py -s seqence.fasta -p peptide.txt -j junction.fasta
<br> ~$ python3 analyze-tool4.py -i circRNA\_info3.txt -g GCF\_000001735.4\_TAIR10.1\_genomic.fna -a Arabidopsis\_thaliana.TAIR10.47.gtf -p peptide2.txt -j junction2.fasta 
<br> ~$ python3 analyze-tool4.py -i circRNA\_info3.txt -g GCF\_000001735.4\_TAIR10.1\_genomic.fna -r GCF\_000001735.4\_TAIR10.1\_transcript.fna -p peptide2.txt -j junction2.fasta


## TWO. Input files format description:
 + peptide.txt is the prelimenary result, generated by mass spectrometry software analyse;<br>
 + junction.fasta is translated by 200nt spanning BSJ sequence by six-frame translation;<br>
 + circRNA_info.txt is downloaded from different online databse such as circBase [here](http://www.circbase.org/) and PlantcircBase [here](http://ibi.zju.edu.cn/plantcircbase/), which demenostrates the detailed information of different circRNAs,such as chromosome number, start position, end position, parental gene and etc;<br>
 + sequence.fasta, circRNAs exon nucleotide sequence file  in .FASTA format.<br>
 

## THREE. Function presention

project21-10-14/
<table>
<tr>
<th>Feature</th><th>Function module Name </th><th> Function description </th><th>Folder</th>                            
</tr>
<tr>
<td>m1</td><td>map_corf </td><td> map MS-based on circular RNA open reading framework (cORF)</td><td> temp_mapp_corf<td>
</tr>
<tr>
<td>m2</td><td>sklearn_codinng</td><td>assess the circRNA-derived peptides</td><td>temp_skcoding</td>
</tr>
<tr>
<td>m3</td><td>map_junct</td><td>screen MS-based peptide spanning back-splice jucntion site (BSJ) </td><td>temp_mapp_junct<td>
</tr>
<tr>
<td>m4</td><td>peptide_merge</td><td>merge the overlapping MS-based peptide</td><td>temp_peptide_merge<td>
</tr>
<tr>
<td>m5</td><td>ires_predict</td><td>predict all the possible Internal Ribosome Entry Site (IRES) elements</td><td>temp_ires_predict</td>
</tr>
<tr>
<td>m6</td><td>path_analysis</td><td>GO enrichment and KEGG pathways</td><td>temp_enrich</td>
</tr>
<tr>
<td>m7</td><td>draw_circ</td><td>draw picture of the circRNAs</td><td>temp_draw_circ</td>
</tr>
<tr>
<td>o1</td><td>rem_peptide</td><td>remove MS-based repeat peptide</td><td>temp_option</td>
</tr>
<tr>
<td>o2</td><td>circ_classify</td><td>classify the circRNAs into six groups</td><td>temp_option</td>
</tr>
<tr>
<td>o3</td><td>circ_annote</td><td>annote the circRNA by host gene function</td><td>temp_option</td>
</tr>
<tr>
<td>o4</td><td>ms_ribo</td><td>add directive evidence Ribo-seq</td><td>temp_option</td>
</tr>
</table>                           
<br>generate 'result' folder
 

## FOUR. Running requirement:
**Mass Spectrometry to Translatable CircRNAs (MStoCIRC)** is available at <u>https://github.com/QUMU00/mstocirc-master</u>, and *it can work well without being installed in operating systerm* after you set the parameters properly, saving the time of installing it in the operating system.

MStoCIRC includes one main program and several small programs and during statistics analysis, those small programs are called to perform their specific functions by main program. Except enrichment analysis part dependent on package 'clusterProfiler' in R environment, the body of other main and small are programmed by Python3 language, meaning that *Python3 installation for your computer is necessary condition to run MStoCIRC*. If you want the tool to take the GO enrichment and KEGG pathways analysis, R is suggested to be installed together. But both `Python3` and `R Language` are the well-used programming languages for science research, therefore it's easy to satisfy the requirements for MStoCIRC. 

MStoCIRC doesn't rely on other tool (e.g. bedtool), but like other application developed by Python3 it needs some dependencies. <u>'matplotlib'</u> and <u>'scikit-learn'</u> module for Python3 and package <u>'clusterProfiler' </u> for R. You can install all the dependencies before running MStoCIRC or tempt to run with files in'/data'. And the tool will tempt firstly to install the dependencies automatically, during this period please keep the Internet connected for MStoCIRC. 
  
We introduce modified <u>IRESfinder [here](https://github.com/xiaofengsong/IRESfinder) </u>with into MStoCIRC for function module " IRES_predict ", and explain here on purpose what we do with their team tool to respect for their copyright.

<u>Please transform your important input file format as same as the test input file in '\data' recognized well by MStoCIRC </u>. And there are not specific requirements for general input file,such as genome sequences file (.fa, .fna format) and genome annotation file (.gtf or .gff format), but we recommend to download them from NCBI websiten [here](https://www.ncbi.nlm.nih.gov/genome/?term= "here"). 

To perform one specific function of small program, here we do not describe them about how to use the small program in detail, and README file in their own folers are available which contain the usage information.

Thank you for your downloading and running this tool. If you have any questions, vcontact with us by email **(glli@snnu.edu.cn)**.

Edited by Zhou Cao on October 4th, 2021.<br>Thank you....


<center>
<img src="https://img0.baidu.com/it/u=4197095029,2338689732&fm=253&fmt=auto&app=138&f=JPEG?w=499&h=359">
</center>

