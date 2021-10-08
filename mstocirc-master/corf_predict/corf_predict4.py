#!/usr/bin/python
#*_coding: utf-8_*
__author__='qumu0000'


import argparse
import os
import re

def proano():
  print('extract the exon positon from the annotion.gtf...')
  fi=open(args.annotion,'r+')
  fo=open('corf_predict/exon_infor.txt','w+')
  fp=open('corf_predict/exon_infor2.txt','w+')
  lock_gene=0
  llock=0
  exon_cstr=[]
  exon_cend=[]
  exon_str2=[]
  exon_end2=[]
  chrm=[]
  chrm2=[]
  std=''
  stdl=[]
  stdl2=[]
  cuunt=0
  #for sp in fi.readlines():
  while(1):
    sp=fi.readline()
    if(sp==''):
      break

    sp=sp.strip()
    spp=sp.split('\t')
    #print(spp)
    if(sp[:2]=='#!' or sp[0]=='#'): 
      continue 
    #print '03'
    fp.write('chr'+spp[0]+'\t'+spp[2]+'\t'+spp[3]+'\t'+spp[4]+'\t%s\n'%spp[6])
    if(spp[2]=='gene'):
      lock_gene=lock_gene+1
    if(spp[2]=='transcript' or spp[2]=='mRNA' or spp[2]=='lnc_RNA'):
      #print('transcript')
      llock=1
    if(lock_gene%2==0):
      llock=1
      lock_gene=lock_gene-1
      
    if(spp[2]=='exon' and llock==0):
      #print(spp[3],spp[4])
      chrm.append('chr'+spp[0])
      chrm2.append('chr'+spp[0])
      exon_str2.append(spp[3])
      exon_end2.append(spp[4])
      #std=spp[6]
      stdl.append(spp[6])
      stdl2.append(spp[6])
      
      
    if(llock==1):
      num=len(exon_str2)
      llock=0
      
      if(num==0):
        continue
      for i in range(0,num):
          
        if(stdl[i]=='+'):
          exon_cstr.append(exon_str2[i])
          exon_cend.append(exon_end2[i])
          #print '01'
          fo.write(chrm[i]+'\texon\t'+exon_str2[i]+'\t'+exon_end2[i]+'\t'+stdl[i]+'\n')
        else:
          exon_cstr.append(exon_str2[num-1-i])
          exon_cend.append(exon_end2[num-1-i])
          fo.write(chrm[num-1-i]+'\texon\t'+exon_str2[num-1-i]+'\t'+exon_end2[num-1-i]+'\t'+stdl[num-1-i]+'\n')
      
      exon_str2=[]
      exon_end2=[]
      chrm=[]
      stdl=[]
      #chrm.append('~')
      chrm2.append('~')
      #stdl.append('~')
      stdl2.append('~')
      exon_cstr.append('~')
      exon_cend.append('~')
      fo.write('~\texon\t~\t~\t~\n')
    
  num=len(exon_str2)
 
  for i in range(0,num):
    if(stdl[i]=='+'):
      exon_cstr.append(exon_str2[i])
      exon_cend.append(exon_end2[i])
          #print '01'
      fo.write(chrm[i]+'\texon\t'+exon_str2[i]+'\t'+exon_end2[i]+'\t'+stdl[i]+'\n')
    else:
      exon_cstr.append(exon_str2[num-1-i])
      exon_cend.append(exon_end2[num-1-i])
      fo.write(chrm[num-1-i]+'\texon\t'+exon_str2[num-1-i]+'\t'+exon_end2[num-1-i]+'\t'+stdl[num-1-i]+'\n')
  #print(len(chrm2),len(exon_cstr),len(exon_cend),len(stdl2))
  #for i in range(0,len(exon_cstr)):
  #  fo.write(chrm[i]+'\t'+'exon'+'\t'+exon_cstr[i]+'\t'+exon_cend[i]+'\t'+stdl[i]+'\n')
  fi.close()
  fo.close()
  fp.close()
  return chrm2,exon_cstr,exon_cend,stdl2
  

def exgain():
  proano()
  print('acquire the exon position of circRNA...')
  fi=open('corf_predict/exon_infor.txt','r+')
  fj=open(args.input,'r+')
  fo=open('corf_predict/__TEMP_exon.txt','w+')
  #fp=open('./__MINOTOR__.txt','w+')
  chck=0
  lock=0
  cnt=0
  std='+'
  cirh=[]
  exoo=[]
  
    #print exon_cstr,exon_cend
  #for st in fj.readlines():
  while(1):
    st=fj.readline()
    if(st==''):
      break
    stt=st.strip().split('\t')
    
    #cnt=cnt+1 
    
    
    #stt3=stt[1].split(' ')
    
    if(len(stt)<5):
      continue
    chrm='chr'+stt[2]
    strand =''
    
    exon_cstr=[]
    exon_cend=[]
    chck=0
    lock=0
    
    exon_str=[]
    exon_end=[]
    fi.seek(0,0)
    stdl=[]
    for sp in fi.readlines():
      spt=sp.strip().split('\t')
      chck=0
      #print(spt[0],chrm)
      if(spt[0]==chrm):
        exon_cstr.append(spt[2])
        exon_cend.append(spt[3])
        stdl.append(spt[4])
    
    for i in range(0,len(exon_cstr)):
      if(stt[3]==exon_cstr[i]):
        strand=stdl[i]
        #print(strand)
        lock=1
      if(lock==1):
        chck=chck+1
        exon_str.append(exon_cstr[i])
        exon_end.append(exon_cend[i])
        if(exon_cend[i]==stt[4]):
          lock=0
          break 
        elif(exon_cend[i]>stt[4] or chck>20):
          exon_str=[]
          exon_end=[]
          lock=0
          break

    fo.write(stt[0]+'\t'+strand+'\t')
    cirh.append(stt[0])
    exos=strand+' '
    for i in range(0,len(exon_str)):
      fo.write(exon_str[i]+',')
      exos=exos+exon_str[i]+','
    fo.write('\t')
    exos=exos+' '
    for i in range(0,len(exon_end)):
      fo.write(exon_end[i]+',')
      exos=exos+exon_end[i]+','
    fo.write('\n')
    #print exos
    exoo.append(exos)
    
    #print 'work%02d Done'%cnt
    #print stt[0],',',stt[1],',DONE'
  fi.close()
  fo.close()
  fj.close()
  return dict(zip(cirh,exoo))


def exgain2():
  #proano()
  chrm,exon_str,exon_end,strd=proano()
  #fi=open('corf_predict/exon_infor.txt','r+')
  print('acquire exon position of the circRNA...')
  fj=open(args.input,'r+')
  
  
  fo=open('corf_predict/__TEMP_exon.txt','w+')
  fq=open('corf_predict/error.txt','w+')
  chrm_exon_dict=dict()
  chrmm=''
  exon_strr=[]
  exon_endd=[]
  strdd=[]
  #chrm,exon_str,exon_end,strd=proano()
  #print(chrm,len(exon_str))
  for i in range(0,len(chrm)):
    if(chrmm!=chrm[i] and chrm[i]!='~'):
      if(chrmm!=''):
        chrm_exon_dict[chrmm+'_exon_str']=exon_strr
        chrm_exon_dict[chrmm+'_exon_end']=exon_endd
        chrm_exon_dict[chrmm+'_strd']=strdd
        exon_strr=[]
        exon_endd=[]
        strdd=[]
      chrmm=chrm[i]
      
    exon_strr.append(exon_str[i])
    exon_endd.append(exon_end[i])
    strdd.append(strd[i])
  
  chrm_exon_dict[chrmm+'_exon_str']=exon_strr
  chrm_exon_dict[chrmm+'_exon_end']=exon_endd
  chrm_exon_dict[chrmm+'_strd']=strdd
  
  exoo=[]
  strddd=''
  fj.readline()
  for sj in fj.readlines():
    sjt=['','','','','','','']
    sjtt=sj[:-1].split('\t')
    #sjt: circID, null,chrm.start,end,gene_symbl,null
    #input: chrom start end strand circRNA_ID transcript gene_symbol
    sjt[0]=sjtt[4]
    sjt[2]=sjtt[0]
    sjt[3]=sjtt[1]
    sjt[4]=sjtt[2]
    sjt[5]=sjtt[6]
     
    print(sjt[0])
    str_ind=0
    exos_str=''
    exos_end=''
    
    
    for i in range(0,len(chrm_exon_dict['chr%s_exon_str'%sjt[2]])):
      if(chrm_exon_dict['chr%s_exon_str'%sjt[2]][i]==sjt[3]):
        str_ind=i
        strddd=chrm_exon_dict['chr%s_strd'%sjt[2]][i]
        fq.write(sjt[0]+'\t'+chrm_exon_dict['chr%s_exon_str'%sjt[2]][i]+'\t'+strddd+'\n')
        print('dfg-1')
        break
    if(str_ind==0):
      continue
    #for i in range(str_ind,len(chrm_exon_dict['chr%s_exon_end'%sjt[2]])):
    for i in range(str_ind,str_ind+50):
      #print(sjt[4])
      if(i>=len(chrm_exon_dict['chr%s_exon_end'%sjt[2]]) or chrm_exon_dict['chr%s_exon_end'%sjt[2]][i]=='~'):
        break
      #if():
      #  break
      if(chrm_exon_dict['chr%s_exon_end'%sjt[2]][i]==sjt[4]):
        print('dfg-2')
        exos_str=exos_str+','+chrm_exon_dict['chr%s_exon_str'%sjt[2]][i]
        exos_end=exos_end+','+sjt[4]
        exoo.append(sjt[0]+' '+sjt[2]+' '+sjt[3]+' '+sjt[4]+' %s '%sjt[5]+strddd+' '+exos_str[1:]+' '+exos_end[1:])
        fo.write(sjt[0]+' '+sjt[2]+' '+sjt[3]+' '+sjt[4]+' %s '%sjt[5]+strddd+' '+exos_str[1:]+' '+exos_end[1:]+'\n')
        break
      
      if(chrm_exon_dict['chr%s_exon_str'%sjt[2]][i]==sjt[3]):
        exos_str=''
        exos_end=''
      exos_str=exos_str+','+chrm_exon_dict['chr%s_exon_str'%sjt[2]][i]
      exos_end=exos_end+','+chrm_exon_dict['chr%s_exon_end'%sjt[2]][i]      
    
  fj.close()
  fo.close()
  return exoo






def recompl(seqn):
  seqn=seqn.upper()
  seqn=seqn.replace('T','a')
  seqn=seqn.replace('A','t')
  seqn=seqn.replace('G','c')
  seqn=seqn.replace('C','g')
  seqn=seqn.upper()
  return seqn[::-1]
  



def fagain():
  exoo=exgain2()
  #fi=open('./__TEMP_exon.txt','r+')
  fj=open(args.genome,'r+')
  #fm=open(args.input,'r+')
  fo=open(args.output+'/circRNA_exon.fasta','w+')
  fp=open('corf_predict/__monitor__.txt','w+')
  fq=open(args.output+'/circ_draw.txt','w+')
  
  exon_str=[]
  exon_end=[]
  
  #chrml=[]
  #seql=[]
  sequen=''
  chrm_num=''
  print ('generate genome dict...')

  #for sj in fj.readlines():
   
  #  if(sj[0]=='>'):
  #    if(sequen!=''):
  #      seql.append(sequen)
  #    chrml.append('chr'+sj[1:].split(' ')[0])
  #    sequen=''
  #    continue
  #  sequen=sequen+sj.strip()
  #seql.append(sequen)
  #for n in range(0,len(chrml)):
  #  fp.write(chrml[n]+'\n'+seql[n]+'\n')
    
  #dit1=dict(zip(chrml,seql))
  
  dit1=dict()
  for sj in fj.readlines():
    if(sj[0]=='>'):
      if(chrm_num!=''):
        dit1[chrm_num]=sequen
      chrm_num=sj.split(' ')[0][1:]
      sequen=''
      continue
    sequen=sequen+sj.strip()
  dit1[chrm_num]=sequen
 
  print ('extract circular RNA sequence fasta...')
   
  for line in exoo:
  #seqs=''
    seqe=''
    exoot=line.split(' ')
    strand=exoot[5]
    chrm=exoot[1] 
    #seqs=dit1[chrm][star-1:end]
    #fo.write('>'+smt[0]+'\n')
    lge=0
    exon_p='0'
    
    #exon_str=ddtt[1].split(', ')[0].split(',')
    #exon_end=ddtt[1].split(', ')[1].strip(',\n').split(',')
    exon_str=exoot[-2].split(',')
    exon_end=exoot[-1].split(',')
      
    for i in range(0,len(exon_str)):
      seqe=seqe+dit1[chrm][int(exon_str[i])-1:int(exon_end[i])]
      exon_p=exon_p+'-'+str(len(seqe))
      lge=lge+len(dit1[chrm][int(exon_str[i])-1:int(exon_end[i])])
    if(len(seqe)<50):
      continue
    seqe=seqe.upper()
    fo.write('>'+exoot[0]+':'+exoot[4]+'\n')
    fq.write(exoot[0]+','+str(len(seqe))+','+exon_p+',')
    if(strand=='-'):
      seqe=recompl(seqe)
    fo.write(seqe+'\n')
    fq.write(seqe+'\n')
      
    
      
      
  #fi.close()
  fj.close()
  #fm.close()
  fo.close()
  fp.close()
  fq.close()    









def fagain2():
  os.system('cp %s %s/circRNA_exon.fasta'%(args.sequence,args.output))
  fi=open(args.sequence,'r+')
  fo=open(args.output+'/circ_draw.txt','w+')
  while(True):
    si=fi.readline()
    if(si==''):
     break
    if(si[0]=='>'):
      fo.write(si[1:].split(':')[0])
      seqn=fi.readline().strip()
      ll=str(len(seqn))
      fo.write(','+ll+',0-'+ll+','+seqn+'\n')
  
  fi.close()
  fo.close()


def fagain3():
  
  #fi=open('GCF_000001405.39_GRCh38.p13_genomic.fna','r+')
  fi=open(args.genome,'r+')
  genome_dict=dict()
  genm_seq=''
  genm_chr=''
  while(1):
    si=fi.readline()
    if(si==''):
      break
    if(si[0]=='>'):
      if(genm_chr!=''):
        genome_dict[genm_chr]=genm_seq
        #print(genm_chr,genm_seq[:20])
      genm_chr=si[1:].split(' ')[0]
      genm_seq=''
      continue
    genm_seq=genm_seq+si.strip().upper()
  genome_dict[genm_chr]=genm_seq

  print('01')
  #fj=open(sys.argv[1],'r+')
  #fo=open(sys.argv[1].replace('ciri','spliced.fna'),'w+')
  fj=open(args.input,'r+')
  fo=open(args.output+'/circRNA_exon.fasta','w+')
  fq=open(args.output+'/circ_draw.txt','w+')

  fj.readline()
  circ_list=[]
  num=0
  print ('seq_extract,wait...')
  for jline in fj.readlines():
    jlinet=jline.split('\t')
    #num=num+1
    #print(jlinet)
    #fi.seek(0,0)
    #genome_seq=''
    #while(1):
    #  si=fi.readline()
    #  if(si==''):
    #    break
    #  if('>'+jlinet[2] in si):
    #    while(1):
    #      si=fi.readline()
    #      if(si=='' or '>' in si):
    #        break
    #      genome_seq=genome_seq+si.strip().upper()
    #    break
    #print(genome_seq[:20])      
    if(jlinet[3]=='-'):
      #circ_list.append('%s,%s,'%(jlinet[0],jlinet[5])+re_compl(genome_seq[int(jlinet[3])-1:int(jlinet[4])]))
      circ_list.append('%s,%s,%s,'%(jlinet[4],jlinet[5],jlinet[6])+recompl(genome_dict[jlinet[0]][int(jlinet[1])-1:int(jlinet[2])]))
    else:
      #circ_list.append('%s,%s,'%(jlinet[0],jlinet[5])+genome_seq[int(jlinet[3])-1:int(jlinet[4])])
      circ_list.append('%s,%s,%s,'%(jlinet[4],jlinet[5],jlinet[6])+genome_dict[jlinet[0]][int(jlinet[1])-1:int(jlinet[2])])
    #print('seq',genome_seq[int(jlinet[2])-1:int(jlinet[3])])
    
    

  rna_dict=dict()
  rna_hed=''
  rna_seq=''
  #fk=open('GCF_000001405.39_GRCh38.p13_rna.fna','r+')
  fk=open(args.transcript,'r+')
  while(1):
    sk=fk.readline()
    if(sk==''):
      break
    if(sk[0]=='>'):
      if(rna_hed!=''):
        rna_dict[rna_hed]=rna_seq  
        #print(rna_hed,rna_seq[:50])
      #rna_hed=sk[sk.find('(')+1:sk.find(')')]
      rna_hed=sk[1:].split(' ')[0]
      rna_seq=''
      continue
    rna_seq=rna_seq+sk.strip().upper()
  rna_dict[rna_hed]=rna_seq
  
  c_seq=''
  for cline in circ_list:
    clinet=cline.split(',')
    
    if(clinet[1] in rna_dict):
      c_seq=rna_dict[clinet[1]][rna_dict[clinet[1]].find(clinet[-1][:20]):rna_dict[clinet[1]].find(clinet[-1][-10:])]+clinet[-1][-20:]
      #print(clinet[0],'02')
    else:
      c_seq=clinet[-1]
    fo.write('>%s:%s%s\n'%(clinet[0],clinet[2],c_seq))
    fq.write('%s,%d,0-%d,%s\n'%(clinet[0],len(c_seq),len(c_seq),c_seq))
  fi.close()
  fj.close()
  fk.close()
  fo.close()
  fq.close()






#---------------------------------------------
#-------------orf predict---------------------
#---------------------------------------------

def RNA_protein(RNA_string):

  #start_code = 'ATG'
  #end_code = ['UAA', 'UAG', 'UGA']
  protein_table = {'TTT': 'F', 'CTT': 'L', 'ATT': 'I', 'GTT': 'V', \
  'TTC': 'F', 'CTC': 'L', 'ATC': 'I', 'GTC': 'V', \
  'TTA': 'L', 'CTA': 'L', 'ATA': 'I', 'GTA': 'V', \
  'TTG': 'L', 'CTG': 'L', 'ATG': 'M', 'GTG': 'V', \
  'TCT': 'S', 'CCT': 'P', 'ACT': 'T', 'GCT': 'A', \
  'TCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A', \
  'TCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A', \
  'TCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A', \
  'TAT': 'Y', 'CAT': 'H', 'AAT': 'N', 'GAT': 'D', \
  'TAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D', \
  'TAA': '*', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E', \
  'TAG': '*', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E', \
  'TGT': 'C', 'CGT': 'R', 'AGT': 'S', 'GGT': 'G', \
  'TGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G', \
  'TGA': '*', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G', \
  'TGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'}
  
  #start_sit = re.search(start_code, RNA_string)
  protein_list=[]
  index_list=[]
  
  RNA_string2=RNA_string+RNA_string[:3]
  RNA_string3=RNA_string*3
  pm=0
  while(1):
    pm=(RNA_string2).find('ATG')
    #print(pm)
    rp='y'
    if(pm<0):
      break
    for ii in range(0,len(index_list)):
      if((pm+1-index_list[ii])%3==0 and pm+1-index_list[ii]<len(protein_list[ii])*3):
        RNA_string2=RNA_string2[:pm]+'@'+RNA_string2[pm+1:]
        rp='n'
        break
    if(rp=='n'):
      continue
    RNA_string2=RNA_string2[:pm]+'@'+RNA_string2[pm+1:]
    pm=pm+3
    protein='M'
    for sit in range(pm,((len(RNA_string3)-pm)//3)*3+pm, 3):
      protein = protein + protein_table[RNA_string3[sit:sit+3]]
      #print(protein)
      if (protein_table[RNA_string3[sit:sit+3]]=='*'):    
        break
    if(len(protein)>20):
      protein_list.append(protein)
      index_list.append(pm-2) 
  
  return protein_list,index_list
  #RNA_string = open('E:\\bioinfo\data\\rosalind_prot\\rosalind_prot.txt', 'r').read().strip()
  #mRNA_protein(RNA_string)




def orf_predict():
  fi=open(args.output+'/circRNA_exon.fasta','r+')
  fo=open(args.output+'/circRNA_corf.txt','w+')
  print('predict circular RNA orf ,,,')
  circn=''
  seq_exon=''
  for si in fi.readlines():
    if(si[0]=='>'):
      circn=si.strip().split(':')
      #print(circn[0])
      continue
    seq_exon=si.strip()
    if(len(seq_exon)>10000):
      continue
    orf,index=RNA_protein(seq_exon)
    #rint(seq_exon,orf)
    for ii in range(0,len(orf)):
      fo.write(circn[0]+'-ORF%d(%d,%d):%s\n'%(ii+1,index[ii],index[ii]+len(orf[ii])*3,circn[1])+orf[ii]+'\n')   

  fi.close()
  fo.close()




#--------------------------------------------
#----------------main function---------------
#-----------------------------------------------
              
parser=argparse.ArgumentParser()
parser.add_argument('--input','-i',help='statics file be treated.')
parser.add_argument('--sequence','-s',help='circRNA spliced sequnce')
parser.add_argument('--output','-o',help='output file name necessary',default='circRNA_exon_result.fasta')
parser.add_argument('--annotion','-a',help='genomic annotion file to extract exon')
parser.add_argument('--genome','-g',help='genome fasta file to extract sequence')
parser.add_argument('--transcript','-t',help='rna.fa')
args=parser.parse_args()

if(args.sequence!=None):
  fagain2()
elif(args.transcript!=None):
  fagain3()
elif(args.annotion!=None):
  fagain()
else:
  exit(0)

orf_predict()
#os.system('rm corf_predict/exon_infor.txt corf_predict/__monitor__.txt corf_predict/__TEMP_exon.txt')

              
      
        
      
      
