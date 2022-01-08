#!/usr/bin/python
#*_coding: utf-8_*
__author__='QUMU00'


import argparse
import os
import re
from multiprocessing import Pool
from multiprocessing import Manager
import time
from gnexon import proano
from transla import RNA_protein
from crexon3 import circRNA_extractp
from listsp import listsplit

curPATH=os.path.split(os.path.realpath(__file__))[0]


 
'''
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
'''

def exgain2():
  #proano()
  chrm,exon_str,exon_end,strd=proano(args.annotion,curPATH)
  #fi=open('corf_predict/exon_infor.txt','r+')
  print('acquire exon position of the circRNA...')
  fj=open(args.input,'r+')
  fo=open(curPATH+'/__TEMP_exon.txt','w+')
  fq=open(curPATH+'error.txt','w+')

  mgr=Manager()
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
  
  exoo=mgr.list()
  fq_temp_list=mgr.list()
  fo_temp_list=mgr.list()
  pool=Pool(processes=16)
  
  fj.readline()
  list_sj=listsplit(fj.readlines(),75)
  for list_list_sj in list_sj:
    pool.apply_async(circRNA_extractp,(exoo,fq_temp_list,fo_temp_list,list_list_sj,chrm_exon_dict))
  pool.close()
  pool.join()     
   
  fo.writelines(fo_temp_list)
  fq.writelines(fq_temp_list) 
  fj.close()
  fo.close()
  fq.close()
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
  fp=open(curPATH+'/__monitor__.txt','w+')
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
t1=time.time()              
parser=argparse.ArgumentParser()
group=parser.add_mutually_exclusive_group()
group.add_argument('--input','-i',help='circRNA detailed information file,eg,the chromosome number, start position, host gene,etc.')
group.add_argument('--sequence','-s',help='circRNA spliced sequnce,fa or fasta, used to straightly predict the circular ORFs')
parser.add_argument('--output','-o',help='output file name necessary',default='.')
parser.add_argument('--annotion','-a',help='genomic annotion file of model species to collect exon absolute position information.')
parser.add_argument('--genome','-g',help='genome file, downloaded from NCBI, fna,fa or fasta, to extract sequence')
parser.add_argument('--transcript','-t',help='rna.fa')
args=parser.parse_args()

mode=0
#mode choice, based  on given files,
if(args.sequence!=None):
  fagain2()
  model=3
elif(args.transcript!=None):
  fagain3()
  model=2
elif(args.annotion!=None):
  fagain()
  model=1
else:
  # cannot satify any conditions above, 
  print('some files lacked makes it fail to run.')
  exit(0)

orf_predict()
if(model==1):
  os.system('rm /exon_infor.txt /__monitor__.txt /__TEMP_exon.txt'.replace('/',curPATH+'/'))
t2=time.time()
print('cost time: %f s.'%(t2-t1))

              
      
        
      
      
