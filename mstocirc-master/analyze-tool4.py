#!/usr/bin/python3
 
__Author__='qumu00'
## edited on July 28th by gedit,
 

import os
import time
import argparse
import re



#refer to link https://blog.csdn.net/weixin_37621790/article/details/86682831
##merge sort for sort the large list efficiently 
def merge(arr,l,m,r): 
  n1=m-l+1
  n2=r-m 
  # 创建临时数组
  L=[0]*(n1)
  R=[0]*(n2)
  # 拷贝数据到临时数组 arrays L[] 和 R[] 
  for i in range(0,n1): 
    L[i]=arr[l+i] 
  for j in range(0,n2):    
    R[j]=arr[m+1+j] 
    # 归并临时数组到 arr[l..r] 
  i=0     # 初始化第一个子数组的索引
  j=0     # 初始化第二个子数组的索引
  k=l     # 初始归并子数组的索引
  while (i<n1 and j<n2): 
    if (L[i][3]<=R[j][3]):
    #if(comp(L[i],R[j])): 
      arr[k]=L[i] 
      i+=1
    else: 
      arr[k] = R[j] 
      j+=1
    k+=1
   # 拷贝 L[] 的保留元素
  while (i<n1): 
    arr[k]=L[i] 
    i+=1
    k+=1
 
  # 拷贝 R[] 的保留元素
  while (j<n2): 
    arr[k]=R[j] 
    j+=1
    k+=1
  
def mergeSort(arr,l,r): 
    if (l < r): 
        m = (l+(r-1))//2
        mergeSort(arr, l, m) 
        mergeSort(arr, m+1, r) 
        merge(arr, l, m, r) 

#------------------------------------------------------
#-----------------remove repeat peptide----------------
#------------------------------------------------------

##functions: when there exist same peptide sequence in the peptide, the func can remove the same ones
def rem_peptide():
  print('begin to remove the same peptide sequence...')
  cc=input('whether to remove the repeat peptides(yse or no):')
  #cc='yes'
  if(cc[0]=='Y' or cc[0]=='y'):
    if not os.path.exists(projectn+'/7option'):
      os.system('mkdir '+projectn+'/7option')
    global file_name
    file_name=file_name.split('.')[0]+'_rep.txt'
    fo=open(projectn+'/7option/'+file_name,'w+')
    global temp_file
    temp_file_rep=[]
    #file format sequence  specificty circRNA q-value

    arr=[]
    for line in temp_file:
      line_t=line[:-1].split('\t')
      aa=line_t[0]
      line_t[0]=line_t[3]
      line_t[3]=aa
      arr.append(line_t)
      
    n=len(arr)
    mergeSort(arr,0,n-1)
    aa_seq=''
    for line_t in arr:
      if (aa_seq!=line_t):
        aa_seq=line_t[3]
        line_t[0],line_t[3]=line_t[3],line_t[0]
        temp_file_rep.append('\t'.join(line_t)+'\n')
        
    temp_file=temp_file_rep  
    fo.writelines(temp_file)
    fo.close()  
  print('------------remove_peptide finished!!------------')
  

#------------------------------------------------------
#----------------map onto the cORF---------------------
#------------------------------------------------------

#function: mapp the peptide sequence on circRNA to obtain the circRNAs whose peptides of the intact coRF
def map_corf():
  print('begin run the map_corf,,,')
  if(not os.path.exists(projectn+'/1mapp_corf')):
    os.system('mkdir '+projectn+'/1mapp_corf')
  
  global file_name
  file_name=file_name.split('.')[0]+'_corf.txt'
  #fj=open(args.corf,'r+')
  print('step1: predict circRNA ORF,,,')
  command_line='-o %s'%(projectn+'/1mapp_corf')
  if(args.sequence!=None):
    command_line=command_line+' -s %s'%args.sequence
  elif(args.transcript!=None):
    command_line=command_line+' -i %s -g %s -t %s'%(args.info,args.genome,args.transcript)
  elif(args.annote!=None):
    command_line=command_line+' -i %s -g %s -a %s'%(args.info,args.genome,args.annote)
  else:
    exit(0)
    print('Error')

  #if(args.sequence!=''):
  #  os.system('python3 corf_predict/corf_predict6.py -s %s -o %s'%(args.sequence,projectn+'/1mapp_corf'))
  #else:
  #  os.system('python3 corf_predict/corf_predict4.py -i '+args.info+' -g corf_predict/has_genome.fa -a corf_predict/hsa_genome.gff -o '+projectn+'/1mapp_corf/')
  os.system('python3 corf_predict/corf_predict4.py %s'%command_line)
  fj=open(projectn+'/1mapp_corf/circRNA_corf.txt','r+')
  fo=open(projectn+'/1mapp_corf/'+file_name,'w+')
  
  
  fo.write('sequence\tspecificity\tcircRNA\thost_gene\tcorf\n')
  
  print('step2: map the peptide onto the corf,,,')
  sjky=[]
  sjval=[]
  circ_n=''
  circ_corf=''
  
  for sj in fj.readlines():
    if(sj[0]=='>'):
      circ_nn=sj[1:].split('-ORF')[0]
      if(circ_n==''):
        circ_n=circ_nn
        circ_corf=sj.strip().split(':')[1]
        continue
      if(circ_nn!=circ_n):
        sjval.append(circ_corf)
        sjky.append(circ_n)
        #print(circ_n)
        circ_corf=sj.strip().split(':')[1]
        circ_n=circ_nn
      continue
    circ_corf=circ_corf+','+sj.strip()
  sjval.append(circ_corf)
  sjky.append(circ_n)  

  dict_sj=dict(zip(sjky,sjval))
  global temp_file
  temp_file_corf=[]  
  for si in temp_file:
    sit=si.split('\t')
    circl=sit[2].strip('/').split('/')
    #print(circl)
    
    for i in range(0,len(circl)):

#      llock=0
#      for dict_sj_ky in dict_sj: 
#        if(circl[i].split('_',1)[1] in dict_sj_ky):
#          if(sit[0] in dict_sj[dict_sj_ky]):
#            temp_file_corf.append(sit[0]+'\t'+sit[1]+'\t'+circl[i]+'\t%s\t'%dict_sj_ky.strip().split(':')[1]+dict_sj[dict_sj_ky]+'\n')
#            llock=1
#        elif(llock==1):
#          break
      if circl[i].split('_',1)[1] not in dict_sj:
        continue
      circ_corf=dict_sj[circl[i].split('_',1)[1]].split(',')
      for ii in range(1,len(circ_corf)):
        if(sit[0] in circ_corf[ii]):
          temp_file_corf.append(sit[0]+'\t'+sit[1]+'\t'+circl[i]+'\t%s\t'%circ_corf[0]+circ_corf[ii]+'\n')         
        
  temp_file=temp_file_corf
  print('write into file...')
  fo.writelines(temp_file)
  
  
  fj.close()
  fo.close()
  print('--------finished map corf-----------')
  


#--------------------------------------------------
#------------------------sequence analysis---------
#---------------------------------------------------

def sequence_composition_analysis():
  print('begin ro run seqComp...')
  print('need created...')
  print('.........SeqComp finsihed!!...............')





#------------------------------------------------
#-------------------map onto the back-splice junction-----
#--------------------------------------------------

# function to map peptide on the BSJ OF circRNA 
def map_junct():
  print('begin to map onto the BSJ,,, ')
  if not os.path.exists(projectn+'/2mapp_junct'):
    os.system('mkdir '+projectn+'/2mapp_junct')
  global file_name
  file_name=file_name.split('.')[0]+'_junct.txt'
  fi=open(args.junction,'r+')
  #fj=open(projectn+'/1mapp_corf/peptide_corf.txt','r+')
  fo=open(projectn+'/2mapp_junct/' +file_name,'w+')
  fo.write('BSJ_span\tpeptide\tspecificity\tcircRNA\thost_gene\tcorf_aa\n')
  reh=[]
  req=[]
  for si in fi.readlines():
    if(si[0]=='>'):
      sih=si[1:].strip()
      reh.append(sih)
      continue
    req.append(si.strip())
  print('create the dict() for junction.fasta...')  
  dict_junt=dict(zip(reh,req))
  global temp_file
  temp_file_junt=[]
  #fj.readline()
  #temp_file=fj.readlines()
  
  print('judge the spanning amino acid amount...')  
  for str1 in temp_file:
    ss=str1.split('\t')
    lgh=len(ss[0])
    
    rfq=dict_junt[ss[2]]
    pp=rfq.find(ss[0])
    
    
    if(pp>=0):
      lgh2=len(rfq)
      if(pp<lgh2/2 and pp+lgh>lgh2/2):
         #if(pp+qq<lgh2/2 and pp+lgh-qq>lgh2/2):
         if(lgh2/2-pp > pp+lgh-lgh2/2):
           temp_file_junt.append('YY%d\t'%(pp+lgh-lgh2/2)+str1)  
         else:
           temp_file_junt.append('YY%d\t'%(lgh2/2-pp)+str1)
           #print ss[1]+' '+pr[j]
  temp_file=temp_file_junt

  print('write into files...')
  fo.writelines(temp_file_junt)
    
  fi.close()
  fo.close()
  print('------------map_junct finished-------------')



#------------------------------------------------------
#------------------------------------------------------
#-----------------------------------------------------
def map_gene():
  print('begin to run map_gene...')
  cc=input('whther map the peptides onto the linear protein(yes or no): ')
  #cc='no'
  if(cc[0]=='y' or cc[0]=='Y'):
    if( not os.path.exists(projectn+'/7option')):
      os.makedirs(projectn+'/7option')
    global temp_file
    temp_file_gen=[]
    global file_name
    #file_nn=input('please input the linear protein file: ')
    file_nn='GCF_000001735.4_TAIR10.1_protein.faa'
    fi=open(file_nn,'r+')
    file_name=file_name.split('.')[0]+'_gen.txt'
    fo=open(projectn+'/7option/'+file_name,'w+')
    protein=[]
    proteinn=''
    print('extract linear proteins...')
    for si in fi.readlines():
      if(si[0]=='>' and protein!=''):
        protein.append(proteinn)
        proteinn=''
        continue
      proteinn=proteinn+si.strip()
    protein.append(proteinn)
    
    print('remove peptides...')
    for line in temp_file:
      lint=line.split('\t')
      wr='y'
      for proteinn in protein:
        if(lint[1] in proteinn):
          wr='n'
          #print(line)
          break
      if(wr=='y'):
        temp_file_gen.append(line)
    
    print('write result into files...')
    fo.writelines(temp_file_gen)
    temp_file=temp_file_gen
    
    fi.close()
    fo.close() 

  print('----------------map_gene finished!!--------------------')   








#-------------------------------------------------------
#---------------------merge the overlappinng peptide--------
#--------------------------------------------------------

def peptide_merge_lr(pep_list):
  #print(pep_list)
  line_mgg=[]
  while(pep_list!=[]):
    pep_list2=[]
    pep_list3=[]
    corf=pep_list[0][-1]
    for jj in range(0,len(pep_list)):
      if (pep_list[jj][-1]==corf):
        pep_list2.append(pep_list[jj])
      else:
        pep_list3.append(pep_list[jj])
    pep_list=pep_list3
    line_mg=['','','','','','']
    for jj in range(3,len(pep_list2[0])):
      line_mg[jj]=pep_list2[0][jj]

    ll=pep_list2[0][-1].find(pep_list2[0][1])
    rr=pep_list2[0][-1].find(pep_list2[0][1])+len(pep_list2[0][1])
    
    for line in pep_list2:
      for jj in range(0,3):
        line_mg[jj]=line_mg[jj]+','+line[jj]
    
      if(line[1] in line[-1][ll:rr]):
         continue
      elif(line[-1][ll:rr] in line[1]):
        ll=line[-1].find(line[1])
        rr=line[-1].find(line[1])+len(line[1])
       
      else:
        if(ll>line[-1].find(line[1])):
          ll=line[-1].find(line[1])
        if(rr<line[-1].find(line[1])+len(line[1])):
          rr=line[-1].find(line[1])+len(line[1])
        
          
    line_mg.insert(2,line_mg[-1][ll:rr])
    line_mgg.append(line_mg)
  return line_mgg

#function: map on the corf and merge the overlaping peptide to the longest and single peptide od spanning the BSJ 
def peptide_merge():
  print('begin run the merge,,,')
  if not os.path.exists(projectn+'/3peptide_merge'):
    os.makedirs(projectn+'/3peptide_merge')
  
  global file_name
  file_name=file_name.split('.')[0]+'_exd.txt'
  fo=open(projectn+'/3peptide_merge/'+file_name,'w+')
  #fi=open(projectn+'/2mapp_junct/peptide_corf_junct.txt','r+')
  fo.write('BSJ_span\tpeptide\tmerge_peptide\tspecificity\tcircRNA\thost_gene\tcorf\n')
  global temp_file
  temp_file_merge=[]
  arr=[]
  #fi.readline()
  #temp_file=fi.readlines()
  for si in temp_file:
    sit=si.split('\t')
    arr.append(sit) 
  n=len(arr)
  print('sort the peptides by circRNA names...')
  mergeSort(arr,0,n-1)
  circn=''
  #print(arr)
  pep_list=[]
  for i in range(0,n):
    if(arr[i][3]!=circn):
      if(circn!=''):  
        line_meg=peptide_merge_lr(pep_list)
        for line_megg in line_meg:
          temp_file_merge.append('\t'.join(line_megg))
      
      pep_list=[]
      circn=arr[i][3]
    pep_list.append(arr[i])
  
  line_meg=peptide_merge_lr(pep_list)
  for line_megg in line_meg:
    temp_file_merge.append('\t'.join(line_megg))
  print('write into file...')
  fo.writelines(temp_file_merge)
  temp_file=temp_file_merge
  #fi.close()
  fo.close()
  
  
  print('------------peptide_merge finished---------')


  
#--------------------------------------------------------
#-------------------predict the IRES elements------------
#--------------------------------------------------------

#function:to predict the ires

try:
  from sklearn.linear_model import LogisticRegressionCV
except ImportError:
  print('sklearn module not exists! mstocirc trys to install it for python. ')
  os.system('pip3 install -U scikit-learn')
  

def ires_predict():
  print('begin to run ires_predict,,,')
  if not os.path.exists(projectn+'/4ires_predict'):
    os.system('mkdir '+projectn+'/4ires_predict')
  
  global file_name
  #file_name='PT6415_18_corf_junct_exd.txt'
  file_name=file_name.split('.')[0]+'_ires.txt'
  global temp_file
  temp_file_ires=[]
  #fi=open(projectn+'/3peptide_merge/PT6415_18_corf_junct_exd.txt','r+')
  #fj=open(args.info,'r+')
  fj=open(projectn+'/1mapp_corf/circRNA_exon.fasta','r+')
  fo=open(projectn+'/4ires_predict/'+file_name,'w+')
  
  fo.write('BSJ_span\tpeptide\tmerge_peptide\tspecificty\tcircRNA\thost_gene\tIRES_element\tcorf\n')
  file_num=1
  fp=open(projectn+'/4ires_predict/circRNA_ires_%d.fasta'%file_num,'w+')
  #fi.readline()
  #temp_file=fi.readlines()
  print('step1:extract sequence...')

  seq_acount=1
  seq_dict=dict()
  circ_name=''
  for sj in fj.readlines():
    if(sj[0]=='>'):
      circ_name=sj[1:].split(':')[0]
      #print(circ_name)
      continue
    seq_dict[circ_name]=sj.strip()
  #print(seq_dict)
  
  for si in temp_file:
    sit=si.split('\t')
    #print(si)
    
    #fj.seek(0,0)
    #wr='n'
    #for sj in fj.readlines():
      #print(sit[4])
    #  if(sj[0]=='>' and sit[4].split('_',1)[1] in sj ):
    #    if(seq_acount>500):
    #      fp.close()
    #      seq_acount=1
    #      file_num=file_num+1
    #      fp=open(projectn+'/4ires_predict/circRNA_ires_%d.fasta'%file_num,'w+')
    #    wr='y'
    #    fp.write(sj)
    #    continue
    #  if(wr=='y'):
    #    fp.write(sj)
    #    seq_acount=seq_acount+1
    #    break

    
    if(seq_acount>200):
      fp.close()
      seq_acount=1
      file_num=file_num+1
      fp=open(projectn+'/4ires_predict/circRNA_ires_%d.fasta'%file_num,'w+')
    else:
      fp.write('>%s\n'%sit[4].split('_',1)[1])
      fp.write(seq_dict[sit[4].split('_',1)[1]][:300]+'\n')
      seq_acount=seq_acount+1
  fp.close()

  #os.system('python3 ./exon_extract/corf_predict4.py -g exon_extract/TAIR10.sequence.fa -a exon_extract/Arabidopsis_thaliana.TAIR10.47.gtf -i exon_extract/circRNA_ires.txt -o ./')
  print('step2:ires predict...')
  #print(file_num)
  for ii in range(0,file_num):
    #print('haha...ires')
    os.system('python3 IRESfinder-master/IRESfinder.py -f '+projectn+'/4ires_predict/circRNA_ires_%d.fasta -o '%(ii+1)+projectn+'/4ires_predict/circRNA_ires_%d.res -m 2 -w 50 -s 50'%(ii+1))
  os.system('cat %s/4ires_predict/circRNA_ires*.res >%s/4ires_predict/circRNA_ires.res'%(projectn,projectn))
  os.system('cat %s/4ires_predict/circRNA_ires*.fasta >%s/4ires_predict/circRNA_ires.fasta'%(projectn,projectn))
  os.system('rm %s/4ires_predict/circRNA_ires_*'%projectn)
  fk=open(projectn+'/4ires_predict/circRNA_ires.res','r+')
  fl=open(projectn+'/1mapp_corf/circ_draw.txt','r+')
  sfl=fl.readlines()
  
  ires_dict=dict()
  circ_n=''
  ires_r=''
  ires_lock=1
  for sk in fk.readlines():
    if('IRES' not in sk):
      continue
    skt=sk.split('\t')
    sktt=skt[0][::-1].split('_',2)
    #print(sktt)
    if(circ_n!=sktt[2][::-1]):
      ires_lock=1
      circ_n=sktt[2][::-1]
    if(skt[1]=='IRES' and ires_lock==1):
      ires_r='%s,%s,%s'%(sktt[1][::-1],sktt[0][::-1],','.join(skt[1:]))
      ires_dict[circ_n]=ires_r
      ires_lock=0

  print('step3: write IRES prediction into file... ')
  #sfk=fk.readlines()
  temp_file_fl=[]
  fl.close()
  fl=open(projectn+'/4ires_predict/circ_draw.txt','w')
  sl_dict=dict()
  for sl in sfl:
    slt=sl.split(',',1)
    sl_dict[slt[0]]=slt[1]
  for si in temp_file:
    sit=si.split('\t')
    ires='non-ires'
    ires_s='0'
    ires_e='0'
    #for sk in sfk:
    #  if(sit[4].split('_',1)[1] in sk and '\tIRES' in sk):
    #    ires='ires'
    #    ires_s=sk.split('\t')[0].split('_')[-2]
    #    ires_e=sk.split('\t')[0].split('_')[-1]
        #print(sk,ires_s)
        #break
    
    #    for sl in sfl:
    #      if sit[4].split('_',1)[1] in sl:
    #        slt=sl.strip().split(',')
    #        slt.append(sit[-1])
    #        slt.insert(3,'NUll')
    #        slt.insert(2,ires_e)
    #        slt.insert(2,ires_s)
            #slt.insert(2,'20,'+str(len(slt[-2])+20))
    #        temp_file_fl.append(','.join(slt))
    #    break

    if(sit[4].split('_',1)[1] in ires_dict):
      ires='ires'
      ires_rt=ires_dict[sit[4].split('_',1)[1]].split(',')
      ires_s=ires_rt[0]
      ires_e=ires_rt[1]
    sl_dt=sl_dict[sit[4].split('_',1)[1]].strip().split(',')
    sl_dt.append(sit[2]+'-'+sit[-1])
    sl_dt.insert(2,'null')
    sl_dt.insert(1,ires_e)
    sl_dt.insert(1,ires_s)
    temp_file_fl.append(sit[4].split('_',1)[1]+','+','.join(sl_dt))
    sit.insert(6,ires) 
    temp_file_ires.append('\t'.join(sit))
    
  fo.writelines(temp_file_ires)
  fl.writelines(temp_file_fl)
  temp_file=temp_file_ires 
  
  print('-----------IRES analyse finished!!----------------')
  #fi.close()
  fj.close()
  fk.close()
  fo.close()
  fp.close()    

#------------------------------------------------------
#-----------------------------------------------------
def m6a_modification_predict():
  print('begin to run m6a predict...')
  print('need to be created!!')
  print('---------------m6a finished!!...')




#------------------------------------------------------
#------------------------------------------------------
def path_analysis():
  print('begin to run path_analysis')
  #cc=input('R must be installed,whether to RUN(yes or no):')
  cc='yes'
  if(cc[0]=='Y' or cc[0]=='y'):
    if not os.path.exists(projectn+'/5enrich'):
      os.makedirs(projectn+'/5enrich')
    fo=open(projectn+'/5enrich/query3.tsv','w+')
    global temp_file
    if(len(temp_file)<300):
       print('circRNA less than <300, skip this part.')
       print('--------------enrichment finished!!-------------')
       return 
    fo.write('BSJ_span\tpeptide\tmerge_peptide\tspecificty\tcircRNA\thost_gene\tIRES_element\n')
    for line in temp_file:
      fo.write('\t'.join(line.split('\t')[:7])+'\n')
    fo.close()
    os.system('Rscript rpathway/pathways.R %s/5enrich'%projectn)
    
  print('--------------enrichment finished!!---------------')



#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------
def ms_ribo():
  print('begin to run ms_ribo...')
  cc=input('need provide ribo evidence in .fasta, whether to run: ')
  #cc='no'
  if(cc[0]=='Y' or cc[0]=='y'):
    if not os.path.exists(projectn+'/7option'):
      os.makedirs(projectn+'/7option')
    global temp_file
    global file_name
    temp_file_ribo=[]
    file_name=file_name.split('.')[0]+'_ribo.txt'
    print('create 2 files for result saving...')
    fo=open(projectn+'/7option/'+file_name,'w+')
    fp=open(projectn+'/7option/'+file_name.replace('_ribo','_ribo2'),'w+')
    #file_ribo=input('please input the ribo ribo-sse evidence file(.faa):')
    file_ribo='ath_pep.fa.faa'
    fj=open(file_ribo,'r+')
    print('loading files over...')
    fm=open(args.junction,'r+')
    
    fo.write('BSJ_span\tpeptide\tmerge_peptide\tspecificty\tcircRNA\thost_gene\tIRES_element\tRibo_evidence\tcorf\n')
    if not os.path.exists(projectn+'/7option'):
      os.makedirs(projectn+'/7option')
    
     
    sjl=fj.readlines()
    temp_ribo=[[],[],[]]
    dict_rb=dict()
    for si in temp_file:
      sit=si.split('\t')
      llk=0
      ribo='N'
      for sj in sjl:
        if(sit[4].split('_',1)[1] in sj):
          temp_ribo[0].append(sit[4].split('_',1)[1]+'\n')
          ribo='Y'
          llk=1
          continue
        if(llk==1 and sj.strip() in sit[-1]):
          temp_ribo[1].append(sit[4].split('_',1)[1]+'\n')
          
      #if(sit[2] in sj.strip()):
      #  temp_3.append(sit[4].split('_',1)[1]+'\n')
          if sit[4].split('_',1)[1] not in dict_rb:
            dict_rb[sit[4].split('_',1)[1]]=sj.strip()
          else:
            dict_rb[sit[4].split('_',1)[1]]=dict_rb[sit[4].split('_',1)[1]]+','+sj.strip()
          llk=0
      sit.insert(-1,ribo)
      temp_file_ribo.append('\t'.join(sit))
      
    for sj_ky in dict_rb:
      #print(dict_rb)
      lllck=0
      fm.seek(0,0)
      while(True):
        sm=fm.readline().strip()
        llllck=1
        if(sm==''):
          break
    
        if(sj_ky in sm):
          lllck=1
          #print(sj_ky,sm)
          sm=fm.readline().strip()
          md=len(sm)//2
      
      #print(sm[md-2:md+2],'@'+dict_rb[sj_ky])
          if(sm[md-2:md+2] in dict_rb[sj_ky]):
            temp_ribo[2].append(sj_ky+'\n')
          else:
            llllck=0
          if(lllck==1):
            if(llllck==0):
              break
    fo.writelines(temp_file_ribo)
    temp_file=temp_file_ribo
    for i in range(0,3):
      fp.writelines(temp_ribo[i])
      fp.write('\n---------\n')

  print('............ms_ribo finished!!........')
  


#------------------------------------------------------
#-----------------------------------------------------
#------------------------------------------------------
def circ_annote():
  print('begin to run circ_annote...')
  cc=input('whether to annote the circRNA(yes or no):')
  #cc='no'
  if(cc.lower()[0]=='y'):
    if not os.path.exists(projectn+'/7option'):
      os.system('mkdir '+projectn+'/7option')
    global file_name
    file_name=file_name.split('.')[0]+'_anno.txt'
    fo=open(projectn+'/7option/'+file_name,'w+')
    #file_p=input('input annotion.gtf name:')
    file_p='Araport11_functional_descriptions_20190930.txt'
    fi=open(file_p,'r+')
    temp_file_si=fi.readlines()
    fo.write('BSJ_span\tpeptide\tmerge_peptide\tspecifity\tcircRNA\thost_gene\tIRES\tcorf\tfunction\n')
    global temp_file
    temp_file_anno=[]
    for sj in temp_file:
      sjt=sj.strip().split('\t')
      function='null\n'
      for si in temp_file_si:
        sit=si.split('\t')
        if(sjt[5] in sit[0]):
          function=sit[2]+'\n'
          break
      sjt.append(function)
      temp_file_anno.append('\t'.join(sjt))    
    
    temp_file=temp_file_anno  
    fo.writelines(temp_file_anno)
    fi.close()
    fo.close()
      
  print('------------circ_annote finished!!-----------------------')





  
def circ_corf_info(strr):
  print('extract the....')
  fi=open(projectn+'/4ires_predict/circ_draw.txt','r+')
  fj=open(projectn+'/1mapp_corf/circRNA_corf.txt','r+')
  fo=open(projectn+'%s/circ_draw.txt'%strr,'w+')
  #temp_file_fj=fj.readlines()
  sj_dict=dict()
  orf_po=''
  circ_n=''
  
  for sj in fj.readlines():
    if(sj[0]=='>'):
      orf_po=re.findall('\d+',sj.split(':',1)[0])[-2:]
      circ_n=sj[1:].split('-')[0]
      continue
    
    sj_dict[circ_n+','+sj]=orf_po
  #print(sj_dict)
  for si in fi.readlines():
    sit=si.split(',')
  
    #dg='n'
    #orf_p=''
    #print('0908')
    #for sj in temp_file_fj:
      #print(sit,sit[0],sj)
    #  if sit[0] in sj:
    #    dg='y'
        #print('0909')
        #print(sj)
    #    orf_p=re.findall('\d+',sj.split(':',1)[0])[-2:]
        
    #    continue
    #  if(dg=='y'):
        #print(sit[-1],sj)
    #    if(sit[-1]==sj):
          
    #      sit.insert(2,orf_p[1])
    #      sit.insert(2,orf_p[0])
    #      break
    #    dg='n' 
    orf_po=sj_dict[sit[0]+','+sit[-1].split('-',1)[-1]]
    sit.insert(2,orf_po[1])
    sit.insert(2,orf_po[0]) 
    fo.write(','.join(sit))
  
  fo.close()
  fi.close()
  fj.close()
  #os.system('drawcirc/drawcirc.py '+projectn+'/6draw_circ/circ_draw.txt '+projectn+'/6draw_circ/')
  
  
  


#-------------------------------------------------------
#-------------------------------------------------------
#-------------------------------------------------------
def circ_classify():
  print('begin to run circ_classify...')
  cc=input('whether to classify circRNA(yes or no):')
  #cc='yes'
  if(cc.lower()[0]=='y'):
    #print('need created...')
    if(not os.path.exists(projectn+'/7option')):
      os.makedirs(projectn+'/7option')
    global file_name
    #file_name='PT6415_18_corf_junct_exd_ires.txt'
    file_name=file_name.split('.')[0]+'_clf.txt'
    circ_corf_info('/7option')
    fi=open(projectn+'/7option/circ_draw.txt','r+')
    #fj=open(projectn+'/4ires_predict/PT6415_18_corf_junct_exd_ires.txt','r+')
    fo=open(projectn+'/7option/circ_draw_clf.txt','w+')
    fp=open(projectn+'/7option/'+file_name,'w+')
    #fq=open(projectn+'/7option/error.txt','w+')
    global temp_file
    #fj.readline()
    #temp_file=fj.readlines()
    temp_file_clf=[[],[],[],[],[],[]]
    temp_file_clf2=[[],[],[],[],[],[]]
    sil=fi.readlines()
    
    for si in sil:
      sit=si.split(',')
      i,j=0,0
      if((int(sit[3])-int(sit[2]))//int(sit[1])<1):
        i=0
      elif((int(sit[3])-int(sit[2]))//int(sit[1])<2):
        i=1
      else:
        i=2
      
      #if(int(sit[2])<180):
      if(sit[5]!='0'):
        j=0
      else:
        j=1
      temp_file_clf[2*i+j].append(si)
      for sj in temp_file:
        if(sit[0] in sj and sit[-1].strip().split('-',1)[-1] in sj):
          #fq.write('%d '%(i*2+j)+' '+sit[0]+' '+sj)
          temp_file_clf2[2*i+j].append(sj)
          print('klkl')
          break
          #print(2*i+j,sj)
    temp_file=[]
    print('write into files...')     
    for ii in range(0,6):
      fo.write('##00%d\n'%(ii+1))
      fo.writelines(temp_file_clf[ii])
      temp_file=temp_file+temp_file_clf2[ii]
      fp.writelines(temp_file_clf2[ii])
   
    fi.close()
    fo.close()
    fp.close()
        
  print('--------circ_classify finished!!-----')




#---------------------------------------------------
#---------------------------------------------------
#----------------------------------------------------

try:
  import matplotlib.pyplot as pltt
except ImportError:
  print ('module matplotlib not exists. mstocirc trys to install it for python3.')
  os.system('pip3 install -U matplotlib')

def draw_circ():
  print('begin to draw circRNAs....')
  if not os.path.exists(projectn+'/6draw_circ'):
    os.makedirs(projectn+'/6draw_circ')
  print('create 6draw_circ file...')
  print('gain the drawing paramters files...')
  if not os.path.exists(projectn+'/7option/circ_draw.txt'):
    circ_corf_info('/6draw_circ')
  else:
    os.system('cp %s/7option/circ_draw_clf.txt %s/6draw_circ/circ_draw.txt'%(projectn,projectn))
  
  os.system('python3 drawcirc/drawcirc.py '+projectn+'/6draw_circ/circ_draw.txt '+projectn+'/6draw_circ/')
  
  print('-----------draw_circ finished!!-------------------------')  





#--------------------------------------------------------------
#---------------------------------------------------------------
#=================main function================================
#-----------------------------------------------------------------

time_str=time.time()
now=time.strftime('%y-%m-%d %H:%M:%S',time.localtime())
print('Data: %s'%now)
ii=1
projectn='project'+now.split()[0]
while(os.path.exists(projectn+'.%d'%ii)):
  ii=ii+1
projectn=projectn+'.%d'%ii
os.system('mkdir '+projectn)




parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
parser.add_argument('-p','--peptide',help="peptide file generated by pfind,amino acide",default=None,required=True)
group.add_argument('-s','--sequence',help="circRNA putitive spliced sequence file",default=None)
#parser.add_argument('-c','--corf',help="circRNA-derived cORF sequqnce,amino acide",default=None)
parser.add_argument('-j','--junction',help="the BSJ reference of circRNA,amino acide",default=None)
group.add_argument('-i','--info',help="info file ,including the name,position etc. of circRNA",default=None)
parser.add_argument("--output", "-o", "-O", help="Output directory", default='./', required=False)
parser.add_argument('--annote','-a',help='genomic.gtf',required=False)
parser.add_argument('--genome','-g',help='genome.fna',required=False)
parser.add_argument('--transcript','-r',help='rna.fna',required=False)
args = parser.parse_args()

#-----------------------------
#judge the files exists
print('files check,wait...')
file_list=[]
if(not os.path.exists(args.peptide)):
  file_list.append('peptide')
if(args.info!=None ):
  if(not os.path.exists(args.info)):
    file_list.append('info')
if(args.sequence!=None):
  if(not os.path.exists(args.sequence)):
    file_list.append('sequence')
if(not os.path.exists(args.junction)):
  file_list.append('junction')
if(file_list!=[]):
  print('file %s not exist,please try again...'%'/'.join(file_list))
  exit(1)
print('analyse will begin soon,,,')
#---------------------------------------------


fii=open(args.peptide,'r+')
fii.readline()
temp_file=fii.readlines()
fii.close()
file_name=args.peptide.split('/')[-1]


#function part
#------------------------------------------- 

#001
rem_peptide()
#002
#map_corf(temp_file)
#why other can change the value of paramater,but I can't do that, 
map_corf()
#003
sequence_composition_analysis()
#004
map_junct()
#
map_gene()
#005
peptide_merge()
#006
ires_predict()
#m6a_modification_predict()
path_analysis()
#007
ms_ribo()
#008
circ_annote()
#009
circ_classify()
#010
draw_circ()

#-----------------------------------------------------------------



print('generate the result files...')
if not os.path.exists('%s/result'%projectn):
  os.system('mkdir %s/result'%projectn)
os.system('cp %s/1mapp_corf/circRNA_corf.txt %s/result'%(projectn,projectn))
#os.system('cp -r /4ires_predict /result'.replace(projectn,'/'))
os.system('cp %s/1mapp_corf/circRNA_exon.fasta %s/result'%(projectn,projectn))
os.system('cp %s/4ires_predict/circRNA_ires.res %s/result'%(projectn,projectn))
os.system('cp %s/6draw_circ/circ_draw.txt %s/result'%(projectn,projectn))
os.system('cp -r %s/5enrich %s/result/'%(projectn,projectn))
os.system('cp -r %s/6draw_circ %s/result/'%(projectn,projectn))
if('_ires.txt' in file_name):
  os.system('cp %s/4ires_predict/%s %s/result'%(projectn,file_name,projectn))
else:
  os.system('cp %s/7option/%s %s/result'%(projectn,file_name,projectn))



tempn=' /1mapp_corf /2mapp_junct /3peptide_merge /4ires_predict /6draw_circ /7option' 
# save the predix of temporal file directory
tempn.replace(projectn+'/','/')
cc=input('delete the temporary file(yes or no??):')
#cc='no'
if(cc[0]=='y' or cc[0]=='Y'):
  os.system('rm -r'+tempn)


time_end=time.time()

print('write log.txt...')
fn=open(projectn+'/log.txt','w+')
fn.write('-----------------------\nrun time: %s\n'%now)
fn.write('\n------------------------------------------------\n')
fn.write('/                                               \\\n')
fn.write('/           analyse-tool                        \\\n')
fn.write('/                                               \\\n')
fn.write('/-----------------------------------------------\\\n')
fn.write('\nmain funtion:\n\n')
fn.write('7option function:\n\n')
fn.write('total cost time: %d seconds\n'%(time_end-time_str))
fn.write('\nparameter sets,\n  peptide file name: %s\n'%args.peptide)
fn.write('  information file name: %s\n'%args.info)
fn.write('\n\n---------------------------------------\n')
fn.close()

os.system('mv ./%s %s'%(projectn,args.output))
print('\ntotal cost time: %dh %dmin %ds.'%((time_end-time_str)//3600,((time_end-time_str)%3600)//60,(time_end-time_str)%60))
print('----------------thankyou!!----------------')
print('===============analyse finished==========')


