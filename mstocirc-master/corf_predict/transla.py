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


