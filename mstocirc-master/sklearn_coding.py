#!/usr/bin/python3
__Author__='QUMU00'
import sys
import random
from sklearn.model_selection import train_test_split
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.naive_bayes import MultinomialNB
from sklearn.metrics import classification_report
from sklearn.metrics import accuracy_score



amino_list=['I','M','T','N','K','S','R','L','P','H','Q','V','A','D','E','G','F','Y','*','C','W']
def pepvitual():
  pep='M'
  while(True):
    nnn=random.randint(0,20)
    #print(nnn)
    pep=pep+amino_list[nnn]
    if(pep[-1]=='*'):
      if(len(pep)<100 or len(pep)>2000):
        pep='M'
      else:
        return pep[:-1]
   


def train_data():
  x = []
  y = []
  fi=open(sys.argv[1],'r+')
  pep_list2=[]
  pep2=''
  for si in fi.readlines():
    if(si[0]=='>'):
      if(pep2!=''):
        pep_list2.append(pep2)
        pep2=''
      continue
    pep2=pep2+si[:-1]
  pep_list2.append(pep2)
  fi.close()

  nn=0
  while(True):
    if(random.randint(0,1)==0):
      x.append(pep_list2[nn])
      y.append(1)
      nn=nn+1
    else:
      x.append(pepvitual())
      y.append(0)
    if(nn>=len(pep_list2)):
      break
    print(x[-1])
  return x,y


def model_train(strx):
  x,y=train_data()
  print('train model...')
  print(y)
  
  #this function is refer tp link,https://blog.csdn.net/weixin_43685844/article/details/88064485
  #we are grate to hin for sharinng thses codes,
  x_train,x_test,y_train,y_test=train_test_split(x, y, test_size = 0.2, random_state = 1)
  vect=CountVectorizer(analyzer='char_wb',ngram_range=(4,4))
  vect.fit(x_train)
  x_train_df=vect.transform(x_train)
  x_test_df=vect.transform(x_test)

  strx_df=vect.transform(strx)
  
  model=MultinomialNB()
  model.fit(x_train_df,y_train)
  ytest_predict=model.predict(x_test_df)
  Accuracyy=accuracy_score(y_test,ytest_predict)
  print('model.selection accuracu: %f%%'%(Accuracyy*100))
  #print(classification_report(y_test,NB_pred,target_names=types))
  
  resy=model.predict(strx_df)
  return resy



#------------------------------------------------------------------
#-------------------main function----------------------------------
#------------------------------------------------------------------

fj=open(sys.argv[2],'r+')
peps=[]
lineh=fj.readline()
temp_file_skl=fj.readlines()
fj.close()

for sj in temp_file_skl:
  sjt=sj[:-1].split('\t')
  if(sjt[-1][-1]=='*'):
    sjt[-1]=sjt[-1][:-1]
  peps.append(sjt[-1])
resyy=model_train(peps)

nn=0
temp_file_skl2=[]
for sj in temp_file_skl:
  sjt=sj.split('\t')
  sjt.insert(-1,str(resyy[nn]))
  nn=nn+1
  temp_file_skl2.append('\t'.join(sjt))

fo=open(sys.argv[3],'w+')
lineh.split('\t').insert(-1,'coding_potential')
fo.write('\t'.join(lineh))
fo.writelines(temp_file_skl2)

















