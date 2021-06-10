#!/usr/bin/env python

import os,sys
from glob import glob

WorkDir=os.getcwd()
os.chdir('/shared/uher/FORBOW/rawdata/')

relSubjs=[]

for t2 in ['B','D','F']:
	r2 = glob('[0-1]??_%s_20*'%(t2))
	print '\n%s:'%(t2), r2
	if t2=='B': t1='A'
	elif t2=='D': t1='C'
	elif t2=='F': t1='E'
	for s in r2:
		id = s[0:3]
		m = glob('%s_%s_20*'%(id,t1))
		if len(m)==1:
			#print s,m
			relSubjs.append(s[0:5])
			relSubjs.append(m[0][0:5])

print 'reliability-IDs:',relSubjs
outF=open('all_reliability_subjects.txt','w')
for S in relSubjs:
	outF.write('%s\n'%(S))
outF.close()



# E=/shared/uher/FORBOW/FSv6_analysis/
# cd $E/
# for S in `cat -v ../rawdata/all_reliability_subjects.txt`; do if [ ! -d "${S}_FLAIR" ]; then echo "$S" >>20190121_new_reliability_subjs.txt; fi ; done
# for S in `cat -v 20190121_new_reliability_subjs.txt`; do if [ -d "/shared/uher/FORBOW/analysis/${S}_FLAIR" ]; then mkdir -p ./${S}_FLAIR/unprocessed/; cp -rpv ../analysis/${S}_FLAIR/unprocessed/T?w1/ ./${S}_FLAIR/unprocessed/; echo "${S}_FLAIR" >>2019_new_reliability_FLAIRs.txt; fi ; done

# E=/shared/uher/FORBOW/NIHPD_analysis/
# cd $E/
# for S in `cat -v ../rawdata/all_reliability_subjects.txt`; do if [ ! -d "${S}_FLAIR" ]; then echo "$S" >>20190121_new_reliability_subjs.txt; fi ; done
# for S in `cat -v 20190121_new_reliability_subjs.txt`; do if [ -d "/shared/uher/FORBOW/analysis/${S}_FLAIR" ]; then mkdir -p ./${S}_FLAIR/unprocessed/; cp -rpv ../analysis/${S}_FLAIR/unprocessed/T?w1/ ./${S}_FLAIR/unprocessed/; echo "${S}_FLAIR" >>2019_new_reliability_FLAIRs.txt; fi ; done


#cd /shared/uher/FORBOW/NIHPD_analysis/
#/bin/ls -d ???_? ???_?_NP ???_?_FLAIR >>20190201_subjlist.txt
#cd /shared/uher/FORBOW/analysis/
#nList=""
#for D in `/bin/ls -d ???_? ???_?_NP ???_?_FLAIR`; do if [ "`grep $D ../NIHPD_analysis/20190201_subjlist.txt`" == "$D" ]; then nList="${nList} $D"; fi; done
#echo "$nList" >>20190201_new_FSv6_datasets.txt
#for S in `cat -v 20190201_new_FSv6_datasets.txt`; do echo "$S"; mkdir -p ../NIHPD_analysis/$D/unprocessed/; cp -rpv $D/unprocessed/T?w1 ../NIHPD_analysis/$D/unprocessed/;  done

