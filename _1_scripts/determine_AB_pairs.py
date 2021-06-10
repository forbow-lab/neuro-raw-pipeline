#!/bin/bash

## - list all subjects with both A and B scan sessions 
## will need to filter these to ensure each session has all necessary scans

## BASH:
for D in `ls -d ???_A_20??????`; do 
	S=${D:0:3}; 
	B=$(find . -maxdepth 1 -type d -name "${S}_B_20??????" -print); 
	if [ -d "$B" ]; then 
		echo "$S"; 
	fi;
done

exit 1


## PYTHON:
SubjList='_3_docs/Subjlist_AB_20210309.txt'
D=open(,'r').read().split()
SSID=[s[0:5] for s in D]
for a in SSID:
	if a[-1:] == 'A':
		for b in SSID:
			if b[-1:] == 'B' and b[0:3] == a[0:3]:
				print(s)

##
013
014
015
019
023
024
025
026
027
030
034
035
036
037
041
042
057
063
065
067
068
069
072
074
075
077
078
079
082
083
085
087
097
098
099
100
103
107
117
128
133
136
138
139
147
151
153
154
157
161
164
169
172
173
176
180
183
187
188
189
190
201
202
212
213