from collections import Counter
import itertools
from collections import OrderedDict
from itertools import combinations
import json,sys

def jvennstat(prodataseries,pepseqdataseries):
	mrmdatabase={'Cptac':['CPTAC','A'],'Panoramaweb':['PanoramaWeb','B'],'PeptideTracker':['PeptideTracker','C'],'SRMAtlas':['SRMAtlas','D'],'Passel':['PASSEL','E']}
	prodatavalues={}
	pepdatavalues={}
	modmrmdatabase={mrmdatabase[k][1]:mrmdatabase[k][0] for k in prodataseries.keys()}

	for i in range(len(prodataseries)):
		for v in combinations(prodataseries.keys(),i+1):
			vsets = [set(prodataseries[x]) for x in v ]
			sets=tuple(sorted(v))
			sets= list(sets)
			mrmdbname=[mrmdatabase[k][1] for k in sets]
			mrmdbname.sort()
			label=list(reduce(lambda x,y: x.intersection(y), vsets))
			prodatavalues[ ''.join(mrmdbname)]=label

	for j in range(len(pepseqdataseries)):
		for v in combinations(pepseqdataseries.keys(),j+1):
			vsets = [set(pepseqdataseries[x]) for x in v ]
			sets=tuple(sorted(v))
			sets= list(sets)
			mrmdbname=[mrmdatabase[k][1] for k in sets]
			mrmdbname.sort()
			label=list(reduce(lambda x,y: x.intersection(y), vsets))
			pepdatavalues[ ''.join(mrmdbname)]=label

	valueskey=prodatavalues.keys()
	sortedvalueskey=sorted(valueskey, key=len,reverse=True)
	newprodatavalues={}
	newpepdatavalues={}
	for vk in sortedvalueskey:
		try:
			tempvalprot=prodatavalues[vk]
			newprodatavalues[vk]=len(tempvalprot)
			tempvalpep=pepdatavalues[vk]
			newpepdatavalues[vk]=len(tempvalpep)
			prodatavalues.pop(vk, None)
			pepdatavalues.pop(vk, None)
			tmpprodatavalues = dict((k, list(set(v)-set(tempvalprot))) for k, v in prodatavalues.items())
			prodatavalues.update(tmpprodatavalues)

			tmppepdatavalues = dict((k, list(set(v)-set(tempvalpep))) for k, v in pepdatavalues.items())
			pepdatavalues.update(tmppepdatavalues)
		except KeyError:
			pass
	return newprodatavalues,newpepdatavalues,modmrmdatabase
