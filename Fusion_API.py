import json
import requests

name1="ESRRA"
name2="Catsperz"

#get number of pages for gene1
URL = str("https://mastermind.genomenon.com/api/v2/articles?api_token={key}&gene="+name1+"&categories[]=fusion&categories[]=breakpoint")
gene1 = requests.get(url=URL)
gene1 = gene1.json()

print(gene1['pages'])
pages1=(gene1['pages'])
y = 1
#list containing PMIDs for gene1
arts1 = []
#iterate through pages
while y < pages1:
    URL = str("https://mastermind.genomenon.com/api/v2/articles?api_token={key}&gene="+name1+"&categories[]=fusion&categories[]=breakpoint" + "&page=" + str(y))
    print(URL)
    gene1 = requests.get(url=URL)
    gene1 = gene1.json()
    j = 0
    #iterate over 5 article limit
    for x in gene1['articles']:
        # print(x)
        arts1.append(gene1['articles'][j]['pmid'])
        j = j + 1
    y=y+1
print(arts1)

#get number of pages for gene2
URL = str("https://mastermind.genomenon.com/api/v2/articles?api_token={key}&gene="+name2+"&categories[]=fusion&categories[]=breakpoint")
gene1 = requests.get(url=URL)
gene1 = gene1.json()

print(gene1['pages'])
pages1=(gene1['pages'])
y = 1
#list containing PMIDs for gene2
arts2 = []
#iterate through pages
while y < pages1:
    URL = str("https://mastermind.genomenon.com/api/v2/articles?api_token={key}&gene="+name2+"&categories[]=fusion&categories[]=breakpoint" + "&page=" + str(y))
    print(URL)
    gene1 = requests.get(url=URL)
    gene1 = gene1.json()
    j = 0
    # iterate over 5 article limit
    for x in gene1['articles']:
        # print(x)
        arts2.append(gene1['articles'][j]['pmid'])
        j = j + 1
    y=y+1
print(arts2)
print(list(set(arts1)&set(arts2)))
#intersect between articles for gene1 and gene2
arts3=list(set(arts1)&set(arts2))
