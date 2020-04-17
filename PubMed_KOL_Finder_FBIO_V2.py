from Bio import Entrez
from Bio import Medline
from tqdm import tqdm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

# Change this email to your email address
Entrez.email = "jwang@fortressbiotech.com"

disease_list = [input('Please enter a search term: ')]

#search and return total number of publications 
def search(x):
    Entrez.email=Entrez.email
    for x in disease_list:
        keyword = x
        handle = Entrez.esearch(db ='pubmed',
                                retmax=10,
                                retmode ='text',
                                term = keyword)
        results= Entrez.read(handle)
        print('Total number of publications that contain the term {}: {}'.format(keyword, results['Count']))    
    return results

def batch(x):

    keyword = disease_list
    handle = Entrez.esearch(db = 'pubmed',
                            retmax = results['Count'],
                            term = keyword,
                            mindate = 2010,
                            maxdate = 2020
                            )
    result = Entrez.read(handle)
    ids = result['IdList']
    batch_size = 100
    batches = [ids[x: x + 100] for x in range(0, len(ids), batch_size)]
    return batches

def fetch(x):
    record_list = []
    for batch in tqdm(batches):
        handle = Entrez.efetch(db="pubmed", 
                               id=batch, 
                               rettype="medline", 
                               retmode="text")
        records = Medline.parse(handle)
        record_list.extend(list(records))
        print('Complete.')
    return record_list
    
if __name__ == '__main__':
    results = search(disease_list)
    batches= batch(results)
    record_list = fetch(batches)


'''Data analysis and visualization section'''

#import the saved CSV
from collections import Counter
import seaborn as sns
import matplotlib.pyplot as plt

sns.set_style("white")
publication_data = pd.DataFrame(record_list)



plt.figure(figsize=(12,12), dpi=100)


# Top authors


n = 20 #Set number of authors

#Count number of authors, and weigh by impact factor 
df = pd.DataFrame(record_list, columns = ['FAU', 'TA', 'AD'])
#This splits up the list into searate rows
df['FAU'] = df['FAU'].apply(pd.Series) \
    .merge(df, left_index = True, right_index = True)
df['TA'] = df['TA'].str.lower()
#df.to_excel('C:/Users/jwang/Desktop/Python/Biopython/CSV/test.xlsx')

mask = {'nature':10,
        'lancet':10,
        'jama': 11, 
        'n engl j med': 20,
        'hepatology': 8,
        r'(?s).*':1
        }

#Basically replicates the journal names column as a new column 
df['weight'] = df['TA'].replace(mask, regex = True)
df['weight'] = df['weight'].fillna(1)
df = df.reindex(df.index.repeat(df.weight))

#df.to_excel('C:/Users/jwang/Desktop/Python/Biopython/CSV/test1.xlsx')

authors_flat = [i for i in list(df['FAU'].values.tolist())]

top_authors = pd.DataFrame.from_records(
    Counter(authors_flat).most_common(n), columns=["Name", "Count"]
)

sns.barplot(x = 'Count', y = 'Name', data = top_authors, palette = 'RdBu_r')


# Publications over Time

publication_data.dropna(subset=['EDAT'], inplace=True)
publication_data["Year"] = (
    publication_data["EDAT"].astype(str).str[0:4].astype(int)
)
yearly = pd.DataFrame(publication_data["Year"].value_counts().reset_index())
yearly.columns = ["Year", "Count"]
sns.lineplot(x="Year", y="Count", data=yearly)
plt.title("Publications over Time")




# TOP 10 Journals


top10journals = pd.DataFrame.from_records(
    Counter(publication_data["TA"]).most_common(10),
    columns=["Journal", "Count"],
)

sns.barplot(x="Count", 
            y="Journal", 
            data=top10journals, 
            palette="RdBu_r")
plt.title("Top 10 Journals")

# Top associated keywords
plt.subplot(2, 2, 4)

flat_kw = [
    _.lower()
    for kws in list(publication_data["OT"].dropna())
    for kw in kws
    for _ in kw.split(" ")
]

top10kw = pd.DataFrame.from_records(
    Counter(flat_kw).most_common(10), columns=["Keyword", "Count"]
)

sns.barplot(x="Count", y="Keyword", data=top10kw, palette="RdBu_r")
plt.title("Top 10 Associated Keywords")
plt.subplots_adjust(top=1, bottom=0, left=0, right=1, hspace=0.3, wspace=0.3)
plt.show()

'''
In the next section we will visualize how often authors are publishing together. 
Each author is a node an the distance betweent wo authors is called an edge. 
'''
from itertools import combinations
import networkx as nx
from nxviz.plots import CircosPlot
'''
# Extract author connections
authors = publication_data["FAU"].dropna()
author_connections = list(
    map(lambda x: list(combinations(x[::-1], 2)), authors)
)
flat_connections = [item for sublist in author_connections for item in sublist]

#create a dataframe with the connections
df = pd.DataFrame(flat_connections, columns = ['From', 'To'])
df_graph = df.groupby(['From', 'To']).size().reset_index()
df_graph.columns = ['From', 'To', 'Count']

#Graph using networkx
G = nx.from_pandas_edgelist(df_graph, source = 'From', target = 'To', edge_attr = 'Count')
#Limit to top 50 authors
top_authors = pd.DataFrame.from_records(Counter(authors_flat).most_common(30), columns = ['Name', 'Count'])
top_nodes = (n for n in list(G.nodes()) if n in list(top_authors['Name']))
G_20 = G.subgraph(top_nodes)

#Give networkx something to rank by color
for n in G_20.nodes():
    G_20.node[n]['publications'] = int(top_authors[top_authors['Name'] == n]['Count'])
#Make the actual graph
c = CircosPlot(G_20, dpi = 600, node_grouping = 'publications', 
               edge_width = 'Count', figsize = (14,14), node_color = 'publications', node_labels = True)
c.draw()
plt.show()
'''
#Get contact info for top 20
'''
Get top author emails
'''
authors = list(top_authors['Name']) #Only 30 top authors in list
emails = pd.DataFrame(record_list, columns = ['FAU', 'AD']).dropna() #Pull out authors and, contact info
#This will sepearate out the names in the list into individual rows
emails['FAU'] = emails.FAU.apply(pd.Series) \
    .merge(emails, left_index = True, right_index = True)
#Checks to see if the name of the top 30 are in the author column
mask = emails.FAU.apply(lambda x: any(item for item in authors if item in x)) 
#Returns only rows that contain at least one of the authors 
emails = emails[mask] 
emails = emails.set_index('FAU')
emails = emails.loc[authors] #Reorder the authors according to their rankings 
#Now reset the index so that FAU can be called again for splitting the first/last name
emails = emails.reset_index()
emails['AD'] = emails['AD'].apply(lambda x: re.findall(r'[\w\.-]+@[\w\.-]+', str(x))) #removes everything except the email(s)
emails['AD'] = emails['AD'].apply(', '.join) #Turns the list of emails into strings
#pd.set_option('display.max_rows', emails.shape[0]+1)
#Rename the columns
emails = emails.rename(columns = {'FAU':'Name', 'AD':'email'})
#Split first and last name in separate columns
splitted = emails['Name'].str.split()
emails['First'] = splitted.str[1]
emails['Last'] = splitted.str[0]
#emails = emails.drop(columns = ['Name'])
#Swap email/last name location
column_titles = ['First', 'Last', 'email', 'Name']
emails = emails.reindex(columns = column_titles)
#Check if last name in email address
emails['Last'] = emails['Last'].str.lower()

#Drop duplicates
emails = emails.drop_duplicates(subset = 'Name',
                                keep = 'first',
                                inplace = False)

print(emails[['First', 'Last', 'email']])




