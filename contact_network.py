__author__ = 'Anindya Guha (aguha@colgate.edu)'

from numpy import *
import networkx as nx
from networkx.algorithms.approximation import * 
import matplotlib.pyplot as plt 
from scipy.stats import rv_discrete

def is_infected(edge_weight, percentage_infected, timestep):
    
    prob = 1 - (1 - edge_weight*percentage_infected/100.0)**timestep

    xk = (0, 1, 2, 3)
    
    pk = ((1-prob), prob/3.0, prob/3.0, prob/3.0)
    
    val = rv_discrete(name='val', values=(xk, pk))

    R = val.rvs(size=1)

    return R

G = nx.Graph() 

all_cities = ['Conakry', 'Gueckdou', 'Macenta', 'Kissidougou', 'Boffa', 'Nzerekore', 'Yomou', 'Lola', 'Dubreka', 'Forecariah', 'Kerouane', 'Coyah']

G.add_nodes_from(all_cities)


G.add_weighted_edges_from([('Nzerekore', 'Lola', 0.51), ('Nzerekore', 'Yomou', 0.15), ('Nzerekore', 'Macenta', 0.709), ('Macenta', 'Gueckdou', 0.1111), ('Macenta', 'Kissidougou', 0.5478), ('Gueckdou', 'Kissidougou', 0.7261), ('Conakry', 'Boffa', 0.6), ('Conakry', 'Kissidougou', 0.9210), ('Gueckdou', 'Yomou', 0.0134)])
G.add_weighted_edges_from([('Dubreka', 'Conakry', 0.6), ('Dubreka', 'Boffa', 0.45), ('Forecariah', 'Conakry', 0.51), ('Kerouane', 'Kissidougou', 0.54), ('Kerouane', 'Macenta', 0.45), ('Coyah','Dubreka', 0.256), ('Coyah', 'Conakry', 0.782), ('Coyah', 'Forecariah', 0.1567)])
G.add_weighted_edges_from([('Coyah', 'Kissidougou', 0.345), ('Forecariah', 'Kissidougou', 0.275), ('Yomou', 'Lola', 0.1)])


very_infected_cities = ['Macenta']
moderately_infected_cities = ['Conakry', 'Gueckdou']
slightly_infected_cities = ['Nzerekore', 'Yomou', 'Dubreka', 'Kerouane']

infected_cities = very_infected_cities+moderately_infected_cities+slightly_infected_cities

for city in (infected_cities):
    G.node[city]['infected'] = True 
    if city in very_infected_cities:
        G.node[city]['percentage_infected'] = 10
    elif city in moderately_infected_cities:
        G.node[city]['percentage_infected'] = 5
    elif city in slightly_infected_cities:
        G.node[city]['percentage_infected'] = 2
    else:
        raise 

uninfected_cities = []

for city in all_cities:
    if city not in infected_cities:
        G.node[city]['infected'] = False 
        uninfected_cities.append(city)

elarge=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight'] >0.5]
esmall=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight'] <=0.5]

pos = {'Conakry': (2,9), 'Boffa':(1,11), 'Gueckdou':(7,4.5), 'Macenta':(10, 6), 'Kissidougou':(7, 6), 'Nzerekore':(11.5, 5), 'Lola':(13.5, 4.5), 'Yomou':(8, 4),
        'Dubreka': (3,10), 'Forecariah': (4, 7), 'Kerouane': (9, 6.5), 'Coyah': (4.5, 8)}


for city in infected_cities: 
    for neighbor in G.neighbors(city):
        if G.node[neighbor]['infected'] == False: 
            infection_amount = is_infected(G.edge[city][neighbor]['weight'], G.node[city]['percentage_infected'], 100)
            if infection_amount > 0:
                G.node[neighbor]['infected'] = True
                if infection_amount == 1:
                    slightly_infected_cities.append(neighbor)
                elif infection_amount == 2:
                    moderately_infected_cities.append(neighbor)
                elif infection_amount == 3:
                    very_infected_cities.append(neighbor)
                else:
                    raise
                uninfected_cities.remove(neighbor)

for (u,v,d) in G.edges(data=True):
    G.edge[u][v]['weight'] = 1.0/G.edge[u][v]['weight']

dominating_set = ['Conakry', 'Kissidougou', 'Nzerekore'] 

remaining_cities = [] 

for city in all_cities:
    if city not in dominating_set:
        remaining_cities.append(city)


nx.draw_networkx_nodes(G, pos, nodelist=very_infected_cities, node_color='#7f0000')

nx.draw_networkx_nodes(G, pos, nodelist=moderately_infected_cities, node_color='r')
nx.draw_networkx_nodes(G, pos, nodelist=slightly_infected_cities, node_color='m')

nx.draw_networkx_nodes(G, pos, nodelist=uninfected_cities, node_color='b')

nx.draw_networkx_edges(G,pos,edgelist=elarge,
                    width=2.5, edge_color='#323234')
nx.draw_networkx_edges(G,pos,edgelist=esmall,
                    width=2,alpha=0.5,edge_color='b')
nx.draw_networkx_labels(G,pos,font_size=10,font_family='sans-serif')

plt.axis('off')
plt.savefig('hundred_days_later_correct.png')

plt.show() 

#print is_infected(0.2, 1, 4)
