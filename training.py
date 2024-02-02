#Bioinformatics of RNA, project G.Postic
# Thémis LEMARCHAND 20222535 M2GENIOMHE
# Training script

import re
import numpy as np
from math import log

#creation of the dictionnary that will contain all computed data
pairs_list = ["AA","AU","AC","AG","UU","UC","UG","CC","CG","GG"]
pairs={}
pairs["all_distances"]=[]
pairs["frequences_ref"]=[]
for p in pairs_list:
     pairs[p]={}
     pairs[p]['distances']=[]
     pairs[p]['frequences']=[]
     pairs[p]['score']=[]

# Creation of a dictionnary containing data from pdb file
f = open("2jyh.pdb","r")
lines=f.readlines()
f.close()

data=[]
for line in lines:
    if re.search("ATOM[\d\s\D]*C3'",line):
        data.append(line)
     
# creation of a dictionnary containing the id of the atom
# each atom is described by a dictionnairy containing the nucleotide, the chain id ,the residue position, x, y and Z positions
# atoms={atom_id={"nucleotide","chain_id","sequence_number",x","y","z"}}
atoms={}
id=0
for l in data:
    field=l.split()
    if(field[4]=='A'):
        atoms[id]={"nucleotide":field[3],"chain_id":field[4],"position":int(field[5]),"x":float(field[6]),"y":float(field[7]),"z":float(field[8])}
        id+=1

# computation of distances between residues
for key in atoms:
    for key2 in atoms :
           if( key2 > key):
                if(atoms[key2]["position"]-atoms[key]["position"]>3):
                    
                    dist = int(np.sqrt((atoms[key2]['x'] - atoms[key]['x'])**2 + (atoms[key2]['y'] - atoms[key]['y'])**2 + (atoms[key2]['z'] - atoms[key]['z'])**2))
                    if dist <= 20 :
                        pair=atoms[key]["nucleotide"]+atoms[key2]["nucleotide"]
                        if pair in pairs:
                            pairs[pair]['distances'].append(dist)
                            pairs["all_distances"].append(dist)
                        else:
                            pair=atoms[key2]["nucleotide"]+atoms[key]["nucleotide"]
                            pairs[pair]['distances'].append(dist) 
                            pairs["all_distances"].append(dist)
 
# computation of frequencies                       
for pair in pairs:
        for i in range(0,20):
            if (pair!="all_distances" and pair!="frequences_ref"):
                Nij=0
                for j in pairs[pair]['distances']:                    
                    if j == i:
                        Nij+=1         
                freq = Nij/len(pairs[pair]['distances'])
                pairs[pair]['frequences'].append(freq)
            elif pair == "all_distances":
                  Nxx=0
                  for j in pairs['all_distances']:
                        if j == i:
                            Nxx+=1         
                  freq_ref = Nxx/len(pairs['all_distances'])
                  pairs['frequences_ref'].append(freq_ref)

# computation of scores
for pair in pairs:
    for i in range (0,20):
        if (pair!="all_distances" and pair!="frequences_ref"):
            if (pairs['frequences_ref'][i] ==0 or pairs[pair]['frequences'][i] ==0):
                u=0
            else:
                u=-log(pairs[pair]['frequences'][i]/pairs['frequences_ref'][i])
            if u <= 10 :
                pairs[pair]["score"].append(u)
            else:
                u=10
                pairs[pair]["score"].append(u)


# saving scoring data into files
for pair in pairs:
     if (pair!="all_distances" and pair!="frequences_ref"):
        name="data/"+pair
        f=open(name,'w')
        for score in pairs[pair]['score']:
            line=str(score)+"\n"
            f.write(line)
        f.close()
            
                
