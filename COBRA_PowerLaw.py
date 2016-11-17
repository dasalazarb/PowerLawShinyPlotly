
# coding: utf-8

# In[1]:

#import all packages
import numpy as np
import re
import networkx as nx
import os
import cobra.test
import matplotlib.pyplot as plt
import graphviz
from igraph import *
#import colorlover as cl


# In[2]:

def is_reversible_recon(expression):
    return "<=>" in str(expression)

def is_irreversible_recon(expression):
    return "-->" in str(expression) or "<--" in str(expression)

def remove_coefficient(expression):
    return re.sub(r"^((\d+\.)?\d+\s)?","", expression)

def split_reaction(expression):
    return expression.replace(" -->", "/").replace("--> ","/").replace(" <=> ", "/").replace(" <-- ", "/").replace("<-- ", "/").split("/")
    
def split_metabolites(expression):
    return expression.replace(" + ", "/").split("/")

def split_metabolites_from_list(expression):
    return expression[0].replace(" + ", "/").split("/"), expression[1].replace(" + ", "/").split("/")


# Red metabolito-metabolito, donde se conectan por metabolito producto con metabolito reactivo. Para las irreversibles se tiene que metabolito producto tambien es reactivo y viceversa.

# In[3]:

for carpeta in os.listdir('C:\Users\Diego Salazar\Google Drive\Project_Power-Law\_modelos_'):
    #nuevo_nombre = 'C:\Users\Diego Salazar\Google Drive\Project_Power-Law\_modelos_\Modelos_faltantes'
    #usar para leer toda las carpetas:
    nuevo_nombre = 'C:\Users\Diego Salazar\Google Drive\Project_Power-Law\_modelos_' + "\\" + carpeta
    if os.path.isdir(nuevo_nombre) == True:
        os.chdir(nuevo_nombre)
        for modelo in os.listdir(nuevo_nombre):
            if modelo.endswith(".xml"):
                nombre_modelo = modelo.replace(".xml", "") + "_metab" + ".txt"
                model = cobra.io.read_sbml_model(modelo)
                irrev_rev = []
                dict_irrev_rev = {}
                for x in xrange(len(model.reactions)):
                    if is_reversible_recon(model.reactions[x].reaction):
                        irrev_rev.append(1)
                    elif is_irreversible_recon(model.reactions[x].reaction):
                        irrev_rev.append(0)

                for x in xrange(len(model.reactions)):
                    dict_irrev_rev[model.reactions[x]] = irrev_rev[x]

                #crear un diccionario: como key -> hmr_index como value -> rxn
                l_hmr = {}
                for x in xrange(len(model.reactions)):
                    l_hmr[model.reactions[x]] = model.reactions[x].reaction.strip()

                #diccionario de E.C. - rxn_izq_der, se dejo [['utp[e]'], ['']] para las reacciones de tipo utp[e] <=>
                l_ec_rxn_final = []
                l_inter = []
                dict_ec_rxn_izq_der = {}

                for key in l_hmr:
                    #poner la depuracion de expresiones regulares
                    para_recon = " " + l_hmr[key] + " "
                    l_ec_rxn = split_reaction(para_recon)
                    try:
                        l_inter = l_ec_rxn[0].replace(" + ", "/").split("/"), l_ec_rxn[1].replace(" + ", "/").split("/")
                    except:
                        l_inter = l_ec_rxn[0].replace(" + ", "/").split("/")
                    for x in xrange(len(l_inter)):
                        for y in xrange(len(l_inter[x])):
                            try:
                                l_inter[x][y] = l_inter[x][y].strip()
                                l_inter[x][y] = remove_coefficient(l_inter[x][y])
                            except:
                                pass
                    dict_ec_rxn_izq_der[key] = list(l_inter)

                #partir las reacciones y meterlas en nuevo diccionario.
                rxn_rev = list()
                rxn_irrev = list()
                rxn_irrev_final = dict()
                for x in xrange(len(model.reactions)):
                    if "<=>" in model.reactions[x].reaction:
                        rxn_rev = model.reactions[x].reaction.replace("<=>","/").split("/")
                    elif "--" in model.reactions[x].reaction:
                        rxn_irrev = model.reactions[x].reaction.replace("--","/").split("/")

                #lista de metabolitos. Se cargan, se quitan string no necesarios y se separan por /
                reaction = [model.reactions[x].reaction for x in xrange(len(model.reactions))]
                reaction = [x.replace(" + ", "/").replace("<=>", "/").replace("-->", "/").replace("<--", "/") for x in reaction]
                reaction = [x.split("/") for x in reaction]
                l = list(reaction[0])
                #a = []

                #Conteo para lista metabolitos. 
                for x in xrange(0, len(reaction)):
                    for j in xrange(0,len(reaction[x])):
                        reaction[x][j] = reaction[x][j].strip()
                        reaction[x][j] = remove_coefficient(reaction[x][j])
                    l+=reaction[x]

                #diccionario de hmr_ con metab en forma de lista.
                d_hmr = {}
                for x in xrange(len(model.reactions)):
                    d_hmr[model.reactions[x]] = reaction[x]
                d_recon_hmr = {}
                for key in d_hmr:
                    if '' in d_hmr[key]:
                        d_hmr[key].sort(reverse=True)
                        d_hmr[key].pop()
                        d_recon_hmr[key] = d_hmr[key]
                    else:
                        d_recon_hmr[key] = d_hmr[key]

                #diccionario de conteo de metabolitos
                counts = dict()
                for metab in l:
                    counts[metab] = counts.get(metab, 0) + 1
                #print counts

                #quitar duplicados (set()) y convertir en lista nuevamente.
                l.sort()
                lista_final_metab = list(set(l))
                #print len(list_final_metab)
                G=nx.DiGraph()
                #red metab_metab
                valor_interact = 1000000
                for i in lista_final_metab:
                    if counts[i] <= valor_interact:
                        for key in dict_ec_rxn_izq_der:
                            if i in dict_ec_rxn_izq_der[key][0]:#re.match('\\b' + pi + '\\b', x)
                                for u in xrange(len(dict_ec_rxn_izq_der[key][1])):
                                    if i == dict_ec_rxn_izq_der[key][1][u]:
                                        continue
                                    else:
                                        try:
                                            if counts[dict_ec_rxn_izq_der[key][1][u]] <= valor_interact:
                                                G.add_edge(i, dict_ec_rxn_izq_der[key][1][u], label=key)                                                
                                            else:
                                                pass
                                        except:
                                            pass
                            else:
                                continue
                            if i in dict_ec_rxn_izq_der[key][1]:
                                for u in xrange(len(dict_ec_rxn_izq_der[key][0])):
                                    if i == dict_ec_rxn_izq_der[key][0][u]:
                                        continue
                                    else:
                                        try:
                                            if counts[dict_ec_rxn_izq_der[key][0][u]] <= valor_interact:
                                                G.add_edge(dict_ec_rxn_izq_der[key][0][u], i, label=key)
                                            else:
                                                pass
                                        except:
                                            pass
                            else:
                                continue
                grados = G.degree()
                grados.pop('', None)
                grados_totales = nombre_modelo.replace('.txt', '') + '_TotalDegree.txt'
                with open (grados_totales, 'w') as fp:
                    for key, valor in sorted(grados.iteritems(), key=lambda (k,v): (v,k), reverse=True):
                        fp.write("%s\t%s\n" % (valor, key))
                in_grados = G.in_degree()
                in_grados.pop('', None)
                grados_in = nombre_modelo.replace('.txt', '') + '_InDegree.txt'
                with open (grados_in, 'w') as fp:
                    for key, valor in sorted(in_grados.iteritems(), key=lambda (k,v): (v,k), reverse=True):
                        fp.write("%s\t%s\n" % (valor, key))
                out_grados = G.out_degree()
                out_grados.pop('', None)
                grados_out = nombre_modelo.replace('.txt', '') + '_OutDegree.txt'
                with open (grados_out, 'w') as fp:
                    for key, valor in sorted(out_grados.iteritems(), key=lambda (k,v): (v,k), reverse=True):
                        fp.write("%s\t%s\n" % (valor, key))


# Red reacion-reaccion

# In[4]:

for carpeta in os.listdir('C:\Users\Diego Salazar\Google Drive\Project_Power-Law\_modelos_'):
    #nuevo_nombre = 'C:\Users\Diego Salazar\Google Drive\Project_Power-Law\_modelos_\Modelos_faltantes'
    #usar para leer toda las carpetas:
    nuevo_nombre = 'C:\Users\Diego Salazar\Google Drive\Project_Power-Law\_modelos_' + "\\" + carpeta
    if os.path.isdir(nuevo_nombre) == True:
        os.chdir(nuevo_nombre)
        for modelo in os.listdir(nuevo_nombre):
            if modelo.endswith(".xml"):
                nombre_modelo = modelo.replace(".xml", "") + "_rxn" + ".txt"
                model = cobra.io.read_sbml_model(modelo)
                irrev_rev = []
                dict_irrev_rev = {}
                for x in xrange(len(model.reactions)):
                    if is_reversible_recon(model.reactions[x].reaction):
                        irrev_rev.append(1)
                    elif is_irreversible_recon(model.reactions[x].reaction):
                        irrev_rev.append(0)

                for x in xrange(len(model.reactions)):
                    dict_irrev_rev[model.reactions[x]] = irrev_rev[x]

                #crear un diccionario: como key -> hmr_index como value -> rxn
                l_hmr = {}
                for x in xrange(len(model.reactions)):
                    l_hmr[model.reactions[x]] = model.reactions[x].reaction.strip()

                #diccionario de E.C. - rxn_izq_der, se dejo [['utp[e]'], ['']] para las reacciones de tipo utp[e] <=>
                l_ec_rxn_final = []
                l_inter = []
                dict_ec_rxn_izq_der = {}

                for key in l_hmr:
                    #poner la depuracion de expresiones regulares
                    para_recon = " " + l_hmr[key] + " "
                    l_ec_rxn = split_reaction(para_recon)
                    try:
                        l_inter = l_ec_rxn[0].replace(" + ", "/").split("/"), l_ec_rxn[1].replace(" + ", "/").split("/")
                    except:
                        l_inter = l_ec_rxn[0].replace(" + ", "/").split("/")
                    for x in xrange(len(l_inter)):
                        for y in xrange(len(l_inter[x])):
                            try:
                                l_inter[x][y] = l_inter[x][y].strip()
                                l_inter[x][y] = remove_coefficient(l_inter[x][y])
                            except:
                                pass
                    dict_ec_rxn_izq_der[key] = list(l_inter)

                #partir las reacciones y meterlas en nuevo diccionario.
                rxn_rev = list()
                rxn_irrev = list()
                rxn_irrev_final = dict()
                for x in xrange(len(model.reactions)):
                    if "<=>" in model.reactions[x].reaction:
                        rxn_rev = model.reactions[x].reaction.replace("<=>","/").split("/")
                    elif "--" in model.reactions[x].reaction:
                        rxn_irrev = model.reactions[x].reaction.replace("--","/").split("/")

                #lista de metabolitos. Se cargan, se quitan string no necesarios y se separan por /
                reaction = [model.reactions[x].reaction for x in xrange(len(model.reactions))]
                reaction = [x.replace(" + ", "/").replace("<=>", "/").replace("-->", "/").replace("<--", "/") for x in reaction]
                reaction = [x.split("/") for x in reaction]
                l = list(reaction[0])
                #a = []

                #Conteo para lista metabolitos. 
                for x in xrange(0, len(reaction)):
                    for j in xrange(0,len(reaction[x])):
                        reaction[x][j] = reaction[x][j].strip()
                        reaction[x][j] = remove_coefficient(reaction[x][j])
                    l+=reaction[x]

                #diccionario de hmr_ con metab en forma de lista.
                d_hmr = {}
                for x in xrange(len(model.reactions)):
                    d_hmr[model.reactions[x]] = reaction[x]
                d_recon_hmr = {}
                for key in d_hmr:
                    if '' in d_hmr[key]:
                        d_hmr[key].sort(reverse=True)
                        d_hmr[key].pop()
                        d_recon_hmr[key] = d_hmr[key]
                    else:
                        d_recon_hmr[key] = d_hmr[key]

                #diccionario de conteo de metabolitos
                counts = dict()
                for metab in l:
                    counts[metab] = counts.get(metab, 0) + 1
                #print counts

                #quitar duplicados (set()) y convertir en lista nuevamente.
                l.sort()
                lista_final_metab = list(set(l))
                #print len(list_final_metab)
                H=nx.DiGraph()
                valor_interact = 100000
                for key in dict_ec_rxn_izq_der:
                    for i in dict_ec_rxn_izq_der[key][0]:
                        for key_2 in dict_ec_rxn_izq_der:
                            try:
                                if counts[i] <= valor_interact:
                                    if i in dict_ec_rxn_izq_der[key_2][1]:
                                        H.add_edge(str(key_2), str(key), label=i)
                                    else:
                                        pass
                                else:
                                    pass
                            except:
                                pass
                    for a in dict_ec_rxn_izq_der[key][1]:
                        for key_2 in dict_ec_rxn_izq_der:
                            try:
                                if counts[i] <= valor_interact:
                                    if j in dict_ec_rxn_izq_der[key_2][0]:
                                        H.add_edge(str(key), str(key_2), label=i)
                                    else:
                                        pass
                                else:
                                    pass
                            except:
                                pass
                #Diccionario de (enzima-enzima): label (metabolito de union)
                label=nx.get_edge_attributes(H,'label')
                grados = H.degree()
                grados.pop('', None)
                grados_totales = nombre_modelo.replace('.txt', '') + '_TotalDegree.txt'
                with open (grados_totales, 'w') as fp:
                    for key, valor in sorted(grados.iteritems(), key=lambda (k,v): (v,k), reverse=True):
                        fp.write("%s\t%s\n" % (valor, key))
                in_grados = H.in_degree()
                in_grados.pop('', None)
                grados_in = nombre_modelo.replace('.txt', '') + '_InDegree.txt'
                with open (grados_in, 'w') as fp:
                    for key, valor in sorted(in_grados.iteritems(), key=lambda (k,v): (v,k), reverse=True):
                        fp.write("%s\t%s\n" % (valor, key))
                out_grados = H.out_degree()
                out_grados.pop('', None)
                grados_out = nombre_modelo.replace('.txt', '') + '_OutDegree.txt'
                with open (grados_out, 'w') as fp:
                    for key, valor in sorted(out_grados.iteritems(), key=lambda (k,v): (v,k), reverse=True):
                        fp.write("%s\t%s\n" % (valor, key))


# In[6]:

G.out_degree()


# In[7]:

G.in_degree()


# In[8]:

G.degree()

