# -*- coding: utf-8 -*-
"""
Genetic Algorithm
for businessman trip trouble
using cities locations only [(x,y)]

@author: RAGAMAR
"""

import csv, importlib, subprocess, sys
from random import sample, random, randint, uniform, gauss, choices
from math import sqrt
from copy import copy
from os.path import exists
import matplotlib.pyplot as plt

def ensureModule(nm:str):
    try:
        importlib.import_module(nm)
    except ImportError:
        print(f"Module {nm} not found. Installing...")
        subprocess.check_call([sys.executable,"-m","pip","install",nm])
        print(f"{nm} installed - OK")


# ========================== DICTIONARIES =============================== #
mutations:dict={
    'random':0, # randomly selects a valid mutation method
    'insert':1, #insert element
    'swap':2, # element swap
    'scramble':3, #scramble methods all-scramble / section-scramble
    'inversion':4 #inverts all values
    }

crossOvers:dict={'random':0, #randomly selects a valid crossover method
                 'tails':1, #mixes parents by swapping tails
                 'section':2, #mixes parents swapping mid sections
                 'uniform':3, #selects 1/n parent's bit at every bit location
                 }

selections:dict = {
    'random':0, #randomly selects a valid selection method
    'roulette':1, #roulette wheel
    'tournament':2, #torunament
    'rank':3,
    'sus':4 #stochastic universal sampling
    }

survivals:dict = {
    'generational':0, #only picks from offspring
    'comma':1, # "μ , λ" method
    'plus':2, # "μ + λ" method
    'ssr':3, #steady-state replacement only few offsrping replace worst parents
    'elitism':4 #best parents remain, offspring fills
    }

conversions:dict = {
    'none':0, #returns as string
    'int':1, #converts data to integer
    'float':2, #converts read data to float type
    }
#--------------------------------------------------------------------------

# ========================== GENERAL CONTROL ================================ #

csvdir:str = 'cities.csv'           #csv location

ppl:list = []                       # actual population
head:list = []                      # headers
useReturn:bool = True               # include returning trace
samples=10                          # population (different routes)
crate = 0.2                         # crossover rate
mrate = 0.2                         # mutation rate
iters = 100                         # iterations
useMax = False                      # max = best fit?
kParents = 2                        # No. of parents to reamin when +

#------------------------------------------------------------------------------

# ========================== FUNCTIONS ====================================== #
                
#read csv and returns population (heads stored in __head)
def read(datafile:str,_convert:int=0,nline='',
         x_indx:int=0,y_indx:int=1)->list:
    '''reads csv and returns coordinates X and Y from it (dots)\n
    "datafile" = csv file location (str)\n
    "convert" = see "conversions" dictionary (default = none)\n
    "nline" = allows newline value modification\n
    "x_indx" = row's "x" location (default = 0)\n
    "y_indx" = row's "y" location (default = 1)\n
    Note: access read data from property "original"
    '''
    ppl = []
    with open(datafile,newline=nline) as f:
        reader = csv.reader(f)
        head = next(reader) #header
        for row in reader:
            if _convert==conversions['none']: 
                x=row[x_indx]; y=row[y_indx]
            elif _convert==conversions['int']: 
                x=int(float(row[x_indx])); y=int(float(row[y_indx]))
            elif _convert==conversions['float']:
                x=float(row[x_indx]); y=float(row[y_indx])
            else:
                raise ValueError(f"conversion type {_convert} isn't supported")
            ppl.append([x,y]) #x,y
    plt.scatter([ct[0]for ct in ppl],
                [ct[1]for ct in ppl],
                marker='o',color='black')
    plt.title('Original Population Read')
    plt.grid(True)
    plt.show()
    return ppl

#gets a mid-point to trace a squared route
def setMidPoint(p0:list, p1:list)->list:
    '''adds a mid point considering least distance first\n
    "p0" = initial ponint (x,y)\n
    "p1" = final point (x,y)\n'''
    d=[p1[0]-p0[0],p1[1]-p0[1]] #x&y distances
    dab=[abs(d[0]),abs(d[1])]
    idx = dab.index(min(dab))
    return[p1[0]if idx>0 else p0[0],p0[1]if idx>0 else p1[1]]

#measures linear distance between 2 points    
def getDistance(p0:list,p1:list)->float:
    '''gets the distance between 2 points\n
    "p0" = initial point (x,y)\n
    "p1" = final point (x,y)'''
    return sqrt((p1[0]-p0[0])**2+(p1[1]-p0[1])**2)

#gets the total route distance (adds initial point)
def routeDistance(route:list)->float:
    d = 0
    #add returning route
    r=route+route[0]
    for _ in range(len(route)-1):
        d+=sqrt(((r[_+1][0]-r[_][0])**2)+((r[_+1][1]-r[_][1])**2))
    return d

#trace (allows start-end or cyclyc routing)
def traceSqr(popl:list,cycle:bool=True)->list:
    '''creates possible chromosomes of a given pupulation\n
    "ppl" = population [(x,y)]\n
    "cycle" = considers a returning distance also 
    -USE WITH MIDPOINTS- (bool)\n
    "plot_num" = number of random chromosomes to plot (int)'''
    #mid points
    mids = [setMidPoint(popl[_],popl[_+1])for _ in range(len(popl)-1)]
    #insert mid points
    for _ in range(len(mids)): popl.insert((_*2)+1,mids[_])
    #create return trace
    if cycle: popl.append(setMidPoint(popl[-1],popl[0]))
    return popl

#shows a route in plot    
def show(popl:list,useMids:bool=False,ttl=None)->None:
    '''plots the indicated chromosome\n
    "chrom" = chromosome to plot ([x,y])\n
    "useMids" = chromosome has mid-poins (bool)
    '''
    fin = len(popl)
    stp = 2 if useMids else 1
    #set title
    if ttl==None: ttl = 'Random Route'
    #Original
    plt.scatter([popl[_][0]for _ in range(0,fin,stp)],
                [popl[_][1]for _ in range(0,fin,stp)],
                marker='o',color='black',
                label='Cities')
    if useMids: #Mid-Points
        plt.scatter([popl[p][0]for p in range(1,fin,2)],
                    [popl[p][1]for p in range(1,fin,2)],
                    marker='o',color='gray',
                    label='Mid-Points')
    col=(random(),random(),random())
    plt.plot([_[0]for _ in popl],[_[1]for _ in popl],
             color=col,label='squared route'if useMids else 'direct route')
    plt.legend()
    plt.title(ttl)
    plt.grid(True)
    plt.show()
    
#selects a determined number of parents
def selection(parents,mode:int=0,quantity:int=0):
    '''performs parent selection\n
    - "parents" = parents list\n
    - "mode" = parent selection method - see Selections dictionary-\n
    - " quantity" = desired selection amount
    '''
    global fit
    #adjust mode
    if mode==0: mode=randint(1,max(selections.values()))
    #adjust quantity
    ls = len(parents)
    #no selection cases
    if quantity<=0 or quantity>=ls:return parents
    #fitness
    fit = [routeDistance(parents[_])for _ in range(len(parents))]
    #total fitness
    tf = sum(fit)
    #selection
    sel=[]

    #probability is prop-to fitnes    
    if mode==selections['roulette']: #0
        #distribute probability
        prob = [_/tf for _ in fit]
        #accumulate probability
        cd = [sum(_ for _ in prob[:i])for i in range(1,len(prob)+1)]
        #spin the wheel
        while len(sel)<quantity:
            for _ in range(len(cd)):
                if cd[_]>=random(): 
                    sel.append(parents[_])
                    break
                
    #fragments parents in small groups and compete
    elif mode==selections['tournament']: #1
        fit = [[fit[_],parents[_]]for _ in range(ls)]
        while len(sel)<quantity:
            gsize = randint(1,ls-1) #group size
            smp = ls//gsize+((ls//gsize)>0) #samples
            #candidates
            cands = [sample(fit,gsize)for _ in range(smp)]
            #pre-selection
            sel+=[max(cands[_])if useMax else min(cands[_])for _ in range(smp)]
            sel=[sel[_][1]for _ in range(len(sel))]
    
    #probability is exponential
    elif mode==selections['rank']: #2
        fit = [[fit[_],parents[_]]for _ in range(ls)]
        fit.sort(reverse=useMax)
        q=0.7
        pi=[(((1-q)*(q**(_-1)))/(1-q**ls))for _ in range(1,ls+1)]
        #cumulative distribution
        cd=[sum([_ for _ in pi[:i]])for i in range(1,ls+1)]
        while len(sel)<quantity:
            for _ in range(len(cd)):
                if cd[_]>=random(): 
                    sel.append(fit[_][1])
                    break
        
    #stochastic universal sampling
    elif mode==selections['sus']: #3
        prob = [_/tf for _ in fit]
        #cumulative
        cd=[sum([_ for _ in prob[:i]])for i in range(1,ls+1)]
        while len(sel)<quantity:
            r=random()/quantity
            pts=[r+_/quantity for _ in range(quantity)]
            for p in pts:
                for i,c in enumerate(cd):
                    if p<=c: sel.append(parents[i]);break
    

    #fail
    else:
        raise ValueError(f"selection mode {mode} is not supported")
    return sel
        
#crossover parents function
def cross(parents:list,mode:int=0,quantity:int=0,location=-1):
    '''performs crossover\n
    - "parents" = parents list\n
    - "mode" = crossover mode - check crossOvers dictionary -\n
    - "quantity" = desired output number
    '''
    #get parents
    p0,p1=parents[0],parents[1]
    #length
    ls=len(p0)
    #adjust quantity
    if quantity<=0 or quantity>2:quantity=2
    #adujst mode
    if mode==0: mode=randint(1,max(crossOvers.values()))
    #offspring
    offs=[]
    
    #print(f'cross mode = {mode}')
    
    #swap tails
    if mode==crossOvers['tails']: #1
        #adjust location
        if location==-1 or location>=ls:location=randint(0,ls-1)
        for _ in range(quantity):
            a,b = (p0,p1) if _==0 else (p1,p0)
            h=a[:location]
            for g in b[location:]:
                if not g in h:h.append(g)
            for g in b[:location]:
                if not g in h:h.append(g)
            offs.append(h)
    
    #swap section
    elif mode==crossOvers['section']: #2
        #adjust location
        if isinstance(location,int):
            location=[randint(0,len(p0)-1)for _ in range(2)]
            while location[0]==location[1]:location[1]=randint(0,len(p0)-1)
        location.sort() #ensure min->max
        for _ in range(quantity):
            a,b = (p0,p1) if _==0 else (p1,p0)
            h = a[location[0]:location[1]]
            for g in b:
                if not g in h:h.append(g)
            offs.append(h)
    
    #uniform distribution
    elif mode==crossOvers['uniform']: #3
        #adjust location
        while len(offs)<quantity:
            h=[]; i=0
            while len(h)<ls:
                g,g1 = (p0[i],p1[i])if random()>0.5 else (p1[i],p0[i])
                if not g in h:
                    h.append(g)
                elif not g1 in h:
                    h.append(g1)
                i=i+1 if i<ls-1 else 0
            offs.append(h)
            
    #fail
    else:raise ValueError(f"crossover mode {mode} not supported")
    return offs
    
#selects individuals
def survive(parents,offs,mode:int=-1,quantity:int=0):
    '''performs survival selection\n
    - "parents" = parents list
    - "offs" = offspring list
    - "mode" = -check survivals dictionary-
    - "quantity" desired survivors number
    '''
    #adjust mode
    if mode==-1:mode = randint(0,max(selections.values()))
    #adjust quantity
    if quantity<=0:quantity=1 #1 child
    #new generation
    new = []
    
    #print(f"survive mode = {mode}")
    
    #take only from offspring
    if mode==survivals['generational']: #0
        ofit = [[routeDistance(route),route]for route in offs]
        ofit.sort(reverse=useMax)
        new= [ofit[_][1]for _ in range(quantity)]
    
    #add parents + offspring
    elif mode==survivals['plus']: #1
        ofit = [[routeDistance(route),route]for route in parents]
        ofit += [[routeDistance(route),route]for route in offs]
        ofit.sort(reverse=useMax)
        new=[ofit[_][1]for _ in range(quantity)]
    
    #save k parents, add remain offspring
    elif mode==survivals['comma']: #2
        ofit = [[routeDistance(route),route]for route in offs]
        ofit.sort(reverse=useMax)
        new=[ofit[_][1]for _ in range(quantity)]
        
    elif mode==survivals['ssr']: #3
        if useMax: #direct
            ofit = [routeDistance(_)for _ in offs]
        else: #inverse
            ofit = [1/routeDistance(_)for _ in offs]
        #total fitness
        tf = sum(ofit)
        #probability
        probs = [_/tf for _ in ofit]
        new = choices(offs,weights=probs,k=quantity)
        
        
    elif mode==survivals['elitism']: #4
        pfit = [[routeDistance(route),route]for route in parents]
        pfit.sort(reverse=useMax)
        new=[pfit[_][1]for _ in range(kParents)]
        ofit = [[routeDistance(route),route]for route in offs]
        ofit.sort(reverse=useMax)
        new+=[ofit[_][1]for _ in range(quantity-kParents)]
        new = new[:quantity]
    
    else:
        raise ValueError(f"survival mode {mode} not supported")
    return new

#mutates chromosomes
def mutate(route,mode:int=0):
    '''performs mutation with given method\n
    - "gene" = gene to mutate (city) = [x,y]\n
    - "mode" = mutation method - check mutations dictionary -
    '''
    #adjust mode
    if mode==0: mode = randint(1,max(mutations.values()))
    ls = len(route)
    
    # remove 1 element and insert it somewhere else
    if mode==mutations['insert']: #1
        #get element by popping it from list
        el = route.pop(randint(0,ls-1))
        #insert it somewhere else
        route.insert(randint(0,ls-1),el)
        
    # swap elements
    elif mode==mutations['swap']: #2
        loc = [randint(0,ls-1)for _ in range(2)]
        while loc[0]==loc[1]: loc[1]=randint(0,ls-1)
        loc.sort() #min->max
        #swap values
        route[loc[0]],route[loc[1]]=route[loc[1]],route[loc[0]]
    
    # scrambles whole elements
    elif mode==mutations['scramble']: #3
        route = sample(route,ls)
    
    # inverts the route
    elif mode==mutations['inversion']: #4
        route = route[::-1]
    
    #fail
    else:
        raise ValueError(f"mutation mode {mode} not supported")
    return route

#performs GA
def evolve(route,pplSize:int,iters:int,
           crossMode:int=0, mutMode:int=0):
    '''performs evolution with given paremeters\n
    - "route" = initial known route\n
    - "pplSize" = desired population size (No. of routes)\n
    - "iters" = No. of iterations to try-off\n
    - "crossMode" = crossover method -check crossOvers dictionary-\n
    - "mutMode" = mutation method -check mutations dictionary-
    '''
    #initial population
    smp = len(route)
    ppl = [sample(route,len(route))for _ in range(pplSize)]
    toff, offs = [],[] #auxiliaries
    #show 1 possible route
    nr = randint(0,pplSize-1)
    show(ppl[nr]+[ppl[nr][0]],False,
         f"Initial Population's Route No. {nr}/{pplSize}")
    
    for itr in range(iters+1):
        #reset offspring
        offs = []
        #parent selection
        sys.stdout.write("\r - {:05.02f} % completed".format(itr/iters*100))
        sys.stdout.flush()
        tppl=[] #temporary
        for _ in range(pplSize):
            #selection|
            toff=selection(ppl,quantity=2)
            #crossover (crossover-rate)
            toff=cross(toff,0)if random()>=crate else toff
            #mutation (mutation-rate)
            offs=[mutate(_,0)if random()>=mrate else _ for _ in toff]
            #survival
            tppl.append(survive(toff,offs,mode=-1,quantity=1)[0])
        plt.plot(min([routeDistance(_)for _ in tppl]))
        ppl=tppl
    
    return ppl
        
def bests(ppl:list,quantity:int):
    bst=[]
    #fitness
    vs=[]
    fit = [routeDistance(_)for _ in ppl]
    while len(vs)<quantity and len(ppl)>0:
        v=max(fit)if useMax else min(fit)
        idx = fit.index(v)
        if not v in vs:
            vs.append(v)
            cfit=fit.pop(idx)
            p=ppl.pop(idx)
            show(p+[p[0]],ttl='top {} result: fitness={:05.02f}m'.format(
                    len(vs),cfit))
            if len(vs)==quantity:break
        else:
            fit.pop(idx)
            ppl.pop(idx)
        
# -----------------------------------------------------------------------------

# ========================== PROGRAM ======================================== #

# read the route from csv file
route = read(csvdir,conversions['int'])
#perform evolution
ppl=evolve(route,pplSize=100,iters=1000)
#
bests(copy(ppl),5)
