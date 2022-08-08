#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 30 16:50:22 2022

@author: aidenlee
"""

from neuron import h, gui
from neuron.units import ms, mV
h.load_file('stdrun.hoc')
import numpy as np
import matplotlib.pyplot as plt
import random
import math
'''

Questions:
    Carla plan: Rely on delay for changing frequency to 50Hz, push the network 
    onto the brink of entering significant inhibition by toying around with
    weight in order to keep realistic tau, then add ChR2 to cause inhibition,
    measuring with LFP to see how ChR2 can cause dips in LFP.
    
    Dr. Bezaire plan: Change tau and weight, but I'm not sure how to do that due
    to the limitations of the refractory period. Would require lower weight, but
    greater tau to compensate, however too high tau keeps synaptic channels in 
    refractory period, and thus can't trigger another action potential even if
    reaching threshold.

'''

class Cell:         #Creation of a generic class called cell
    def __init__(self, gid, x, y, z, cell_type):  
        self._cell_type = cell_type
        self._gid = gid
        if self._cell_type == 'e':
            self._cell_type = "Excitatory"
        else:
            self._cell_type = 'Inhibitory'
        self._name = 'BallAndStick, {}'.format(self._cell_type)
        self._setup_morphology()
        self.all = self.soma.wholetree()    
            #will return a list of the whole tree of sections (somas and dendrites)
            #self.all is all the sections of the neuron
            #more future proof than creating a list of only soma and dendrite
        self._setup_biophysics(cell_type)       #Runs biophysics definition but with input of cell_type to determine biophysics
        self._setup_IClamp(cell_type)
        self.x = self.y = self.z = 0   
        h.define_shape()                     
        self._set_position(x, y, z)
        self.vrec_soma = None
        self.vrec_dend = None
        self._spike_detector = h.NetCon(self.soma(0.5)._ref_v, None, sec=self.soma)
        self.spike_times = h.Vector()
        self._spike_detector.record(self.spike_times)
        
    def __repr__(self):
        return '{} [{}]'.format(self._name, self._gid)
            #Formatted so it will return name[gid(aka the neuron number)]
    def _set_position(self, x, y, z):
        self.x, self.y, self.z = x, y, z
            
            
            
class BallAndStick(Cell):
    #Creation of subclass under cell so BallAndStick will have methods defined under both Cell and BallAndStick
    #Both init and repr removed because they are already defined in Cell
    def _setup_IClamp(self, cell_type):
        if cell_type == "e":
            self.stimClampE = h.IClamp(self.soma(0.5))
            self.stimClampE.delay = 30
            self.stimClampE.dur = 400
            self.stimClampE.amp = 0.036
        else:
            self.stimClampI = h.IClamp(self.soma(0.5))
            self.stimClampI.delay = 30
            self.stimClampI.dur = 400
            self.stimClampI.amp = 0.036
            #Amp which neuron fires on its own is 0.036nA and greater
            #Stim set to 30 ms to avoid all firing immediately
    def _setup_morphology(self):
        self.soma = h.Section(name='soma', cell=self)
        self.dend = h.Section(name='dend', cell=self)
        self.dend.connect(self.soma)
        self.soma.L = self.soma.diam = 12.6157
        self.dend.L = 200
        self.dend.diam = 1
    def _setup_biophysics(self, cell_type):
        for sec in self.all:
            sec.Ra = 100    # Axial resistance in Ohm * cm
            sec.cm = 1      # Membrane capacitance in micro Farads / cm^2
        self.soma.insert('hh')  
        self.dend.insert('hh')                                        
        for seg in self.soma:
            seg.hh.gnabar = 0.12  # Sodium conductance in S/cm2
            seg.hh.gkbar = 0.036  # Potassium conductance in S/cm2
            seg.hh.gl = 0.0003    # Leak conductance in S/cm2
            seg.hh.el = -54.3     # Reversal potential in mV
        for seg in self.dend:
            seg.hh.gnabar = 0.12  # Sodium conductance in S/cm2
            seg.hh.gkbar = 0.036  # Potassium conductance in S/cm2
            seg.hh.gl = 0.0003    # Leak conductance in S/cm2
            seg.hh.el = -54.3     # Reversal potential in mV
        # Insert passive current in the dendrite
        self.dend.insert('pas')                 
        for seg in self.dend:
            seg.pas.g = 0.001  # Passive conductance in S/cm2
            seg.pas.e = -65    # Leak reversal potential mV
        
        #if cell_type == 'e':
        self.esyn = h.ExpSyn(self.dend(0.5))
        self.esyn.tau = 2 * ms
        self.esyn.e = 0 * mV
        #else:
        self.isyn = h.Exp2Syn(self.dend(0.5))
        self.isyn.tau1 = 0.2 * ms
        self.isyn.tau2 = 10 * ms
        self.isyn.e = -80 * mV
            
        

            
            
            
def create_n_BallAndStick(n):
    """n = number of cells; r = radius of circle"""
    cells = []
    for i in range(int(0.80*n)):
        cells.append(BallAndStick(i, 50*i, 0, 0, 'e'))
    for i in range(int(0.20*n)):
        cells.append(BallAndStick(i, 50*i, -50, 0, 'i'))
    return cells


#Setup connection matrix for a ixi and number of cells

# E -> E = 1
# E -> I = 1
# I -> E = 0.5 of excitatories
# All self connections = 0

my_cells = create_n_BallAndStick(100)



    

connection_matrix = np.ones((len(my_cells),len(my_cells)))
for i in range (0,len(my_cells)):
    connection_matrix[i][i] = 0         #Self connections
    
''' Cut out matrix generation
#for each inhibitory neuron, removes connections with the first half of excitatory neurons
for i, cell in enumerate(my_cells):
    if cell._cell_type == "Inhibitory":
        inhibitory_con = []
        for j, value in enumerate(connection_matrix[i]):
            if value == 1 and my_cells[j]._cell_type == "Excitatory":
                inhibitory_con.append(j)
            if value == 1 and my_cells[j]._cell_type == "Inhibitory":
                connection_matrix[i][j] = 0
        for index in random.sample(inhibitory_con, len(inhibitory_con)//2):
            connection_matrix[i][index] = 0

'''

#Hard code in the data
import csv

with open('connection_matrix.csv', 'r', newline='') as file:
  myreader = csv.reader(file)
  for i, rows in enumerate(myreader):
      for j, item in enumerate(rows):
          connection_matrix[i][j] = int(item[0])
      



newlist = []
with open('delayI.csv', 'r') as f:
    myreader = csv.reader(f)
    for list in myreader:
        newlist.append(list)
    

    delayI = newlist[0]
    for i, value in enumerate(delayI):
        delayI[i] = float(value)


newlist = []
with open('delayE.csv', 'r') as f:
    myreader = csv.reader(f)
    for list in myreader:
        newlist.append(list)
    

    delayE = newlist[0]
    for i, value in enumerate(delayE):
        delayE[i] = float(value)





'''
connection_matrix = np.ones((len(my_cells),len(my_cells)))
for i in range (0,len(my_cells)):
    connection_matrix[i][i] = 0   
connection_matrix[8] = [1., 0., 1., 0., 0., 1., 0., 1., 0., 0.]
connection_matrix[9] = [1., 0., 0., 0., 1., 1., 1., 0., 0., 0.]
'''






#Set up stimulations

# TODO Set an iclamp to each cell

    
# set a synaptic input onto my_cells[0]
# from an artificial stim cell


syn_ = h.ExpSyn(my_cells[0].dend(0.5))

stim = h.NetStim()

stim.number = 1     #Assigned number of stimulation
stim.start = 9      #start of stimulus in ms
ncstim = h.NetCon(stim, syn_)   #Looks at potential of stim, and once that hits threshold, triggers syn_
ncstim.delay =  1 * ms           #1 ms delay to model time for NTs or signal propogation
ncstim.weight[0] = 0.04    #how strong the transmittions is from pre to post synapse
                                #Removed stim for now

syn_.tau = 2 * ms




#Connecting the cells


con_loc = []
for count_row, row in enumerate(connection_matrix):
    for count_column, value in enumerate(row):
        if value == 1:
            con_loc.append([count_row, count_column])
'''
Basically goes through the each row and saves the index of the connections (both the row(source) and column (target))
'''
        
        
syns = []
netcons = []

'''
delaysE = [10.037896272691983, 9.327329873811465, 9.936889354537374, 8.932833109642614, 9.328916471457706, 13.288206248617328, 18.279080704931317, 5.609313831481031, 12.43208810896164, 6.6185309959502865, 6.018205412096546, 9.47962326661358, 7.606579944415668, 7.270614057281771, 14.92692690923485, 13.902215186457884, 4.616788961644186, 4.208612975144392, 6.289016366332285, 5.0303570238423045, 7.307372576505347, 5.536223827074861, 13.681536610146793, 10.101942441773678, 12.091710474693269, 7.819977705461372, 5.381679270596206, 4.007839467474556, 9.355026268191592, 10.698713108464391, 9.424330611368909, 13.266900996477975, 9.120172122795607, 8.21446249930859, 15.240910762962375, 7.864246285149289, 9.993021352629063, 6.884273933656798, 7.850460182486554, 4.381035147221762, 9.419281649571172, 13.540730377429592, 9.27186962780612, 8.602280645738894, 7.325025057181831, 9.86710317390649, 7.562870653951009, 6.774783148375746, 6.555252079831622, 6.89903358612259, 4.546223222669831, 8.369821472505192, 7.365597221350233, 11.054733474288007, 13.887433444435228, 14.007677254984614, 9.51660209930858, 12.838313097958455, 3.748143671837105, 5.059239273442815, 12.70316947777893, 14.878545032704448, 13.449596601150962, 10.274272153881494, 12.394288540551761, 12.810714878319382, 13.694397320777473, 10.62230006425422, 11.27039212767095, 10.23729002902207, 15.283470858666542, 11.137311561332412]
delayECount = 0
delaysI = [9.64217834565535, 13.150316784882234, 8.61879927069436, 11.475344500981558, 10.329778837595299, 10.228503624102556, 9.929274995451383, 9.843321615639276]
delayICount = 0
'''
delayECount = 0
delayICount = 0

for row in con_loc:
    source_index, target_index = row
    source = my_cells[source_index]
    target = my_cells[target_index]
    #distance = math.sqrt((source.x - target.x)**2 + (source.y - target.y)**2)
   # syn = h.ExpSyn(target.dend(0.5))   #syn is the synapse and later appended to list syns
    if source._cell_type == "Excitatory":
        nc = h.NetCon(source.soma(0.5)._ref_v, target.esyn, sec=source.soma)    
        nc.weight[0] = random.gauss(0.1,0.0001)
        nc.delay = delayE[delayECount] 
        delayECount += 1
            #random.gauss(abs((distance)/500)+6,2)  
        
        ''' nc.delay = delaysE[delayECount]
        delayECount +=1                            #random.gauss(abs((distance)/50)+6,3)      
        delaysE.append(nc.delay)
        '''
    else: # Inhibitory
        nc = h.NetCon(source.soma(0.5)._ref_v, target.isyn, sec=source.soma) 
        nc.weight[0] = random.gauss(0.8,0.0001)  #strength of connection between pre and post cells (maximum conductance)
        nc.delay = delayI[delayICount] 
        delayICount += 1        #random.gauss(abs((distance)/500)+6,2)  
        
        '''  nc.delay = delaysI[delayICount]
        delayICount += 1                             #random.gauss(abs((distance)/50)+6,3)        #ms - delay between presynaptic cell reaching threshold and postsynaptic cell triggering (time to propogate down axon and NT)
        delaysI.append(nc.delay)
        '''

    netcons.append(nc)
  #  syns.append(syn)




'''
colors = ['g', 'b', 'r', 'y']    
fig, ax = plt.subplots()    #Sets up figure with legend and certain number of subplots (in this case 0)
plt.xlabel('t (ms)')        #Add labels to plot
plt.ylabel('v (mV)')
for cell, colour in zip(my_cells, colors):
    recording_cell = cell
    soma_v = h.Vector().record(recording_cell.soma(0.5)._ref_v)
    dend_v = h.Vector().record(recording_cell.dend(0.5)._ref_v)
    t = h.Vector().record(h._ref_t)
    h.finitialize(-65 * mV)
    h.continuerun(100 * ms)
    ax.plot(t, soma_v, f'-{colour}', label='soma(0.5)' if recording_cell == my_cells[0] else None)
    ax.plot(t, dend_v, f'--{colour}', label='dend(0.5)' if recording_cell == my_cells[0] else None)

leg = ax.legend();
'''

'''
with open('delayE', 'f') as f:
    writer = csv.writer(f)
    writer.writerows(delaysE)
'''

#Recording values

for cell in my_cells:
    #recording_cell = my_cells[0]
    cell.vrec_soma = h.Vector().record(cell.soma(0.5)._ref_v)
    cell.vrec_dend = h.Vector().record(cell.dend(0.5)._ref_v)
    
t = h.Vector().record(h._ref_t)
        

#%% Plotting
runtime = 400
soma_Vms = np.zeros(16001)
#%% Plotting
h.finitialize(-65 * mV)
h.continuerun(400 * ms)
fig = plt.figure(figsize= (7,11))
for i, cell in enumerate(my_cells):
    ax = fig.add_subplot(len(my_cells),1,i+1)
    plt.plot(t, cell.vrec_soma, label='soma(0.5)')
    plt.plot(t, cell.vrec_dend, label='dend(0.5)')
    plt.legend()
    plt.title(cell._cell_type + ' ' + str(cell._gid))
    if i<(len(my_cells) - 1):
        plt.tick_params(axis='x', labelbottom=False) 
    soma_Vms = soma_Vms + np.array(cell.vrec_soma)


soma_Vms = soma_Vms/100
fig2 = plt.figure()
plt.plot(t, soma_Vms)
plt.show()
#%%

#%%
#Showing neuron firing and times

fig = plt.figure() 
inh = 0
exc = 0
for cell in my_cells:
    if cell._cell_type == "Excitatory":
        exc += 1
    else:
        inh += 1
    print(cell._gid if cell._cell_type == "Excitatory" else -cell._gid)
    gidvec = [cell._gid if cell._cell_type == "Excitatory" else -1-cell._gid]*len(cell.spike_times)  
    plt.plot(cell.spike_times, gidvec, marker = 'o', linestyle="None", color = 'b' if cell._cell_type == "Excitatory" else 'r')    

plt.ylabel("Cell #")
plt.title("Spike Raster")
plt.xlabel("Time (ms)")
plt.yticks(range(-inh,exc))
plt.show()

#%%

soma_v = h.Vector().record(my_cells[-1].soma(0.5)._ref_v)
dend_v = h.Vector().record(my_cells[-1].dend(0.5)._ref_v)
t = h.Vector().record(h._ref_t)
h.finitialize(-65 * mV)
h.continuerun(400 * ms)

fig3= plt.figure()
plt.title('Inhibitory cell [-1]')
plt.ylabel('Potential (mV)')
plt.xlabel('Time (ms)')
plt.plot(t, soma_v, label = 'soma(0.5)')
plt.plot(t, dend_v, label = 'dend(0.5)')
plt.show()

soma_v = h.Vector().record(my_cells[0].soma(0.5)._ref_v)
dend_v = h.Vector().record(my_cells[0].dend(0.5)._ref_v)
t = h.Vector().record(h._ref_t)
h.finitialize(-65 * mV)
h.continuerun(400 * ms)

fig4= plt.figure()
plt.title('Excitatory Cell [0]')
plt.ylabel('Potential (mV)')
plt.xlabel('Time (ms)')
plt.plot(t, soma_v, label = 'soma(0.5)')
plt.plot(t, dend_v, label = 'dend(0.5)')
plt.show()








