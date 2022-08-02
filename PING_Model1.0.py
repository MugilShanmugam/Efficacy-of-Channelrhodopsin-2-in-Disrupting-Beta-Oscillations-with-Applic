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
        for seg in self.soma:
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
        cells.append(BallAndStick(i, 50*i, 0, 0, 'i'))
    return cells


#Setup connection matrix for a 4x4 and number of cells
connection_matrix = np.ones((5,5))
# E -> E = 1
# E -> I = 1
# I -> E = 0.5 of excitatories
# All self connections = 0

my_cells = create_n_BallAndStick(5)
h.topology()

for i in range (0,5):
    connection_matrix[i][i] = 0         #Self connections
    

#for each inhibitory neuron, removes connections with the first half of excitatory neurons
for i, cell in enumerate(my_cells):
    if cell._cell_type == "Inhibitory":
        inhibitory_con = []
        for j, value in enumerate(connection_matrix[i]):
            if value == 1:
                inhibitory_con.append(j)
                print (j)
        for index in random.sample(inhibitory_con, len(inhibitory_con)//2):
                 connection_matrix[i][index] = 0
                
            

print(connection_matrix)
    

#Set up stimulations

# TODO Set an iclamp to each cell

# set a synaptic input onto my_cells[0]
# from an artificial stim cell
syn_ = h.ExpSyn(my_cells[0].dend(0.5))

stim = h.NetStim()
stim.number = 1     #Assigned number of stimulation
stim.start = 9      #start of stimulus in ms
ncstim = h.NetCon(stim, syn_)   #Looks at potential of stim, and once that hits threshold, triggers syn_
ncstim.delay = 1 * ms           #1 ms delay to model time for NTs or signal propogation
ncstim.weight[0] = 0.04       #how strong the transmittions is from pre to post synapse

syn_.tau = 2 * ms




#Connecting the cells


con_loc = []
for count_row, row in enumerate(connection_matrix):
    for count_column, value in enumerate(row):
        if value == 1:
            con_loc.append([count_row, count_column])
'''
Basically goes through the each row and 
'''
        
        
syns = []
netcons = []
for row in con_loc:
    source_index, target_index = row
    source = my_cells[source_index]
    target = my_cells[target_index]
   # syn = h.ExpSyn(target.dend(0.5))   #syn is the synapse and later appended to list syns
    if source._cell_type == "Excitatory":
        nc = h.NetCon(source.soma(0.5)._ref_v, target.esyn, sec=source.soma)    #Still need help understanding netcons
        nc.weight[0] = 0.055    #strength of connection between pre and post cells (maximum conductance)
        nc.delay = 5            #ms - delay between presynaptic cell reaching threshold and postsynaptic cell triggering (time to propogate down axon and NT)
                                #if planning to try to shorten frequency of gamma waves in ping model, probably mess around with tau (decay factor of synapse) rather than delay
    else: # Inhibitory
        nc = h.NetCon(source.soma(0.5)._ref_v, target.isyn, sec=source.soma)    #Still need help understanding netcons
        nc.weight[0] = 0.055    #strength of connection between pre and post cells (maximum conductance)
        nc.delay = 5            #ms - delay between presynaptic cell reaching threshold and postsynaptic cell triggering (time to propogate down axon and NT)
                                #if planning to try to shorten frequency of gamma waves in ping model, probably mess around with tau (decay factor of synapse) rather than delay

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



#Recording values
print (my_cells)
for cell in my_cells:
    #recording_cell = my_cells[0]
    cell.vrec_soma = h.Vector().record(cell.soma(0.5)._ref_v)
    cell.vrec_dend = h.Vector().record(cell.dend(0.5)._ref_v)
    
t = h.Vector().record(h._ref_t)
        

#%% Plotting
h.finitialize(-65 * mV)
h.continuerun(100 * ms)
fig = plt.figure(figsize= (7,11))
for i, cell in enumerate(my_cells):
    ax = fig.add_subplot(len(my_cells),1,i+1)
    plt.plot(t, cell.vrec_soma, label='soma(0.5)')
    plt.plot(t, cell.vrec_dend, label='dend(0.5)')
    plt.legend()
    plt.title(cell._cell_type + ' ' + str(cell._gid))
    if i<(len(my_cells) - 1):
        plt.tick_params(axis='x', labelbottom=False) 

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











