#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 30 16:50:22 2022

@author: Aidan Lee, Jeffrey Wen, Mugil Shanmugam
"""
from neuron import h, gui
from neuron.units import ms, mV
h.load_file('stdrun.hoc')
import numpy as np
import matplotlib.pyplot as plt
import random
import math
import time
import pandas as pd

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
            self.stimClampE.delay = 5
            self.stimClampE.dur = 400
            self.stimClampE.amp = 0.032
        else:
            self.stimClampI = h.IClamp(self.soma(0.5))
            self.stimClampI.delay = 5
            self.stimClampI.dur = 400
            self.stimClampI.amp = 0.032
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

num_cells = 100
my_cells = create_n_BallAndStick(num_cells)

connection_matrix = np.ones((num_cells, num_cells))
for i in range (0,len(my_cells)):
    connection_matrix[i][i] = 0         #Self connections
    
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

syn_ = h.ExpSyn(my_cells[0].dend(0.5))
stim = h.NetStim()
stim.number = 1     #Assigned number of stimulation
stim.start = 9      #start of stimulus in ms
ncstim = h.NetCon(stim, syn_)   #Looks at potential of stim, and once that hits threshold, triggers syn_
ncstim.delay =  1 * ms           #1 ms delay to model time for NTs or signal propogation
ncstim.weight[0] = 0.04    #how strong the transmittions is from pre to post synapse
syn_.tau = 2 * ms

con_loc = []
for count_row, row in enumerate(connection_matrix):
    for count_column, value in enumerate(row):
        if value == 1:
            con_loc.append([count_row, count_column])
              
syns = []
netcons = []
delays_list = []
for row in con_loc:
    source_index, target_index = row
    source = my_cells[source_index]
    target = my_cells[target_index]
    distance = math.sqrt((source.x - target.x)**2 + (source.y - target.y)**2)
   # syn = h.ExpSyn(target.dend(0.5))   #syn is the synapse and later appended to list syns
    if source._cell_type == "Excitatory":
        nc = h.NetCon(source.soma(0.5)._ref_v, target.esyn, sec=source.soma)    
        nc.weight[0] = 0.1  
        nc.delay = random.gauss(abs((distance)/100)+6,1)
        delays_list.append(nc.delay)
    else: # Inhibitory
        nc = h.NetCon(source.soma(0.5)._ref_v, target.isyn, sec=source.soma) 
        nc.weight[0] = 0.8  #strength of connection between pre and post cells (maximum conductance)
        nc.delay = random.gauss(abs((distance)/100)+6,1)  
        delays_list.append(nc.delay)
    netcons.append(nc)
    
#Recording values
for cell in my_cells:
    cell.vrec_soma = h.Vector().record(cell.soma(0.5)._ref_v)
    cell.vrec_dend = h.Vector().record(cell.dend(0.5)._ref_v) 
t = h.Vector().record(h._ref_t)
        
runtime = 400
soma_Vms = np.zeros(16001)
dend_Vms = np.zeros(16001)
#%% 
# Plotting

h.finitialize(-65 * mV)
h.continuerun(400 * ms)
#fig = plt.figure(figsize= (7,11))
for i, cell in enumerate(my_cells):
    #ax = fig.add_subplot(len(my_cells),1,i+1)
    #plt.plot(t, cell.vrec_soma, label='soma(0.5)')
    #plt.plot(t, cell.vrec_dend, label='dend(0.5)')
    #plt.legend()
    #plt.title(cell._cell_type + ' ' + str(cell._gid))
    #if i<(len(my_cells) - 1):
        #plt.tick_params(axis='x', labelbottom=False) 
    soma_Vms = soma_Vms + np.array(cell.vrec_soma)
    dend_Vms = dend_Vms + np.array(cell.vrec_dend)


soma_Vms = soma_Vms/100
dend_Vms = dend_Vms/100
fig2 = plt.figure()
plt.plot(t, soma_Vms)
plt.plot(t, dend_Vms)
plt.show()
#%%
# FFT time

from numpy.fft import fft, ifft
sr = 1000/h.dt
X = fft(soma_Vms - np.mean(soma_Vms))
N = len(X)
n = np.arange(N)
T = N/sr
soma_freq = n/T 

plt.figure(figsize = (12, 6))
plt.subplot(121)

plt.plot(soma_freq, np.abs(X), 'b')
plt.title("Soma FFT")
plt.xlabel('Freq (Hz)')
plt.ylabel('FFT Amplitude |X(freq)|')
plt.xlim(5, 100)
plt.ylim(0,20000)


sr = 1000/h.dt
X = fft(dend_Vms - np.mean(dend_Vms))
N = len(X)
n = np.arange(N)
T = N/sr
dend_freq = n/T 

plt.figure(figsize = (12, 6))
plt.subplot(121)

plt.plot(dend_freq, np.abs(X), 'b')
plt.title("Dendrite FFT")
plt.xlabel('Freq (Hz)')
plt.ylabel('FFT Amplitude |X(freq)|')
plt.xlim(5, 100)
plt.ylim(0,20000)
#%%
# Implementing bandpower

beta_range = []
for frequency in soma_freq:
    if frequency >= 19 and frequency <= 21:
        beta_range.append(frequency)

bandpower = sum(fft(beta_range)) / len(beta_range)
print(bandpower)
#%%
# Showing neuron firing and times

fig = plt.figure() 
inh = 0
exc = 0
for cell in my_cells:
    if cell._cell_type == "Excitatory":
        exc += 1
    else:
        inh += 1
    #print(cell._gid if cell._cell_type == "Excitatory" else -cell._gid)
    gidvec = [cell._gid if cell._cell_type == "Excitatory" else -1-cell._gid]*len(cell.spike_times)  
    plt.plot(cell.spike_times, gidvec, marker = 'o', linestyle="None", color = 'b' if cell._cell_type == "Excitatory" else 'r')    

plt.ylabel("Cell #")
plt.title("Spike Raster")
plt.xlabel("Time (ms)")
plt.yticks(range(-inh,exc))
plt.show()
#%%
#Saving parameters and results

answer = input("Do you want to save the data?[y/n] ")
print("Note that you need to save the plots manually. Find the saved data where your spyder files are saved or check your downloads.")
if answer == "y":
    con_loc_start = []
    con_loc_end = []
    for i in range(0, len(con_loc)):
        con_loc_start.append(con_loc[i][0])
        con_loc_end.append(con_loc[i][1])
    
    df = pd.DataFrame({'DELAYS': delays_list + [""] * 7281,
                       'CM_STA': con_loc_start + [""] * 7281,
                       'CM_END': con_loc_end + [""] * 7281,
                       'E0 S': my_cells[0].vrec_soma,
                       'E0 D': my_cells[0].vrec_dend,
                       'I0 S': my_cells[80].vrec_soma,
                       'I0 D': my_cells[80].vrec_dend})

    curr = time.ctime(time.time())
    file_name = "rundata_" + curr[4:7] + curr[9:10] + "_" + curr[11:13] + curr[14:16] + ".csv" 
        # saves to a file that includes the date and time
    df.to_csv(file_name, index=False)
