# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 10:55:06 2022

@author: maria
"""

from neuron import h
h.load_file("stdlib.hoc")

class OptoLight:
    
    def __init__(self,cellref):
        self.cell = cellref
        
    def intensity(self,sec):
        return sec.irradiance_chanrhod * sec.Tx_chanrhod
    def photon_flux(self, sec):
        """ Determine the light intensity at a given section (photons/ms cm2)
        """
        return h.photons_chanrhod * sec.Tx_chanrhod #  photons/ms cm2
    def photons(self,sec):
        section_area = h.area(0.5,sec=sec) # um2
        section_intensity = self.intensity(sec) #  photons/ms cm2
    def root_section(self):
        return self.cell.soma #h.SectionRef().root
    # def build_tree(self, func,segfunc=False):
    #     """ func must act on a neuron section
    #     """
    #     from numpy import array
    #     print("-"*100)
    #     def append_data(sec, xyzdv, parent_id, connections,func,segfunc):
    #         """ Append data to xyzdv
    #         """
    #         if not segfunc: v=func(sec)
    #         n = int(h.n3d(sec=sec))
    #         for ii in range(1, n):
    #             x = h.x3d(ii,sec=sec)
    #             y = h.y3d(ii,sec=sec)
    #             z = h.z3d(ii,sec=sec)
    #             d = h.diam3d(ii,sec=sec)
    #             if segfunc:
    #                 if n==1:v=func(sec(0.5))
    #                 else:v = func(sec(ii/float(n-1)))
    #             xyzdv.append([x,y,z,d,v])
    #             child_id = len(xyzdv)-1
    #             if len(xyzdv)>1:
    #                 connections.append([child_id, parent_id])
    #             parent_id = child_id
    #         return xyzdv, connections

    #     def append_children_data(parent, parent_id, xyzdv, connections, func, segfunc):
    #         sref = h.SectionRef(sec=parent)
    #         if sref.child:
    #             for child in sref.child:
    #                 xyzdv, connections = append_data(child, xyzdv, parent_id, connections, func, segfunc)
    #                 xyzdv, connections = append_children_data(parent = child,
    #                                                           parent_id = len(xyzdv)-1,
    #                                                           xyzdv = xyzdv,
    #                                                           connections = connections,
    #                                                           func = func,
    #                                                           segfunc = segfunc)
    #         return xyzdv, connections

    #     # Find data and connections
    #     root_section = self.root_section()
    #     if segfunc:
    #         if root_section.nseg==1:
    #             v = func(root_section(0.5))
    #         else:
    #             v = func(root_section(0.0))
    #     else:
    #         v=func(root_section)
    #     xyzdv = [[h.x3d(0,sec=root_section),h.y3d(0,sec=root_section),h.z3d(0,sec=root_section),h.diam3d(0,sec=root_section),v]]
    #     xyzdv, connections = append_data(root_section, xyzdv, 0, [],func,segfunc)
    #     xyzdv, connections = append_children_data(root_section,len(xyzdv)-1,xyzdv,connections,func,segfunc)
    #     self.xyzdv = array(xyzdv)
    #     self.connections = array(connections)
    # def move(self, xyz, move_mlab=False):
    #     """ Move visualization and cell """
    #     from neuron import h
    #     if move_mlab:
    #         if self.mlab_cell:
    #             self.mlab_cell.mlab_source.x = self.mlab_cell.mlab_source.x + xyz[0]
    #             self.mlab_cell.mlab_source.y = self.mlab_cell.mlab_source.y + xyz[1]
    #             self.mlab_cell.mlab_source.z = self.mlab_cell.mlab_source.z + xyz[2]
    #     tree = h.SectionList()
    #     tree.wholetree(sec=self.root)
    #     for sec in tree:
    #         for ii in range(h.n3d(sec=sec).__int__()):
    #             x=h.x3d(ii,sec=sec)
    #             y=h.y3d(ii,sec=sec)
    #             z=h.z3d(ii,sec=sec)
    #             d=h.diam3d(ii,sec=sec)
    #             h.pt3dchange(ii,x+float(xyz[0]),y+float(xyz[1]),z+float(xyz[2]),d)
    # def retrieve_coordinates(self, sec):
    #     xyzds = []
    #     for ii in range(int(h.n3d(sec=sec))):
    #         xyzds.append([h.x3d(ii,sec=sec),
    #                       h.y3d(ii,sec=sec),
    #                       h.z3d(ii,sec=sec),
    #                       h.diam3d(ii,sec=sec)])
    #     return xyzds
    # def plot(self, func, scaling = 1, segfunc=False, clim=None,cmap=None):
    #     """ plot cell in matplotlib line plot collection
    #     """
    #     from numpy import array, linspace
    #     from matplotlib.collections import LineCollection
    #     from matplotlib import pyplot
    #     self.build_tree(func,segfunc)
    #     pts   = self.xyzdv[:,:2]
    #     edges = self.connections
    #     diam  = self.xyzdv[:,3]
    #     data  = self.xyzdv[:,4]
    #     print("DATA RANGE: ",data.min(),data.max())
    #     # Define colors
    #     if not cmap:
    #         from matplotlib.cm import jet as cmap
    #     if not clim:
    #         clim=[data.min(),data.max()]
    #     a = (data - clim[0])/(clim[1]-clim[0])
    #     # Define line segments
    #     segments = []
    #     for edge in edges:
    #         segments.append([pts[edge[0],:], pts[edge[1],:]])
    #     # Build Line Collection
    #     collection = LineCollection(segments = array(segments),
    #                                 linewidths = diam*scaling,
    #                                 colors=cmap(a))
    #     collection.set_array(data)
    #     collection.set_clim(clim[0], clim[1])
    #     pyplot.gca().add_collection(collection,autolim=True)
    #     pyplot.axis('equal')
    #     return collection
    # def channels_in_list(self,seclist):
    #     channels = 0
    #     for sec in seclist:
    #         if h.ismembrane('chanrhod',sec = sec):
    #             for seg in sec:
    #                 rho  = seg.channel_density_chanrhod/1e8 # 1/cm2 --> 1/um2
    #                 area = h.area(seg.x, sec=sec) # um2
    #                 n = rho * area
    #                 channels += n
    #     return channels
    def area_in_list(self,seclist):
        area = 0
        for sec in seclist:
            if h.ismembrane('chanrhod',sec = sec):
                for seg in sec:
                    area += h.area(seg.x, sec=sec) # um2
        return area
    def illuminated_area_in_list(self,seclist,Tx_threshold = 0.001):
        area = 0
        for sec in seclist:
            if h.ismembrane('chanrhod',sec = sec):
                for seg in sec:
                    if seg.Tx_chanrhod>Tx_threshold:
                        area += h.area(seg.x, sec=sec) # um2
        return area
    def open_channels_in_list(self,seclist):
        open_channels = 0
        for sec in seclist:
            if h.ismembrane('chanrhod', sec = sec):
                for seg in sec:
                    rho  = seg.channel_density_chanrhod/1e8 # 1/cm2 --> 1/um2
                    area = h.area(seg.x, sec=sec) # um2
                    try:
                        f_open = seg.o2_chanrhod + seg.o1_chanrhod # open fraction # 4 state model
                    except:
                        f_open = seg.o1_chanrhod # open fraction # 3 state model
                    n = f_open * rho * area
                    open_channels += n
        return open_channels
    def get_axonal(self):
        """ Additional iseg compartment
        """
        secs = [h.hill]
        secs.extend([sec for sec in h.ais])
        secs.extend([h.nakeaxon])
        secs.extend([sec for sec in h.myelin])
        secs.extend([sec for sec in h.node])
        return secs
    def get_axonal_channels(self):
        return self.channels_in_list(self.axonal)
    def get_open_axonal_channels(self):
        return self.open_channels_in_list(self.axonal)
    def get_dendritic(self):
        secs = [sec for sec in h.somatodendritic]
        #secs.extend([h.hill])
        return secs
    def get_dendritic_channels(self):
        return self.channels_in_list(self.dendritic)
    def get_open_dendritic_channels(self):
        return self.open_channels_in_list(self.dendritic)
    def get_apical_tuft(self):
        """ Return a list of all sections which make up the apical tuft, starting
        at the branch point
        """
        secs=[]
        for ii in range(23,len(h.dend11)):
            secs.append(h.dend11[ii])
        return secs
    def get_apical_shaft(self):
        """ Return the sections which compose the apical shaft
        """
        secs=[]
        for ii in [0,4,10,16,18,20,22]:
            secs.append(h.dend11[ii])
        return secs
    def get_basilar_tuft(self):
        """ Return the dendritic sections which compose the basilar tuft
        """
        secs=[]
        for ii,dendrite in enumerate((h.dend1,h.dend2,h.dend3,h.dend4,h.dend5,h.dend6,
                                      h.dend7,h.dend8,h.dend9,h.dend10,h.dend11)):
            if ii==10:
                for jj,sec in enumerate(dendrite):
                    if jj < 23: # Apical Tuft
                        if jj not in [0,4,10,16,18,20,22]: # Apical Shaft
                            secs.append(sec)
            else:
                for sec in dendrite:
                    secs.append(sec)
        return secs
    def get_somatic(self):
        """ Return the sections which compose the Soma
        """
        return [h.soma]
    def set_density_distribution(self, distribution=0.5, n_channels = 1e7):
        """ Set density in dendritic compartments
        distribution: 0.0 - Higher Somatic density
                      0.5 - Uniform distribution
                      1.0 - Higher Apical density
        """
        # Find the maximal distance between the soma and all dendrities
        max_distance = 0
        for sec in self.dendritic:
            for seg in sec:
                max_distance = max([max_distance, self.seg_section_distance(seg)])

        for sec in self.dendritic:
            for seg in sec:
                distance=self.seg_section_distance(seg)
                s = (max_distance-distance)/(max_distance) # soma centric weighting
                a = (distance)/(max_distance)              # apical centric weighting
                W = distribution
                seg.channel_density_chanrhod  =  s*(1-W) + a*W
        scale = n_channels/self.dendritic_channels
        for sec in self.dendritic:
            for seg in sec:
                seg.channel_density_chanrhod = scale*seg.channel_density_chanrhod
        assert 0.001 > (n_channels - self.dendritic_channels)
    def get_open_channels(self):
        return self.open_channels_in_list(h.allsec())
    def get_total_channels(self):
        return self.channels_in_list(h.allsec())
    def get_icat(self):
        """ Determine the total amount of channelrhodopsin current in the cell
        """
        icat = 0
        for sec in h.allsec():
            if h.ismembrane('chanrhod',
                            sec = sec):
                for seg in sec:
                    i = seg.icat_chanrhod # (mA/cm2)
                    area = h.area(seg.x, sec=sec)/1e8      # cm2
                    icat += area * i # mA
        return icat