from group import *
import AutF as aut
import networkx as nx


class FPGraphOfGroups(nx.MultiDiGraph):
    """
    A directed graph with groups assigned to each edge and vertex and inclusions of each edge group into its vertex groups.
    """
    def __init__(self, data=None, **attr):
        nx.MultiDiGraph.__init__(self, data, **attr)

        self.edgekeys={}  # dictionary that takes name of edge and gives corresponding vertices
        for edge in self.edges(keys=True):
            self.edgekeys[edge[2]]=(edge[0],edge[1])

    def __repr__(self):
        return "vertices:"+str(self.nodes(data=True))+" edges:"+str(self.edges(keys=True, data=True))+" edgekeys:"+repr(self.edgekeys)

    def __str__(self):
        nontrivialedges=[e for e in self.edges(keys=True,data=True) if e[3]['group'].gens]
        thestring="Graph:\n"
        thestring+="Vertices: "+str(self.nodes())+"\n"
        theedges=[]
        for edge in self.edges(keys=True,data=True):
            if 'label' in edge[3]:
                if edge[3]['label']:
                    theedges.append("("+str(edge[0])+", "+str(edge[1])+", "+str(edge[2])+", label="+str(edge[3]['label'])+")")
                else:
                    theedges.append("("+str(edge[0])+", "+str(edge[1])+", "+str(edge[2])+")")
            else:
                theedges.append("("+str(edge[0])+", "+str(edge[1])+", "+str(edge[2])+")")
        thestring+="Edges: "+", ".join(theedges)+"\n"
        thestring+="Vertex Groups:\n"
        thestring+="\n".join([str(v)+": "+str(self.node[v]['group']) for v in self])+"\n"
        thestring+="Nontrivial Edge Maps:\n"
        thestring+="\n".join([str(e[:3])+" origin\n"+str(self[e[0]][e[1]][e[2]]['omap'])+"\n"+str(e[:3])+" terminus\n"+str(self[e[0]][e[1]][e[2]]['tmap']) for e in nontrivialedges])
        return thestring
    
    def origin(self,edge):
        try:
            return self.edgekeys[edge][0]
        except (KeyError, TypeError):
            return edge[0]

    def terminus(self,edge):
        try:
            return self.edgekeys[edge][1]
        except (KeyError, TypeError):
            return edge[1]

    def ekey(self,edge):
        try:
            if edge in self.edgekeys:
                return edge
            else:
                return edge[2]
        except TypeError:
            return edge[2]

    def incident_edges(self, vert,**kwargs):
        return [edge for edge in self.in_edges(vert,**kwargs)+self.out_edges(vert,**kwargs)]

    def incident_edge_keys(self, vert):
        return [edge[2] for edge in self.in_edges(vert,keys=True)+self.out_edges(vert,keys=True)]
        
    def add_vertex(self, vert, vertgroup=FPGroup()):
        self.add_node(vert,{'group':vertgroup})

    def gettmap(self,edge):
        return self[self.origin(edge)][self.terminus(edge)][self.ekey(edge)]['tmap']

    def getomap(self,edge):
        return self[self.origin(edge)][self.terminus(edge)][self.ekey(edge)]['omap']

    def localgroup(self,vertoredge):
        try: # maybe vertoredge is a vertex
            return self.node[vertoredge]['group']
        except (KeyError, TypeError): # otherwise it's an edge tuple
            return self[vertoredge[0]][vertoredge[1]][vertoredge[2]]['group']

    def add_edge(self, u, v,  group=None, omap=None, tmap=None ,**kwargs):
        """
        Add an edge. If key is None a unique key will be generated. If group, omap, and tmap are None then the edge group is trivial.
        """
        if not hasattr(group, 'gens'):
            edgegroup=FPGroup()  # the trivial group
        elif group.gens==[]:
            edgegroup=group  # still the trivial group
        else:  # if the edge group is nontrivial the inclusion maps can not be trivial
            edgegroup=group
            assert(edgegroup.gens)
            assert(omap.variant_generators())
            assert(tmap.variant_generators())
        if u not in self:  # if u is not already in the graph add it with group isomorphic to the edge group
            self.add_vertex(u, edgegroup)
            vertgroup=self.node[u]['group']
            omap=Homomorphism(edgegroup,vertgroup,[vertgroup.word([i]) for i in range(1,1+len(edgegroup.gens))])
        if v not in self:
            self.add_vertex(v, edgegroup)
            vertgroup=self.node[v]['group']
            tmap=Homomorphism(edgegroup,vertgroup,[vertgroup.word([i]) for i in range(1,1+len(edgegroup.gens))])
        if omap is None:
            assert(edgegroup.gens==[])
            omap=Homomorphism(edgegroup,self.node[u]['group'],[])
        if tmap is None:
            assert(edgegroup.gens==[])
            tmap=Homomorphism(edgegroup,self.node[v]['group'],[])

        try:
            key=kwargs['key']
        except KeyError:
            counter=0
            key='e'+str(counter)
            while key in self.edgekeys:
                counter+=1
                key='e'+str(counter)

        if type(key)==list:
            key=tuple(key)
        datadict={'group': group,'omap': omap, 'tmap': tmap}
        for kw in kwargs:
            if kw!='key':
                datadict[kw]=kwargs[kw]
        nx.MultiDiGraph.add_edge(self,u,v,key,datadict)
        self.edgekeys[key]=(u,v)

    def neighbors(self, vert):
        nbs=set([])
        for e in self.out_edges(vert):
            nbs.add(e[1])
        for e in self.in_edges(vert):
            nbs.add(e[0])
        return nbs

    def remove_vertex(self,vert):
        for edge in self.edges_iter(vert, keys=True):
            del self.edgekeys[edge[2]]
        self.remove_node(vert)

    def remove_edge(self, edge):
        del self.edgekeys[edgekey]
        self.remove_edge(edge)
        assert(nx.is_connected(self))

    def change_vertex_group_by_automorphism(self, vert, aut):
        """
        Change all the incident edgemaps by composing with aut.
        """
        for edge in self.incident_edges(vert, keys=True):
	    ekey=edge[2]
            origin=self.origin(ekey)
            terminus=self.terminus(ekey)
            if vert==terminus:
                self[origin][terminus][ekey]['tmap']=compose(aut,self.gettmap(ekey))
            if vert==origin:
                self[origin][terminus][ekey]['omap']=compose(aut,self.getomap(ekey))

    def change_edge_map(self, edge, end, vertexconjugator):
        """
        Change an edge map by conjugating the image.
        """
        # End should be 0 for origin or 1 for terminus
        if end:
            self[self.origin(edge)][self.terminus(edge)][edge]['tmap']=compose(group.InnerAutomorphism(self.terminus(edge).group,vertexconjugator),self.gettmap(edge))
        else:
            self[self.origin(edge)][self.terminus(edge)][edge]['omap']=compose(group.InnerAutomorphism(self.terminus(edge).group,vertexconjugator),self.getomap(edge))

    def refine_vertex(self, vert, newgog, edgeupdate):
        """
        Refine graph of groups by replacing vertex vert with graph of groups newgog.
        """
        # edgeupdate={edge: (newvert ,newmap)}
        # e is an edge incident to vert
        # if e is a loop at vert then newvert and newmap are tuples

        # add a copy of newgog to self, with vertex and edgekeys prefixed
        for v in newgog:
            self.add_vertex(v,newgog.localgroup(v))
        for e in newgog.edges(keys=True, data=True):
            # if the edge has some other data attached we shoud remember that
            edgedata=dict([(k,e[3][k]) for k in e[3] if k!='tmap' and k!='omap' and k!='group' and k!='key'])
            edgedata['key']=e[2]
            self.add_edge(newgog.origin(e),newgog.terminus(e),newgog.localgroup(e),newgog.getomap(e),newgog.gettmap(e), **edgedata)
        # for each edge incident to vert add an edge connecting to the new subgraph
	edgestobeadded=[]
        for e in edgeupdate:
            if e[0]==vert:
                if e[1]==vert:  # e is a loop at vert
		    edgestobeadded.append(dict({"u":edgeupdate[e][0][0],"v":edgeupdate[e][0][1],"group":self.localgroup(e),"omap": edgeupdate[e][1][0],"tmap": edgeupdate[e][1][1], "key": self.ekey(e)}))
                else:  # e has only origin==vert
		    edgestobeadded.append(dict({"u": edgeupdate[e][0],"v": self.terminus(e),"group": self.localgroup(e),"omap": edgeupdate[e][1],"tmap": self.gettmap(e), "key": self.ekey(e)}))
            else:  # e has only terminus==vert
                assert(e[1]==vert)
                edgestobeadded.append(dict({"u": self.origin(e),"v": edgeupdate[e][0],"group":self.localgroup(e),"omap": self.getomap(e),"tmap": edgeupdate[e][1], "key": self.ekey(e)}))
        # remove vert and all incident edges
        self.remove_vertex(vert)
	# add the new edges
	for newedge in edgestobeadded:
	    self.add_edge(**newedge)
