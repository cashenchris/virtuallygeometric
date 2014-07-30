import networkx as nx

class OMGError(Exception):
    """Exception raised for errors in the input.

    Attributes:
        expr -- expression generating the error
        msg  -- explanation of the error
    """

    def __init__(self, expr, msg):
        self.expr = expr
        self.msg = msg

# class OrderedMultiGraph is a multigraph that also fixes an ordering of the incident edges at each vertex.
# This allows for combinig graphs by "splicing" by specifying graphs G1 and G2, vertices v1 in G1 and v2 in G2
# of the same valence b, and a splicemap which is list of length b containing the numbers 0 thru b-1 specifying
# which edge of G1 incident to v1 should be spliced to which edge of G2 incident to v2, there the numbers in
# the splicemap refer to the index of the edge in the ordered lists of edges incident to v1 and v2

class OrderedMultiGraph(nx.MultiGraph):
    """ 
    A networkx.MultiGraph that remembers an ordering of the incident edges at each vertex.
    """
    # Each vertex has data an ordered lsit of incident edges.
    # Additionally add a dictionary edgedict={edgekey:(origin, terminus),...}
    


    def __init__(self, data=None, edgeorders=None, **attr):
        self.edgekeys={}
        nx.MultiGraph.__init__(self, data, **attr)
        for edge in self.edges(keys=True):
            self.edgekeys[edge[2]]=(edge[0],edge[1])
        if edgeorders==None:
            for vert in self.nodes_iter():
                self.node[vert]['edgeorder']=[edge[2] for edge in self.edges_iter(vert, keys=True)]
        else:
            for vert in self.nodes_iter():
                if vert in edgeorders:
                    self.node[vert]['edgeorder']=edgeorders[vert]
                else:
                    self.node[vert]['edgeorder']=[edge[2] for edge in self.edges_iter(vert, keys=True)]
            
    def __repr__(self):
        return "vertices:"+str(self.nodes(data=True))+" edges:"+str(self.edges(keys=True, data=True))+" edgekeys:"+repr(self.edgekeys)      
        
    def origin(self,edge):
        return self.edgekeys[edge][0]

    def terminus(self,edge):
        return self.edgekeys[edge][1]

    def incident_edges(self, vert):
        return self.node[vert]['edgeorder']

    def edge_order(self, vert):
        return self.node[vert]['edgeorder']
        
    def add_vertex(self,vert):
        self.add_node(vert)
        self.node[vert]['edgeorder']=[]

    def remove_vertex(self,vert):
        for edge in self.edges_iter(vert, keys=True):
            opvert=edge[1]
            edgekey=edge[2]
            self.node[opvert]['edgeorder']=self.edge_order(opvert)[:self.edge_order(opvert).index(edgekey)]+self.edge_order(opvert)[1+self.edge_order(opvert).index(edgekey):]
            del self.edgekeys[edgekey]
        self.remove_node(vert)


    def add_edge(self, u, v, key=None, uposition=-1, vposition=-1, attr_dict=None, **attr):
        # Adds an edge from u to v with given key. If no key given a unique one will
        # be generated. The edge is also inserted into the lists of incident edges to u and v
        # before the current entries in uposition and vposition, respectively.
        # Note that if u=v and then the greater of uposition and vposition will
        # end up one spot higher that given.

        # If no key is given we need to generate one and pass the information
        # to networkx.MultiGraph.add_edge
        # Start with key=0 and if necessary count up until find unused key
        # Note: we want key to be unique among all edge keys, not just u to v edges
        # this is more restrictive than in networkx.MultiGraph.add_edge
        if u in self.nodes():
            uisnew=False
        else:
            uisnew=True
        if v in self.nodes():
            visnew=False
        else:
            visnew=True
            
        if key==None:
            key=0
            while key in self.edgekeys:
                key+=1
        nx.MultiGraph.add_edge(self, u, v, key, attr_dict, **attr)
        self.edgekeys[key]=(u,v)
        if uisnew:
            self.node[u]['edgeorder']=[key]
        else:
            if uposition!=-1:
                self.node[u]['edgeorder'].insert(uposition,key)
            else:
                self.node[u]['edgeorder']+=[key]
                
        if visnew:
            self.node[v]['edgeorder']=[key]
        else:
            if vposition==-1:
                self.node[v]['edgeorder']+=[key]
            elif u==v and uposition<=vposition:
                self.node[v]['edgeorder'].insert(vposition+1,key)
            else:
                self.node[v]['edgeorder'].insert(vposition,key)
        
    def remove_edge(self,edge):
        u=self.origin(edge)
        ulist=self.edge_order(u)
        uposition=ulist.index(edge)
        self.node[u]['edgeorder']=ulist[:uposition]+ulist[1+uposition:]
        
        v=self.terminus(edge)
        vlist=self.edge_order(v)
        vposition=vlist.index(edge)
        self.node[v]['edgeorder']=vlist[:vposition]+vlist[1+vposition:]
        
        del self.edgekeys[edge]
        self.remove_edge(u, v, edge)

    def change_edge_key(self,oldkey, newkey):
        u=self.origin(oldkey)
        uposition=self.edge_order(u).index(oldkey)
        
        v=self.terminus(oldkey)
        vposition=self.edge_order(v).index(oldkey)
        
        del self.edgekeys[oldkey]
        self.remove_edge(u, v, oldkey)
        self.add_edge(u,v,newkey,uposition,vposition)

    def opposite_end(self,edge,vert):
        if vert==self.terminus(edge):
            return self.origin(edge)
        elif vert==self.origin(edge):
            return self.terminus(edge)
        else:
            try:
                raise OMGError(self.__repr__()+'.opposite_end('+str(edge)+','+str(vert)+')',str(vert)+' is not a vertex of '+str(edge))
            except OMGError as e:
                print 'Error in', e.expr,': ',e.msg

    def valence(self, vertex):
        return self.degree(vertex)
    
    def is_connected(self):
        if len(self)==0:
            return False 
        else:
            return nx.is_connected(self)
    
    def connected_component(self,vertex):
        return nx.node_connected_component(self,vertex)

    def connected_components(self):
        return [concom for concom in nx.connected_components(self)] # nx.connected_components changed to return generator
    
    def connected_component_minus_a_vertex(self,vertex1, vertex2):
        """
        find the connected component of vertex1 in the graph-vertex2
        """
        ondeck = set([vertex1])
        seen = set()
        while ondeck!=set():
            nextvert = ondeck.pop()
            nextneighbors = set(self.neighbors(nextvert))-set([vertex2])
            ondeck |= nextneighbors
            ondeck -= seen
            seen.add(nextvert)
        return seen
    
    def connected_components_minus_two_vertices(self,vertex1, vertex2):
        """
        Connected components of W-{v1,v2}
        """
        G=self.copy()
        G.remove_nodes_from([vertex1,vertex2])
        return [concom for concom in nx.connected_components(G)]  # nx.connected_components changed to return generator
      
    def is_cut_vertex(self,vertex):
        """
        check if given vertex is a cut vertex
        """
        if self.incident_edges(vertex)==[]:
            return True
        else:
            aneighbor=self.neighbors(vertex)[0]
            return self.connected_component_minus_a_vertex(aneighbor, vertex)!=(set(self.nodes())-set([vertex]))
            
    def find_cut_vertex(self):
        for vert in self.nodes_iter():
            if self.is_cut_vertex(vert):
                return vert
        return None                
            
    def is_circle(self):
        if any(self.valence(vertex)!=2 for vertex in self.nodes_iter()):
            return False
        elif not self.is_connected():
            return False
        else:
            return True
      
def splice(G1, G2, v1, v2, splicemap,G1prefix=(),G2prefix=(),lookforisolatedvertices=False):
    """
    Splice together two ordered multigraphs at vertices vi in Gi according to splicemap, prefixing names from Gi with Giprefix
    """
    # splicemap maps edge incident to v1 to edge incident to v2
    assert(G1.valence(v1)==G2.valence(v2))
    def rename(prefix,base):
        return prefix+base if type(base)==tuple else prefix+(base,)
    G1edges=[(rename(G1prefix,edge[0]), rename(G1prefix,edge[1]), rename(G1prefix,edge[2]), edge[3]) for edge in G1.edges_iter(keys=True, data=True)]
    G2edges=[(rename(G2prefix,edge[0]), rename(G2prefix,edge[1]), rename(G2prefix,edge[2]), edge[3]) for edge in G2.edges_iter(keys=True, data=True)]
    G1edgeorders= dict((rename(G1prefix,vert),[rename(G1prefix,edge) for edge in G1.node[vert]['edgeorder']]) for vert in G1.nodes_iter()) # {rename(G1prefix,vert):[rename(G1prefix,edge) for edge in G1.node[vert]['edgeorder']] for vert in G1.nodes_iter()}
    G2edgeorders=dict((rename(G2prefix,vert),[rename(G2prefix,edge) for edge in G2.node[vert]['edgeorder']]) for vert in G2.nodes_iter())# {rename(G2prefix,vert):[rename(G2prefix,edge) for edge in G2.node[vert]['edgeorder']] for vert in G2.nodes_iter()}
    newedgeorders={}
    newedgeorders.update(G2edgeorders)
    newedgeorders.update(G1edgeorders)
    for i in range(G1.valence(v1)):
        G1edge=G1.incident_edges(v1)[i] # an edge in G1 incident to v1
        G2edge=G2.incident_edges(v2)[splicemap[i]] # the edge of G2 incident to v2 to which we splice G1edge
        newedgename=rename(G1prefix,G1edge) # new edge will take the name from G1 with the G1prefix
        neworigin=rename(G1prefix,G1.opposite_end(G1edge, v1))
        newterminus=rename(G2prefix,G2.opposite_end(G2edge, v2))
        #k={k for k in range(len(G1edges)) if G1edges[k][2]==newedgename}.pop() # position of the edge named newedgename in the edgelist G1edges
        for k in range(1+len(G1edges)):
            if G1edges[k][2]==newedgename:
                break
        G1edges[k]=(neworigin,newterminus,newedgename,G1edges[k][3]) # IndexError here would mean we didn't find newedgename in G1edges in previous for loop, which should never happen
        # G1edgeorders don't need to change
        G2edgeorders[newterminus][G2edgeorders[newterminus].index(rename(G2prefix,G2edge))]=newedgename # G2edgeorders we need to replace G2edge with newedgename
        # remove the unneeded edge from G2edges
        #k={k for k in range(len(G2edges)) if G2edges[k][2]==rename(G2prefix,G2edge)}.pop() # position of the edge G2edge in the edgelist G2edges
        for k in range(1+len(G2edges)):
            if G2edges[k][2]==rename(G2prefix,G2edge):
                break
        del G2edges[k] # IndexError here would mean we didn't find G2edge in G2edges in previous for loop, which should never happen
    del newedgeorders[rename(G1prefix,v1)]
    del newedgeorders[rename(G2prefix,v2)]
    
    newgraph=OrderedMultiGraph(G1edges+G2edges,newedgeorders)

    # We have missed any isolated vertices
    if lookforisolatedvertices:
        for vert in G1.nodes_iter():
            if G1.valence(vert)==0:
                newgraph.add_vertex(rename(G1prefix,vert))
        for vert in G2.nodes_iter():
            if G2.valence(vert)==0:
                newgraph.add_vertex(rename(G2prefix,vert))
    
    return newgraph
