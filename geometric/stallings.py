
class xgraph(object):
    """ 
    Representation is as a dictionary { g: [(gs,s),(gt,t),...],...}
    except everything is an integer.
    """
    # TODO: add switch_nodes method to switch two nodes.  When
    # folding, if we're about to delete node 0, we need to switch
    # first, so result always has a node 0.
    def __init__(self,wordlist_or_edgelist=[],as_edges=False):
        if not as_edges: #construct from words
            self.dict={0:[]}
            wordlist=wordlist_or_edgelist
            for word in wordlist:
                # Attach a loop to the base vertex:
                offset = len(self.dict.keys())-1
                vertexpath = [0]+range(1+offset,len(word)+offset)+[0]
                for i in range(len(vertexpath)-1):
                    letter = word.letters[i]
                    self.dict[vertexpath[i]].append( (vertexpath[i+1],letter) )
                    if vertexpath[i+1] not in self.dict.keys():
                        self.dict[vertexpath[i+1]]=[]
                    self.dict[vertexpath[i+1]].append( (vertexpath[i],-letter))
        else: #construct from ***positive*** edgelist
            edgelist = wordlist_or_edgelist
            self.dict={}
            for edge in edgelist:
                self.add_node(edge[0])
                self.add_node(edge[1])
                self.add_edge(edge[0],edge[1],edge[2])
    
    def __copy__(self):  # This doesn't work properly if there are isolated vertices.
        return xgraph(self.posedges(),as_edges=True)
    
    def __repr__(self):
        return str(self.dict)

    def add_node(self,vertex): # if vertex already there, do nothing
        if vertex not in self.nodes():
            self.dict[vertex]=[]

    def del_node(self,vertex):
        while self.dict[vertex]!=[]:
            target,letter = self.dict[vertex][0]
            self.del_edge(vertex,target,letter)
        del self.dict[vertex]

    def add_edge(self,source,target,label):
        self.dict[source].append((target,label))
        self.dict[target].append((source,-label))

    def del_edge(self,source,target,label):
        j = self.dict[source].index((target,label))
        del self.dict[source][j]
        k = self.dict[target].index((source,-label))
        del self.dict[target][k]
        
    def clone_node(self,v,new_node = None):
        if new_node==None:
            new_node = max(self.dict.keys())+1
        self.add_node(new_node)
        for edge in self.dict[v]:
            if edge[0]==v:
                if edge[1]>0:
                    self.add_edge(new_node,new_node,edge[1])
            else:
                self.add_edge(new_node,edge[0],edge[1])
        return new_node

    def switch_nodes(self,one,two):
        """
        switch two nodes, carrying edges along.
        """
        new = self.clone_node(one)
        self.del_node(one)
        self.clone_node(two,one)
        self.del_node(two)
        self.clone_node(new,two)
        self.del_node(new)
        


    def neighbors(self,vertex):
        """
        return as list
        """
        nlist = map(lambda x: x[0], self.dict[vertex])
        return nlist

    def distinct_neighbors(self,vertex):
        """
        return as set
        """
        nlist = map(lambda x: x[0], self.dict[vertex])
        nset = set(nlist)-set([vertex])
        return nset

    def nodes(self):
        return self.dict.keys()

    def posedges(self):
        elist = []
        for node in self.nodes():
            for edge in self.dict[node]:
                if edge[1]>0:
                    elist.append( (node,edge[0],edge[1]) )
        return elist

    def edges(self):
        elist = []
        for node in self.nodes():
            for edge in self.dict[node]:
                elist.append( (node,edge[0],edge[1]) )
        return elist

    def is_connected(self):
        previsited = set([0])
        visited = set([])
        while previsited!=set([]):
            next = previsited.pop()
            previsited = (previsited | self.distinct_neighbors(next))-visited
            visited.add(next)
        return visited==set(self.nodes())


    def fold(self):
        while foldonce(self)!=None:
            pass

def unfolded(g):
    "Return key for vertex which is not folded "
    for vertex,edges in g.dict.iteritems():
        for i in range(len(edges)-1):
            for j in range(i+1,len(edges)):
                if edges[i][1]==edges[j][1]:
        # We're going to delete one of the targets if they differ.  Make sure not to delete the source.
                    target1 = edges[i][0]; target2 = edges[j][0]
                    if target1!=target2 and target2==vertex:
                        return (vertex,j,i)
                    else:
                        return (vertex,i,j)
    return None

def foldonce(g):
    "Do a single fold (identify two edges) if possible.  Returns 1 if something happened. Returns 'None' otherwise."
    place_to_fold = unfolded(g)
    if place_to_fold==None:
        return None
    else:
        source,i,j = place_to_fold
        target1 = g.dict[source][i][0]
        target2 = g.dict[source][j][0]
        foldletter = g.dict[source][i][1]
        # First delete extra edge coming from source
        g.del_edge(source,target2,foldletter)
        # If target1=target2, we are done.
        if target1==target2:
            return g
        # We are going to delete target2, but we don't want to delete the zero node.
        if target2==0:
            g.switch_nodes(target1,target2)
            dummy = target2
            target2 = target1
            target1 = dummy
        # Now redirect edges from target2
        for vertex,letter in g.dict[target2]:
            if vertex==target2:
                # move the self loops over to target1
                g.add_edge(target1,target1,letter)
            else:
                # now move the rest of the edges
                g.add_edge(target1,vertex,letter)
        #finally, remove target2 from dictionary
        g.del_node(target2)
        return 1
    
def max_forest(g,rootvertex=0):
    "return a subgraph of g which is a maximal forest"
    vertexset = set(g.nodes())
    tree = xgraph()
    while set(tree.nodes())!=vertexset:
        diff = vertexset - set(tree.nodes())
        newvertex = diff.pop()
        possible_edges = g.dict[newvertex][:]
        while newvertex not in tree.nodes():
            newedge=possible_edges.pop()
            if newedge[0]!=newvertex:
                tree.add_node(newvertex)
                if newedge[0] not in tree.nodes():
                    tree.add_node(newedge[0])
                tree.add_edge(newvertex,newedge[0],newedge[1])
    return tree
                
def max_tree(g,rootvertex=0):
    """Return a subgraph of g which is a maximal tree."""
    vertexset = set(g.nodes())
    tree = xgraph()
    while set(tree.nodes())!=vertexset:
        targets = vertexset - set(tree.nodes())
        sources = tree.nodes()[:]
        done = False
        while not done:
            source = sources.pop()
            for edge in g.dict[source]:
                if edge[0] in targets:
                    if source not in tree.nodes():
                        tree.add_node(source)
                    if edge[0] not in tree.nodes():
                        tree.add_node(edge[0])
                    tree.add_edge(source,edge[0],edge[1])
                    done = True
                    break
    return tree
        



