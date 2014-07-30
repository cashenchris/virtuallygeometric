#!/usr/bin/env python

## usage:  ./extract_and_plot2.py rank gvgfilename fullwordfilename
# outputs graphs in files "rank'rank'gvgfull.pdf" and "rank'rank'loggvg.pdf" 


from scipy.stats import linregress
import matplotlib.pyplot as plt
import numpy as np
from scipy import linalg, log, exp, sqrt
import sys

####

rank=int(sys.argv[1])
gvgfilename=sys.argv[2]
try:
    fullwordfilename=sys.argv[3]
except:
    fullwordfilename=None
####



def extractgvgdata(filename):
    f=open(filename)
    lines = list(f)
    f.close()
    data = []
    for line in lines:
        sline=line.split()
        if len(sline)==4:
            intline=[int(entry) for entry in sline]
            data.append(intline)

    lengths = [line[0] for line in data]
    g_file = [line[1] for line in data]
    vg_file = [line[2] for line in data]
    nvg_file = [line[3] for line in data]

    gs=[]
    gvgs=[]
    trials=[]
    for i in range(len(lengths)):
        trials.append(g_file[i]+vg_file[i]+nvg_file[i])
        gs.append(float(g_file[i])/trials[i])
        gvgs.append(float(g_file[i]+vg_file[i])/trials[i])
    return lengths,gs,gvgs,trials

def extractfullworddata(filename):
    f=open(filename)
    lines = list(f)
    f.close()
    data = []
    for line in lines:
        sline=line.split()
        if len(sline)==3:
            intline=[int(entry) for entry in sline]
            data.append(intline)

    lengths = [line[0] for line in data]
    full_file = [line[1] for line in data]
    notfull_file = [line[2] for line in data]

    fullrate=[]
    trials=[]
    for i in range(len(lengths)):
        trials.append(full_file[i]+notfull_file[i])
        fullrate.append(float(full_file[i])/trials[i])
    return lengths,fullrate,trials







def weightedexpapproximation(xs,ys,trials,outputdomain=None):
    """
    Returns best approxiation of ys as y=exp(mx+b)
    """
    W=np.diag(np.array([log(ys[i]/(ys[i]-sqrt(ys[i]*(1-ys[i])/float(trials[i])))) for i in range(len(xs))]))
    X=np.array([[1,x] for x in xs])
    Y=np.array([[log(y)] for y in ys])
    theapprox=linalg.solve(X.T.dot(linalg.inv(W).dot(X)),X.T.dot(linalg.inv(W).dot(Y)))
    intercept=theapprox[0,0]
    slope=theapprox[1,0]
    if outputdomain is not None:
        zs=[min(1,max(0,exp(slope*x+intercept))) for x in outputdomain]
    else:
        zs=[min(1,max(0,exp(slope*x+intercept))) for x in xs]
    return zs, slope, intercept

def weightedapproximatebottomhalf(xs,ys,trials,outputdomain=None):
    for halfindex in range(len(ys)):
        if ys[halfindex]<.5:
            break
    else:
        raise RuntimeError("data does not go below .5")
    if outputdomain is not None:
        return weightedexpapproximation(xs[halfindex:],ys[halfindex:],trials[halfindex:],outputdomain)
    else:
        return weightedexpapproximation(xs[halfindex:],ys[halfindex:],trials[halfindex:],xs)



def getWeightedDataPlot(rank,gvgfile,fullfile=None):
    """
    Returns fig with plot of data.
    gvgfile is output of VGexperiment.py
    optional fullfile is output of fullwordexperiment.py
    Also plots exponential approximation to the portion of each data series that comes after the perentage had dipped below 50% for the first time.
    """
    gvglengths,gs,gvgs,gvgtrials=extractgvgdata(gvgfile)
    if fullfile:
        fulllengths,fullrate,fulltrials=extractfullworddata(fullfile)
        notfullrate=[1.0-x for x in fullrate]
    else:
        fulllengths=[0]
    lengths=[x for x in range(1,1+max(fulllengths+gvglengths))]
    approxg,gslope,gintercept=weightedapproximatebottomhalf(gvglengths,gs,gvgtrials,lengths)
    approxgvg,gvgslope,gvgintercept=weightedapproximatebottomhalf(gvglengths,gvgs,gvgtrials,lengths)
    if fullfile:
        approxfull,fullslope,fullintercept=weightedapproximatebottomhalf(fulllengths,notfullrate,fulltrials,lengths)
    fig,ax = plt.subplots()
    plt.xlabel('Word length')
    plt.ylabel('Proportion')
    plt.title('Random rank '+str(rank)+' words')
    plt.axis([0,lengths[-1],0,1])
    if fullfile:
        plt.scatter(fulllengths,notfullrate,c='g')
        plt.plot(lengths,approxfull,color='green',linestyle="--",label="not full          $\sim \exp("+"%.3f"%(fullslope)+"x+"+"%.3f"%(fullintercept)+")$")
    plt.scatter(gvglengths,gvgs,c='r')
    plt.plot(lengths,approxgvg,color='r',linestyle="dotted",label="v. geometric $\sim \exp("+"%.3f"%(gvgslope)+"x+"+"%.3f"%(gvgintercept)+")$")
    plt.scatter(gvglengths,gs,c='b')
    plt.plot(lengths,approxg,c='b',label="geometic      $\sim \exp("+"%.3f"%(gslope)+"x+"+"%.3f"%(gintercept)+")$")
    ax.legend()
    return fig

def weightedlogapproximation(xs,ys,trials,outputdomain=None):
    """
    Returns best approxiation of ys as y=exp(mx+b)
    """
    W=np.diag(np.array([log(ys[i]/(ys[i]-sqrt(ys[i]*(1-ys[i])/float(trials[i])))) for i in range(len(xs))]))
    X=np.array([[1,x] for x in xs])
    Y=np.array([[log(y)] for y in ys])
    theapprox=linalg.solve(X.T.dot(linalg.inv(W).dot(X)),X.T.dot(linalg.inv(W).dot(Y)))
    intercept=theapprox[0,0]
    slope=theapprox[1,0]
    if outputdomain is not None:
        zs=[slope*x+intercept for x in outputdomain]
    else:
        zs=[slope*x+intercept for x in xs]
    return zs, slope, intercept

def getWeightedLogDataPlot(rank,gvgfile):
    """
    Returns fig with plot of data.
    gvgfile is output of VGexperiment.py
    Also plots exponential approximation to the portion of each data series that comes after the perentage had dipped below 50% for the first time.
    """
    gvglengths,gs,gvgs,gvgtrials=extractgvgdata(gvgfile)
    lengths=np.array([i for i in range(1+max(gvglengths))])
    gserrors=np.array([sqrt(gs[i]*(1-gs[i])/float(gvgtrials[i])) for i in range(len(gvgtrials))])
    gvgserrors=np.array([sqrt(gvgs[i]*(1-gvgs[i])/float(gvgtrials[i])) for i in range(len(gvgtrials))])
    approxg,gslope,gintercept=weightedlogapproximatebottomhalf(gvglengths,gs,gvgtrials,lengths)
    approxgvg,gvgslope,gvgintercept=weightedlogapproximatebottomhalf(gvglengths,gvgs,gvgtrials,lengths)
    loggs=log(np.array(gs))
    loggvgs=log(np.array(gvgs))
    glogerrors=[loggs-log(np.array(gs)-gserrors),log(np.array(gs)+gserrors)-loggs]
    gvglogerrors=[loggvgs-log(np.array(gvgs)-gvgserrors),log(np.array(gvgs)+gvgserrors)-loggvgs]
    y0=min(loggs-glogerrors[0])-.1
    y1=0
    fig,ax = plt.subplots()
    plt.xlabel('Word length')
    plt.ylabel('$\log($Proportion$)$')
    plt.title('Random rank '+str(rank)+' words')
    plt.axis([0,.3+lengths[-1],y0,y1])
    ax.errorbar(gvglengths,loggvgs,yerr=gvglogerrors,fmt='o',c='r')
    plt.plot(lengths,approxgvg,color='r',linestyle="dotted",label="$\log($v. geometric $) \sim "+"%.3f"%(gvgslope)+"x+"+"%.3f"%(gvgintercept)+"$")
    ax.errorbar(gvglengths,loggs,yerr=glogerrors,fmt='o',c='b')
    plt.plot(lengths,approxg,c='b',label="$\log($geometic$)$      $\sim "+"%.3f"%(gslope)+"x+"+"%.3f"%(gintercept)+"$")
    ax.legend(loc=3)
    return fig

def weightedlogapproximatebottomhalf(xs,ys,trials,outputdomain=None):
    for halfindex in range(len(ys)):
        if ys[halfindex]<.5:
            break
    else:
        raise RuntimeError("data does not go below .5")
    if outputdomain is not None:
        return weightedlogapproximation(xs[halfindex:],ys[halfindex:],trials[halfindex:],outputdomain)
    else:
        return weightedlogapproximation(xs[halfindex:],ys[halfindex:],trials[halfindex:],xs)



fig=getWeightedDataPlot(rank,gvgfilename,fullwordfilename)
fig.set_size_inches(12,2.9)
fig.savefig("rank"+str(rank)+"gvgfull.pdf",bbox_inches='tight')
fig=getWeightedLogDataPlot(rank,gvgfilename)
fig.set_size_inches(12,2.9)
fig.savefig("rank"+str(rank)+"loggvg.pdf",bbox_inches='tight')
