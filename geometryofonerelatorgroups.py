import networkx as nx
import automaticgroups.automaticgroups as ag
import smallcancellation.smallcancellation as sc
import freegroups.freegroup as fg
import freegroups.whiteheadgraph as wg
import pexpect
from fractions import Fraction
import itertools
import os
from collections import deque
#import freegroups.stallings_folding as sf
import freegroups.enumeratefreegroupwords as enum
import glob
import subprocess32


# One of the checks that certify_hyperbolicity does is to use the walrus package in GAP. The GAP startup takes a while. If multiple checks are to be run it is faster to spawn GAP only once and reuse it as follows:
# words = list of words in free group of rank r, each defining a 1=relator quotient
# >>> gap=spawngapforwalrus(r)
# >>> for word in words:
#             certify_hyperbolicity(word, gap=gap)
def certify_hyperbolicity(relator,reportreason=False,**kwargs):
    """
    Check if one-relator group defined by input relator is hyperbolic.
    Returns 1 if group is known hyperbolic, -1 if known not hyperbolic, 0 if inconclusive.

    If reportreason=True output additionally contains string with name of the method certifying hyperbolicity, or None is all tests are inconclusive.
    tryhard=2 will make kbmag retry several times with random orderings of the generators if the first attempt fails.


    >>> certify_hyperbolicity([1,2,3],reportreason=True)
    (True, 'free')
    >>> certify_hyperbolicity('bbb',reportreason=True)
    (True, 'torsion')
    >>> certify_hyperbolicity([1,2,-1,-2,3,4,-3,-4],reportreason=True)
    (True, 'cyclically pinched')
    >>> certify_hyperbolicity('abcabcdeDEdeDE',reportreason=True)
    (False, 'cyclically pinched')
    >>> certify_hyperbolicity('CCBBCAAbbcaa',reportreason=True)
    (True, 'walrus')
    """
    def format_return(result,reason):
        if reportreason:
            return result,reason
        else:
            return result
    if not relator:
        return format_return(True,'free')
    F,r1=fg.parseinputword(relator)
    if 'is_minimal' in kwargs and kwargs["is_minimal"]: # input relator guaranteed to be Whitehead minimal
        r2=r1
    else: # find a Whitehead minimal element in the same automorphism orbit of the relator
        r2=wg.whitehead_minimal_representative(r1)
    if len(r2)<=1:
        return format_return(True,'free')
    if F.degree(r2)>1:# one-relator groups with torsion are hyperbolic
        return format_return(True,'torsion')
    # Check if relarot is cyclically pinched
    (m,n)=is_cyclically_pinched(r2,reportpowers=True)
    if m is not None:
        if m>1 and n>1:
            return format_return(False,'cyclically pinched')
        else:
            return format_return(True, 'cyclically pinched')
    # Try Ivanov-Scupp criteria, always gives an answer if some generator appears at most 3 times
    if 'verbose' in kwargs and kwargs["verbose"]:
        print "Checking Ivanov-Schupp criteria."
    ISresult=IvanovSchupp(r2,reportreason=False)
    if ISresult is not None:
        return format_return(ISresult,'Ivanov Schupp')
    if 'verbose' in kwargs and kwargs["verbose"]:
        print "Ivanonv-Schupp does not apply."
    # check if relator defines small cancellation presentation
    Cprime=sc.Cprime([r2])
    if sc.smallcancellation([r2],Cprime):
        return format_return(True,'small cancellation')
    # check Blufstein-Minian condition
    if BlufsteinMinian(r2,Cprime):
        return format_return(True,'Blufstein Minian')
    if 'verbose' in kwargs and kwargs["verbose"]:
        print "Not small cancellation."
    # try GAP with walrus
    if 'walrus' in kwargs and kwargs['walrus']==False:
        pass
    else:
        if 'verbose' in kwargs and kwargs['verbose']:
            print "Trying GAP with walrus."
        if 'gapprompt' in kwargs:
            gapprompt=kwargs["gapprompt"]
        else:
            gapprompt='gap>'
        if 'gapfreegroupname' in kwargs:
            gapfreegroupname=kwargs["gapfreegroupname"]
        else:
            gapfreegroupname='f'
        if 'walrusparameter' in kwargs:
            walrusparameter=kwargs["walrusparamter"]
        else:
            walrusparameter='1/100'
        if 'gap' in kwargs:
            walrus=checkhyperbolicitywithwalrus(r2.letters,gap=kwargs["gap"],gapfreegroupname=gapfreegroupname,gapprompt=gapprompt,theparameter=walrusparameter)
        else:
            walrus=checkhyperbolicitywithwalrus(r2.letters,gap=None,gapfreegroupname=gapfreegroupname,gapprompt=gapprompt,theparameter=walrusparameter)
        if walrus:
            return format_return(True,'walrus')
        if 'verbose' in kwargs and kwargs["verbose"]:
            print "walrus failed."
    # try kbmag
    if 'kb' in kwargs and kwargs['kb']==False:
        return format_return(None,None)
    if 'verbose' in kwargs and kwargs["verbose"]:
        print "Trying kbmag"
    autandhyp=ag.certify_hyperbolicity(r2(),**kwargs)
    if autandhyp:
        return format_return(True,'kbmag')
    else: # all tests inconclusive
        return format_return(None, None)


def is_cyclically_pinched(relator,reportwords=False,reportpowers=False):
    F,rels=fg.parseinputwords([relator])
    r=rels[0].letters
    rr=r+r
    for (startingindex,L) in itertools.product(range(len(r)),range(1,len(r))):
        prefix=rr[startingindex:startingindex+L]
        prefixletters={abs(x) for x in prefix}
        suffix=rr[startingindex+L:startingindex+len(r)]
        suffixletters={abs(x) for x in suffix}
        if not prefixletters&suffixletters:
            result=True
            words=(F.word(prefix),F.word(suffix))
            if reportpowers:
                powers=(F.degree(words[0]),F.degree(words[1]))
            break
    else:
        result=False
        words=(None,None)
        powers=(None,None)
    if reportwords:
        return words
    if reportpowers:
        return powers
    return result
            
            

    


def BlufsteinMinian(relator,Cprime=None):
    """
    Decide if the one-relator group defined by the given relator satisfies the hyperbolicity condition of Blufstein Minian that it is Cprime(1/4) and satisfies a condition they call Tprime.
    """
    if Cprime is None:
        Cprime=sc.Cprime([relator])
    return Cprime<Fraction(1,4) and BlufsteinMinianTprime(relator)

def BlufsteinMinianTprime(relator):
    """
    Decide if the one-relator group defined with given relator satisfies the Blufstein Minian T-prime condition.

    >>> BlufsteinMinianTprime('aaaababaBAbbbAb')
    True
    >>> BlufsteinMinianTprime('DCabFEcdBAef')
    False
    """
    F,rels=fg.parseinputwords([relator])
    R=rels[0].letters
    Rinv=(rels[0]**(-1)).letters
    def longestcommonprefix(list1,list2):
        matchingprefixlength=0
        for i in range(min(len(list1),len(list2))):
            if list1[i]!=list2[i]:
                return matchingprefixlength
            else:
                matchingprefixlength+=1
        return matchingprefixlength
    def overlap(a,b):
        if a[1]==1 and b[1]==1:
            return longestcommonprefix((Rinv+Rinv)[len(R)-a[0]-1:],(R+R)[b[0]+1:])
        elif a[1]==-1 and b[1]==1:
            if len(R)-a[0]==b[0]:
                return float('inf')
            else:
                return longestcommonprefix((R+R)[len(R)-a[0]-1:],(R+R)[b[0]+1:])
        elif a[1]==1 and b[1]==-1:
            if len(R)-a[0]==b[0]:
                return float('inf')
            else:
                return longestcommonprefix((Rinv+Rinv)[len(R)-a[0]-1:],(Rinv+Rinv)[b[0]+1:])
        elif a[1]==-1 and b[1]==-1:
            return longestcommonprefix((R+R)[len(R)-a[0]-1:],(Rinv+Rinv)[b[0]+1:])
        else:
            raise InputError
    G=nx.Graph(wg.WGraph(rels))
    A=nx.adjacency_matrix(G)
    B=A**3
    nodelist=list(G.nodes())
    threecycles={(nodelist[i],nodelist[j],nodelist[k]) for (i,j,k) in itertools.combinations([h for h in range(len(nodelist)) if B[h,h]!=0],3) if A[i,j] and A[j,k] and A[k,i]}
    for (i,j,k) in threecycles:
        possiblefirstindex=[(h,1) for h in range(len(R)) if (R+R)[h:h+2]==[-j,i]]+[(h,-1) for h in range(len(R)) if (Rinv+Rinv)[h:h+2]==[-j,i]]
        possiblesecondindex=[(h,1) for h in range(len(R)) if (R+R)[h:h+2]==[-k,j]]+[(h,-1) for h in range(len(R)) if (Rinv+Rinv)[h:h+2]==[-k,j]]
        possiblethirdindex=[(h,1) for h in range(len(R)) if (R+R)[h:h+2]==[-i,k]]+[(h,-1) for h in range(len(R)) if (Rinv+Rinv)[h:h+2]==[-i,k]]
        for (f,s,t) in itertools.product(possiblefirstindex,possiblesecondindex,possiblethirdindex):
            tripodlength=overlap(f,s)+overlap(s,t)+overlap(t,f)
            if tripodlength<float('inf') and 2*tripodlength>=len(R):
                return False
    return True
                   
    
def IvanovSchupp(relator,reportreason=False):
    """
    Check if one-relator group defined by input relator is hyperbolic according to Ivanov and Schupp criteria.
    Returns True if group is known hyperbolic, False if known not hyperbolic, None if inconclusive

    >>> IvanovSchupp('',reportreason=True)
    (True, 'free')
    >>> IvanovSchupp('a',reportreason=True)
    (True, 'free')
    >>> IvanovSchupp('aa',reportreason=True)
    (True, 'torsion')
    >>> IvanovSchupp('abaCBC',reportreason=True)
    (False, 'Thm3(1)')
    >>> IvanovSchupp('ababcbc',reportreason=True)
    (True, 'Thm3(1)')
    >>> IvanovSchupp('abbAbcBC',reportreason=True)
    (True, 'Thm3(2)')
    >>> IvanovSchupp('bcBCCBAcbCBcbCBabc',reportreason=True)
    (False, 'Thm3(2a)')
    >>> IvanovSchupp('acaacacaacacaacbAAB',reportreason=True)
    (False, 'Thm3(2b)')
    >>> IvanovSchupp('ababbacb',reportreason=True)
    (True, 'Thm3(3)')
    >>> IvanovSchupp('ababaccBCbccBCbb',reportreason=True)
    (False, 'Thm3(3a)')
    >>> IvanovSchupp('abaBcbCCBcbCCbaccBCbccBCbb',reportreason=True)
    (False, 'Thm3(3b)')
    >>> IvanovSchupp('abaccBCbbaBcbCCb',reportreason=True)
    (False, 'Thm3(3c)')
    >>> IvanovSchupp('abaBcbCCbaBcbCCBcbCCb',reportreason=True)
    (False, 'Thm3(3d)')
    >>> IvanovSchupp('ababbcBCbcBCbcBCBAcbCB',reportreason=True)
    (True, 'Thm3(4)')
    >>> IvanovSchupp('ababbcBCBAcbCB',reportreason=True)
    (False, 'Thm3(4a)')
    >>> IvanovSchupp('ababbcBCbcBCBAcbCB',reportreason=True)
    (False, 'Thm3(4b)')
    >>> IvanovSchupp('ababcabccabcbcbcbcabcbcbcbcbc',reportreason=True)
    (True, 'Thm4')
    >>> IvanovSchupp('abacabcabCbc',reportreason=True)
    (False, 'Thm4(3)')
    >>> IvanovSchupp('ababcBCababccc',reportreason=True)
    (None, None)
    >>> IvanovSchupp('aaaaabbbbbccccc',reportreason=True)
    (None, None)
    """
    def format_return(result,reason):
        if reportreason:
            return result,reason
        else:
            return result
    if not relator:
        return format_return(True,'free')
    F,r1=fg.parseinputword(relator)
    r2=F.cyclic_reduce(r1)
    if F.degree(r2)>1:# one-relator groups with torsion are hyperbolic
        return format_return(True,'torsion')
    # Try Ivanov-Schupp criteria
    for a in range(1,1+F.rank):
        acount=[abs(x) for x in r2.letters].count(a)
        if acount==0: # neither theorem applies, check next a
            continue
        elif acount==1: # relator is primitive
            return format_return(True, 'free')
        # Ivanov-Schupp Theorem 3 applies when acount=2 or 3
        elif acount==2: 
            if r2.letters.count(a)==0:
                therel=(r2**(-1)).letters
            else:
                therel=r2.letters
            firsta=therel.index(a)
            therel=therel[firsta:]+therel[:firsta]
            if therel[1:].count(a)>0: # Ivanov-Schupp Thm 3 case (1)
                nexta=1+therel[1:].index(a)
                B=F.word(therel[1:nexta])
                if nexta<len(therel)-1:
                    C=F.word(therel[nexta+1:])
                else:
                    C=F.word([])
                if F.degree(B*C**(-1))>1:
                    return format_return(False,'Thm3(1)')
                else:
                    return format_return(True,'Thm3(1)')
            else: # Ivanov-Schupp Thm 3 case (2)
                nexta=1+therel[1:].index(-a)
                B=F.word(therel[1:nexta])
                C=F.word(therel[nexta+1:]) # no index error here since word is cyclically reduced
                if F.degree(B)>1 and F.degree(C)>1:
                    return format_return(False,'Thm3(2b)')
                elif F.is_conjugate_into(B,C) or F.is_conjugate_into(C,B):# group is not hyperbolic by case (2a), simplified because by previous case we know either B or C is not a proper power
                    return format_return(False,'Thm3(2a)')
                else:
                    return format_return(True,'Thm3(2)')
        elif acount==3:
            if r2.letters.count(-a)>r2.letters.count(a):
                therel=(r2**(-1)).letters
            else:
                therel=r2.letters
            if therel.count(-a)>0: # Theorem 3 case (4)
                firsta=therel.index(a)
                seconda=1+firsta+therel[1+firsta:].index(a)
                nega=therel.index(-a)
                if nega>seconda or nega<firsta:
                    therel=therel[firsta:]+therel[:firsta]
                else:
                    therel=therel[seconda:]+therel[:seconda]
                firsta=therel.index(a)
                seconda=1+firsta+therel[1+firsta:].index(a)
                nega=therel.index(-a)
                B=F.word(therel[1+firsta:seconda])
                C=F.word(therel[1+seconda:nega])
                D=F.word(therel[1+nega:])
                Z1,n1=F.max_root(B**(-1)*C*B)
                Z2,n2=F.max_root(D)
                if Z2==Z1**(-1):
                    Z2=Z1
                    n2=-n2
                if Z1==Z2 or Z1==F.word([]) or Z2==F.word([]):
                    if abs(n1)==abs(n2):
                        return format_return(False,'Thm3(4a)')
                    elif n1==-2*n2 or n2==-2*n1:
                        return format_return(False,'Thm3(4b)')
                    else:
                        return format_return(True,'Thm3(4)')
                else:
                    return format_return(True,'Thm3(4)')
            else: # Theorem 3 case (3)
                therel=therel[therel.index(a):]+therel[:therel.index(a)]
                firsta=0
                seconda=1+therel[1:].index(a)
                thirda=1+seconda+therel[1+seconda:].index(a)
                B=F.word(therel[1+firsta:seconda])
                C=F.word(therel[1+seconda:thirda])
                try:
                    D=F.word(therel[1+thirda:])
                except IndexError:
                    D=F.word([])
                Z1,n1=F.max_root(C*B**(-1))
                Z2,n2=F.max_root(D*B**(-1))
                if Z2==Z1**(-1):
                    Z2=Z1
                    n2=-n2
                if Z1==Z2 or Z1==F.word([]) or Z2==F.word([]):
                    if (n1==0 and abs(n2)>1) or (n2==0 and abs(n1)>1):
                        return format_return(False,'Thm3(3a)')
                    elif abs(n1)==abs(n2) and abs(n1)>1:
                        return format_return(False,'Thm3(3b)')
                    elif n1!=0 and n1==-n2:
                        return format_return(False,'Thm3(3c)')
                    elif n1!=0 and (n1==2*n2 or n2==2*n1):
                        return format_return(False,'Thm3(3d)')
                    else:
                        return format_return(True,'Thm3(3)')
                else:
                    return format_return(True,'Thm3(3)')
        else: # Theorem 3 does not apply. Theorem 4 applies to some words with higher acount. Need all a's same sign plus some other conditions.
            therel=r2.letters
            if therel.count(a)==0: # there are no a's, but we know there are some a's or a^-1's, so the inverse must have some a's
                therel=(r2**(-1)).letters
            if therel.count(a)>0 and therel.count(-a)==0: # only positive a's appear. Now check words between a's are all distinct.
                firstaindex=therel.index(a)
                R=F.word(therel[firstaindex:]+therel[:firstaindex]) #cylically permute r so that a is first letter
                T=[]
                currentword=[]
                for i in range(1,len(R)):
                    if R.letters[i]==a:
                        T.append(F.word(currentword))
                        currentword=[]
                    else:
                        currentword.append(R.letters[i])
                    if i==len(R)-1:
                        T.append(F.word(currentword))
                if len(T)==len(set(T)): #intermediate words are distinct, so theorem 4 applies
                    if R.letters.count(a)>4:
                        return format_return(True,'Thm4')
                    else: # check Ivanov-Schupp Thm 4 case (3)
                        assert(R.letters.count(a)==4) # if it were fewer we should have caught this in the theorem 3 part above
                        for i in range(4):
                            if T[i%4]*T[(i+1)%4]**(-1)*T[(i+2)%4]*T[(i+3)%4]**(-1)==F.word([]): 
                                return format_return(False,'Thm4(3)')
                        else: 
                            return format_return(True,'Thm4')
            else: # Theorem 4 does not apply with this a.
                pass
    # Ivanov-Schupp does not apply
    return format_return(None,None)


   


def checkhyperbolicitywithwalrus(theword,theparameter='1/100',gap=None,gapfreegroupname='f',gapprompt='gap>',fulloutput=False):
    """
    Check hyperbolicity of the one relator group with relator defined by theword using the walrus package in GAP.

    Input 'gap' should be a pexpect process of GAP with walrus loaded, a free group defined, ready for input. If None, process will be spawned.
    """ 
    if gap is None:
        gap=spawngapforwalrus(max(abs(x) for x in theword),gapfreegroupname=gapfreegroupname,gapprompt=gapprompt)
    gap.send('IsHyperbolic(PregroupPresentationFromFp('+gapfreegroupname+',[],['+converttogapword(theword,gapfreegroupname)+']),'+theparameter+');\r')
    gap.expect(gapprompt)
    output=gap.before
    if 'true' in output:
        result=True
    elif 'fail' in output:
        result=False
    else:
        raise ValueError(output)
    if fulloutput:
        return result,output
    else:
        return result

def spawngapforwalrus(rank,gapfreegroupname='f',gapprompt='gap>'):
    """
    Sapwns a GAP process, loads package walrus, and defines free group of specified rank.
    """
    gap=pexpect.spawn('gap -b')
    gap.expect(gapprompt)
    gap.send('LoadPackage("walrus");;\r')
    gap.expect(gapprompt)
    gap.send(gapfreegroupname+':=FreeGroup('+str(rank)+');;\r')
    gap.expect(gapprompt)
    return gap

def converttogapword(listsofints,gapfreegroupname='f'):
    """
    Convert list of nonzero intergers representing a word in a free group to a string readable by GAP.

    >>> converttogapword([1,2,-1,-2,])
    'f.1*f.2*f.1^-1*f.2^-1'
    """
    outputstring=list()
    for i in range(len(listsofints)):
        x=listsofints[i]
        if x>0:
            outputstring.append(gapfreegroupname+'.'+str(x))
        else:
            outputstring.append(gapfreegroupname+'.'+str(abs(x))+"^-1")
        if i<len(listsofints)-1:
            outputstring.append("*")
    return ''.join(outputstring)



#----------------------------------

def girth(relator,verbose=False,**kwargs):
    """
    Given a word defining a shortlex automatic one-relator group, return the girth of the Cayley graph for the corresponding one-relator presentation and a pair of equivalent words witnessing shortest relation.

    Default generator ordering is [...,'B','A','a','b',...]. Use keyword 'generators' to specify different ordering of generators and inverses.
    Use filename=kbmagbasefilename to specify existing kbmag files defining the automatric structure.
    Use cleanup=True to delete all auxilliary kbmag files or cleanup=False to keep them. Default is to cleanup if this function creates the files and not if they already exist.
    Set timeout=n to set time limit of n seconds on attempt to generate automatric group structure.
    """
    if not relator:
        return float('inf')
    F,r=fg.parseinputword(relator)
    relatorasstring=r()
    relatoraslist=r.letters
    if 'filename' in kwargs:
        thefilename=kwargs['filename']
    else:
        thefilename="OneRelatorGroup-"+relatorasstring
    if 'cleanup' in kwargs:
        cleanup=kwargs['cleanup']
    else:
        cleanup=not os.path.isfile(thefilename) # default is to cleanup if the files are not already existing
    if not os.path.isfile(thefilename):
        if 'generators' in kwargs:
            generators=kwargs['generators']
        else:
            generators=[x.upper() for x in F.gens[::-1]]+F.gens
        ag.writetokbmagfile(thefilename,generators,[relatorasstring])
    if 'timeout' in kwargs:
        timeout=kwargs['timeout']
    else:
        timeout=10 # default timeout 10s       
    if not os.path.isfile(thefilename+'.diff1'):
        try:
            subprocess32.run(['autgroup','-silent',thefilename],check=True,timeout=timeout)
        except subprocess32.TimeoutExpired, subprocess32.CalledProcessError:
            if cleanup:
                files = glob.glob('./'+thefilename+"*")
                for file in files:
                    os.remove(file)
            raise
    foundrelation=False
    currentlength=0
    while not foundrelation:
        currentlength+=1
        if 2*currentlength==len(relatorasstring) or 2*currentlength-1==len(relatorasstring):
            return len(relatorasstring),relatorasstring,''
        if verbose:
            print "Searching words of length "+str(currentlength)
        wordsofcurrentlength=enum.generate_words(F.rank,currentlength,currentlength)
        for theintlist in wordsofcurrentlength:
            thestring=F.word(theintlist)()
            reducedstring=ag.groupelement(thestring,thefilename).string
            if len(reducedstring)<len(thestring):
                assert(len(reducedstring)+1==len(thestring))
                if 'cleanup' in kwargs and kwargs['cleanup']==True:
                    files = glob.glob('./'+thefilename+"*")
                    for file in files:
                        try:
                            os.remove(file)
                        except:
                            pass
                return 2*currentlength-1,thestring,reducedstring
            if thestring!=reducedstring: # found a relation of lengh currentlength*2, but keep searching for one of length 2*currentlength-1
                foundrelation=True
                u=thestring
                v=reducedstring
    if 'cleanup' in kwargs and kwargs['cleanup']==True:
        files = glob.glob('./'+thefilename+"*")
        for file in files:
            try:
                os.remove(file)
            except:
                pass
    return currentlength*2,u,v
    

def ball(F,relator,thefilename,verbose=False,**kwargs):
    NotImplemented
    try:
        d=deque(relator.letters)
    except AttributeError:
        d=deque(relator)
    cyclicpermutations=[]
    for i in range(len(d)):
        cyclicpermutations.append(list(d))
        d.rotate()
    B=nx.MultiDiGraph()
    B.add_node((0,0))
    for i in range(len(relator)):
        if cyclicpermutations[i][0]>0:
            B.add_edge((0,0),(i,1),label=cyclicpermutations[i][0])
        else:
            B.add_edge((i,1),(0,0),label=-cyclicpermutations[i][0])
        if cyclicpermutations[i][-1]>0:
            B.add_edge((i,len(relator)-1),(0,0),label=cyclicpermutations[i][-1])
        else:
            B.add_edge((0,0),(i,len(relator)-1),label=-cyclicpermutations[i][-1])
        for j in range(1,len(relator)-1):
            if cyclicpermutations[i][j]>0:
                B.add_edge((i,j),(i,j+1),label=cyclicpermutations[i][j])
            else:
                B.add_edge((i,j+1),(i,j),label=-cyclicpermutations[i][j])
    sf.fold(B,(0,0))
    # now use automatric structure to check that all vertices are really distinct
    T=nx.maximum_spanning_tree(nx.MultiGraph(B))
    for n in B:
        theword=F.word(edgeword(B,T,(0,0),n))
        minword=ag.groupelement(theword(),thefilename).string
        B.nodes[n]['representative']=minword
    del T
    originalnodelist=[n for n in B.nodes() if n!=(0,0)]
    for n in originalnodelist:
        if n not in B:
            continue
        for m in [x for x in B if x!=n]:
            if B.nodes[n]['representative']==B.nodes[m]['representative']:
                sf.wedge(B,n,m)
                break
    sf.fold(B)
    return B

def edgeword(B,T,v,w):
    path=nx.shortest_path(T,v,w)
    theword=[]
    for i in range(len(path)-1):
        for e in T.edges(path[i],keys=True,data=True):
            if e[1]==path[i+1]:
                break
        else:
            raise RuntimeError
        if e in B.edges(data=True,keys=True):
            theword.append(e[3]['label'])
        elif (e[1],e[0],e[2],e[3]) in B.edges(data=True,keys=True):
            theword.append(-e[3]['label'])
        else:
            raise RuntimeError
    return theword
            













# in terminal, do
# python geometryofonerelatorgroups.py
# to run doctests
if __name__ == "__main__":
    import doctest
    doctest.testmod()
