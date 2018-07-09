# Personalized-Page-Rank
Serial and MPI parallel PPR.

PageRank is a famous algorithm from Google to rank the website page we query about. Personalized PageRank adds personal preference to pages thus it's useful in recommendation and discovery. The concept is in below figure. In this repository, I implement the serial and parallel version of PPR with MPI due to the perfect parallel feature of this problem.
<div align=center><img src="https://github.com/GitBoSun/Personalized-Page-Rank/blob/master/images/concept.png" width="50%">

## Algorithm

We regard page relationship as a directed graph in which pages are nodes and an edge appears when a page points to another. We take random walk to travel on the graph and finally, we'll get a list of number that records the probability of going there which can be considered to be importance or rank of that node. Moreover, we allow random jumping to other nodes to avoid self-locking. We follow two rules to evaluate the relative importance(i.e. rank) of pages:
  
1. Pages that are pointed to by many pages are important.
2. pages that are pointed to by important pages are also important. 

Note that the importance that one node pointed out is divided by its all outing nodes. 

As for Personalized PageRank, we allow jumping to pages that a person prefers such as the origin nodes(source). 
Detailing algorithm is in [Algorithm.pdf](https://github.com/GitBoSun/Personalized-Page-Rank/blob/master/Algorithm.pdf).

## Results

We use a toy graph to show our PPR results. The origin graph is in Fig.2. Fig.3(a) and Fig.3(b) is results of ragarding B and F respectively. Here source means we start from that node and we see it as person's prefering node. We can find that the source will culmulate importance in the whole random walk process. 

<div align=center><img src="https://github.com/GitBoSun/Personalized-Page-Rank/blob/master/images/graph.png" width="70%">
                                                
<div align=center><img src="https://github.com/GitBoSun/Personalized-Page-Rank/blob/master/images/results.png" width="70%">

