COMP428/4 Assignment Submission #3

	@author Stephen Young
	@email st_youn@encs.concordia.ca
	@studentid 9736247
	
	
QUESTION 1
----------

For this assignment I implemented Floyd's all-pairs shortest path algorithm (option b)

Partitioning Strategy
---------------------

My partitioning strategy is row-wise, and used the approach described here:
http://www.mcs.anl.gov/~itf/dbpp/text/node35.html
with the addition of all the necessary MPI logic. 

As K is processed, if k is within the block space of a given processor, the processor
will send its row to the others, to ensure that the kth row is always the same while processing.

Caveats
-------

This program is designed to work best for square matrices whose width/height are evenly
divisible by the number of processors implemented. I had some issues with the processors
that, due to downtime on cirrus, i didn't get a chance to look over. Supplied are four such
matrices of different sizes

Execution
---------

Sample input is in /input/, and the program can be run as <mpirun> ./parallel ../input/adjN.mtx
where N is one of 4,16,32,64.

Benchmark
---------

Please see the included benchmark.png file for the efficiency benchmark graph

QUESTION 2
----------

Disadvantages:
	*   For clusters with access to a large set of processors, the initial node that is 
	    accessed may be beneath the node being searched for, and a significant part of the
	  	tree must be accessed before this node may be reached
	*   a poorly balanced tree will yield a lot of wasted processing time
	
Advantages:
	*	In the case of a well balanced tree, processors will be well used and busy
	*	In data sets where the searched for node is deep within the tree, this approach
		allows us to find the result a lot faster than in other topdown search methods
		
QUESTION 3
----------

Proof of modified Diskstra’s token based termination detection algorithm

First we assume that there is a case that the initiator processor can receive a white token
without all processes being terminated.

This implies that the last processor in the ring will pass a white token. Tokens are passed
in two ways, either the token is passed as-is, or it converts to black if the process is black.

Since we know that the token originates as white from the initiator processor, then if at any point
the token encounters a non-idle processor, it will be turned black, and cannot be turned white again.

So either the token will encounter only white processors, meaning the program is complete, or the token
will arrive black. This implies therefor that it is impossible to receive a white token unless
the processes have all successfully terminated, proving this algorithm by negation.



