# Parallel-Distributed-Programming

A parallel implementation of a solution to the Aggressive Queens problem (AQ). This problem is similar to the classic N-Queens problem, where you try to find a possible mapping of N queens on an N × N board, such that no Queen can take another.

Aggressive Queens, are ones which are attacking other Queens. In AQ, you are given two parameters, N > 3 and k ≥ 0, and you are asked to place on an N × N board a maximum number of Queens such that each of them attacks exactly k other Queens. Note that when k = 0, the problem reduces to the standard N-Queens problem. Your program should find ALL maximum solutions.

# Notes

Attack here means direct attack. For example, if two Queens A and B are on the same line with (inbetween) another Queen, then A and B do not attack each other.
There are two other variables of interest, l and w. The program will return and display two types of results, either (if l = 0) returning just the maximum value (in this case 8) as in this line:
```4,3:8:```

or, if l = 1, the locations of the Queens as in this line (note that this might return multiple results): ```4,3:8:1,4,5,6,11,12,13,16```
In addition, your program will operate in two ways, either (if w = 0) using a normal board as in the above example, or if ```w = 1```, using a wrap-around board, where (for example) cells 8 and 5 are adjacent, and cells 15 and 3.
The program should be called with four parameters, like this:

    ``` ./findAQ N k l w```
    
• N is the size of a side of the board
• k is the numer of pieces that every queen must attack
• ```l=0``` is short format N,k:number: and l=1 is long format N,k:number:pc1,pc2,pc3...
• ```w=0``` is no wraparound, ```w=1``` is wraparound

So:
should return
     ```4,3:8:1,4,5,6,11,12,13,16```
as in the example above.
