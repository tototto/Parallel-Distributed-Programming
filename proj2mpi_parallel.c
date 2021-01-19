#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#define LEFT 1
#define RIGHT 2
#define TOP 3
#define BOTTOM 4

#define LEFT_TOP_DIAG 5
#define LEFT_BTM_DIAG 6
#define RIGHT_TOP_DIAG 7
#define RIGHT_BTM_DIAG 8

#define TOP_LEFT_CORNER 9
#define BTM_LEFT_CORNER 10
#define TOP_RIGHT_CORNER 11
#define BTM_RIGHT_CORNER 12

#define NON_CORNER_CASE_TOP 13
#define NON_CORNER_CASE_BTM 14
#define NON_CORNER_CASE_LEFT 15
#define NON_CORNER_CASE_RIGHT 16

#define NO_CORNER_CASE 0

#define OPP_IS_SAFE 0

typedef struct Points {
   int x;
   int y;
} Point ;

/** 
 * Determines the current time
 *
 **/
long long wall_clock_time()
{
#ifdef LINUX
    struct timespec tp;
    clock_gettime(CLOCK_REALTIME, &tp);
    return (long long)(tp.tv_nsec + (long long)tp.tv_sec * 1000000000ll);
#else
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (long long)(tv.tv_usec * 1000 + (long long)tv.tv_sec * 1000000000ll);
#endif
}


int* decToBinary(int n, int size, int binaryNum[size*size])
{

    // counter for binary array
    int i = 0;
    while (n > 0) {

        // storing remainder in binary array
        binaryNum[i] = n % 2;
        n = n / 2;
        i++;
    }

    // return binary array
    return binaryNum;
}

// Function to print the output
void printQueenLocation(int n, int **arr)
{
    int pos = 1;
    int row = 0;
    for ( row = 0; row < n; row++) {
        int col = 0; 
        for( col = 0; col < n; col++, pos++)
            if(arr[row][col] == 1)
                printf("%d,", pos);
    }
    printf("\n");
}

int numOfQueen(int n, int **arr) {
    int count = 0;
    int i = 0;
    for(i = 0; i < n; i++){
        int j = 0;
        for(j = 0; j < n; j++) {
            if(arr[i][j]==1)
                count++;
        }
    }
    return count;
}

// Function to generate all binary strings
void generateAllBoardPermutation(int n, int result[(int) pow(2, n*n)][n*n])
{
    int totalPermutationSize = pow(2, n*n);
    int currNum = 0;
    for(currNum = 0; currNum < totalPermutationSize; currNum++) {

        int oneBoard[n*n];
        memset(oneBoard, 0, (n*n)*sizeof(int));

        int *eachFilledUpBoard = decToBinary(currNum, n, oneBoard);
        int eachDigit = 0;

        for( eachDigit = 0; eachDigit < n*n; eachDigit++) {
            result[currNum][eachDigit] = eachFilledUpBoard[eachDigit];
        }

    }
}

// convert each Board to a 2D Array
int** convertTo2DArray(int size, int arr[size]) {

    int newArraySize = sqrt(size);;

    // Create 2D array using Malloc
    int **convertedArray = malloc(sizeof(int*) * newArraySize);
    int i = 0;

    for(i = 0; i < newArraySize; i++) {
        convertedArray[i] = malloc(sizeof(int*) * newArraySize);
    }

    // Begin 1D to 2D conversion
    int element = 0;
    int row = 0;
    for( row = 0; row < newArraySize; row++) {
        int col = 0;
        for( col = 0; col < newArraySize; col++, element++) {
            convertedArray[row][col] = arr[element];
        }
    }

    return convertedArray;
}

// Function that returns true if the queen
// can attack the opponent
bool canQueenAttack(int qR, int qC, int oR, int oC)
{
    // If queen and the opponent are in the same row
    if (qR == oR)
        return true;

    // If queen and the opponent are in the same column
    if (qC == oC)
        return true;

    // If queen can attack diagonally
    if (abs(qR - oR) == abs(qC - oC))
        return true;

    // Opponent is safe
    return false;
}

// Get the Point coordinate for all the Queen on ChessBoard //layout[size][size]
Point* getAllQueenPosition(int size, int** layout) {
    Point *result = malloc(sizeof(Point) * (size*size));
    int resultArrIndx = 0;
    int row = 0;
    for( row = 0; row < size; row++) {
        int col = 0;
        for( col = 0; col < size; col++) {
            if(layout[row][col] == 1) {
                result[resultArrIndx].x = row;
                result[resultArrIndx].y = col;
                resultArrIndx++;
            }
        }
    }
    return result;
}

// Check if all array elements are equal
int checkAllEqual(const int a[], int n)
{
    int firstElement = a[0];
    int i = 0;
    for(i = 1; i < n; i++) {
        if(a[i] != firstElement)
            return false;
    }

    return true;
}

// Find the direction the Queen attacks the Opp in
int dirQueenAttack(int qR, int qC, int oR, int oC)
{
    // If queen and the opponent are in the same row
    if (qR == oR) {
        if(qC < oC) {
            // Queen will atk Opp on Queen's Right
            return RIGHT;
        }
        if(qC > oC) {
            // Queen will atk Opp on Queen's Left
            return LEFT;
        }
    }

    // If queen and the opponent are in the same column
    if (qC == oC) {
        if(qR < oR) {
            // Queen will atk Opp on Queen's Bottom
            return BOTTOM;
        }
        if(qR > oR) {
            // Queen will atk Opp on Queen's Top
            return TOP;
        }
    }

    // If queen can attack diagonally
    if (abs(qR - oR) == abs(qC - oC)) {

        if(qR > oR && qC < oC) {
            return RIGHT_TOP_DIAG;
        }
        if(qR < oR && qC < oC) {
            return RIGHT_BTM_DIAG;
        }
        if(qR > oR && qC > oC) {
            return LEFT_TOP_DIAG;
        }
        if(qR < oR && qC > oC) {
            return LEFT_BTM_DIAG;
        }
    }

    // Opponent is safe
    return OPP_IS_SAFE;
}


// Takes in an array of Queen Point determine if the 2 POINTS supplied is blocked
// Returns true if blocked
bool checkIfBlock(int qR, int qC, int oR, int oC, int size, Point* queenPos) {

    // Check direction of Queen Attk
    int attkDir = dirQueenAttack(qR, qC, oR, oC);

    // Check if exists any queen that block the Dir of Attk
    switch(attkDir) {
        int other = 0;
        case TOP:  // If there is another queen in the upper rows && same column
                    for( other = 0; other < size; other++)
                        if(queenPos[other].x < qR && qC == queenPos[other].y && queenPos[other].x > oR && queenPos[other].y == oC)
                            return true;
                    break;

        case BOTTOM: // If there is another queen in the lower rows && same column 
                    for( other = 0; other < size; other++)
                        if(queenPos[other].x > qR && qC == queenPos[other].y && queenPos[other].x < oR && queenPos[other].y == oC)
                            return true;
                    break;

        case LEFT: // If there is another queen in the same rows && left column
                    for( other = 0; other < size; other++)
                        if(queenPos[other].x == qR && queenPos[other].y < qC && queenPos[other].x == oR && queenPos[other].y > oC)
                            return true;
                    break;

        case RIGHT: // If there is another queen in the same rows && right column
                    for( other = 0; other < size; other++)
                        if(queenPos[other].x == qR && queenPos[other].y > qC && queenPos[other].x == oR && queenPos[other].y < oC)
                            return true;
                    break;

        case RIGHT_BTM_DIAG: // If there is another queen in the left diagonal btm before Opp
                              while(qR != oR && qC != oC) {
                                 ++qR;
                                 ++qC;

                                 if(qR == oR && qC == oC)
                                    break;

                                 for( other = 0; other < size; other++)
                                     if(queenPos[other].x == qR &&  queenPos[other].y == qC)
                                        return true;
                              }
                              break;

        case RIGHT_TOP_DIAG: // If there is another queen in the right diagonal top before Opp
                              while(qR != oR && qC != oC) {
                                --qR;
                                ++qC;

                                if(qR == oR && qC == oC)
                                    break;


                                 for( other = 0; other < size; other++)
                                     if(queenPos[other].x == qR &&  queenPos[other].y == qC)
                                        return true;
                              }
                              break;

        case LEFT_TOP_DIAG: // If there is another queen in the left diagonal top before Opp
                            while(qR != oR && qC != oC) {
                                --qR;
                                --qC;

                                 if(qR == oR && qC == oC)
                                    break;

                                 for( other = 0; other < size; other++)
                                     if(queenPos[other].x == qR &&  queenPos[other].y == qC)
                                        return true;
                             }
                             break;

        case LEFT_BTM_DIAG: // If there is another queen in the left diagonal top before Opp
                            while(qR != oR && qC != oC) {
                                ++qR;
                                --qC;

                                if(qR == oR && qC == oC)
                                    break;

                                 for( other = 0; other < size; other++)
                                     if(queenPos[other].x == qR &&  queenPos[other].y == qC)
                                        return true;
                             }
                             break;
    }

    return false;
}

// Check if Queen is at corner
// return the corner it is at, otherwise return border its at
int queenAtWhichCorner(int qR, int qC, int N) {

    if(qR == 0 && qC == 0)
        return TOP_LEFT_CORNER;

    if(qR == 0 && qC == N-1)
        return TOP_RIGHT_CORNER;

    if(qR == N-1 && qC == 0)
        return BTM_LEFT_CORNER;

    if(qR == N-1 && qC == N-1)
        return BTM_RIGHT_CORNER;

    if(qR == 0)
        return NON_CORNER_CASE_TOP;

    if(qR == N-1)
        return NON_CORNER_CASE_BTM;

    if(qC == 0)
        return NON_CORNER_CASE_LEFT;

    if(qC == 0)
        return NON_CORNER_CASE_RIGHT;

    return NO_CORNER_CASE; // if Not at any corner or border
}

// Check if Queen is at ChessBoard border
// return TRUE if so, FALSE otherwise
bool isQueenAtBorder(int qR, int qC, int N) {

    // Check if Queen @ Top Row
    if(qR == 0)
        return true;

    // Check if Queen @ Bottom Row
    if(qR == N-1)
        return true;

    // Check if Queen @ LEFT MOST column
    if(qC == 0)
        return true;

    // Check if Queen @ RIGHT MOST column
    if(qC == N-1)
        return true;

    return false;
}

// Calculate the num of Other Queen that Curr Queen can Attk in Wrap Around
int wrapArndAttkCount(int qR, int qC, int N, Point* posOfQueenList, int queenPos) {

    int queenListSize = sizeof(posOfQueenList);
    int attkCount = 0;
    int eachQueen = 0;

    switch(queenPos) {

        case NON_CORNER_CASE_TOP:
                        for( eachQueen = 0; eachQueen < queenListSize; eachQueen++) {
                            if(posOfQueenList[eachQueen].x == N-1 && posOfQueenList[eachQueen].y == qC)
                                attkCount++;
                            if(posOfQueenList[eachQueen].x == N-1 && posOfQueenList[eachQueen].y == qC-1)
                                attkCount++;
                            if(posOfQueenList[eachQueen].x == N-1 && posOfQueenList[eachQueen].y == qC+1)
                                attkCount++;
                        }
                        break;

        case NON_CORNER_CASE_BTM:
                        for( eachQueen = 0; eachQueen < queenListSize; eachQueen++) {
                            if(posOfQueenList[eachQueen].x == 0 && posOfQueenList[eachQueen].y == qC)
                                attkCount++;
                            if(posOfQueenList[eachQueen].x == 0 && posOfQueenList[eachQueen].y == qC-1)
                                attkCount++;
                            if(posOfQueenList[eachQueen].x == 0 && posOfQueenList[eachQueen].y == qC+1)
                                attkCount++;
                        }
                        break;

        case NON_CORNER_CASE_LEFT:
                        for( eachQueen = 0; eachQueen < queenListSize; eachQueen++) {
                            if(posOfQueenList[eachQueen].x == qR && posOfQueenList[eachQueen].y == N-1)
                                attkCount++;
                            if(posOfQueenList[eachQueen].x == qR+1 && posOfQueenList[eachQueen].y == N-1)
                                attkCount++;
                            if(posOfQueenList[eachQueen].x == qR-1 && posOfQueenList[eachQueen].y == N-1)
                                attkCount++;
                        }
                        break;

        case NON_CORNER_CASE_RIGHT:
                        for( eachQueen = 0; eachQueen < queenListSize; eachQueen++) {
                            if(posOfQueenList[eachQueen].x == qR && posOfQueenList[eachQueen].y == 0)
                                attkCount++;
                            if(posOfQueenList[eachQueen].x == qR+1 && posOfQueenList[eachQueen].y == 0)
                                attkCount++;
                            if(posOfQueenList[eachQueen].x == qR-1 && posOfQueenList[eachQueen].y == 0)
                                attkCount++;
                        }
                        break;

        case TOP_LEFT_CORNER:
                        for( eachQueen = 0; eachQueen < queenListSize; eachQueen++) {
                            if(posOfQueenList[eachQueen].x == qR && posOfQueenList[eachQueen].y == N-1)
                                attkCount++;
                            if(posOfQueenList[eachQueen].x == qR+1 && posOfQueenList[eachQueen].y == N-1)
                                attkCount++;
                            if(posOfQueenList[eachQueen].x == N-1 && posOfQueenList[eachQueen].y == qC)
                                attkCount++;
                        }
                        break;

        case TOP_RIGHT_CORNER:
                        for( eachQueen = 0; eachQueen < queenListSize; eachQueen++) {
                            if(posOfQueenList[eachQueen].x == qR && posOfQueenList[eachQueen].y == 0)
                                attkCount++;
                            if(posOfQueenList[eachQueen].x == qR+1 && posOfQueenList[eachQueen].y == 0)
                                attkCount++;
                            if(posOfQueenList[eachQueen].x == N-1 && posOfQueenList[eachQueen].y == qC)
                                attkCount++;
                        }
                        break;

        case BTM_LEFT_CORNER:
                        for( eachQueen = 0; eachQueen < queenListSize; eachQueen++) {
                            if(posOfQueenList[eachQueen].x == 0 && posOfQueenList[eachQueen].y == qC)
                                attkCount++;
                            if(posOfQueenList[eachQueen].x == qR && posOfQueenList[eachQueen].y == N-1)
                                attkCount++;
                            if(posOfQueenList[eachQueen].x == qR-1 && posOfQueenList[eachQueen].y == N-1)
                                attkCount++;
                        }
                        break;

        case BTM_RIGHT_CORNER:
                        for( eachQueen = 0; eachQueen < queenListSize; eachQueen++) {
                            if(posOfQueenList[eachQueen].x == 0 && posOfQueenList[eachQueen].y == qC)
                                attkCount++;
                            if(posOfQueenList[eachQueen].x == qR && posOfQueenList[eachQueen].y == 0)
                                attkCount++;
                            if(posOfQueenList[eachQueen].x == qR-1 && posOfQueenList[eachQueen].y == 0)
                                attkCount++;
                        }
                        break;

        default: return 0;
    }

    return attkCount;
}

// Check each Queen can attack how many other Queen
// return K value if each Queen attk K times
// return -1 otherwise
int checkAllQueen_K_value(Point* posOfQueenList, int wrapArnd, int N) {
    int arrSize = sizeof(posOfQueenList); // Number of Queens
    int eachQueen_K_Value[arrSize];
    memset(eachQueen_K_Value, 0, arrSize*sizeof(int));

    int curr = 0;
    for( curr = 0; curr < arrSize; curr++) {
        int attkTarget = 0;
        int other = 0;
        for( other = 0; other < arrSize; other++) {

            // Check if Curr Queen can attk Other Queen AND yet is not blocked
            // capture the number of such possible attks
            if(canQueenAttack(posOfQueenList[curr].x, posOfQueenList[curr].y, posOfQueenList[other].x, posOfQueenList[other].y) && curr != other) {
                if(checkIfBlock(posOfQueenList[curr].x, posOfQueenList[curr].y, posOfQueenList[other].x, posOfQueenList[other].y, arrSize, posOfQueenList) == false)
                    attkTarget++;
            }

            // Check if curr position of queen is a border case AND is there another Queen at corresponding WrapAround?
            // & finish calculating the K-value of each Queen if we allow wrap around
            if(wrapArnd && isQueenAtBorder(posOfQueenList[curr].x, posOfQueenList[curr].y, N)) {
                int queenPos = queenAtWhichCorner(posOfQueenList[curr].x, posOfQueenList[curr].y, N);
                // calculate Num of other Queen it can attak in its curr position if wrap around
                attkTarget = attkTarget + wrapArndAttkCount(posOfQueenList[curr].x, posOfQueenList[curr].y, N, posOfQueenList, queenPos);
            }

        }
        eachQueen_K_Value[curr] = attkTarget;
    }

    if(checkAllEqual(eachQueen_K_Value,arrSize)) {
        return eachQueen_K_Value[0]; // return K-Value
    }

    return -1; // if each Queen don't share same K-Value
}

// Driver Code
int main(int argc, char *argv[]) {

    int N= 0, K= 0, L= 0, W= 0;

    if(argc == 3) {
        N = atoi(argv[1]);
        K = atoi(argv[2]);
    }

    if(argc == 4) {
        N = atoi(argv[1]);
        K = atoi(argv[2]);
        L = atoi(argv[3]);
    }

    if(argc == 5) {
        N = atoi(argv[1]);
        K = atoi(argv[2]);
        L = atoi(argv[3]);
        W = atoi(argv[4]);
    }

    MPI_Init(0, NULL);

    // Find out rank, size
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Calculate Board size using N
    int boardSize = N*N;

    // Total permutation possible for the board
    int totalPermutationSize = pow(2,boardSize);

    // Container to store all board permutation e.g. [2^16][16]
    int result[totalPermutationSize][boardSize];
    memset(result, 0, totalPermutationSize*(boardSize) * sizeof(int));

    // Create all permutation for ChessBoard - e.g. permute [2^16][16]
    generateAllBoardPermutation(N, result);

    // Create Array to capture K-value for each Board
    int indx = 0;
    int K_valueBox[totalPermutationSize];
    memset(K_valueBox, -1, totalPermutationSize * sizeof(int));
    int K_valueResult = 0;
    long long before = 0;
    long long after = 0;

    // Convert each 1D array to 2D array
    int eachPermutation = 0;
    long long start = wall_clock_time();

    // slave 1
    for( eachPermutation = 0; eachPermutation < totalPermutationSize/2; eachPermutation++) {

        if(world_rank == 1) {
            before = wall_clock_time();
            int permuNum = eachPermutation;
            int *eachPermu = result[eachPermutation];
            MPI_Send(eachPermu, boardSize, MPI_INT, (world_rank + 1) % world_size, eachPermutation, MPI_COMM_WORLD);
            MPI_Send(&permuNum, 1, MPI_INT, (world_rank + 3) % world_size, eachPermutation, MPI_COMM_WORLD);
            after = wall_clock_time();
            long long comm_time = after - before;
            //fprintf(stderr, " --- Process %d: communication_time=%6.2f seconds; \n", world_rank, comm_time / 1000000000.0);

        }   

        if(world_rank == 2) {
            before = wall_clock_time();
            int *permutation = malloc(sizeof(int) * boardSize);
            MPI_Recv(permutation, boardSize, MPI_INT, world_rank - 1, eachPermutation, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Perform Conversion
            int **eachBoard_2D = convertTo2DArray(boardSize, permutation);
            Point* posOfQueen = getAllQueenPosition(N, eachBoard_2D);
            MPI_Send(posOfQueen, boardSize*4*2, MPI_BYTE, (world_rank + 1) % world_size, eachPermutation, MPI_COMM_WORLD);
            after = wall_clock_time();
            long long comm_time = after - before;
            //fprintf(stderr, " --- Process %d: communication_time=%6.2f seconds; \n", world_rank, comm_time / 1000000000.0);
        }

        // Check this Board's Queen K-Value
        if(world_rank == 3) {
            before = wall_clock_time();
            Point* posOfQueenRecv = malloc(sizeof(Point) * boardSize); // need to create mem to do the recv
            MPI_Recv(posOfQueenRecv, boardSize*4*2, MPI_BYTE, world_rank - 1, eachPermutation, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //fprintf(stderr, "process 3 : (%d,%d)\n", posOfQueenRecv[0].x, posOfQueenRecv[0].y);
            K_valueResult = checkAllQueen_K_value(posOfQueenRecv, W, N);

            MPI_Send(&K_valueResult, 1, MPI_INT, (world_rank + 1) % world_size, eachPermutation, MPI_COMM_WORLD);
            after = wall_clock_time();
            long long comm_time = after - before;
            //fprintf(stderr, " --- Process %d: communication_time=%6.2f seconds; \n", world_rank, comm_time / 1000000000.0);
        }

        // Save the Permutation Number (if current permutation matches user input K)
        if(world_rank == 4) {
            before = wall_clock_time();
            int *permNumRecv = malloc(sizeof(int));
            MPI_Recv(&K_valueResult, 1, MPI_INT, world_rank - 1, eachPermutation, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(permNumRecv, boardSize, MPI_INT, world_rank - 3, eachPermutation, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            if(K_valueResult == K) {
                //K_valueBox[indx] = eachPermutation;
                K_valueBox[indx] = *permNumRecv;
                indx++;
                //fprintf(stderr, " --- permutation recv is %d and its k value result: %d; \n", permNumRecv[0], K_valueResult);
            }

            if(eachPermutation == (totalPermutationSize/2)-1) {
                fprintf(stderr, " --- k value box sent: %d; \n", K_valueBox[0]);
                MPI_Send(K_valueBox, totalPermutationSize, MPI_INT, (world_rank - 4) % world_size, eachPermutation, MPI_COMM_WORLD);
            }

            after = wall_clock_time();
            long long comm_time = after - before;
            //fprintf(stderr, " --- Process %d: communication_time=%6.2f seconds; \n", world_rank, comm_time / 1000000000.0);
        }
    }

    // slave 2
    for( eachPermutation = totalPermutationSize/2; eachPermutation < totalPermutationSize; eachPermutation++) {

        if(world_rank == 5) {
            before = wall_clock_time();
            int permuNum = eachPermutation;
            int *eachPermu = result[eachPermutation];
            MPI_Send(eachPermu, boardSize, MPI_INT, (world_rank + 1) % world_size, eachPermutation, MPI_COMM_WORLD);
            MPI_Send(&permuNum, 1, MPI_INT, (world_rank + 3) % world_size, eachPermutation, MPI_COMM_WORLD);
            after = wall_clock_time();
            long long comm_time = after - before;
            //fprintf(stderr, " --- Process %d: communication_time=%6.2f seconds; \n", world_rank, comm_time / 1000000000.0);

        }   

        if(world_rank == 6) {
            before = wall_clock_time();
            int *permutation = malloc(sizeof(int) * boardSize);
            MPI_Recv(permutation, boardSize, MPI_INT, world_rank - 1, eachPermutation, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Perform Conversion
            int **eachBoard_2D = convertTo2DArray(boardSize, permutation);
            Point* posOfQueen = getAllQueenPosition(N, eachBoard_2D);
            MPI_Send(posOfQueen, boardSize*4*2, MPI_BYTE, (world_rank + 1) % world_size, eachPermutation, MPI_COMM_WORLD);
            after = wall_clock_time();
            long long comm_time = after - before;
            //fprintf(stderr, " --- Process %d: communication_time=%6.2f seconds; \n", world_rank, comm_time / 1000000000.0);
        }

        // Check this Board's Queen K-Value
        if(world_rank == 7) {
            before = wall_clock_time();
            Point* posOfQueenRecv = malloc(sizeof(Point) * boardSize); // need to create mem to do the recv
            MPI_Recv(posOfQueenRecv, boardSize*4*2, MPI_BYTE, world_rank - 1, eachPermutation, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //fprintf(stderr, "process 3 : (%d,%d)\n", posOfQueenRecv[0].x, posOfQueenRecv[0].y);
            K_valueResult = checkAllQueen_K_value(posOfQueenRecv, W, N);

            MPI_Send(&K_valueResult, 1, MPI_INT, (world_rank + 1) % world_size, eachPermutation, MPI_COMM_WORLD);
            after = wall_clock_time();
            long long comm_time = after - before;
            //fprintf(stderr, " --- Process %d: communication_time=%6.2f seconds; \n", world_rank, comm_time / 1000000000.0);
        }

        // Save the Permutation Number (if current permutation matches user input K)
        if(world_rank == 8) {
            before = wall_clock_time();
            int *permNumRecv = malloc(sizeof(int));
            MPI_Recv(&K_valueResult, 1, MPI_INT, world_rank - 1, eachPermutation, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(permNumRecv, boardSize, MPI_INT, world_rank - 3, eachPermutation, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            if(K_valueResult == K) {
                //K_valueBox[indx] = eachPermutation;
                K_valueBox[indx] = *permNumRecv;
                indx++;
                //fprintf(stderr, " --- permutation recv is %d and its k value result: %d; \n", permNumRecv[0], K_valueResult);
            }

            if(eachPermutation == totalPermutationSize-1) {
                fprintf(stderr, " --- k value box sent: %d; \n", K_valueBox[0]);
                MPI_Send(K_valueBox, totalPermutationSize, MPI_INT, (world_rank - 8) % world_size, eachPermutation, MPI_COMM_WORLD);
            }

            after = wall_clock_time();
            long long comm_time = after - before;
            //fprintf(stderr, " --- Process %d: communication_time=%6.2f seconds; \n", world_rank, comm_time / 1000000000.0);
        }
    }

    // Merge Slave 1 & Slave 2 result together
    // Process Calculate MAX K-Value & print result
    if(world_rank == 0) {

        before = wall_clock_time();

        int p = 0;
        for(p = 4; p <= 8; p = p + 4) {
                // Recv Operation here
                if(p== 4)
                    MPI_Recv(K_valueBox, totalPermutationSize, MPI_INT, (world_rank + p), (totalPermutationSize/2)-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                else
                    MPI_Recv(K_valueBox, totalPermutationSize, MPI_INT, (world_rank + p), totalPermutationSize-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                //fprintf(stderr, " --- k value box recv is %d \n", K_valueBox[0]);

                int max = 0;
                int i = 0; 

                // Loop through permutation to find Max num Of Queen
                for( i = 0; i < totalPermutationSize; i++) {
                    if(K_valueBox[i] != -1) {
                        int indxNum = K_valueBox[i];
                        int **eachBoard_2D = convertTo2DArray(boardSize, result[indxNum]);
                        int queenCount = numOfQueen(N, eachBoard_2D);
                        if(queenCount > max)
                            max = queenCount;
                    }
                }

                // Loop entire K_value Box to Display Result
                for( i = 0; i < totalPermutationSize; i++) {

                    if(K_valueBox[i] != -1) {
                        int indxNum = K_valueBox[i];
                        int **eachBoard_2D = convertTo2DArray(boardSize, result[indxNum]);
                        int queenCount = numOfQueen(N, eachBoard_2D);

                        if(queenCount == max) {
                            // Control L display
                            printf("%d,%d:%d:", N,K,queenCount);
                            if(!L) printf("\n");
                            if(L) printQueenLocation(N,eachBoard_2D);
                        }
                    }
                }
        }

        after = wall_clock_time();
        long long comm_time = after - before;
        fprintf(stderr, " --- Process %d: communication_time=%6.2f seconds; \n", world_rank, comm_time / 1000000000.0);
    }
    

    long long end = wall_clock_time();
    long long comp_time = end - start;
    fprintf(stderr, " --- Total computation_time=%6.2f seconds\n", comp_time / 1000000000.0);

    MPI_Finalize();

    return 0;
}