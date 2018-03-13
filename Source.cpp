/*
	Source.cpp
	Author:		Jake Wagner
	About:		This project is designed to test the time complexity of three different methods
				for matrix multiplication: Naive/Traditional (theoretical O(n^3)), Divide and 
				Conquer (theoretical O(n^3)), and Strassen's (theoretical O(n^2.8)). It assumes
				that all matrices are square and of size 2^k, where integer k >= 1.
*/
#include <iostream>
#include <time.h>
using namespace std;


//CONTROLS ========================================================================================================
const int n = 4; //size of matrices
const int t = 1000; //numbers of trials of each algorithm for average time
const bool v = false; //selection for displaying matrices or not

//STRUCTS =========================================================================================================
struct corners { int topRow, btmRow, topCol, btmCol; };//contains corners of matrices/submatrices

//FUNCTION PROTOTYPES =============================================================================================
void dispMtx(int Z[n][n], int N);
void naiveMM(int A[n][n], int B[n][n], int C[n][n], int N);
void daqMM(int A[n][n], int B[n][n], int C[n][n], int N);
void daqMul(int A[n][n], int B[n][n], int C[n][n], corners crnA, corners crnB, corners crnC);
void setCorners(corners crnZ, int row, int col, corners *Z_ii);
void strassenMM(int A[n][n], int B[n][n], int C[n][n], int N);
void strassenMul(int A[n][n], int B[n][n], int C[n][n], int N, corners crnA, corners crnB, corners crnC);
void mtxAdd(int A[n][n], int B[n][n], int C[n][n], int N, corners crnA, corners crnB, corners crnC);
void mtxSub(int A[n][n], int B[n][n], int C[n][n], int N, corners crnA, corners crnB, corners crnC);
void reset(int Z[n][n], corners crnZ);


//main driver
int main() {
	//COMMON VARIABLES
	int mtxA[n][n]; //first input matrix
	int mtxB[n][n]; //second input matrix
	int mtxC[n][n]; //resultant matrix
	clock_t temp;
	clock_t timeResult = 0;

    //First instatiate the matrices
    srand(time(0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            mtxA[i][j] = rand() % 10; //range of 0-10;
            mtxB[i][j] = rand() % 10; //range of 0-10
   	        mtxC[i][j] = 0; //initiate to zero
        }
    }
    
   	if (v) {
        cout << "Matrix A ===============\n";
		dispMtx(mtxA, n);
		cout << "Matrix B ===============\n";
		dispMtx(mtxB, n);
    }
    
    //Display Run Info
	cout << "For square matrices size N = " << n << endl << "Number of trials ran: " << t << endl;

    for(int i = 0; i < n; ++i) {
        //call naive
        temp = clock();	
        naiveMM(mtxA, mtxB, mtxC, n);
        temp = clock() - temp;
        timeResult += temp;
    }
    cout << "Naive/Traditional Matrix Multiplication average time: " << ((float)timeResult) / CLOCKS_PER_SEC / t << " seconds.\n";
   	//display?
	if (v) {
		cout << "Matrix C ===============\n";
		dispMtx(mtxC, n);
	}
    
    //reset resultatnt matrix and time
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
   	        mtxC[i][j] = 0; //initiate to zero
        }
    }
    timeResult = 0;
    
    
    for(int i = 0; i < n; ++i) {
        //call divide and conquer
        temp = clock();
        daqMM(mtxA, mtxB, mtxC, n);
        temp = clock() - temp;
        timeResult += temp;
    }
    cout << "Divide and Conquer Matrix Multiplication average time: " << ((float)timeResult) / CLOCKS_PER_SEC / t << " seconds.\n";
   	//display?
	if (v) {
		cout << "Matrix C ===============\n";
		dispMtx(mtxC, n);
	}
     
   	//reset resultatnt matrix and time
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
   	        mtxC[i][j] = 0; //initiate to zero
        }
    }
    timeResult = 0;
    
    
    for(int i = 0; i < n; ++i) {
        //call strassen
        temp = clock();
	    strassenMM(mtxA, mtxB, mtxC, n);
        cout << "Strassen's Matrix Multiplication average time: " << ((float)timeResult) / CLOCKS_PER_SEC / t << " seconds.\n";
	}
	//display?
	if (v) {
		cout << "Matrix C ===============\n";
		dispMtx(mtxC, n);
	}
	
	system("pause");

return 0;
}//end main



//FUNCTION DECLARATIONS ===========================================================================================
/*
	dispMtx
	- displays the given matrix
	@param Z the matrix
	@param N the size of the square matrix
*/
void dispMtx(int Z[n][n], int N) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			cout << Z[i][j] << " ";
		}
		cout << endl;
	}
}//end dispMtx


/*
	naiveMM
	- matrix multiplication in the traditional definition, three forloops
	@param A first input matrix
	@param B second input matrix
	@param C the resultant matrix
	@param N size of the square matrices
*/
void naiveMM(int A[n][n], int B[n][n], int C[n][n], int N) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}//end naiveMM


/*
	daqMM
	- implements matrix multiplication using a divide and conquer method
	@param A first input matrix
	@param B second input matrix
	@param C the resultant matrix
	@param N size of the square matrices
*/
void daqMM(int A[n][n], int B[n][n], int C[n][n], int N) {
	//set original corners, the whole of starting matrices
	corners crnA = { 0, N, 0, N };
	corners crnB = { 0, N, 0, N };
	corners crnC = { 0, N, 0, N };
	//call recursive engine
	daqMul(A, B, C, crnA, crnB, crnC);
}//end daqMM


/*
	daqMul
	- the recursive engine for the daqMM
	@param A first input matrix
	@param B second input matrix
	@param C the resultant matrix
	@param crnA contains the corners for mtx A
	@param crnB contains the corners for mtx B
	@param crnC contains the corners for mtx C
*/
void daqMul(int A[n][n], int B[n][n], int C[n][n], corners crnA, corners crnB, corners crnC) {
	//12 submatrices
	corners A_ii[2][2], B_ii[2][2], C_ii[2][2];

	//check for base cases
	if ((crnA.btmRow - crnA.topRow) < 0)
		return;
	if ((crnA.btmCol - crnA.topCol) == 1) {
		C[crnC.topRow][crnC.topCol] += (A[crnA.topRow][crnA.topCol] * B[crnB.topRow][crnB.topCol]);
		return;
	}

	//set corners for submatrices
	for (int i = 0; i<2; i++) {
		for (int j = 0; j<2; j++) {
			setCorners(crnA, i, j, &A_ii[i][j]);
			setCorners(crnB, i, j, &B_ii[i][j]);
			setCorners(crnC, i, j, &C_ii[i][j]);
		}
	}

	//call the 8 submatrices recursivvely
	//
	//    A			B				C
	// --  ---    --  ---    -------  -------
	// | a b |  x | e f | =  | ae+bg  af+bh |
	// | c d |    | g h |    | ce+dg  cf+dh |
	// --  ---    --  ---    -------  -------
	//
	//C_11 = 
	daqMul(A, B, C, A_ii[0][0], B_ii[0][0], C_ii[0][0]);//ae
	daqMul(A, B, C, A_ii[0][1], B_ii[1][0], C_ii[0][0]);//+ bg
	//C_12 =
	daqMul(A, B, C, A_ii[0][0], B_ii[0][1], C_ii[0][1]);//af
	daqMul(A, B, C, A_ii[0][1], B_ii[1][1], C_ii[0][1]);//+ bh
	//C_21 =
	daqMul(A, B, C, A_ii[1][0], B_ii[0][0], C_ii[1][0]);//ce
	daqMul(A, B, C, A_ii[1][1], B_ii[1][0], C_ii[1][0]);//+ dg
	//C_22 =
	daqMul(A, B, C, A_ii[1][0], B_ii[0][1], C_ii[1][1]);//cf
	daqMul(A, B, C, A_ii[1][1], B_ii[1][1], C_ii[1][1]);//+ dh
}//end daqMul


/*
	setCorners
	- sets the corners for the 4 submatrices from the given matrix's corners
	@param crnZ the corners of a given matrix
	@param row index for the four submatices matrix
	@param col index for the matrix of 4 submatrices
	@param Z_ii a reference to the matrix of 4 corners objects
*/
void setCorners(corners crnZ, int row, int col, corners *Z_ii) {
	//get new dimmensions N/2
	int newRow = crnZ.topRow + ((crnZ.btmRow - crnZ.topRow) / 2);
	int newCol = crnZ.topCol + ((crnZ.btmCol - crnZ.topCol) / 2);
	//give submatrix bounds of matrix
	*Z_ii = crnZ;
	//assign new values
	if (row == 0)
		Z_ii->btmRow = newRow;// top row
	else
		Z_ii->topRow = newRow;// bottom row
	if (col == 0)
		Z_ii->btmCol = newCol;// left col
	else
		Z_ii->topCol = newCol;// right col
}//end setCorners


/*
	strassenMM
	- implements matrix multiplication using Strassen's Method
	@param A first input matrix
	@param B second input matrix
	@param C the resultant matrix
	@param N size of the square matrices
*/
void strassenMM(int A[n][n], int B[n][n], int C[n][n], int N) {
	//set original corners, the whole of starting matrices
	corners crnA = { 0, N, 0, N };
	corners crnB = { 0, N, 0, N };
	corners crnC = { 0, N, 0, N };
	//call recursive engine
	strassenMul(A, B, C, N, crnA, crnB, crnC);
}//end strassenMM


/*
	strassemMul
	- thje recursive engine for strassenMM
	@param A first input matrix
	@param B second input matrix
	@param C the resultant matrix
	@param N size of the square matrices
	@param crnA contains the corners for mtx A
	@param crnB contains the corners for mtx B
	@param crnC contains the corners for mtx C
*/
void strassenMul(int A[n][n], int B[n][n], int C[n][n], int N, corners crnA, corners crnB, corners crnC) {
	//check for base case
	if ((crnC.btmRow - crnC.topRow) == 1) {
		C[crnC.topRow][crnC.topCol] += A[crnA.topRow][crnA.topCol] * B[crnB.topRow][crnB.topCol];
		return;
	}

	//otherwise calculate recursively...
	//local variables
	corners A_ii[2][2], B_ii[2][2], C_ii[2][2], crnM; //12 submatrices plus temp corners for Ms
	int m1[n][n], m2[n][n], m3[n][n], m4[n][n], m5[n][n], m6[n][n], m7[n][n], AA[n][n], BB[n][n];//temp matrices

	//set corners of submatrices
	for (int i = 0; i<2; i++) {
		for (int j = 0; j<2; j++) {
			setCorners(crnA, i, j, &A_ii[i][j]);
			setCorners(crnB, i, j, &B_ii[i][j]);
			setCorners(crnC, i, j, &C_ii[i][j]);
		}
	}
	int m = (crnA.btmRow - crnA.topRow) / 2;
	crnM = { 0, m, 0, m };

	//clear new matrices
	for (int i = 0; i < 7; i++) {
		reset(m1, crnM);
		reset(m2, crnM);
		reset(m3, crnM);
		reset(m4, crnM);
		reset(m5, crnM);
		reset(m6, crnM);
		reset(m7, crnM);
	}

	//calculate
	reset(AA, crnM);
	reset(BB, crnB);
	//m1 = (a_11 + a_22)(b_11 + b_22)
	mtxAdd(A, A, AA, N/2, A_ii[0][0], A_ii[1][1], crnM);
	mtxAdd(B, B, BB, N/2, B_ii[0][0], B_ii[1][1], crnM);
	strassenMul(AA, BB, m1, N/2, crnM, crnM, crnM);

	reset(AA, crnM);
	reset(BB, crnB);
	//m2 = (a_21 + a_22)b_11
	mtxAdd(A, A, AA, N/2, A_ii[1][0], A_ii[1][1], crnM);
	strassenMul(AA, B, m2, N/2, crnM, B_ii[0][0], crnM);

	reset(AA, crnM);
	reset(BB, crnB);
	//m3 = a_11(b_12 - b_22)
	mtxSub(B, B, BB, N / 2, B_ii[0][1], B_ii[1][1], crnM);
	strassenMul(A, BB, m3, N / 2, A_ii[0][0], crnM, crnM);

	reset(AA, crnM);
	reset(BB, crnB);
	//m4 = a_22(b_21 - b_11)
	mtxSub(B, B, BB, N / 2, B_ii[1][0], B_ii[0][0], crnM);
	strassenMul(A, BB, m4, N / 2, A_ii[1][1], crnM, crnM);

	reset(AA, crnM);
	reset(BB, crnB);
	//m5 = (a_11 + a_12)b_22
	mtxAdd(A, A, AA, N / 2, A_ii[0][0], A_ii[0][1], crnM);
	strassenMul(AA, B, m5, N / 2, crnM, B_ii[1][1], crnM);

	reset(AA, crnM);
	reset(BB, crnB);
	//m6 = (a_21 - a_11)(b_11 + b_12)
	mtxSub(A, A, AA, N / 2, A_ii[1][0], A_ii[0][0], crnM);
	mtxAdd(B, B, BB, N / 2, B_ii[0][0], B_ii[0][1], crnM);
	strassenMul(AA, BB, m6, N / 2, crnM, crnM, crnM);

	reset(AA, crnM);
	reset(BB, crnB);
	//m7 = (a_12 - a_22)(b_21 + b_22)
	mtxSub(A, A, AA, N / 2, A_ii[0][1], A_ii[1][1], crnM);
	mtxAdd(B, B, BB, N / 2, B_ii[1][0], B_ii[1][1], crnM);
	strassenMul(AA, BB, m7, N / 2, crnM, crnM, crnM);

	//load to C
	//c_11 = m1 + m4 - m5 + m7
	mtxAdd(m1, m4, AA, N / 2, crnM, crnM, crnM);//m1+m4
	mtxSub(AA, m5, BB, N / 2, crnM, crnM, crnM);//m1+m4-m5
	mtxAdd(BB, m7, C, N / 2, crnM, crnM, C_ii[0][0]);//m1+m4-m5+m7

	//c_12 = m3 + m5
	mtxAdd(m3, m5, C, N / 2, crnM, crnM, C_ii[0][1]);//m3+m5

	//c_21 = m2 + m4
	mtxAdd(m2, m4, C, N / 2, crnM, crnM, C_ii[1][0]);//m2+m4

	//c_22 = m1 + m3 - m2 + m6
	mtxAdd(m1, m3, AA, N / 2, crnM, crnM, crnM);//m1+m3
	mtxSub(AA, m2, BB, N / 2, crnM, crnM, crnM);//m1+m3-m2
	mtxAdd(BB, m6, C, N / 2, crnM, crnM, C_ii[1][1]);//m1+m3-m2+m6
}//end strassenMul


/*
	reset
	- sets the matrix values to 0
	@param Z the matrix to be reset
	@param crnZ the corner object containing the bounds for what is to be reset
*/
void reset(int Z[n][n],  corners crnZ) {
	for (int i = crnZ.topRow; i < crnZ.btmRow; i++) {
		for (int j = crnZ.topCol; j < crnZ.btmCol; j++) {
			Z[i][j] = 0;
		}
	}
}//end reset


/*
	mtxAdd
	- adds two matrices together and stores them
	@param A first input matrix
	@param B second input matrix
	@param C the resultant matrix
	@param crnA contains the corners for mtx A
	@param crnB contains the corners for mtx B
	@param crnC contains the corners for mtx C
*/
void mtxAdd(int A[n][n], int B[n][n], int C[n][n], int N, corners crnA, corners crnB, corners crnC) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			C[i + crnC.topRow][j + crnC.topCol] = A[i + crnA.topRow][j + crnA.topCol] + B[i + crnB.topRow][j + crnB.topCol];
		}
	}
}//end mtxAdd


 /*
 mtxSub
 - subtracts one matrix from another and stores the result
 @param A first input matrix
 @param B second input matrix
 @param C the resultant matrix
 @param crnA contains the corners for mtx A
 @param crnB contains the corners for mtx B
 @param crnC contains the corners for mtx C
 */
void mtxSub(int A[n][n], int B[n][n], int C[n][n], int N, corners crnA, corners crnB, corners crnC) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			C[i+crnC.topRow][j+crnC.topCol] = A[i+crnA.topRow][j+crnA.topCol] - B[i+crnB.topRow][j+crnB.topCol];
		}
	}
}//end mtxSub
