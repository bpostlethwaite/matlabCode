////////////////////////////////////////////////////////////////////////////////
// File: get_column.c                                                         //
// Routine(s):                                                                //
//    Get_Column                                                              //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Get_Column(double v[], double *A, int col, int nrows, int ncols)     //
//                                                                            //
//  Description:                                                              //
//     Copy the column 'col' from the nrows x ncols matrix A to the vector v. //
//     Note that v should be declared "double v[N]", with N >= nrows in the   //
//     calling routine.                                                       //
//                                                                            //
//  Arguments:                                                                //
//     double v[]   Destination address of the column col of the matrix A.    //
//     double *A    Pointer to the first element of the matrix A.             //
//     int    col   The column of matrix A to copy to the vector v,           //
//                  0 <= col < ncols.                                         //
//     int    nrows The number of rows matrix A.                              //
//     int    ncols The number of columns of the matrix A.                    //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N],  v[N];                                                 //
//                                                                            //
//     (your code to set the matrix A and the column number i)                //
//                                                                            //
//     if ( (i >= 0) && (i < N) ) Get_Column(v, &A[0][0], i, M, N);           //
//     printf("The vector v is \n"); ...                                      //
////////////////////////////////////////////////////////////////////////////////
void Get_Column(double v[], double *A, int col, int nrows, int ncols)
{
   int i;

   for (A += col, i = 0; i < nrows; A += ncols, i++) v[i] = *A;
}
