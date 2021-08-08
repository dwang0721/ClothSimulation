using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class MathLib : MonoBehaviour
{

    static public void MatrixAdd2D(ref double[,] A, ref double[,] B)
    {
        int ay = A.GetLength(0); 
        int ax = A.GetLength(1);

        int by = B.GetLength(0);
        int bx = B.GetLength(1);

        Debug.Assert(ax == bx);
        Debug.Assert(ay == by);

        for (int row = 0; row < by; row++)
        {
            for (int col = 0; col < bx; col++)
            {
                B[row, col] += A[row, col];
            }
        }
    }

    static public void MatrixAdd3D(ref double[,,] A, ref double[,,] B)
    {
        int ay = A.GetLength(0);
        int ax = A.GetLength(1);
        int az = A.GetLength(2);

        int by = B.GetLength(0);
        int bx = B.GetLength(1);
        int bz = A.GetLength(2);

        Debug.Assert(ax == bx);
        Debug.Assert(ay == by);
        Debug.Assert(az == bz);

        for (int row = 0; row < by; row++)
        {
            for (int col = 0; col < bx; col++)
            {
                for (int i = 0; i < bz; i++)
                {
                    B[row, col, i] += A[row, col, i];
                }
                
            }
        }
    }

    static public void printMatrix2D(ref double[,] M, string info)
    {
        int y = M.GetLength(0);
        int x = M.GetLength(1);
        string msg = info;
        for (int row = 0; row < y; row++)
        {
            msg += "\n|";
            for (int col = 0; col < x; col++)
            {
                msg += M[row, col] + ", ";
            }
            msg += "|";
        }
        msg += "\n";
        Debug.Log(msg);
    }

    static public void printMatrix1D(ref double[] M, string info)
    {
        int length = M.Length;
        string msg = info;

        msg += "\n|";
        for (int col = 0; col < length; col++)
        {
            msg += M[col] + ", ";
        }
        msg += "|";

        msg += "\n";
        Debug.Log(msg);
    }


    static public void testAlglib()
    {
        // https://www.alglib.net/translator/man/manual.csharp.html#sub_spdmatrixcholesky
        // https://comnuan.com/cmnn01014/cmnn01014.php
        // https://stackoverflow.com/questions/5607631/matrix-multiplication-alglib

        int n = 2;
        double[,] a, b;
        alglib.matinvreport rep;
        int info;


        Debug.Log("============== Matrix Inverse Example ==========");
        // input
        a = new double[,] { { 2, 1 }, { 1, 2 } };
        printMatrix2D(ref a, "a value");

        alglib.spdmatrixinverse(ref a, out info, out rep);

        // expected result:
        // |-1/3, 2/3|
        // | 2/3,-1/3|
        printMatrix2D(ref a, "Inverse a:");



        Debug.Log("==============Cheolesky Decompose Example==========");
        // input
        b = new double[,] { { 2, 1 }, { 1, 2 } };
        printMatrix2D(ref b, "b Value:");

        // expected L result:
        // | 1.4142   0.0000  |
        // | 0.7071   1.2247  |
        bool result = alglib.spdmatrixcholesky(ref b, n, false);
        Debug.Log(" b is Symmetric Positive Definite:  " + result);
        printMatrix2D(ref b, "Decompose Result L:");



        Debug.Log("============== Matrix Solver Example ==========");
        // input
        double[,] A = new double[,] { { 5, 8, -4 }, { 6, 9, -5 }, { 4, 7, -2 } };
        double[] B = new double[] { -18, -20, -15 };
        printMatrix2D(ref A, "A Value:");
        printMatrix1D(ref B, "B Value:");

        // expected L result:
        // x = (2, -3, 1)
        alglib.rmatrixsolvefast(A, 3, ref B, out info);
        Debug.Log(" Info:   " + info);
        printMatrix1D(ref B, "Solve Ax=B Result x:");


        Debug.Log("============== Matrix Addition Example ==========");
        // input 
        double[,] C = new double[,] { { 4, 8 }, { 6, 9 } };
        double[,] D = new double[,] { { 1, 2 }, { 2, 3 } };
        
        printMatrix2D(ref C, "C Value:");
        printMatrix2D(ref D, "D Value:");

        // expected L result:
        // | 6   10  |
        // | 8   12  |
        // double[,] I = new double[,] { { 1, 0 }, { 0, 1 } };
        // alglib.rmatrixgemm(2, 2, 2, 1.0, C, 0, 0, 0, I, 0, 0, 0, 1, ref D, 0, 0);
        MatrixAdd2D(ref C, ref D);
        printMatrix2D(ref D, "Compute C + D Result:");

        double[,,] E = new double[,,] {
                                        { { 1, 1, 1 }, { 2, 2, 2 } } ,
                                        { { 3, 3, 3 }, { 4, 4, 4 } } 
                                       };
        double[,,] F = new double[,,] {
                                        { { 1, 2, 3 }, { -2, -1, 0 } } ,
                                        { { 2, 1, 0 }, { 4, 3, 2 } }
                                       };
        MatrixAdd3D(ref E, ref F);

        int m = 3;
        Debug.Log(m);
    }
}
