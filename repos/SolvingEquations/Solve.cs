using System;
using System.Collections;
namespace SolvingEquations
{
    public class Solve
    {
        static void matrixTransport(double[,] a, ref double[,] b, int row, int column)
        {//row和column均为a的行和列 
            int i, j;
            //int row, column; 
            for (i = 0; i < column; i++)
            {
                for (j = 0; j < row; j++)
                {
                    b[i, j] = a[j, i];
                }
            }
        }
        /* 后继 */
        static int getNext(int i, int m, int n)
        {
            return (i % n) * m + i / n;
        }

        /* 前驱 */
        static int getPre(int i, int m, int n)
        {
            return (i % m) * n + i / m;
        }

        /* 处理以下标i为起点的环 */
        static void movedata(double[,] mtx, int i, int m, int n)
        {
            double temp = mtx[i / n, i % n]; // 暂存
            int cur = i;    // 当前下标
            int pre = getPre(cur, m, n);
            while (pre != i)
            {
                mtx[cur / n, cur % n] = mtx[pre / n, pre % n];
                cur = pre;
                pre = getPre(cur, m, n);
            }
            mtx[cur / n, cur % n] = temp;
        }

        /* 转置，即循环处理所有环 */
        static void transpose(ref double[,] mtx, int m, int n)
        {
            for (int i = 0; i < m * n; ++i)
            {
                int next = getNext(i, m, n);
                while (next > i) // 若存在后继小于i说明重复,就不进行下去了（只有不重复时进入while循环）
                    next = getNext(next, m, n);
                if (next == i)  // 处理当前环
                    movedata(mtx, i, m, n);
            }
        }

        //LUP分解
        static void LUP_Descomposition(double[,] A, ref double[,] L, ref double[,] U, ref int[] P, int N)
        {
            int row = 0;
            for (int i = 0; i < N; i++)
            {
                P[i] = i;
            }
            for (int i = 0; i < N - 1; i++)
            {
                double p = 0.0;
                for (int j = i; j < N; j++)
                {
                    if (System.Math.Abs(A[j, i]) > p)
                    {
                        p = System.Math.Abs(A[j, i]);
                        row = j;
                    }
                }
                if (0 == p)
                {
                    Console.WriteLine("矩阵奇异，无法计算逆");
                    return;
                }

                //交换P[i]和P[row]
                int tmp = P[i];
                P[i] = P[row];
                P[row] = tmp;

                double tmp2 = 0.0;
                for (int j = 0; j < N; j++)
                {
                    //交换A[i][j]和 A[row][j]
                    tmp2 = A[i, j];
                    A[i, j] = A[row, j];
                    A[row, j] = tmp2;
                }

                //以下同LU分解
                double u = A[i, i], l = 0.0;
                for (int j = i + 1; j < N; j++)
                {
                    l = A[j, i] / u;
                    A[j, i] = l;
                    for (int k = i + 1; k < N; k++)
                    {
                        A[j, k] = A[j, k] - A[i, k] * l;
                    }
                }

            }

            //构造L和U
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    if (i != j)
                    {
                        L[i, j] = A[i, j];
                    }
                    else
                    {
                        L[i, j] = 1;
                    }
                }
                for (int k = i; k < N; k++)
                {
                    U[i, k] = A[i, k];
                }
            }

        }

        //LUP求解方程
        static double[] LUP_Solve(double[,] L, double[,] U, int[] P, double[] b, int N)
        {
            double[] x = new double[20];
            double[] y = new double[20];
            //double* x = new double[N]();
            //double* y = new double[N]();

            //正向替换
            for (int i = 0; i < N; i++)
            {
                y[i] = b[P[i]];
                for (int j = 0; j < i; j++)
                {
                    y[i] = y[i] - L[i, j] * y[j];
                }
            }
            //反向替换
            for (int i = N - 1; i >= 0; i--)
            {
                x[i] = y[i];
                for (int j = N - 1; j > i; j--)
                {
                    x[i] = x[i] - U[i, j] * x[j];
                }
                x[i] /= U[i, i];
            }
            return x;
        }

        //LUP求逆(将每列b求出的各列x进行组装)
        static void LUP_solve_inverse(double[,] A, ref double[,] inv_A, int N)
        {
            //创建矩阵A的副本，注意不能直接用A计算，因为LUP分解算法已将其改变
            double[,] A_mirror = new double[20, 20];
            double[] inv_A_each = new double[20];//矩阵逆的各列
            double[] b = new double[20];//b阵为B阵的列矩阵分量
  
            for (int i = 0; i < N; i++)
            {
                double[,] L = new double[20, 20];//N*N
                double[,] U = new double[20, 20];//N*N
                int[] P = new int[20];//N
 
                //构造单位阵的每一列
                for (int i0 = 0; i0 < N; i0++)
                {
                    b[i0] = 0;
                }
                b[i] = 1;

                //每次都需要重新将A复制一份
                for (int i0 = 0; i0 < N; i0++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        A_mirror[i0, j] = A[i0, j];
                    }
                }
                LUP_Descomposition(A_mirror, ref L, ref U, ref P, N);

                inv_A_each = LUP_Solve(L, U, P, b, N); 
                for (int i0 = 0; i0 < N; i0++)
                {
                    inv_A[i, i0] = inv_A_each[i0];
                }
            }
            transpose(ref inv_A, N, N);//由于现在根据每列b算出的x按行存储，因此需转置
        }

        /********************************************
        A:矩阵A；B:矩阵B；C:相乘结果矩阵；rowA：A的行数；columnB：B的列数；columnA：A的列数
        ********************************************/
        //必须要保证两个数组能够进行乘法运算，这里没有加入判断条件，默认能够进行，望读者能够明白这一点。 
        static void matrixMul(double[,] A, double[,] B, ref double[,] C, int rowA, int columnB, int columnA)
        {
            for (int i = 0; i < rowA; i++)
            {
                for (int j = 0; j < columnB; j++)
                {
                    C[i, j] = 0;
                    for (int k = 0; k < columnA; k++)
                    {
                        C[i, j] += A[i, k] * B[k, j];
                    }
                }
            }
        }

        /// <summary>
        /// 最小二乘法求解，
        /// </summary>
        /// <param name="A">方程组A矩阵，m行n列</param>
        /// <param name="b">方程组b矩阵，m行1列</param>
        /// <param name="x">方程组结果x矩阵，n行1列</param>
        /// <param name="m">表示A矩阵行数</param>
        /// <param name="n">表示A矩阵列数</param>
        /// <param name="info">错误信息</param>
        /// <returns>是否正常执行</returns>
        /// <remarks>x = inv(A'*A)*A'*b; A'表示A的转置，inv()表示求逆</remarks>
        public static bool Maqr(double[,] A, double[,] b, ref double[,] x, int m, int n, ref string info)
        {
            try
            {
                double[,] At = new double[20, 20];//n*m 
                matrixTransport(A, ref At, m, n);

                double[,] mul = new double[20, 20];// n*n 
                matrixMul(At, A, ref mul, n, n, m);

                double[,] invOfA = new double[20, 20];//n*n 
                LUP_solve_inverse(mul, ref invOfA, n);

                double[,] mul1 = new double[20, 20];//n*m 
                matrixMul(invOfA, At, ref mul1, n, m, n);
                matrixMul(mul1, b, ref x, n, 1, m); 
            }
            catch(Exception ex)
            {
                info = ex.ToString();
                return false;
            }
            return true;
        }
    }
}
