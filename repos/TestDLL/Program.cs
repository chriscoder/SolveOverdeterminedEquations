using System; 
using System.Collections;
 
using AddDLL;
using SolvingEquations;
namespace TestDLL
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("Hello World!");
            int res = (new MyAdd()).AddInt(10, 30);
            Console.WriteLine("测试DLL！10 + 30 = " + res);
            Console.WriteLine("\n");

            //A*x=b; A表示m行n列矩阵，b表示m行1列矩阵
            int m = 4, n = 3;//m行数，n列数
            double[,] A = new double[20, 20]; 
            A[0,0] = 1; A[0,1] = 2; A[0,2] = 4;
            A[1,0] = 2; A[1,1] = 1; A[1,2] = 1;
            A[2,0] = 1; A[2,1] = 1; A[2,2] = 2;
            A[3,0] = 1; A[3,1] = -1; A[3,2] = -2;
             
            for (int i = 0; i < m; i++)
            {
                string str = "矩阵A的第" + (i + 1) + "行元素：";
                for (int j = 0; j < n; j++)
                    str += " " + A[i,j];
                Console.WriteLine(str);
            }
            Console.WriteLine("\n");

            double[,] b = new double[20, 1]; 
            b[0, 0] = -1;
            b[1, 0] = 4;
            b[2, 0] = 2;
            b[3, 0] = 1;
            for (int i = 0; i < m; i++)
            {
                string str = "矩阵b的第" + (i + 1) + "行元素：";
                for (int j = 0; j < 1; j++)
                    str += " " + b[i,j];
                Console.WriteLine(str);
            }
            Console.WriteLine("\n");

            string info = "";
            double[,] x = new double[20, 1];//n*1
            bool isOK = Solve.Maqr(A, b, ref x, m, n, ref info);
            if(isOK)
            {
                Console.WriteLine("超定方程解算完毕！解如下:");
                for (int i = 0; i < n; i++)
                {
                    Console.WriteLine("x[" + (i + 1) + "]：" + (x[i, 0]));
                } 
            }
            else
            {
                Console.WriteLine("超定方程解算错误！错误信息:" + info);
            } 
        }
    }
}
