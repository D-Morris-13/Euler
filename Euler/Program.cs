using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Euler
{
    class Program
    {
        static void Main(string[] args)
        {
/*            double[] throughpoints = { 1, 8 ,27,64};
//           double[,] matrixin = new double[3, 4] { { 1, 2, 3, 4 }, { 4, 5, 6, 7 }, { 7, 8, 9, 10 } };
//            Matrix testmatrix = new Matrix(matrixin);
//            testmatrix.PrintMatrix();

//            testmatrix.RecudedForm();
//            testmatrix.PrintMatrix();
            Polynomial testpoly = new Polynomial(throughpoints);
            for (int n = 0; n < 10; n++)
            {
                Console.WriteLine("n = " + (double) n + " P(n) = " + testpoly.Evaluate((double)n));
            }
            Console.WriteLine();

            double[] testthrough = testpoly.Throughpoints(5);
            for (int n = 0; n < testthrough.Length; n++)
            {
                Console.WriteLine(testthrough[n]);
            }

            double[] subpoints;
            for (int n = 1; n < throughpoints.Length; n++)
            {
                subpoints = throughpoints.Take(n).ToArray();
                testpoly = new Polynomial(subpoints);
                Console.WriteLine("n = " + n + " Eval = " + testpoly.Evaluate((double)n+1));
            }*/

            double[] coeffs = { 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1 };

            Polynomial qpoly = new Polynomial(coeffs,true);

            double eval;

            for (int i = 1; i < 11; i++)
            {
                eval = qpoly.Evaluate((double)i);
                Console.WriteLine(eval);
            }

            int n = 0;
            double[] throughpoints;
            Polynomial subpoly;
            double sum = 0;
       
            for (n = 1; n < 11; n++)
            {
                throughpoints = qpoly.Throughpoints(n);
                subpoly = new Polynomial(throughpoints);
                eval = subpoly.Evaluate((double)n + 1);
                Console.WriteLine("Eval = " + eval);
                sum += eval;
            }

            Console.WriteLine();
            Console.WriteLine(sum);
        }

    }

    class Polynomial
    {
        double[] coefficients;

        public Polynomial(double[] throughpoints)
        {
            int length = throughpoints.Length;
            double[,] coefficientarray = new double[length, length+1];
            for (int i = 0; i < length; i++)
            {
                coefficientarray[i, 0] = 1;
                coefficientarray[i, length] = throughpoints[i];
                if (length > 1)
                {
                    for (int j = 1; j < length; j++)
                    {
                        coefficientarray[i, j] = coefficientarray[i, j - 1] * (i+1);
                    }
                }
            }
            Matrix coefficientmatrix = new Matrix(coefficientarray);
            coefficientmatrix.RecudedForm();
//            coefficientmatrix.PrintMatrix();

            coefficients = coefficientmatrix.GetColumn(length);
            Console.WriteLine("n = " + throughpoints.Length);
            for (int i = 0; i < coefficients.Length; i++) Console.Write(" a_" + i + " = " + coefficients[i]);
            Console.WriteLine();

        }

        public Polynomial(double[] coeffs,bool blank)
        {
            coefficients = (double[])coeffs.Clone();
        }

        public double Evaluate(double point)
        {
            double[] monomials = (double[]) coefficients.Clone();
            for (int i = 0; i < coefficients.Length; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    monomials[i] = monomials[i] * point;
                }
            }
            double eval = 0;
            for (int i = 0; i < coefficients.Length; i++)
            {
                eval += monomials[i];
//                Console.WriteLine("i = " + i + " monomials[i] = " + monomials[i] + " eval = " + eval);
            }
            return eval;
        }

        public double[] Throughpoints(int n)
        {
            double[] throughpoints = new double[n];
            for (int i = 1; i <= n; i++)
            {
                throughpoints[i-1] = this.Evaluate((double)i);
            }
            return throughpoints;
        }


    }

    class Matrix
    {
        double[,] matrix;
        int ilength;
        int jlength;

        public Matrix(double[,] input)
        {
            matrix = input;
            ilength = input.GetLength(0);
            jlength = input.GetLength(1);

        }

        public void RecudedForm()
        {
//            this.PrintMatrix();
            for (int i = 0; i < ilength && i < jlength; i++)
            {
                if (Math.Abs(matrix[i, i]) < 0.00001)
                {
                    int i2;
                    for (i2 = i+1; i2 < ilength; i2++)
                    {
                        if (Math.Abs(matrix[i2,i]) > 0.00001)
                        {
//                            Console.WriteLine("Exchanged rows " + i + " and " + i2);
                            this.Interchange(i, i2);
//                            this.PrintMatrix();
                            break;
                        }
                    }
                    if (i2 == ilength) continue;
                }
//                Console.WriteLine("Scaled row " + i + " by " + 1 / matrix[i, i]);
                this.Scale(i, 1/matrix[i, i]);
//                this.PrintMatrix();
                for (int i2 = 0; i2 < ilength; i2++)
                {
                    if (i2 != i)
                    {
//                        Console.WriteLine("Offset row " + i2 + " by " + -matrix[i2, i] + " times row " + i);
                        this.Offset(i, i2, -matrix[i2, i]);
 //                       this.PrintMatrix();
                    }
                }
            }
        }

        public void PrintMatrix()
        {
            String output = "";
            for (int i = 0; i < ilength; i++)
            {
                if (i == 0) output += "/ ";
                else if (i == ilength - 1) output += "\\ ";
                else output += "| ";
                for (int j = 0; j < jlength; j++)
                {
                    output += String.Format("{0,5} ", matrix[i, j]);
                }
                if (i == 0) output += "\\\n";
                else if (i == ilength - 1) output += "/";
                else output += "|\n";
            }
            output += "\n";
            Console.WriteLine(output);
        }

        void Interchange(int i1, int i2)
        {
            if (i1 >= ilength || i2 >= ilength)
            {

            } else
            {
                double[] row1 = this.GetRow(i1);
                double[] row2 = this.GetRow(i2);
                this.SetRow(i1, row2);
                this.SetRow(i2, row1);
            }
        }

        void Scale(int i, double scale)
        {
            double[] row = this.GetRow(i);
            for (int j = 0; j < jlength; j++)
            {
                row[j] = row[j] * scale;
            }
            this.SetRow(i, row);
        }

        void Offset(int i1, int i2, double scale)
        {
            double[] row1 = this.GetRow(i1);
            double[] row2 = this.GetRow(i2);
            for (int j = 0; j < jlength; j++)
            {
                row2[j] = row2[j] +  row1[j] * scale;
            }
            this.SetRow(i2, row2);
        }

        public double[] GetRow(int i)
        {
            double[] row = new double[jlength];
            for (int j = 0; j < jlength; j++)
            {
                row[j] = matrix[i, j];
            }
            return row;
        }

        public void SetRow(int i, double[] row)
        {
            for (int j = 0; j < jlength; j++)
            {
                matrix[i, j] = row[j];
            }
        }

        public double[] GetColumn(int j)
        {
            double[] column = new double[ilength];
            for (int i = 0; i < ilength; i++)
            {
                column[i] = matrix[i, j];
            }
            return column;
        }

        public void SetColumn(int j, double[] column)
        {
            for (int i = 0; i < ilength; i++)
            {
                matrix[i, j] = column[i];
            }
        }
    }
}
