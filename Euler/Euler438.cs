﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace Euler
{
    public static class Globals
    {
        public const int N = 7;
        public static readonly int[] PRIMELIST = { 2, 3, 5, 7, 11, 13, 17, 19 , 23, 29, 31, 37, 41, 43, 47, 53 , 59, 61, 67, 71};
    }

    class Program
    {

        static void Main(string[] args)
        {
            double[,] matrix = new double[Globals.N + 1, Globals.N + 1];
            for(int i = 0; i <= Globals.N; i++)
            {
                for (int j = 0; j <= Globals.N; j++)
                {
                    matrix[i, j] = 1;
                    for (int k = 1; k <= j; k++) matrix[i, j] *= (i + 1);
                }
                matrix[i, Globals.N] *= -1; // The cooefficient on x^n is fixed at 1, so the requirement that f(x) ~ 0 requires negating this term
            }
            double modifier = 0;
            for (int i = Globals.N; i >= 0; i--)
            {
                matrix[i, Globals.N] += modifier;
                modifier *= -1;
            }
            Matrix polysystem = new Matrix(matrix);
            Simplex fullsimplex = polysystem.RegionalSimplex();
            fullsimplex.PrintSimplex();
            HashSet<Polynomial> polys = new HashSet<Polynomial>();
            fullsimplex.LatticePoints(new int[] { }, polys);
            Console.WriteLine("Found Lattice Points");

            long sumint = 0;
            int passcount = 0;

            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"C:\Users\Demo\source\repos\Euler\Euler\Output.txt"))
            {
                foreach (Polynomial poly in polys)
                {
                    bool pass = true;
                    for (int i = 1; i < Globals.N + 1; i++)
                    {
                        if (!poly.RootCheck(i, i + 1)) pass = false;
                    }
                    if (pass)
                    {
                        string output = "Found passing polynomial: " + poly.ToString();
                        //Console.WriteLine(output);
                        file.WriteLine(output);
                        sumint += poly.AbsCoeffSum();
                        passcount++;
                    }
                    else
                    {
                        string output = "Found failing polynomial: " + poly.ToString();
                        //Console.WriteLine(output);
                        //file.WriteLine(output);
                    }

                }
            }

            Console.WriteLine("passcount = " + passcount);
            Console.WriteLine("sumint = " + sumint);
            
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

        // Perform elementary row operations to reduce the matrix
        // Returns the determinant if square, and 0 if not.
        public double RecudedForm()
        {
            double det = 1;
            for (int i = 0; i < ilength && i < jlength; i++)
            {
                // Exchange the ith row with the first row with a non-zero ith term
                if (Math.Abs(matrix[i, i]) < 0.00001)
                {
                    int i2;
                    for (i2 = i + 1; i2 < ilength; i2++)
                    {
                        if (Math.Abs(matrix[i2, i]) > 0.00001)
                        {
                            this.Interchange(i, i2);
                            break;
                        }
                    }
                    if (i2 == ilength) continue;
                }
                // Scale the ith row so that a_{i,i} = 1
                det *= matrix[i, i];
                this.Scale(i, 1 / matrix[i, i]);
                // Add a multiple of the ith row to each other row so that a_{i',i} = 0 for all i' \neq i
                for (int i2 = 0; i2 < ilength; i2++)
                {
                    if (i2 != i)
                    {
                        this.Offset(i, i2, -matrix[i2, i]);
                    }
                }
            }
            if (ilength == jlength)
            {
                return det;
            }
            else
            {
                return 0;
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

        // Perform elementary row operations on the matrix

        // Exchange rows i1 and i2
        void Interchange(int i1, int i2)
        {
            if (i1 >= ilength || i2 >= ilength)
            {

            }
            else
            {
                double[] row1 = this.GetRow(i1);
                double[] row2 = this.GetRow(i2);
                this.SetRow(i1, row2);
                this.SetRow(i2, row1);
            }
        }

        // Multiply row i by scale
        void Scale(int i, double scale)
        {
            double[] row = this.GetRow(i);
            for (int j = 0; j < jlength; j++)
            {
                row[j] = row[j] * scale;
            }
            this.SetRow(i, row);
        }

        // Add scale times row i1 to row i2
        void Offset(int i1, int i2, double scale)
        {
            double[] row1 = this.GetRow(i1);
            double[] row2 = this.GetRow(i2);
            for (int j = 0; j < jlength; j++)
            {
                row2[j] = row2[j] + row1[j] * scale;
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

        public int GetiLength()
        {
            return ilength;
        }

        public int GetjLength()
        {
            return jlength;
        }

        /* Given a square matrix of the form
         * 
         * /                                             \
         * |  a_{0,0}    a_{0,1} ... a_{0,d-1}    b_{0}  |
         * |  a_{1,0}    a_{1,1} ... a_{1,d-1}    b_{1}  |
         * |    .           .           .          .     | = A|B
         * |    .           .           .          .     |
         * |  a_{d,0}    a_{d,1} ... a_{d,d-1}    b_{d}  |
         * \                                             /
         * 
         * returns the simplex bounded by the hyperplanes
         * A_{!= i} v = B_{!= i}
         * (i.e. deleting the ith row from the above system of equations)
         *
         */
        public Simplex RegionalSimplex()
        {
            int dimension = this.jlength - 1;
            double[,] verts = new double[dimension + 1, dimension];

            double[,] subequations = new double[dimension, dimension + 1];

            Matrix submatrix = new Matrix(subequations);

            for (int noti = 0; noti < dimension + 1; noti++)
            {
                for(int i = 0; i < dimension + 1; i++)
                {
                    // Populate the submatrix with the entries of the full matrix less
                    // the deleted row
                    if(i < noti)
                    {
                        submatrix.SetRow(i, this.GetRow(i));
                    } else if (i > noti)
                    {
                        submatrix.SetRow(i - 1, this.GetRow(i));
                    }
                }
                // Solve the system of equations defined by the submatrix and set
                // the notith vertex of verts to be the solution
                submatrix.RecudedForm();
                double[] irow = submatrix.GetColumn(dimension);
                for (int j = 0; j < dimension; j++) verts[noti, j] = irow[j];
            }

            Simplex solution = new Simplex(verts);

            return solution;
        }
    }

    class Simplex
    /*  An array to define a simplex in "dimension"-dimensional-space.
     *  "dimension"+1 verticies, indexed by the first index
     *  Each has "dimension" coodinates, indexed by the second index
     */
    {
        int dimension;
        double[,] verticies;
        bool[] covered;

        public Simplex(double[,] verts)
        {
            this.verticies = (double[,])verts.Clone();
            this.dimension = verts.GetLength(0) - 1;
            this.covered = new bool[dimension + 1];
        }

        /*  Returns the "dimension"-1 simplex that is the intersection
         *  of this simplex with the plane with "dimention"th coordinate xn
         */
        public Simplex[] Cull(double xn)
        {
            int dimind = this.dimension - 1;
            double epsilon = 0.00000000001;

            double[,] subverts = new double[(this.dimension*this.dimension+1)/2, this.dimension - 1];

            int ind = 0;
            double xi1n, xi2n, t;
            for (int i1 = 0; i1 < dimension; i1++)
            {
                xi1n = this.verticies[i1, dimind];
                for (int i2 = i1 + 1; i2 < dimension + 1; i2++)
                {
                    xi2n = this.verticies[i2, dimind];
                    if ((xn - xi1n) * (xn - xi2n) < 0 && Math.Abs(xn - xi1n) > epsilon && Math.Abs(xn - xi2n) > epsilon) // Check if the "i1"th and "i2"th vertex lie on opposite sides of the plane
                    {
                        t = (xn - xi2n) / (xi1n - xi2n);
                        for (int j = 0; j < dimind; j++)
                        {
                            subverts[ind, j] = Math.Round(t * this.verticies[i1, j] + (1 - t) * this.verticies[i2, j]);
                        }
                        ind++;
                    } else if (Math.Abs(xn - xi1n) < epsilon && i2 == dimension) // Each equality case only gets stored once
                    {
                        for (int j = 0; j < dimind; j++)
                        {
                            subverts[ind, j] = this.verticies[i1, j];
                        }
                        ind++;
                    } else if (Math.Abs(xn - xi2n) < epsilon && i1 == dimension - 1)
                    {
                        for (int j = 0; j < dimind; j++)
                        {
                            subverts[ind, j] = this.verticies[i2, j];
                        }
                        ind++;
                    }
                }
            }

            int vertcount = 0;
            for (int i = 0; i < subverts.GetLength(0); i++)
            {
                if (subverts[i, 0] != 0) vertcount++;  // No vertex will ever have a 0 coordinate (at least in this problem)
            }

            if (vertcount - dimension > 2)
            {
                Console.WriteLine("Many things!");
            }


            Simplex[] culled;
            double[,] cullverts = new double[dimind + 1, dimind];

            if (vertcount > dimension)
            {
                // Triangulate Convex Hull
                int culledsize = 1;
                for (int i = 1; i <= vertcount; i++) culledsize *= i;
                for (int i = 1; i <= vertcount - dimension; i++) culledsize /= i;
                for (int i = 1; i <= dimension; i++) culledsize /= i;
                culled = new Simplex[culledsize];
                int[] index = new int[dimension];
                bool newvert = true;
                for (int i = 1; i < dimension; i++)
                {
                    while (index[i] <= index[i - 1]) index[i]++;
                }
                int i0 = 0;
                while (newvert)
                {
                    newvert = false;
                    for (int i = 0; i < dimension; i++)
                    {
                        for (int j = 0; j < dimind; j++)
                        {
                            cullverts[i, j] = subverts[index[i], j];
                        }
                    }
                    culled[i0] = new Simplex(cullverts);
                    i0++;
                    for (int i = dimind; i >= 0; i--)
                    {
                        if (index[i] < vertcount - 1 && (i == dimind || index[i] < index[i + 1] - 1))
                        {
                            index[i]++;
                            newvert = true;
                            break;
                        }
                        else if (i > 0 && index[i] > index[i - 1] + 1)
                        {
                            index[i - 1]++;
                            for (int j = i; j <= dimind; j++)
                            {
                                index[j] = index[j - 1] + 1;
                            }
                            newvert = true;
                            break;
                        }
                    }
                }
            } else if (vertcount < dimension)
            {
                culled = new Simplex[1];
                for (int i = 0; i < dimension; i++)
                {
                    // If too few verticies have been found, duplicate the last vertex
                    // No vertex will ever have a 0 coordinate (at least in this problem)
                    int fill = 0;
                    while (Math.Abs(subverts[i - fill, 0]) < epsilon)
                    {
                        fill++;
                    }
                    for (int j = 0; j < dimind; j++)
                    {
                        cullverts[i, j] = subverts[i - fill, j];
                    }
                }
                culled[0] = new Simplex(cullverts);
            }
            else
            {
                culled = new Simplex[1];
                for (int i = 0; i < dimension; i++)
                {
                    for (int j = 0; j < dimind; j++)
                    {
                        cullverts[i, j] = subverts[i, j];
                    }
                }
                culled[0] = new Simplex(cullverts);
            }
            return culled;
        }

        double[] Extremes(int index)
        {
            double[] extremes = { double.MaxValue, double.MinValue };
            for (int i = 0; i < dimension + 1; i++)
            {
                if (this.verticies[i, index] < extremes[0]) extremes[0] = this.verticies[i, index];
                if (this.verticies[i, index] > extremes[1]) extremes[1] = this.verticies[i, index];
            }
            return extremes;
        }

        /* Recursively adds the coordinates of all integer points contained within the simplex
         * The argument "partial" is used to pass the sum of all high coordinate dimensions down
         * to the current dimension
         */
        public void LatticePoints(int[] partial, HashSet<Polynomial> coords)
        {
            double[] extremes = this.Extremes(this.dimension - 1);
            int start = (int)Math.Ceiling(extremes[0]);
            int end = (int)Math.Floor(extremes[1]);
            int[] prepartial = new int[partial.Length + 1];
            partial.CopyTo(prepartial, 1);
            if (this.dimension == 1)
            {

                for (int i = start; i <= end; i++)
                {
                    prepartial[0] = i;
                    coords.Add(new Polynomial(prepartial,true));
                }
            }
            else
            {
                for (int i = start; i <= end; i++)
                {
                    Simplex[] culled = this.Cull(i);
                    prepartial[0] = i;
                    for (int j = 0; j < culled.Length; j++)
                    {
                        culled[j].LatticePoints(prepartial,coords);
                    }
                }
            }
        }

        public void PrintSimplex()
        {
            Console.WriteLine("The simplex has the following verticies: ");
            string vertex = "";
            for (int i = 0; i <= this.dimension; i++)
            {
                vertex = "( ";
                for (int j = 0; j < this.dimension- 1; j++)
                {
                    vertex += this.verticies[i, j] + " , ";
                }
                vertex += this.verticies[i, this.dimension - 1] + " )";
                Console.WriteLine(vertex);
            }
        }

    }

    class Polynomial : IEquatable<Polynomial>
    {
        public int[] coefficients { get; }
        public int degree { get; }

        public Polynomial(int[] coeffs, bool monic)
        {
            degree = coeffs.Length;
            if (!monic) degree--;
            coefficients = new int[degree + 1];
            coeffs.CopyTo(coefficients,0);
            if (monic) coefficients[degree] = 1;
        }

        public int Evaluate(int m)
        {
            int output = this.coefficients[degree];
            for(int i = this.degree - 1; i >= 0; i--)
            {
                output *= m;
                output += this.coefficients[i];
            }
            return output;
        }

        public int DerivEval(int m)
        {
            int output = 0;
            int monoderiv = 0;
            for (int i = 1; i <= this.degree; i++)
            {
                monoderiv = i * this.coefficients[i];
                for (int j = 1; j < i; j++)
                {
                    monoderiv *= m;
                }
                output += monoderiv;
            }
            return output;
        }

        public bool RootCheck(int from, int to)
            /* Checks sufficient conditions for this to have exactly one root in the interval [from,to).
             * It is possible for the polynomial to have multiple roots in the interval and this to return
             * false.
             */
        {
            int polyfrom = this.Evaluate(from);
            if (polyfrom == 0) return true;

            int polyto = this.Evaluate(to);
            if (polyfrom * polyto != 0) return (polyfrom*polyto < 0);

            int derivto = this.DerivEval(to);
            return ( polyfrom * derivto > 0);
        }

        public int AbsCoeffSum()
        {
            int sum = 0;
            for (int i = 0; i < degree; i++)
            {
                sum += Math.Abs(this.coefficients[i]);
            }
            return sum;
        }

        override public string ToString()
        {
            string output = "";
            for (int i = 0; i <= this.degree; i++)
            {
                output += this.coefficients[i].ToString();
                if (i > 0)
                {
                    output += " x";
                    if (i > 1)
                    {
                        output += "^" + i.ToString();
                    }
                }
                output += " ";
                if (i < this.degree && this.coefficients[i + 1] > 0) output += "+";
            }
            return output;
        }
        public bool Equals(Polynomial poly)
        {
            if (this.degree != poly.degree) return false;
            return (this.GetHashCode() == poly.GetHashCode());
        }

        public override int GetHashCode()
        {
            string polystring = this.ToString();
            return polystring.GetHashCode();
        }
    }

    class Vertex
    {
        int[] coord;
        int dimension;

        public Vertex(int[] coords)
        {
            dimension = coords.Length + 1;
            coord = new int[dimension - 1];
            coords.CopyTo(coord, 0);
        }

        // Assumes everthing has the same dimension
        public void PointToHull(List<Simplex> hull, Vertex lastpoint)
        {

        }
    }
}