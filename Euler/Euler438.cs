using System;
using System.Collections.Generic;
using System.Collections.Concurrent;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace Euler
{
    class Program
    {

        static void Main(string[] args)
        {
            const int N = 7;
            double[,] matrix = new double[N + 1, N + 1];
            for(int i = 0; i <= N; i++)
            {
                for (int j = 0; j <= N; j++)
                {
                    matrix[i, j] = 1;
                    for (int k = 1; k <= j; k++) matrix[i, j] *= (i + 1);
                }
                matrix[i, N] *= -1; // The cooefficient on x^n is fixed at 1, so the requirement that f(x) ~ 0 requires negating this term
            }
            double modifier = 0;
            for (int i = N; i >= 0; i--)
            {
                matrix[i, N] += modifier;
                modifier *= -1;
            }
            Matrix polysystem = new Matrix(matrix);
            Simplex fullsimplex = polysystem.RegionalSimplex();
            fullsimplex.PrintSimplex();
            double[,] simplexvectors = new double[N, N];
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    simplexvectors[i, j] = fullsimplex.verticies[i].coord[j] - fullsimplex.verticies[N].coord[j];
                }
            }
            Matrix simplexmatrix = new Matrix(simplexvectors);
            Console.WriteLine("Approx volume = " + simplexmatrix.RecudedForm().ToString());
            ConcurrentDictionary<int,int> polyshash = new ConcurrentDictionary<int, int>();
            fullsimplex.LatticePoints(new int[] { }, polyshash);
            Console.WriteLine("Found " + polyshash.Count.ToString() + " Lattice Points");

            long sumint = 0;

            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"C:\Users\Demo\source\repos\Euler\Euler\Output.txt"))
            {
                foreach (KeyValuePair<int,int> enumerator in polyshash)
                {
                    sumint += enumerator.Value;
                }
                
                file.WriteLine("sumint = " + sumint);

            }

            //Console.WriteLine("passcount = " + polyshash.Count);
            Console.WriteLine("sumint = " + sumint);
            Console.Read();
            
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
        public int dimension;
        public List<Vertex> verticies;
        public bool[] covered;

        public Simplex(double[,] verts)
        {
            this.dimension = verts.GetLength(0) - 1;
            double[] coord = new double[dimension + 1];
            this.verticies = new List<Vertex>();
            for(int i = 0; i <= dimension; i++)
            {
                for(int j = 0; j < dimension; j++)
                {
                    coord[j] = verts[i, j];
                }
                this.verticies.Add(new Vertex(coord));
            }           
            this.covered = new bool[dimension + 1];
        }

        public Simplex(List<Vertex> verts)
        {
            this.dimension = verts.Count - 1;
            this.verticies = new List<Vertex>();
            foreach (Vertex vert in verts)
            {
                this.verticies.Add(vert.Clone());
            }
            this.covered = new bool[dimension + 1];
        }

        /*  Returns the "dimension"-1 simplex that is the intersection
         *  of this simplex with the plane with "dimention"th coordinate xn
         */
        public List<Simplex> Cull(double xn)
        {
            int dimind = this.dimension - 1;
            double epsilon = 0.00000000001;

            List<Vertex> subverts = new List<Vertex>();

            double xi1n, xi2n, t;
            for (int i1 = 0; i1 < dimension; i1++)
            {
                xi1n = this.verticies[i1].coord[dimind];
                for (int i2 = i1 + 1; i2 < dimension + 1; i2++)
                {
                    xi2n = this.verticies[i2].coord[dimind];
                    if ((xn - xi1n) * (xn - xi2n) < 0 && Math.Abs(xn - xi1n) > epsilon && Math.Abs(xn - xi2n) > epsilon) // Check if the "i1"th and "i2"th vertex lie on opposite sides of the plane
                    {
                        t = (xn - xi2n) / (xi1n - xi2n);
                        double[] subvert = new double[dimension];
                        for (int j = 0; j < dimind; j++)
                        {
                            subvert[j] = Math.Round(t * this.verticies[i1].coord[j] + (1 - t) * this.verticies[i2].coord[j]);
                        }
                        subverts.Add(new Vertex(subvert));
                    } else if (Math.Abs(xn - xi1n) < epsilon && i2 == dimension) // Each equality case only gets stored once
                    {
                        double[] subvert = new double[dimension];
                        for (int j = 0; j < dimind; j++)
                        {
                            subvert[j] = this.verticies[i1].coord[j];
                        }
                        subverts.Add(new Vertex(subvert));
                    } else if (Math.Abs(xn - xi2n) < epsilon && i1 == dimension - 1)
                    {
                        double[] subvert = new double[dimension];
                        for (int j = 0; j < dimind; j++)
                        {
                            subvert[j] = this.verticies[i2].coord[j];
                        }
                        subverts.Add(new Vertex(subvert));
                    }
                }
            }

            int vertcount = subverts.Count;

            List<Simplex> culled = new List<Simplex>();
 
            if (vertcount > dimension)
            {
                // Triangulate Convex Hull
                subverts.Sort();
                List<Vertex> seedsimplex = new List<Vertex>();
                for(int i = 0; i < dimension; i++)
                {
                    seedsimplex.Add(subverts[i]);
                }
                subverts.RemoveRange(0, dimension);

                culled.Add(new Simplex(seedsimplex));

                foreach (Vertex vert in subverts)
                {
                    vert.PointToHull(culled);
                }
            } else if (subverts.Count < dimension)
            {
                for (int i = vertcount; i < dimension; i++)
                {
                    // If too few verticies have been found, duplicate the last vertex
                    // No vertex will ever have a 0 coordinate (at least in this problem)
                    subverts.Add(subverts[vertcount - 1].Clone());
                }
                culled.Add(new Simplex(subverts));
            }
            else
            {
                culled.Add(new Simplex(subverts));
            }
            return culled;
        }

        double[] Extremes(int index)
        {
            double[] extremes = { double.MaxValue, double.MinValue };
            for (int i = 0; i < dimension + 1; i++)
            {
                if (this.verticies[i].coord[index] < extremes[0]) extremes[0] = this.verticies[i].coord[index];
                if (this.verticies[i].coord[index] > extremes[1]) extremes[1] = this.verticies[i].coord[index];
            }
            return extremes;
        }

        /* Recursively adds the coordinates of all integer points contained within the simplex
         * The argument "partial" is used to pass the sum of all high coordinate dimensions down
         * to the current dimension
         */
        public void LatticePoints(int[] partial, ConcurrentDictionary<int,int> coords)
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
                    Polynomial poly = new Polynomial(prepartial,true);
                    int hash = poly.GetHashCode();
                    if (coords.TryAdd(hash, 0))
                    {
                        bool pass = true;
                        for (int j = 1; j < poly.degree + 1; j++)
                        {
                            if (!poly.RootCheck(j, j + 1)) pass = false;
                        }
                        if (pass)
                        {
                            coords.TryUpdate(hash,poly.AbsCoeffSum(),0);
                        }
                    }
                }
            }
            else
            {
                List<Task> badidea = new List<Task>();
                for (int i = start; i <= end; i++)
                {
                    var threadi = i;
                    Task task = new Task(() =>
                    {
                        int[] threadpartial = (int[])prepartial.Clone();
                        threadpartial[0] = threadi;
                        List<Simplex> culled = this.Cull(threadi);
                        foreach (Simplex slice in culled)
                        {
                            slice.LatticePoints(threadpartial, coords);
                        }
                    });
                    badidea.Add(task);
                }
                foreach (Task task in badidea)
                {
                    task.Start();
                }
                Task.WaitAll(badidea.ToArray());
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
                    vertex += this.verticies[i].coord[j] + " , ";
                }
                vertex += this.verticies[i].coord[this.dimension - 1] + " )";
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
            int hash = this.coefficients.Length;
            for (int i = 0; i < this.coefficients.Length; i++)
            {
                hash = unchecked(hash * 314159 + this.coefficients[i]);
            }
            return hash;
        }
    }

    class Vertex : IComparable<Vertex>, IEquatable<Vertex>
    {
        public double[] coord;
        int dimension;

        public Vertex(double[] coords)
        {
            dimension = coords.Length;
            coord = new double[dimension];
            coords.CopyTo(coord, 0);
        }

        public int CompareTo(Vertex vert)
        {
            if (this.dimension != vert.dimension) return -1;
            int i = 0;
            while(i < dimension && vert.coord[i] == this.coord[i])
            {
                i++;
            }
            if (i == dimension)
            {
                return 0;
            }
            else if (this.coord[i] < vert.coord[i])
            {
                return -1;
            }
            else return 1;
        }

        public bool Equals(Vertex vert)
        {
            return (this.CompareTo(vert) == 0);
        }

        public override int GetHashCode()
        {
            return this.ToString().GetHashCode();
        }

        public Vertex Clone()
        {
            Vertex clone = new Vertex(this.coord);
            return clone;
        }

        public override string ToString()
        {
            string output = "( ";
            for (int i = 0; i < dimension - 1; i++)
            {
                output += this.coord[i].ToString() + " , ";
            }
            output += this.coord[dimension - 1].ToString() + " )";
            return output;
        }

        // Assumes everthing has the same dimension
        public void PointToHull(List<Simplex> hull)
        {
            List<Simplex> added = new List<Simplex>();
            foreach(Simplex simplex in hull)
            {
                for (int i0 = 0; i0 < simplex.verticies.Count; i0++)
                {
                    if (!simplex.covered[i0]) 
                    {
                        /* Check where this point lies on the opposite side of the face of simplex not containing
                         * simplex[i0] by calculating the signed volume of each simplex and comparing sign.
                         */
                        double[,] simplexsign = new double[simplex.dimension, simplex.dimension];
                        double[,] pointsign = new double[simplex.dimension, simplex.dimension];
                        int i = 0;
                        while (i < i0)
                        {
                            for (int j = 0; j < simplex.dimension; j++)
                            {
                                simplexsign[i, j] = simplex.verticies[i].coord[j] - simplex.verticies[i0].coord[j];
                                pointsign[i, j] = simplex.verticies[i].coord[j] - this.coord[j];
                            }
                            i++;
                        }
                        i++;
                        while (i < simplex.dimension + 1)
                        {
                            for (int j = 0; j < simplex.dimension; j++)
                            {
                                simplexsign[i-1, j] = simplex.verticies[i].coord[j] - simplex.verticies[i0].coord[j];
                                pointsign[i-1, j] = simplex.verticies[i].coord[j] - this.coord[j];
                            }
                            i++;
                        }
                        if ((new Matrix(simplexsign).RecudedForm())*(new Matrix(pointsign).RecudedForm()) <= 0)
                        {
                            List<Vertex> newsimplexverts = new List<Vertex>();
                            for (int j = 0; j < simplex.dimension + 1; j++)
                            {
                                if (j != i0) newsimplexverts.Add(simplex.verticies[j]);
                            }
                            newsimplexverts.Add(this);
                            Simplex newsimplex = new Simplex(newsimplexverts);
                            newsimplex.covered[simplex.dimension] = true;

                            // Check if the newly added simplex shares a d-2 face with another added simplex
                            // We then know that the d-1 face with that as a subface and 'this' as the last
                            // vertex is a common face between them.
                            foreach(Simplex preadded in added)
                            {
                                bool shared = true;
                                for (int j0 = 0; j0 < simplex.dimension; j0++)
                                {
                                    HashSet<int> indicies = new HashSet<int>();
                                    for (int j = 0; j < simplex.dimension; j++)
                                    {
                                        indicies.Add(j);
                                    }
                                    for (int j = 0; j < simplex.dimension; j++)
                                    {
                                        if (j != j0 && j != i0)
                                        {
                                            if (!indicies.Remove(preadded.verticies.IndexOf(newsimplex.verticies[j]))) // Eliminates the index of the found vertex from the list of indicies
                                            {
                                                shared = false;
                                                break;
                                            }
                                        }
                                    }
                                    if (shared)
                                    {
                                        newsimplex.covered[j0] = true;
                                        foreach (int index in indicies)
                                        {
                                            preadded.covered[index] = true;
                                        }
                                    }
                                }
                            }
                            added.Add(newsimplex);
                            simplex.covered[i0] = true;
                        }
                    }
                }
            }
            hull.AddRange(added);
        }
    }
}
