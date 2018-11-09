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
            ConcurrentDictionary<int, int> polyshash = new ConcurrentDictionary<int, int>();
            Hull fullhull = new Hull(fullsimplex);
            long sumint = fullhull.LatticePoints(new int[] { });

            Console.WriteLine("sumint = " + sumint);
            Console.Read();
            
        }       
    }

    class Matrix
    {
        double[,] matrix;
        int ilength;
        int jlength;
        const double epsilon = 0.000001;

        public Matrix(double[,] input)
        {
            matrix = input;
            ilength = input.GetLength(0);
            jlength = input.GetLength(1);

        }

        // Perform elementary row operations to reduce the matrix
        // Returns the determinant if square, and 0 if not.
        public double ReducedForm()
        {
            double det = 1;
            for (int i = 0; i < ilength && i < jlength; i++)
            {
                // Exchange the ith row with the first row with a non-zero ith term
                if (Math.Abs(matrix[i, i]) < epsilon)
                {
                    int i2;
                    for (i2 = i + 1; i2 < ilength; i2++)
                    {
                        if (Math.Abs(matrix[i2, i]) > epsilon)
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

        // Performs elementary row operations to reduce the square matrix
        // to an upper triangular one and finds the determinant whilst doing so.
        public double Determinant()
        {
            if (ilength != jlength) return 0;
            double det = 1;
            for (int i = 0; i < ilength; i++)
            {
                // Exchange the ith row with the first row with a non-zero ith term
                if (Math.Abs(matrix[i, i]) < epsilon)
                {
                    int i2;
                    for (i2 = i + 1; i2 < ilength; i2++)
                    {
                        if (Math.Abs(matrix[i2, i]) > epsilon)
                        {
                            this.Interchange(i, i2);
                            break;
                        }
                    }
                    if (i2 == ilength) return 0; // Reduced matrix has a 0 column
                }
                // Scale the ith row so that a_{i,i} = 1
                det *= matrix[i, i];
                this.Scale(i, 1 / matrix[i, i]);
                // Add a multiple of the ith row to each other row so that a_{i',i} = 0 for all i' \neq i
                for (int i2 = i + 1; i2 < ilength; i2++)
                {
                    this.Offset(i, i2, -matrix[i2, i]);
                }
            }
            return det;
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
                double temp;
                for (int j = 0; j < jlength; j++)
                {
                    temp = this.matrix[i1, j];
                    this.matrix[i1, j] = this.matrix[i2, j];
                    this.matrix[i2, j] = temp;
                }
            }
        }

        // Multiply row i by scale
        void Scale(int i, double scale)
        {
            for (int j = 0; j < jlength; j++)
            {
                this.matrix[i, j] *= scale;
            }
        }

        // Add scale times row i1 to row i2
        void Offset(int i1, int i2, double scale)
        {
            for (int j = 0; j < jlength; j++)
            {
                this.matrix[i2, j] += this.matrix[i1, j] * scale;
            }
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
                submatrix.ReducedForm();
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
            double epsilon = 0.000001;

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

    class Hull
    /*  Defines the complex hull of a collection of points in
     *  "dimension"-dimensional space.
     *  There is no assumption on the number of points; a Hull
     *  may lie in a codimension > 0 subspace and/or not be a
     *  simplex.
     */
    {
        public List<Vertex> points;
        public int dimension;
        public bool reduced;
        static double epsilon = 0.000000000001;
        //static Vector nullvector = new Vector();

        public Hull (List<Vertex> pointlist, int dim, bool prereduced, bool toreduce = false)
        {
            this.dimension = dim;
            this.points = pointlist;
            if (prereduced)
            {
                this.reduced = true;
            }
            else
            {
                this.reduced = false;
                if (toreduce)
                {
                    this.Reduce();
                }
            }
        }

        public Hull (Simplex simplex)
        {
            this.dimension = simplex.dimension;
            this.points = simplex.verticies;
            this.reduced = true;
        }

        double[] Extremes(int index)
        {
            double[] extremes = { double.MaxValue, double.MinValue };
            foreach (Vertex point in this.points)
            {
                if (point.coord[index] < extremes[0]) extremes[0] = point.coord[index];
                if (point.coord[index] > extremes[1]) extremes[1] = point.coord[index];
            }
            extremes[0] = Math.Floor(extremes[0]);
            extremes[1] = Math.Ceiling(extremes[1]);
            return extremes;
        }

        // Removes verticies from this.points until a minimal set
        // with the same convex hull is found.
        public void Reduce()
        {
            // Deduplicate
            for (int i = 0; i < this.points.Count; i++)
            {
                for (int j = i + 1; j < this.points.Count; j++)
                {
                    if (this.points[i].Equals(this.points[j]))
                    {
                        this.points.RemoveAt(j);
                        j--;
                    }
                }
            }

            if (this.points.Count <= 1)
            {
                this.reduced = true;
                return;
            }

            if (this.dimension == 1 && this.points[0].coord[0] == 3584)
            {
                Console.WriteLine("It begins again");
            }
            // Finds a starting point with minimal 0th coordinate.
            // Such a point must lie on the boundary of the hull.
            double x0min = int.MaxValue;
            int index = 0;

            for(int i = 0; i < this.points.Count; i ++)
            {
                if (this.points[i].coord[0] < x0min)
                {
                    index = i;
                    x0min = this.points[i].coord[0];
                }
            }

            List<Vector> span = new List<Vector> { };
            
            List<Vertex> reducedpoints = new List<Vertex> { this.points[index] };
            this.points.RemoveAt(index);

            bool addpoint = true;
            double maxdistance = epsilon;

            // Finds the verticies of a simplex of maximal dimension contained
            // in the hull.
            // Iteratively defines the kth vertex to be the one with the
            // maximal perpendicular distance to the space containing the
            // first k-1.
            while (addpoint)
            {
                List<Vector> orthvectors = new List<Vector> { };
                addpoint = false;
                maxdistance = epsilon;
                for (int i = 0; i < this.points.Count; i++)
                {
                    Vector orthvector = this.points[i] - reducedpoints[0];
                    foreach (Vector vector in span)
                    {
                        orthvector -= ((orthvector.Dot(vector)) / (vector.Dot(vector))) * vector;
                    }
                    orthvectors.Add(orthvector);
                    if (orthvector.Dot(orthvector) > maxdistance)
                    {
                        index = i;
                        maxdistance = orthvector.Dot(orthvector);
                        addpoint = true;
                    }
                }
                if (addpoint)
                {
                    span.Add(/*(1 / Math.Sqrt(orthvectors[index].Dot(orthvectors[index]))) */ orthvectors[index]);
                    reducedpoints.Add(this.points[index]);
                    this.points.RemoveAt(index);
                }
            }

            // Finds the faces of the simplex and the normal to
            // each face.
            List<List<Vertex>> faces = new List<List<Vertex>>();
            List<Vector> normals = new List<Vector>();
            

            for (int i = 0; i < reducedpoints.Count; i++)
            {
                List<Vertex> face = new List<Vertex>();
                for (int j = 0; j < reducedpoints.Count; j++)
                {
                    if (j != i)
                    {
                        face.Add(reducedpoints[j]);
                    }
                }
                if (face.Count > 1)
                {
                    for (int j = 1; j < face.Count; j++)
                    {
                        (face[j] - face[0]).AddToOrthList(span);
                    }
                }

                Vector normal;
                if (face.Count > 0)
                {
                    normal = reducedpoints[i] - face[0];
                    for (int j = 0; j < span.Count; j++)
                    {
                        normal -= normal.Dot(span[j]) * span[j];
                    }
                    normals.Add(normal);
                }
                else
                {
                    normals.Add(new Vector());
                }
                faces.Add(face);
                span.Clear();
            }

            // Discards all points inside the hull of the found points.
            // For each face, add the point of greates perpendicular distance
            // to that face to the hull, and then generate all the new faces
            // between that new point and the edges of the original face.
            List<double> maxdistances = new List<double>();
            List<Vertex> outside = new List<Vertex>();
            Dictionary<int, int> maxdistindicies = new Dictionary<int, int>();
            while (this.points.Count > 0)
            {
                maxdistances.Clear();
                outside.Clear();
                foreach(Vector normal in normals)
                {
                    maxdistances.Add(epsilon);
                }
                int ind = 0;
                foreach (Vertex point in this.points)
                {
                    for (int i = 0; i < normals.Count; i++)
                    {
                        double dist = 0;
                        if (faces[i].Count > 0)
                        {
                            dist = -(point - faces[i][0]).Dot(normals[i]);

                            if (dist > epsilon)
                            {
                                outside.Add(point);
                                if (dist > maxdistances[i])
                                {
                                    if (maxdistindicies.ContainsKey(i))
                                    {
                                        maxdistindicies[i] = ind;
                                    }
                                    else
                                    {
                                        maxdistindicies.Add(i, ind);
                                    }
                                    maxdistances[i] = dist;
                                }
                                break;
                            }
                        }
                    }
                    ind++;
                }

                List<List<Vertex>> newfaces = new List<List<Vertex>>();
                List<Vector> newnormals = new List<Vector>();
                for (int i = 0; i < faces.Count; i ++)
                {
                    if (maxdistindicies.ContainsKey(i))
                    {
                        reducedpoints.Add(this.points[maxdistindicies[i]]);
                        faces[i].Add(this.points[maxdistindicies[i]]);
                        for (int j = 0; j < faces[i].Count - 1; j++)
                        {
                            List<Vertex> face = new List<Vertex>();
                            for (int k = 1; k < faces[i].Count; k++)
                            {
                                if (k != j)
                                {
                                    face.Add(faces[i][k]);
                                }
                            }
                            for (int k = 1; k < face.Count; k++)
                            {
                                (face[k] - face[0]).AddToOrthList(span);
                            }
                            Vector normal = faces[i][j] - face[0];
                            for (int k = 0; k < span.Count; k++)
                            {
                                normal -= normal.Dot(span[k]) * span[k];
                            }
                            newfaces.Add(face);
                            newnormals.Add(normal);
                            span.Clear();
                        }
                    }
                    else
                    {
                        newfaces.Add(faces[i]);
                        newnormals.Add(normals[i]);
                    }
                }

                maxdistindicies.Clear();
                faces = newfaces;
                normals = newnormals;
                this.points = outside;
            }


            this.points = reducedpoints;
            this.reduced = true;
        }

        public Hull Cull(double xn)
        {
            int dimind = this.dimension - 1;

            List<Vertex> subverts = new List<Vertex>();

            double xi1n, xi2n, t;
            for (int i1 = 0; i1 < this.points.Count; i1++)
            {
                xi1n = this.points[i1].coord[dimind];
                for (int i2 = i1; i2 < this.points.Count; i2++)
                {
                    xi2n = this.points[i2].coord[dimind];
                    if ((xn - xi1n) * (xn - xi2n) < 0 && Math.Abs(xn - xi1n) > epsilon && Math.Abs(xn - xi2n) > epsilon) // Check if the "i1"th and "i2"th vertex lie on opposite sides of the plane
                    {
                        t = (xn - xi2n) / (xi1n - xi2n);
                        double[] subvert = new double[dimind];
                        for (int j = 0; j < dimind; j++)
                        {
                            subvert[j] = t * this.points[i1].coord[j] + (1 - t) * this.points[i2].coord[j];
                        }
                        subverts.Add(new Vertex(subvert));
                    }
                    else if (Math.Abs(xn - xi1n) < epsilon && i2 == this.points.Count - 1) // Each equality case only gets stored once
                    {
                        double[] subvert = new double[dimind];
                        for (int j = 0; j < dimind; j++)
                        {
                            subvert[j] = this.points[i1].coord[j];
                        }
                        subverts.Add(new Vertex(subvert));
                    }
                }
            }

            return new Hull(subverts, dimind, false, true);
        }

        public long LatticePoints(int[] partial)
        {
            long sum = 0;
            const bool multithread = false;

            if (this.points.Count == 0) return 0;

            double[] extremes = this.Extremes(this.dimension - 1);
            int start = (int)Math.Ceiling(extremes[0]);
            int end = (int)Math.Floor(extremes[1]);
            int[] prepartial = new int[partial.Length + 1];
            partial.CopyTo(prepartial, 1);
            if (this.dimension == 1)
            {
                int partsum = 0;
                for (int i = 0; i < partial.Length; i++)
                {
                    partsum += Math.Abs(partial[i]);
                }
                for (int i = start; i <= end; i++)
                {
                    prepartial[0] = i;
                    Polynomial poly = new Polynomial(prepartial, true);
                    int hash = poly.GetHashCode();
                    bool pass = true;
                    for (int j = 1; j < poly.degree + 1; j++)
                    {
                        if (!poly.RootCheck(j, j + 1)) pass = false;
                    }
                    if (pass)
                    {
                        sum += partsum + Math.Abs(i);
                    }
                    
                }
            }
            else
            {
                if (multithread)
                {
                    List<Task> tasklist = new List<Task>();
                    ConcurrentBag<long> partialsums = new ConcurrentBag<long>();
                    for (int i = start; i <= end; i++)
                    {
                        var threadi = i;
                        Task task = new Task(() =>
                        {
                            int[] threadpartial = (int[])prepartial.Clone();
                            threadpartial[0] = threadi;
                            Hull culled = this.Cull(threadi);
                            partialsums.Add(culled.LatticePoints(threadpartial));
                        });
                        tasklist.Add(task);
                    }
                    foreach (Task task in tasklist)
                    {
                        task.Start();
                    }
                    Task.WaitAll(tasklist.ToArray());
                    foreach (long partialsum in partialsums)
                    {
                        sum += partialsum;
                    }
                }
                else
                {
                    for (int i = start; i <= end; i++)
                    {
                        prepartial[0] = i;
                        Hull culled = this.Cull(i);
                        sum += culled.LatticePoints(prepartial);
                    }
                }
            }
            return sum;
        }
    }
       

    class Polynomial : IEquatable<Polynomial>
    {
        public int[] coefficients { get; }
        public int degree { get; }
        const double epsilon = 0.000001;

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
        const double epsilon = 0.000001;

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
            while(i < dimension && Math.Abs(vert.coord[i] - this.coord[i]) < epsilon)
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

        public static Vector operator- (Vertex vertex1, Vertex vertex2)
        {
            int dim = vertex1.dimension;
            double[] comps = new double[dim];
            for (int i = 0; i < dim; i++)
            {
                comps[i] = vertex1.coord[i] - vertex2.coord[i];
            }
            return new Vector(comps);
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
                        if ((new Matrix(simplexsign).Determinant())*(new Matrix(pointsign).Determinant()) <= 0)
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

    class Vector
    {
        public double[] components;
        public int dimension;
        const double epsilon = 0.000001;

        public Vector (double[] comps)
        {
            this.components = comps;
            this.dimension = comps.Length;
        }

        public Vector()
        {
            this.components = new double[0];
            this.dimension = -1;
        }

        // Initialise a dimension dimension 0 vector.
        public Vector(int dimension)
        {
            this.dimension = dimension;
            this.components = new double[dimension];
        }

        public double Dot(Vector vector)
        {
            if (this.dimension != vector.dimension) return Double.NaN;
            double dot = 0;
            for (int i = 0; i < this.dimension; i++)
            {
                dot += this.components[i] * vector.components[i];
            }
            return dot;
        }

        public void Scale(double lambda)
        {
            for (int i = 0; i < this.dimension; i++)
            {
                this.components[i] *= lambda;
            }
        }

        public double Abs()
        {
            double abssquared = this.Dot(this);
            return Math.Sqrt(abssquared);
        }

        public void Negate()
        {
            for (int i = 0; i < this.dimension; i++)
            {
                this.components[i] *= -1;
            }
        }

        public static Vector operator+ (Vector vector1, Vector vector2)
        {
            int dim = vector1.dimension;
            double[] comps = new double[dim];
            for (int i = 0; i < dim; i++)
            {
                comps[i] = vector1.components[i] + vector2.components[i];
            }
            return new Vector(comps);
        }

        public static Vector operator- (Vector vector1, Vector vector2)
        {
            int dim = vector1.dimension;
            double[] comps = new double[dim];
            for (int i = 0; i < dim; i++)
            {
                comps[i] = vector1.components[i] - vector2.components[i];
            }
            return new Vector(comps);
        }

        public static Vector operator* (double lambda, Vector vector)
        {
            int dim = vector.dimension;
            double[] comps = new double[dim];
            for (int i = 0; i < dim; i++)
            {
                comps[i] = lambda * vector.components[i];
            }
            return new Vector(comps);
        }

        // Finds the orthogonal component of this vector to a
        // list of orthogonal vectors and adds it to the list.
        // Assumes that all vectors in list are mutually orthogonal.
        public void AddToOrthList(List<Vector> list)
        {
            if(list.Count < 1)
            {
                list.Add(this);
            }
            else
            {
                Vector thisvector = this;
                foreach(Vector vector in list)
                {
                    thisvector -= ((thisvector.Dot(vector)) / (vector.Dot(vector))) * vector;
                }
                list.Add(/*(1/Math.Sqrt(thisvector.Dot(thisvector)))*/thisvector);
            }
        }

        // Takes a list of Vectors and turns it into a list of
        // orthonormal Vectors which span the same subspace.
        public void OrthoLise(List<Vector> list)
        {
            double scale = 1;
            while (list.Count > 0)
            {
                scale = list[0].Abs();
                if (Math.Abs(scale) < epsilon)
                {
                    list.RemoveAt(0);
                }
                else break;
            }
            
            if (list.Count > 0)
            {
                list[0].Scale(1 / scale);
                int place = 1;
                while (place < list.Count)
                {
                    for (int i = 0; i < place; i++)
                    {
                        list[place] -= list[place].Dot(list[i]) * list[place];
                    }
                    if (list[place].Dot(list[place]) < epsilon)
                    {
                        list.RemoveAt(place);
                    }
                    else
                    {
                        list[place].Scale(1 / list[place].Abs());
                        place++;
                    }
                }
            }
        }
    }
}
