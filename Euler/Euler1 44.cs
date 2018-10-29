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
            Ellipse mirror = new Ellipse(4, 1, 10);
            double[] source = new double[2] { 0, 10.1 };
            double[] target = new double[2] { 1.4,-9.6 };
            Beam firstbeam = new Beam(source,target,0,mirror);
            firstbeam.PrintBeam();
            Beam nextbeam = firstbeam.NextBeam();
            nextbeam.PrintBeam();
            while (!nextbeam.CheckExit() && nextbeam.GetCount() < 5000)
            {
                nextbeam = nextbeam.NextBeam();
//                nextbeam.PrintBeam();
                nextbeam.NextBeam();
            }
            nextbeam.PrintBeam();

        }
        
        
    }

    class Beam
    {
        double[] source;
        double[] target;
        long count;
        Ellipse mirror;

        public Beam(double[] src, double[] tgt, long n, Ellipse mirrors)
        {
            this.source = (double[]) src.Clone();
            this.target = (double[]) tgt.Clone();
            this.count = n;
            this.mirror = mirrors;
        }

        public Beam NextBeam()
        {
            double[] src = (double[]) this.target.Clone();
            double[] tgt = this.NextBounce();
            long nextcount = this.count + 1;
            Beam nextbeam = new Beam(src, tgt, nextcount, this.mirror);
            return nextbeam;
        }

        double[] NextBounce()
        {
            double Mx, My, mx, my;
            mx = this.target[1] - this.source[1];
            my = this.source[0] - this.target[0];
            Mx = this.target[0] * this.mirror.a;
            My = this.target[1] * this.mirror.b;

            double bmx, bmy;
            bmx = (My * My - Mx * Mx) * mx - 2 * my * Mx * My;
            bmy = (Mx * Mx - My * My) * my - 2 * mx * Mx * My;
            
//             Console.WriteLine("mx = " + mx + "  my = " + my + "  Mx = " + Mx + "  My = " + My + "  bmx = " + bmx + "  bmy = " + bmy);

            double[] xy = new double[2];
            xy[0] = (mirror.b * bmx * bmx * target[0] - mirror.a * bmy * bmy * target[0] + 2 * mirror.b * bmx * bmy * target[1]) / (mirror.b * bmx * bmx + mirror.a * bmy * bmy);
            xy[1] = (mirror.a * bmy * bmy * target[1] - mirror.b * bmx * bmx * target[1] + 2 * mirror.a * bmx * bmy * target[0]) / (mirror.b * bmx * bmx + mirror.a * bmy * bmy);

            return xy;
        }

        public bool CheckExit()
        {
            if (Math.Abs(target[0]) <= 0.01 && target[1] > 0)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        public void PrintBeam()
        {
            string output = "";
            output += "Source    : ( " + this.source[0] + " , " + this.source[1] + "  ) \n";
            output += "Target    : ( " + this.target[0] + " , " + this.target[1] + "  ) \n";
            output += "Count     : " + this.count;
            Console.WriteLine(output);
        }

        public long GetCount()
        {
            return this.count;
        }

    }

    class Ellipse
    {
        // a x^2 + b y^2 = r^2
        public double a;
        public double b;
        public double r;

        public Ellipse(double ai, double bi, double ri)
        {
            this.a = ai;
            this.b = bi;
            this.r = ri;
        }
    }

}
