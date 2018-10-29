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
            long[,] startcounts = new long[10, 3] { { 0, 0, 0 },{ 0, 1, 0 }, { 0, 1, 0 }, { 0, 1, 0 }, { 0, 1, 0 }, { 0, 1, 0 }, { 0, 1, 0 }, { 0, 1, 0 }, { 0, 1, 0 }, { 0, 1, 0 } };
            long totalcount = SumCounts(startcounts);
            int maxlength = 100;
            long[,] nextcounts = (long[,])startcounts.Clone();
            for (int length = 2; length <= maxlength; length++)
            {
                nextcounts = (long[,])NextLength(nextcounts).Clone();
                totalcount += SumCounts(nextcounts);
            }
            PrintCounts(nextcounts);
            Console.WriteLine(totalcount);
        }
        
        static long[,] NextLength(long[,] currentcounts)
        {
            long[,] nextcounts = new long[10, 3];
            // Count Increasing Numbers.
            for(int nextdigit = 0; nextdigit < 10; nextdigit++)
            {
                long increasecount = 0;
                for(int lastdigit = 0; lastdigit < nextdigit; lastdigit++)
                {
                    increasecount += currentcounts[lastdigit, 0] + currentcounts[lastdigit, 1];
                }
                increasecount += currentcounts[nextdigit, 0];
                nextcounts[nextdigit, 0] = increasecount;
            }

            // Count Decreasing Numbers.
            for (int nextdigit = 0; nextdigit < 10; nextdigit++)
            {
                long increasecount = 0;
                for(int lastdigit = nextdigit + 1; lastdigit < 10; lastdigit++)
                {
                    increasecount += currentcounts[lastdigit, 2] + currentcounts[lastdigit, 1];
                }
                increasecount += currentcounts[nextdigit, 2];
                nextcounts[nextdigit, 2] = increasecount;
            }

            // Count Constant Numbers.
            for (int nextdigit = 0; nextdigit < 10; nextdigit++)
            {
                nextcounts[nextdigit, 1] = currentcounts[nextdigit, 1];
            }

            return nextcounts;
        }

        static long SumCounts(long[,] currentcounts)
        {
            long count = 0;
            for (int digit = 0; digit < 10; digit++)
            {
                for (int index = 0; index < 3; index++)
                {
                    count += currentcounts[digit, index];
                }
            }
            return count;
        }

        static void PrintCounts(long[,] currentcounts)
        {
            string output = "/";
            for (int i = 0; i < 10; i++)
            {
                output += "|  "+ i + "  " + String.Format("{0,5} ", currentcounts[i, 0]) + " |";
            }
            output += "\\ \n";
            output += "|";
            for (int i = 0; i < 10; i++)
            {
                output += "|  " + i + "  " + String.Format("{0,5} ", currentcounts[i, 1]) + " |";
            }
            output += "| \n";
            output += "\\";
            for (int i = 0; i < 10; i++)
            {
                output += "|  " + i + "  " + String.Format("{0,5} ", currentcounts[i, 2]) + " |";
            }
            output += "/ \n";
            Console.Write(output);
        }
    }

}
