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
            int maxlength = 100;
            long count = 0;
            long subcount;
            for (int length = 1; length <= maxlength; length++)
            {
                List<int[]> subintervals = SubIntervals(length);
                subcount = CountNumbers(subintervals, length, maxlength);
                count += subcount;
                Console.WriteLine("Length = " + length + " Monic = " + subcount + " Total count = " + count);
            }
            Console.WriteLine("Count = " + count);

        }

        static bool IsBouncy(long number)
        {
            bool increasing = true;
            bool decreasing = true;
            int index = 0;

            string numberstring = number.ToString();

            while (index < numberstring.Length - 1 && (increasing || decreasing))
            {
                if (increasing && numberstring[index] > numberstring[index + 1]) increasing = false;
                if (decreasing && numberstring[index] < numberstring[index + 1]) decreasing = false;
                index++;
            }

            return !(increasing || decreasing);
        }

        static List<int[]> SubIntervals(int length)
        {
            int[] divisions = new int[8];
            List<int[]> subintervals = new List<int[]>();
            if (length <= 0) return subintervals;           
            while(true)
            {
                subintervals.Add(divisions);
                divisions = (int[]) LexicographicIncrease(divisions, length).Clone();
                if (divisions[7] == 0) break;
//                string output = DivisionsToString(divisions, length);
//                Console.WriteLine(output);
            }

            return subintervals;
        }

        static long CountNumbers(List<int[]> intervals, int length, int maxlength)
        {
            long count = 0;
            foreach (int[] division in intervals)
            {
                count += 2 + maxlength - length;
                if (division[0] == length)
                {
                    count--;
                    continue;
                }
                for(int digit = 1; digit < division.Length; digit++)
                {
                    if (division[digit-1] == 0 && division[digit] == length)
                    {
                        count--;
                        break;
                    }
                }
                if (division[division.Length - 1] == 0) count--;
            }
            return count;
        }

        static int[] LexicographicIncrease(int[] division, int length)
        {
            int index = 0;
            int[] divisions = (int[])division.Clone();
            while (index < divisions.Length -1 && divisions[index] < length) index++;
            if (divisions[index] < length)
            {
                divisions[index]++;
            }
            else if (index > 0)
            {
                divisions[index - 1]++;
                while (index < divisions.Length) {
                    divisions[index] = divisions[index - 1];
                    index++;
                }
            }
            else divisions = new int[division.Length];

              return divisions;
        }

        static string DivisionsToString(int[] divisions, int length)
        {
            int digit = 0;
            string output = "";
            for (int pos = 0; pos < length; pos++)
            {
                while (digit < 8 && divisions[digit] <= pos) digit++;
                output += (digit+1);
            }
            return output;
        }
    }

}
