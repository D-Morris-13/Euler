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
            List<string> numbers = ExtendNumber("");
            for (int length = 1; length < 3; length++)
            {
                numbers = ExtendNumberList(numbers,length);
            }
            
            Console.WriteLine("Count = " + numbers.Count);
//                        int maxlength = 20;
  //                      List<string> numberlist = NumberList(maxlength);
    //                    foreach (string number in numberlist) Console.WriteLine(number);
      //                  Console.WriteLine("Count = " + CountNumbers(numberlist, maxlength));

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

        static List<string> ConvertList(List<int[]> subintervals, int length)
        {
            List<string> numberlist = new List<string>();
            foreach(int[] division in subintervals)
            {
                numberlist.Add(DivisionsToString(division, length));
            }
            return numberlist;
        }

        static List<string> NumberList(int maxlength)
        {
            List<string> numberlist = new List<string>();
            for (int length = maxlength; length > 0; length--)
            {
                numberlist.AddRange(ConvertList(SubIntervals(length), length));
            }
            return numberlist;
        }

        static int CountNumbers(List<string> numbers, int length)
        {
            int count = 0;
            foreach (string number in numbers)
            {
                count += (number[0] == number[number.Length - 1] ? 1 : 2) + length - number.Length;
            }
            return count;
        }

        static int[] LexicographicIncrease(int[] division, int length)
        {
            int index = 0;
            int[] divisions = (int[])division.Clone();
            while (index < divisions.Length && divisions[index] < length) index++;
            if (index == divisions.Length) index--;
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

        static List<string> ExtendNumber(string number)
        {
            List<string> validextensions = new List<string>();
            int end = number.Length - 1;
            if (end != -1)
            {
                for (int digit = 0; digit < 10; digit++)
                {
                    if (digit <= (number[end] - '0') && number[end] <= number[0])
                    {
                        validextensions.Add(number + digit);
                    }
                    else if (digit > (number[end] - '0') && number[end] >= number[0])
                    {
                        validextensions.Add(number + digit);
                    }
                }
            }
            else
            {
                for (int digit = 1; digit < 10; digit++)
                {
                    validextensions.Add(digit.ToString());
                }
            }
            foreach (string enumber in validextensions) Console.WriteLine(enumber);
            return validextensions;
        }

        static List<string> ExtendNumberList(List<string> numbers, int maxlength)
        {
            List<string> extendablenumbers = numbers.FindAll( delegate(string number) { return number.Length == maxlength; });
            foreach (string number in extendablenumbers)
            {
                numbers.AddRange(ExtendNumber(number));
            }
            return numbers;

        }
    }

}
