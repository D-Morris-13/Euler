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
            bool bouncy;
            long test = 5;
            bouncy = IsBouncy(test);

            long[] testinputs = {1,100,3543654,123456,654321 };

            for (int i = 0; i < testinputs.Length; i++)
            {
//                Console.WriteLine(testinputs[i] + " is bouncy? " + IsBouncy(testinputs[i])); 
            }

            long bouncecount = 0;
            long testnumber = 1;

            while(true)
            {
//                Console.Write(testnumber + " is ");
                if (IsBouncy(testnumber))
                {
                    bouncecount++;
//                    Console.Write(" bouncy.");
                } else {
//                    Console.Write(" not bouncy.");
                }
                
  //              Console.WriteLine();
//                Console.WriteLine("Bouncecount = " + bouncecount);
                if (bouncecount * 100 >= testnumber * 99) break;
                testnumber++;
            }

            Console.WriteLine(testnumber);
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
    }

}
