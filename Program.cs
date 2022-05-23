using System;

class Program
{
    static void Main(string[] args)
    {
        const float Ru = 8314.462f; //J⋅(K⋅kmol)
        const float g0 = 9.806f;

        float Mw = 2.016f;
        float T0 = 2200f;
        float k = 1.41f;

        float De = 0.25f; //m
        float Dt = 0.05f; //m
        float Ae;
        float At;
        float Aratio;
        
        float p0 = 50000000f;
        float pe = 1000000f;
        float pa = 0f;

        float R = Ru / Mw;
        float ρ0 = p0 / (R * T0); //kg/m3

        float engineMass = 250f;

        float mdot;
        float ve;
        float veq;
        float isp;
        float thrust;
        float twr;
        float mach;
        float a;

        void Inital()
        {
            Console.Title = "rnoz";

            Console.WriteLine("Type (help) for more comands \n");
        }

        Inital();

        while (true)
        {
            string comand = Console.ReadLine();

            if (comand == "help")
            {
                Console.WriteLine(HelpOutput());
            }
            else if (comand == "clear")
            {
                Console.Clear();
            }
            else if (comand == "run")
            {
                Run();
                Console.WriteLine(RunOutput());
            }
            else if (comand == "temp")
            {
                T0 = float.Parse(Console.ReadLine());
            }
            else if (comand == "gamma")
            {
                k = float.Parse(Console.ReadLine());
            }
            else if (comand == "molar-mass")
            {
                Mw = float.Parse(Console.ReadLine());
            }
            else if (comand == "engine-mass")
            {
                engineMass = float.Parse(Console.ReadLine());
            }
            else if (comand == "chamber-pressure")
            {
                p0 = float.Parse(Console.ReadLine()) * 1000000;
            }
            else if (comand == "ambient-pressure")
            {
                pa = float.Parse(Console.ReadLine()) * 1000000;
            }
            else if (comand == "exit-pressure")
            {
                pe = float.Parse(Console.ReadLine()) * 1000000;
            }
            else if (comand == "exit-diameter")
            {
                De = float.Parse(Console.ReadLine());
            }
            else if (comand == "#throat-diameter")
            {
                Dt = float.Parse(Console.ReadLine());
            }
        }

        void Run()
        {
            Ae = DiamToRadi(De);
            At = DiamToRadi(Dt);
            Aratio = Ae / At;

            mdot = (float)(At * p0 / Math.Sqrt(R * T0) * Math.Sqrt(k * Math.Pow((1 + k) / 2, (1 + k) / (1 - k))));

            ve = (float)Math.Sqrt(R * T0 * (2 * k / (k - 1)) * (1 - Math.Pow(pe / p0, (k - 1) / k)));

            //calculate performance
            veq = ve + (pe - pa) / mdot * Ae;
            isp = veq / g0;

            thrust = mdot * veq;
            twr = thrust / (engineMass * g0);

            //calculate mach number in exhaust
            a = (float)Math.Sqrt(k * R * T0);
            mach = ve / a;
        }

        string RunOutput()
        {
            string text = string.Format(
                "\n" +
                "Molar mass         {13} kmol \n" +
                "Spicific gas       {14} J⋅(K⋅kmol) \n" +
                "Gamma              {15} \n\n" +
                "Temp               {9} K \n" +
                "Isp                {5} s \n" +
                "Thrust             {1} kN \n" +
                "Mass flow rate     {2} kg/s \n" +
                "Exit velocity      {3} m/s \n" +
                "Equivalent         {4} m/s \n" +
                "                   {10} mach \n\n" +
                "Thrust to Weight   {0} \n" +
                "Engine mass        {11} kg \n\n" +
                "Chamber pressure   {6} MPa \n" +
                "Ambient pressure   {16} MPa \n" +
                "Exit pressure      {17} MPa \n\n" +
                "Exit diameter      {7} m \n" +
                "Throat diameter    {8} m \n" +
                "Area ratio         {12} m^2 \n", 
                twr, thrust / 1000, mdot, ve, veq, isp, p0 / 1000000, 
                RadiToDiam(Ae), RadiToDiam(At), T0, mach, 
                engineMass, Aratio, Mw, R, k, pa / 1000000, pe / 1000000);
            return text;
        }

        string HelpOutput()
        {
            string text = string.Format(
                "\n" +
                "(clear)            -clears the console \n" +
                "(run)              -runs rocket performance calculation \n" +
                "(temp)             -changes the chamber tempature \n" +
                "(gamma)            -changes the Heat capacity ratio \n" +
                "(molar-mass)       -changes the fuels molar mass \n" +
                "(engine-mass)      -changes the total engine mass \n" +
                "(chamber-pressure) -changes the chamber pressure \n" +
                "(ambient-pressure) -changes the atmoshperes pressure \n" +
                "(exit-pressure)    -changes the exit pressure \n" +
                "(exit-diameter)    -changes the exit diameter \n" +
                "(throat-diameter)  -changes the throat diameter \n");
                
            return text;
        }

        float RadiToDiam(float value)
        {
            value = ((float)Math.Sqrt((value / Math.PI) * 4));
            return value;
        }

        float DiamToRadi(float value)
        {
            value = ((float)(Math.PI * Math.Pow(value, 2) / 4));
            return value;
        }
    }
}