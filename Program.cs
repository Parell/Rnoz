using System;

internal class Program
{
    static void Main(string[] args)
    {
        const float Ru = 8314.462f; //J⋅(K⋅kmol)
        const float g0 = 9.806f;

        float Mw = 2.016f;
        float T0 = 2200f;
        float k = 1.41f;

        float Aratio = 20f;
        float De = 0.25f;
        float Ae = (float)(Math.PI * Math.Pow(De, 2) / 4);
        float At = Ae / Aratio;

        float p0 = 50000000f;
        float pe = 101325f;
        float pa = 0f;

        float R = Ru / Mw;
        float ρ0 = p0 / (R * T0); //kg/m3

        var engineMass = 250; // kg

        float mdot;
        float ve;
        float veq;
        float isp;
        float thrust;
        float twr;

        void Run()
        {
            mdot = (float)(At * p0 / Math.Sqrt(R * T0) * Math.Sqrt(k * Math.Pow((1 + k) / 2, (1 + k) / (1 - k))));

            ve = (float)Math.Sqrt(R * T0 * (2 * k / (k - 1)) * (1 - Math.Pow(pe / p0, (k - 1) / k)));

            veq = ve + (pe - pa) / mdot * Ae;
            isp = veq / g0;

            thrust = mdot * veq;
            twr = thrust / (engineMass * g0);
        }

        Run();

        string Output()
        {
            string text = string.Format("\n" +
                                     "Temp             {9} K \n" +
                                     "Thrust           {1} kN \n" +
                                     "Mass flow rate   {2} kg/s \n" +
                                     "Exit velocity    {3} m/s \n" +
                                     "                 {4} m/s \n\n" +
                                     "Thrust to Weight {0} \n" +
                                     "Chamber Pressure {6} MPa \n" +
                                     "Isp              {5} s \n\n" +
                                     "Exit diameter    {7} cm \n" +
                                     "Throat diameter  {8} cm \n", twr, thrust / 1000, mdot, ve, veq, isp, p0 / 1000000,
                                     (float)Math.Sqrt((Ae / Math.PI) * 4) * 100, (float)Math.Sqrt((At / Math.PI) * 4) * 100, T0);
            return text;
        }

        Console.WriteLine(Output());
    }
}