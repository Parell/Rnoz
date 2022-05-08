internal class Program
{
    static void Main(string[] args)
    {
        const float Ru = 8314.462f; //J⋅(K⋅kmol)
        const float g0 = 9.806f;

        float Mw = 2.016f;
        float T0 = 3500f;
        float k = 1.405f;

        float Aratio = 54f;
        float De = 1f;
        float Ae = (float)(Math.PI * Math.Pow(De, 2) / 4);
        float At = Ae / Aratio;

        float p0 = 20000000f;
        float pe = 101325f;
        float pa = 101325f;

        float R = Ru / Mw;
        float ρ0 = p0 / (R * T0); //kg/m3

        var mass = 600; // kg

        float Γ = (float)Math.Sqrt(k * Math.Pow((1 + k) / 2, (1 + k) / (1 - k)));
        float mdot = (float)(At * p0 / Math.Sqrt(R * T0) * Γ);

        float ve = (float)Math.Sqrt(R * T0 * (2 * k / (k - 1)) * (1 - Math.Pow(pe / p0, (k - 1) / k)));

        float veq = ve + (pe - pa) / mdot * Ae;
        float isp = veq / g0;

        float thrust = mdot * veq;


        float twr = thrust / (mass * g0);

        string s = string.Format("\nTemp             {9} K \n" +
                                 "Thrust           {1} kN \n" +
                                 "Mass flow rate   {2} kg/s \n" +
                                 "Exit velocity    {3} m/s \n" +
                                 "                 {4} m/s \n\n" +
                                 "Thrust to Weight {0} \n" +
                                 "Chamber Pressure {6} MPa \n" +
                                 "Isp              {5} s \n\n" +
                                 "Exit area        {7} m2 \n" +
                                 "Throat area      {8} m2 \n", twr, thrust / 1000, mdot, ve, veq, isp, p0 / 1000000, Ae, At, T0);

        Console.WriteLine(s);
    }
}