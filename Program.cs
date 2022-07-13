using System;

class Program
{
    static void Main(string[] args)
    {
        const double Ru = 8314.462; //J⋅(K⋅kmol)
        const double g0 = 9.806;

        double Mw = 2.016;
        double T0 = 3200; //K
        double k = 1.66;

        double De = 3; //m
        double Dt = 0.25; //m
        double Ae;
        double At;
        double Aratio;

        double p0 = 5000000;
        double pe = 1000;
        double pa = 10000;

        double R = Ru / Mw;
        //double ρ0 = p0 / (R * T0); //kg/m3

        double engineMass = 1200;

        double mdot;
        double ve;
        double veq;
        double isp;
        double thrust;
        double twr;
        double mach;
        double a;

        void Inital()
        {
            Run();
            Console.WriteLine(RunOutput());
        }

        Inital();

        void Run()
        {
            Ae = DiamToRadi(De);
            At = DiamToRadi(Dt);
            Aratio = Ae / At;

            mdot = (At * p0 / Math.Sqrt(R * T0) * Math.Sqrt(k * Math.Pow((1 + k) / 2, (1 + k) / (1 - k))));

            ve = Math.Sqrt(R * T0 * (2 * k / (k - 1)) * (1 - Math.Pow(pe / p0, (k - 1) / k)));

            veq = ve + (pe - pa) / mdot * Ae;
            isp = veq / g0;

            thrust = mdot * veq;
            twr = thrust / (engineMass * g0);

            a = Math.Sqrt(k * R * T0);
            mach = ve / a;
        }

        string RunOutput()
        {
            string text = string.Format(
                "\n" +
                "Molar mass         {13} kmol \n" +
                "Specific gas       {14} J⋅(K⋅kmol) \n" +
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
                "Chamber pressure   {6} bar \n" +
                "Ambient pressure   {16} bar \n" +
                "Exit pressure      {17} bar \n\n" +
                "Exit diameter      {7} m \n" +
                "Throat diameter    {8} m \n" +
                "Area ratio         {12} m^2 \n",
                twr, thrust / 1000, mdot, ve, veq, isp, p0 / 100000,
                RadiToDiam(Ae), RadiToDiam(At), T0, mach,
                engineMass, Aratio, Mw, R, k, pa / 100000, pe / 100000);
            return text;
        }

        double RadiToDiam(double value)
        {
            value = Math.Sqrt((value / Math.PI) * 4);
            return value;
        }

        double DiamToRadi(double value)
        {
            value = (Math.PI * Math.Pow(value, 2)) / 4;
            return value;
        }
    }
}