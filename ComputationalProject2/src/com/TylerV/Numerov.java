package com.TylerV;

public class Numerov {



    public void numerov(double xmax, double n, int iterations){
        double h = xmax /iterations;
        double eMax;
        double eMin;
        double energy;
        double energyRange;
        double nodes;
        boolean cont = true;

        double[] yn = new double[iterations];
        double[] fn = new double[iterations];
        double[] gn = new double[iterations];
        double[] psi = new double[yn.length * 2];
        double[] v = new double[iterations];

        for (int i = 0; i < v.length; i++){
            double x = i * h;
            v[i] = potential(x);
        }

        eMax = findEmax(v);
        eMin = findEmin(v);


        while (cont) {
            energy = (eMax + eMin) / 2;
            System.out.println("Energy is " + energy);
            for (int i = 0; i < v.length; i++) {
                gn[i] = gn(energy, v[i]);
                fn[i] = fn(gn[i], h);
            }

            if (n % 2 == 0) {
                yn[0] = h;
                yn[1] = ((12 - 10 * fn[0]) * yn[0]) / 2 * fn[1];
                for (int i = 2; i < yn.length; i++) {
                    yn[i] = (((12 - 10 * fn[i - 1]) * yn[i - 1]) - fn[i - 2] * yn[i - 2]) / fn[i];
                }
                for (int i = 0; i<psi.length/2; i++){
                    psi[i] = yn[yn.length -(1+i) ];
                }
                for (int i = psi.length/2; i<psi.length; i++){
                    psi[i] = yn[i - psi.length/2];
                }
            } else if (n % 2 != 0) {
                yn[0] = 0;
                yn[1] = h;
                for (int i = 2; i < yn.length; i++) {
                    yn[i] = (((12 - 10 * fn[i - 1]) * yn[i - 1]) - fn[i - 2] * yn[i - 2]) / fn[i];
                }
                for (int i = 0; i<psi.length/2; i++){
                    psi[i] = -yn[yn.length -(1 + i)];
                }
                for (int i = psi.length/2; i<psi.length; i++){
                    psi[i] = yn[i - psi.length/2];
                }
            }

            nodes = checkNodes(psi);
            System.out.println(nodes);
            if (nodes > n){
                eMax = energy;
            } else if (nodes <= n){
                eMin = energy;
            }

            energyRange = eMax - eMin;
            if (energyRange < .000000001){
                cont = false;
            }

//            for (int i = 0; i < yn.length; i++){
//                System.out.println(yn[i]);
//            }
        }

        for (int i = 0; i < psi.length; i++){
            System.out.println("x = " + (-10 + i*h) + "y = " + psi[i]);
        }
    }

    public double potential (double x){
        return (x*x)/2;
    }

    public double gn(double e, double v){
        return 2 * (e - v);
    }

    public double fn(double g, double h){
        return 1 + (g * (h*h) /12);
    }



    public int checkNodes(double[] array){
        int nodes=0;
        for (int i = 1; i < array.length; i++){
            if ((array[i] > 0 && array[i-1] < 0) || (array[i] < 0 && array[i-1] > 0))
                nodes++;
        }
        return nodes;
    }

    public double findEmax(double[] array){
        double eMax = -1000000;
        for (int i = 0; i < array.length; i++){
            double e = array[i];
            if (e > eMax){
                eMax = e;
            }
        }
        return eMax;
    }

    public double findEmin(double[] array){
        double eMin = 10000000;
        for (int i = 0; i < array.length; i++){
            double e = array[i];
            if (e < eMin){
                eMin = e;
            }
        }
        return eMin;
    }

/*
go from -10 to 10
300 equally spaced intervals
create array V(i)
Emax = max of V(i)
Emin = min of V(i)
E = (Emax + Emin)/2
integrate from x = 0 in positive direction, count nodes. if nnodes > n, Emin = E, do again
if nnodes < n, E max = E, do again
when Emax - Emin < some given value, convergence ( 1 * 10^-10)

yn1 = [2*yn * (1 - 5*h^2*gn/12)- y(n-1) * (1 + h^2*g(n-1)/12) ] / (1 + h^2*g(n+1)/12)
g(x) = 2(E-V(x)

even better: fn = 1 + gn*h^2/12
with gn = 2*(E-V(xn))
so y(n+1) = ( (12 - 10 *fn)yn - f(n-1)*y(n-1) ) / f(n+1)

for n odd, y0 = 0 and y1 = an arbitrary finite value
for n even, y0 is arbitrary and finite, y1 is determined by numerov's formula f1 = f(-1) and y1 = y(-1)
n even y1 = (12-10f0)*y0/2f1
*/
}
