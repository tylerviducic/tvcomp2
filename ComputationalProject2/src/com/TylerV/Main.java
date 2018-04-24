package com.TylerV;

public class Main {

    public static void main(String[] args) {
//        Look for a numerical solution of the Schr¨odinger equation equation for the
//        quantum harmonic oscillator using Numerov’s algorithm and the shooting
//        method as outlined in lectures. Your code should return the energy eigenvalue
//        En and eigenfunction ψn(x) with a pre-determined number n of nodes.
//        Compare your results with the analytic results for the lowest 3 energies. Your
//        write-up should contain

//       V(x) = m w^2 x^2 / 2
//       xc -> critical point, where V = E
//       hbar = w = 1
//       use n = 3 for test??
//       expected E = (n + 1/2)

        Numerov test1 = new Numerov();
        test1.numerov(10, 2, 300);






    }
}
