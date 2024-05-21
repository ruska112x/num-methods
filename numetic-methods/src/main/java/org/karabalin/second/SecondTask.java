package org.karabalin.second;

import org.karabalin.Generator;

import java.util.Arrays;

public class SecondTask {
    public static void main(String[] args) {
        Generator generator = new Generator();

        int n = 3;
        double alpha; // alpha
        double beta; // beta
        double aNorm; // norm of the matrix a
        double aInvNorm; // norm of the inverse matrix a
        double vA; // conditioning factor
        double errorNorm; // norm of the error matrix
        double dzeta; // dzeta (division of norm of the error matrix and norm of the inverse matrix a)
        double rNorm; // norm of the incoherence matrix
        double rho; // relative nevyazka

        int cycles = 12;

        System.out.println("-".repeat(136));
        System.out.println("|    alpha     |     beta     |   norm A     |  norm A^-1   |    obusl     |   norm err   |    dzeta     |   nevyazka   |  ro nevyazka |");
        System.out.println("-".repeat(136));
        System.out.println("|    alpha     |     beta     |  ||A||inf    | ||A^-1||inf  |   Vinf(A)    |   ||Z||inf   |    dzeta     |   ||R||inf   |       ro     |");
        System.out.println("-".repeat(136));
        for (int i = 0; i < cycles; i++) {
            alpha = Math.pow(10, -(i + 1));
            beta = 1;

            double[][] a = new double[n][n];
            double[][] aInv = new double[n][n];
            generator.myGen(a, aInv, n, alpha, beta, 1, 2, 1, 1);
            double[] xActual = new double[n];
            Arrays.fill(xActual, 1);

            double[] b = MinimumResiduals.multiplyMatrixAndVector(a, xActual);

//            generator.printMatrix(a, n);

            double[] x0 = new double[n];
            Arrays.fill(x0, Math.PI);
            double[] x = MinimumResiduals.solve(a, b, x0);

//            MinimumResiduals.printVector(x);

            aNorm = generator.matrixInfNorm(a, n);
            aInvNorm = generator.matrixInfNorm(aInv, n);
            vA = aNorm * aInvNorm;
            errorNorm = MinimumResiduals.vectorInfNorm(MinimumResiduals.vectorDifference(x, xActual));
            dzeta = errorNorm / MinimumResiduals.vectorInfNorm(xActual);
            rNorm = MinimumResiduals.vectorInfNorm(MinimumResiduals.vectorDifference(MinimumResiduals.multiplyMatrixAndVector(a, x), b));
            rho = rNorm / MinimumResiduals.vectorInfNorm(b);

            System.out.printf("| %e | %e | %e | %e | %e | %e | %e | %e | %e |\n", alpha, beta, aNorm, aInvNorm, vA, errorNorm, dzeta, rNorm, rho);
        }
        System.out.println("-".repeat(136));
        for (int i = 0; i < cycles; i++) {
            alpha = 1;
            beta = Math.pow(10, (i + 1));

            double[][] a = new double[n][n];
            double[][] aInv = new double[n][n];
            generator.myGen(a, aInv, n, alpha, beta, 1, 2, 1, 1);
            double[] xActual = new double[n];
            Arrays.fill(xActual, 1);

            double[] b = MinimumResiduals.multiplyMatrixAndVector(a, xActual);

//            generator.printMatrix(a, n);

            double[] x0 = new double[n];
            Arrays.fill(x0, Math.PI);
            double[] x = MinimumResiduals.solve(a, b, x0);

//            MinimumResiduals.printVector(x);

            aNorm = generator.matrixInfNorm(a, n);
            aInvNorm = generator.matrixInfNorm(aInv, n);
            vA = aNorm * aInvNorm;
            errorNorm = MinimumResiduals.vectorInfNorm(MinimumResiduals.vectorDifference(x, xActual));
            dzeta = errorNorm / MinimumResiduals.vectorInfNorm(xActual);
            rNorm = MinimumResiduals.vectorInfNorm(MinimumResiduals.vectorDifference(MinimumResiduals.multiplyMatrixAndVector(a, x), b));
            rho = rNorm / MinimumResiduals.vectorInfNorm(b);

            System.out.printf("| %e | %e | %e | %e | %e | %e | %e | %e | %e |\n", alpha, beta, aNorm, aInvNorm, vA, errorNorm, dzeta, rNorm, rho);
        }
        System.out.println("-".repeat(136));
    }
}
