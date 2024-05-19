package org.karabalin;

public class Main {
    public static void main(String[] args) {
        Generator generator = new Generator();

        int n = 100;
        double alpha; // alpha
        double beta; // beta
        double aNorm; // norm of the matrix a
        double aInvNorm; // norm of the inverse matrix a
        double vA; // conditioning factor
        double errorNorm; // norm of the error matrix
        double dzeta; // dzeta (division of norm of the error matrix and norm of the inverse matrix a)
        double rNorm; // norm of the incoherence matrix

        int cycles = 8;

        System.out.println("-".repeat(121));
        System.out.println("|    alpha     |     beta     |   norm A     |  norm A^-1   |    obusl     |   norm err   |    dzeta     |   nevyazka   |");
        System.out.println("-".repeat(121));
        System.out.println("|    alpha     |     beta     |  ||A||inf    | ||A^-1||inf  |   Vinf(A)    |   ||Z||inf   |    dzeta     |   ||R||inf   |");
        System.out.println("-".repeat(121));
        for (int i = 0; i < cycles; i++) {
            alpha = Math.pow(10, -(i + 1));
            beta = 1;

            double[][] a = new double[n][n];
            double[][] aInv = new double[n][n];
            generator.myGen(a, aInv, n, alpha, beta, 1, 2, 0, 1);
            double[][] b = Gauss.inverse(a);

            aNorm = generator.matrixInfNorm(a, n);
            aInvNorm = generator.matrixInfNorm(aInv, n);
            vA = aNorm * aInvNorm;
            errorNorm = generator.matrixInfNorm(Gauss.matrixDifference(aInv, b), n);
            dzeta = errorNorm / aInvNorm;
            rNorm = generator.matrixInfNorm(Gauss.matrixDifference(Gauss.matrixMultiplication(b, a), Gauss.createIdentityMatrix(n)), n);

            System.out.printf("| %e | %2.10f | %2.10f | %e | %e | %e | %e | %e |\n", alpha, beta, aNorm, aInvNorm, vA, errorNorm, dzeta, rNorm);
        }
        System.out.println("-".repeat(121));
        for (int i = 0; i < cycles; i++) {
            alpha = 1;
            beta = Math.pow(10, -(i + 1));

            double[][] a = new double[n][n];
            double[][] aInv = new double[n][n];
            generator.myGen(a, aInv, n, alpha, beta, 1, 2, 0, 1);
            double[][] b = Gauss.inverse(a);

            aNorm = generator.matrixInfNorm(a, n);
            aInvNorm = generator.matrixInfNorm(aInv, n);
            vA = aNorm * aInvNorm;
            errorNorm = generator.matrixInfNorm(Gauss.matrixDifference(aInv, b), n);
            dzeta = errorNorm / aInvNorm;
            rNorm = generator.matrixInfNorm(Gauss.matrixDifference(Gauss.matrixMultiplication(b, a), Gauss.createIdentityMatrix(n)), n);

            System.out.printf("| %2.10f | %e | %2.10f | %e | %e | %e | %e | %e |\n", alpha, beta, aNorm, aInvNorm, vA, errorNorm, dzeta, rNorm);
        }
        System.out.println("-".repeat(121));
    }
}