package org.karabalin.second;

import java.util.Arrays;
import java.util.OptionalDouble;

public class MinimumResiduals {
    public static double[] multiplyMatrixAndVector(double[][] m, double[] v) {
        double[] r = new double[v.length];
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < v.length; j++) {
                r[i] += m[i][j] * v[j];
            }
        }
        return r;
    }

    public static double scalarProduct(double[] a, double[] b) {
        double r = 0;
        for (int i = 0; i < a.length; i++) {
            r += a[i] * b[i];
        }
        return r;
    }

    public static double vectorInfNorm(double[] v) {
        OptionalDouble optionalDouble = Arrays.stream(v).map(Math::abs).max();
        if (optionalDouble.isPresent()) {
            return optionalDouble.getAsDouble();
        } else {
            return -1;
        }
    }

    public static double[] vectorDifference(double[] a, double[] b) {
        double[] r = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            r[i] = a[i] - b[i];
        }
        return r;
    }

    public static double matrixInfNorm(double[][] m) {
        double norm = 0.0;
        for (double[] doubles : m) {
            double s = 0.0;
            for (int j = 0; j < m.length; ++j) {
                s += Math.abs(doubles[j]);
            }
            if (s > norm) {
                norm = s;
            }
        }
        return norm;
    }

    public static void printVector(double[] v) {
        for (int i = 0; i < v.length; i++) {
            if (i + 1 != v.length) {
                System.out.printf("%f, ", v[i]);
            } else {
                System.out.printf("%f", v[i]);
            }
        }
        System.out.println();
    }

    public static double[] solve(double[][] a, double[] b, double[] x0) {
        double epsilon = 0.000000001;
        double maxIterations = 1000;

        int n = b.length;
        double[] x = Arrays.copyOf(x0, n);
        double[] r = new double[n];
        double[] Ar = new double[n];
        double[] Ax = new double[n];

        for (int i = 0; i < maxIterations; i++) {
            // r = b - A * x
            Ax = multiplyMatrixAndVector(a, x);
            for (int j = 0; j < n; j++) {
                r[j] = b[j] - Ax[j];
            }

            // A * r
            Ar = multiplyMatrixAndVector(a, r);

            // (r, r) и (Ar, Ar)
            double rr = scalarProduct(r, r);
            double ArAr = scalarProduct(Ar, Ar);

            // alpha
            double alpha = rr / ArAr;

            for (int j = 0; j < n; j++) {
                x[j] = x[j] + alpha * r[j];
            }

            double residualNorm = vectorInfNorm(r);
            if (residualNorm < epsilon) {
                break;
            }
        }
        return x;
    }
}
