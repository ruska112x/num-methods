package org.karabalin.second;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class MinimumResidualsTest {

    @Test
    void testSolve() {
        double[][] a = new double[][]{
                {5, 1, 2},
                {1, 3, 4},
                {2, 4, 10}
        };
        double[] b = new double[]{12, 13, 14};
        double[] r = MinimumResiduals.solve(a, b, new double[]{1, 1, 1});
        MinimumResiduals.printVector(r);
    }
}