package org.karabalin.first;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class GaussTest {
    @Test
    public void testInverse() {
        double[][] a = new double[][]{
                {4, 3},
                {3, 2}
        };
        Gauss.printMatrix(Gauss.inverse(a));
    }
}