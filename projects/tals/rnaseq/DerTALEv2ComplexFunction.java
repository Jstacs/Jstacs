
package projects.tals.rnaseq;

import java.util.Arrays;
import de.jstacs.algorithms.optimization.EvaluationException;
import de.jstacs.algorithms.optimization.DimensionException;
import de.jstacs.algorithms.optimization.DifferentiableFunction;

public class DerTALEv2ComplexFunction extends DifferentiableFunction
{
    protected double[] yCounts;
    protected int[][] xIndicators;
    
    public DerTALEv2ComplexFunction( double[] yCounts,  int[][] xIndicators) throws IllegalArgumentException {
        this.yCounts = yCounts;
        this.xIndicators = xIndicators;
    }
    
    @Override
    public double evaluateFunction( double[] x) throws DimensionException, EvaluationException {
        double value = 0.0;
        for (int j = 0; j < this.xIndicators.length; ++j) {
            double tempBetaSum = 0.0;
            for (int c = 0; c < this.xIndicators[0].length; ++c) {
                tempBetaSum += this.xIndicators[j][c] * x[1 + c];
            }
            value += Math.pow(this.yCounts[j] - (x[0] + tempBetaSum), 2.0);
        }
        return value;
    }
    
    @Override
    public int getDimensionOfScope() {
        return 1 + this.xIndicators[0].length;
    }
    
    @Override
    public double[] evaluateGradientOfFunction( double[] x) throws DimensionException, EvaluationException {
        double dalpha = 0.0;
         double alpha = x[0];
        for (int j = 0; j < this.xIndicators.length; ++j) {
            double tempBetaSum = 0.0;
            for (int c = 0; c < this.xIndicators[0].length; ++c) {
                tempBetaSum += this.xIndicators[j][c] * x[1 + c];
            }
            dalpha += 2.0 * (alpha - this.yCounts[j] + tempBetaSum);
        }
         double[] dbeta = new double[this.xIndicators[0].length];
        Arrays.fill(dbeta, 0.0);
        for (int ct = 0; ct < this.xIndicators[0].length; ++ct) {
            for (int i = 0; i < this.xIndicators.length; ++i) {
                double tempBetaSum2 = 0.0;
                for (int c2 = 0; c2 < this.xIndicators[0].length; ++c2) {
                    tempBetaSum2 += this.xIndicators[i][c2] * x[1 + c2];
                }
                 double[] array = dbeta;
                 int n = ct;
                array[n] += -this.xIndicators[i][ct] * 2 * (this.yCounts[i] - alpha - tempBetaSum2);
            }
        }
         double[] dArray = new double[x.length];
        dArray[0] = dalpha;
        for (int c3 = 0; c3 < this.xIndicators[0].length; ++c3) {
            dArray[1 + c3] = dbeta[c3];
        }
        return dArray;
    }
}