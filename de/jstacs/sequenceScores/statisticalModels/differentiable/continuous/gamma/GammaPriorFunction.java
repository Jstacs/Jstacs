package de.jstacs.sequenceScores.statisticalModels.differentiable.continuous.gamma;

import umontreal.iro.lecuyer.util.Num;
import de.jstacs.algorithms.optimization.DifferentiableFunction;
import de.jstacs.algorithms.optimization.DimensionException;
import de.jstacs.algorithms.optimization.EvaluationException;
import de.jtem.numericalMethods.calculus.specialFunctions.Gamma;


public class GammaPriorFunction extends DifferentiableFunction implements IntegrableFunction {

	private double chi1;
	private double chi2;
	private double nu;
	
	public GammaPriorFunction(double mua, double mug, double ess){
		this.chi1 = Math.log( mug );
		this.chi2 = -mua;
		this.nu = ess;
	}
	

	public double evaluateSecondOrderDerivativeAt( double[] x, int dimension ) throws DimensionException, EvaluationException {
		
		double temp = Math.exp( nu*x[0]*Math.log( x[1] ) - nu*Gamma.logOfGamma( x[0] ) + nu*(chi1*x[0] + chi2*x[1]) );
		
		if(dimension == 0){
			return temp * nu * ( nu*Math.pow( chi1 + Math.log( x[1] ) - Num.digamma( x[0] ), 2.0) - Num.trigamma( x[0] ) );
		}else if(dimension == 1){
			return temp / x[1] / x[1] * nu * ( x[0]*(nu*x[0] - 1.0) + 2.0*nu*chi2*x[0]*x[1] + nu*chi2*chi2*x[1]*x[1] );
		}else{
			throw new EvaluationException("wrong dimension");
		}
	}


	public double[][] getLimits() {
		return new double[][]{{1E-10,Double.POSITIVE_INFINITY},{1E-10,Double.POSITIVE_INFINITY}};
	}


	public double[] evaluateGradientOfFunction( double[] x ) throws DimensionException, EvaluationException {
		
		double temp = Math.exp( nu*x[0]*Math.log( x[1] ) - nu*Gamma.logOfGamma( x[0] ) + nu*(chi1*x[0] + chi2*x[1]) );
		
		double[] grad = new double[2];
		
		grad[0] = temp*nu*(chi1 + Math.log( x[1] ) - Num.digamma( x[0] ));
		grad[1] = temp/x[1]*nu*(x[0] + chi2*x[1]);
		
		return grad;
	}

	public double evaluateFunction( double[] x ) throws DimensionException, EvaluationException {
		double currVal = -Gamma.logOfGamma( x[0] )
								+ Math.log( x[1] )*x[0] 
								+chi2*x[1] 
								+ x[0]*chi1;
		currVal *= nu;
		currVal = Math.exp( currVal );
		return currVal;
	}

	public int getDimensionOfScope() {
		return 2;
	}


	public boolean isUnimodal() {
		return true;
	}

}
