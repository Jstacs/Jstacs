package de.jstacs.sequenceScores.statisticalModels.differentiable.continuous.gamma;

import java.util.LinkedList;

import de.jstacs.algorithms.optimization.DimensionException;
import de.jstacs.algorithms.optimization.EvaluationException;
import de.jstacs.algorithms.optimization.Function;


/**
 * Class for numerical integration of {@link IntegrableFunction}s using either the trapezoid rule (slow but accurate for small delta), or the adaptive
 * trapezoid rule using nested intervals (fast, accuracy depends on accepted error).
 * @author Jan Grau
 *
 */
public class NumericalIntegration {

	/**
	 * Computes the integral under an {@link IntegrableFunction} <code>fun</code> within the limits given by {@link IntegrableFunction#getLimits()}
	 * using an adaptive trapezoid rules. It initially starts with a number of intervals for which the integral is computed using the trapezoid rule. Each of these
	 * intervals is divided into two intervals of half the size. If the sum of the integrals of these sub-intervals does not differ more than
	 * <code>acceptedEps</code> percent from the integral computed for the enclosing interval, the division is stopped. Otherwise, the two sub-intervals
	 * are again split into two new intervals each. If no intervals remain to be split, the evaluation terminates and return the computed integral.
	 * 
	 * 
	 * @param fun the function
	 * @param minValue the minimum function value that is considered important for unimodal functions when searching for a bracketing interval for integration
	 * @param initDelta the initial width of the trapezoids for searching a bracketing interval for unimodal functions
	 * @param acceptedEps the accepted error between two steps of nesting of a sub-integral, if the error falls below this value, nesting is stopped
	 * @return the value of the integral
	 * @throw DimensionException is thrown if the function does not accept parameter values of dimension {@link Function#getDimensionOfScope()}
	 * @throws EvaluationException is thrown if any specified limit of a non-unimodal function is infinite or both limits are infinite for a unimodal function
	 */
	public static double getIntegralByNestedIntervals(IntegrableFunction fun, double minValue, double initDelta, double acceptedEps) throws DimensionException, EvaluationException{
		
		int dim = fun.getDimensionOfScope();
		
		double[][] limits = fun.getLimits();
		double[] parValues = new double[dim];
		for(int i=0;i<parValues.length;i++){
			parValues[i] = limits[i][0];
		}
		
		return getOneDimIntegralByNestedIntervals( fun, minValue, initDelta, acceptedEps, parValues, limits, 0 );
		
	}
	
	private static double getOneDimIntegralByNestedIntervals(IntegrableFunction fun, double minValue, double initDelta, double acceptedEps, double[] parValues, double[][] limits, int dimension) throws DimensionException, EvaluationException{
		//find bracketing interval
		double currVal = 0;
		double lastVal = 0;
		double lastBorder = 0;
		double currDelta = initDelta;
		
		LinkedList<Interval> list = new LinkedList<Interval>();

		if(Double.isInfinite( limits[dimension][0] ) && Double.isInfinite( limits[dimension][0] )){
			throw new EvaluationException("You must specify at least one finite border for the domain of dimension "+dimension+".");
		}
		if(!fun.isUnimodal() && (Double.isInfinite( limits[dimension][0] ) || Double.isInfinite( limits[dimension][0] ))){
			throw new EvaluationException("You may specify infinite borders only for unimodal functions.");
		}
		
		if(fun.isUnimodal()){
			boolean leftFin = !Double.isInfinite( limits[dimension][0] );
			if(leftFin){
				parValues[dimension] = limits[dimension][0];
			}else{
				parValues[dimension] = limits[dimension][1];
			}
			if(dimension == fun.getDimensionOfScope() - 1){
				currVal = fun.evaluateFunction( parValues );
			}else{
				currVal = getOneDimIntegralByNestedIntervals( fun, minValue, initDelta, acceptedEps, parValues, limits, dimension+1 );
			}
			boolean desc = false;
			while( currVal > minValue || !desc ){
				lastBorder = parValues[dimension];
				lastVal = currVal;
				if(leftFin){
					parValues[dimension] += currDelta;
				}else{
					parValues[dimension] -= currDelta;
				}
				currDelta *= 1.1;
				
				if(dimension == fun.getDimensionOfScope() - 1){
					currVal = fun.evaluateFunction( parValues );
				}else{
					currVal = getOneDimIntegralByNestedIntervals( fun, minValue, initDelta, acceptedEps, parValues, limits, dimension+1 );
				}
				if(leftFin){
					list.add( new Interval(lastBorder,parValues[dimension],lastVal,currVal) );
				}else{
					list.add( new Interval(parValues[dimension],lastBorder,currVal,lastVal) );
				}
				if(currVal < lastVal){
					desc = true;
				}
			}
		}else{
			parValues[dimension] = limits[dimension][0];
			double left = fun.evaluateFunction( parValues );
			parValues[dimension] = limits[dimension][0];
			double right = fun.evaluateFunction( parValues );
			list.add( new Interval(limits[dimension][0],limits[dimension][1],left,right) );
		}

		//System.out.print(lastBorder+" ");
		//if(dimension == 0){System.out.println();}
		//refine borders
		double integral = 0.0;
		while(list.size() > 0){
			Interval curr = list.removeFirst();
			double currInt = curr.getIntegral();
			double midpoint = curr.getMidpoint();
			parValues[dimension] = midpoint;
			double midVal = 0;
			if(dimension == fun.getDimensionOfScope()-1){
				midVal = fun.evaluateFunction( parValues );
			}else{
				midVal = getOneDimIntegralByNestedIntervals( fun, minValue, initDelta, acceptedEps, parValues, limits, dimension+1 );
			}
			
			Interval i1 = new Interval(curr.leftBorder,midpoint,curr.leftValue,midVal);
			Interval i2 = new Interval(midpoint,curr.rightBorder,midVal,curr.rightValue);
			double newInt = i1.getIntegral()+i2.getIntegral();
			if( newInt < acceptedEps || Math.abs( newInt-currInt ) /Math.abs( newInt+currInt )*2.0 < acceptedEps ){
				integral += newInt;
			}else{
				list.addFirst( i2 );
				list.addFirst( i1 );
			}
		}
		return integral;

	}
	
	private static class Interval{
		double leftBorder;
		double rightBorder;
		double leftValue;
		double rightValue;

		Interval( double leftBorder, double rightBorder, double leftValue, double rightValue ) {
			this.leftBorder = leftBorder;
			this.rightBorder = rightBorder;
			this.leftValue = leftValue;
			this.rightValue = rightValue;
		}

		double getIntegral(){
			return (leftValue + rightValue)/2.0*(rightBorder-leftBorder);
		}
		
		double getMidpoint(){
			return (rightBorder+leftBorder)/2.0;
		}
		
	}
	
	/**
	 * Computes the integral under an {@link IntegrableFunction} <code>fun</code> within the limits given by {@link IntegrableFunction#getLimits()} using
	 * the trapezoid rule with a fixed width <code>delta</code> of the trapezoids.
	 * @param fun the function
	 * @param minValue minValue the minimum function value that is considered important for unimodal functions
	 * @param delta the width of the trapezoids
	 * @return the integral
	 * @throws DimensionException is thrown if the function does not accept parameter values of dimension {@link Function#getDimensionOfScope()}
	 * @throws EvaluationException  is thrown if any specified limit of a non-unimodal function is infinite or both limits are infinite for a unimodal function
	 */
	public static double getIntegralByTrapezoidRule(IntegrableFunction fun, double minValue, double delta) throws DimensionException, EvaluationException{
		int dim = fun.getDimensionOfScope();
		
		double[][] limits = fun.getLimits();
		double[] parValues = new double[dim];
		for(int i=0;i<parValues.length;i++){
			parValues[i] = limits[i][0];
		}
		
		return getOneDimIntegral( fun, minValue, delta, 0, limits, parValues );
	}
	
	
	private static double getOneDimIntegral(IntegrableFunction fun, double minValue, double delta, int dimension, double[][] limits, double[] parValues) throws DimensionException, EvaluationException{
		
		if(Double.isInfinite( limits[dimension][0] ) && Double.isInfinite( limits[dimension][0] )){
			throw new EvaluationException("You must specify at least one finite border for the domain of dimension "+dimension+".");
		}
		if(!fun.isUnimodal() && (Double.isInfinite( limits[dimension][0] ) || Double.isInfinite( limits[dimension][0] ))){
			throw new EvaluationException("You may specify infinite borders only for unimodal functions.");
		}
		
		boolean leftFin = !Double.isInfinite( limits[dimension][0] );
		
		if(leftFin){
			parValues[dimension] = limits[dimension][0];
		}else{
			parValues[dimension] = limits[dimension][1];
		}
		double funValue = 0;
		double lastValue = 0;
		
		double val = 0;
		
		boolean unimodal = fun.isUnimodal();
		boolean desc = false;

		if(dimension == fun.getDimensionOfScope() - 1){
			funValue = fun.evaluateFunction( parValues );
		}else{
			funValue = getOneDimIntegral( fun, minValue, delta, dimension+1, limits, parValues );
		}
		
		do{
			
			if(leftFin){
				parValues[dimension] += delta;
			}else{
				parValues[dimension] -= delta;
			}

			lastValue = funValue;
			if(dimension == fun.getDimensionOfScope() - 1){
				funValue = fun.evaluateFunction( parValues );
			}else{
				funValue = getOneDimIntegral( fun, minValue, delta, dimension+1, limits, parValues );
			}

			if(unimodal && (!desc) && funValue < lastValue){
				desc = true;
			}

			val += (funValue + lastValue) / 2.0 * delta;

		}while( (funValue > minValue || (!desc) ) && parValues[dimension] < limits[dimension][1] );
		
		return val;
	}
	
}
