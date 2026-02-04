package de.jstacs.classifiers.differentiableSequenceScoreBased;

import java.io.BufferedWriter;
import java.io.IOException;

import de.jstacs.algorithms.optimization.DimensionException;
import de.jstacs.algorithms.optimization.EvaluationException;
import de.jstacs.data.DataSet;

/**
 * This function allows to save the parameters in a tab-seperated file while a numerical optimization.
 * Each time the method {@link ParameterSavingFunction#evaluateGradientOfFunction(double[])} is called a line is added to the file.
 * 
 * @author Jens Keilwagen
 */
public class ParameterSavingFunction extends OptimizableFunction {

	OptimizableFunction f;
	BufferedWriter w;
	
	/**
	 * Constructor for a ParameterSavingFunction.
	 * 
	 * @param f the function to be optimized
	 * @param w a BufferedWriter to add the parameters
	 */
	public ParameterSavingFunction( OptimizableFunction f, BufferedWriter w ) {
		this.f = f;
		this.w = w;
	}
	
	@Override
	public double evaluateFunction(double[] x) throws DimensionException, EvaluationException {
		return f.evaluateFunction(x);
	}

	@Override
	public int getDimensionOfScope() {
		return f.getDimensionOfScope();
	}

	@Override
	public double[] getParameters(KindOfParameter kind) throws Exception {
		return f.getParameters(kind);
	}

	@Override
	public void setParams(double[] current) throws DimensionException {
		f.setParams(current);
	}

	@Override
	public void reset() throws Exception {
		f.reset();
	}

	@Override
	public DataSet[] getData() {
		return f.getData();
	}

	@Override
	public double[][] getSequenceWeights() {
		return f.getSequenceWeights();
	}

	@Override
	public void setDataAndWeights(DataSet[] data, double[][] weights) throws IllegalArgumentException {
		f.setDataAndWeights(data, weights);
	}

	@Override
	public double[] evaluateGradientOfFunction(double[] x) throws DimensionException, EvaluationException {
		try {
			w.append(""+x[0]);
			for( int i = 1; i < x.length; i++ ) {
				w.append("\t" + x[i]);
			}
			w.newLine();
			w.flush();
		} catch( IOException io ) {
			io.printStackTrace();
		}
		return f.evaluateGradientOfFunction(x);
	}
}