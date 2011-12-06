package de.jstacs.classifier.performanceMeasures;

import de.jstacs.NonParsableException;
import de.jstacs.results.ResultSet;

/**
 * This class implements the performance measure confusion matrix.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class ConfusionMatrix extends AbstractPerformanceMeasure {

	/**
	 * Constructs a new instance of the performance measure {@link ConfusionMatrix}.
	 */
	public ConfusionMatrix() {
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link ConfusionMatrix} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link ConfusionMatrix} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public ConfusionMatrix(StringBuffer xml) throws NonParsableException {
		super(xml);
	}

	@Override
	public ResultSet compute(double[] sortedScoresClass0, double[] sortedScoresClass1) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public ResultSet compute(double[][][] classSpecificScores) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int getAllowedNumberOfClasses() {
		return 0;
	}

	@Override
	public String getName() {
		return "Confusion matrix";
	}

	@Override
	protected void loadParameters() throws Exception {
		initParameterList();
	}
}