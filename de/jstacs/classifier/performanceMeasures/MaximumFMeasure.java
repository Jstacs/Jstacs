package de.jstacs.classifier.performanceMeasures;


import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.NumberValidator;

//http://en.wikipedia.org/wiki/F1_score

/**
 * 
 * @author Jens Keilwagen
 */
public class MaximumFMeasure extends MaximumNumericalTwoClassMeasure {

	/**
	 * Constructs a new instance of the performance measure {@link MaximumFMeasure} with empty parameters.
	 * 
	 * @throws Exception if the internal parameters can not be created
	 */
	public MaximumFMeasure() throws Exception {
		loadParameters();
	}
	
	/**
	 * Constructs a new instance of the performance measure {@link MaximumFMeasure} with given <code>beta</code>.
	 * 
	 * @param beta the beta for which the maximum F-measure should be computed
	 * 
	 * @throws Exception if the internal parameters can not be created or the value can not be set
	 */
	public MaximumFMeasure( double beta ) throws Exception {
		this();
		parameters.get(0).setValue( beta );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link MaximumFMeasure} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link MaximumFMeasure} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public MaximumFMeasure( StringBuffer xml ) throws NonParsableException {
		super(xml);
	}

	@Override
	protected String getMeasureName() {
		return "F-Measure";
	}
	
	@Override
	protected String getSpecificName() {
		if (parameters == null) {
			try {
				loadParameters();
				if (ranged) {
					replaceParametersWithRangedInstance();
				}
			} catch (Exception e) {
				e.printStackTrace();
				return null;
			}
		}
		return getMeasureName() + " with beta=" + parameters.get("beta").getValue();
	}

	@Override
	protected void loadParameters() throws Exception {
		initParameterList( 1 );
		parameters.add( new SimpleParameter( DataType.DOUBLE, "beta", "the beta defining the F measure", true, new NumberValidator<Double>(0d,Double.POSITIVE_INFINITY), 1d ) );
	}

	@Override
	protected double getMeasure(double tp, double fp, double fn, double tn) {
		double b = (Double) parameters.get( "beta").getValue();
		b *= b;
		double precision = tp / (tp+fp);
		double recall = tp / (tp+fn);
		return (1+b) * precision*recall / (b*precision + recall);
	}
}