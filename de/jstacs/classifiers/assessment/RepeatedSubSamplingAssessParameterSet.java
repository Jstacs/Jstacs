/*
 * This file is part of Jstacs.
 * 
 * Jstacs is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 * 
 * Jstacs is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.classifiers.assessment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;

import de.jstacs.DataType;
import de.jstacs.io.NonParsableException;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.Result;

/**
 * This class implements a {@link ClassifierAssessmentAssessParameterSet} that
 * must be used to call method <code>assess( ... )</code> of a
 * {@link RepeatedSubSamplingExperiment}. It contains user specific parameters
 * necessary for a run of a {@link RepeatedSubSamplingExperiment}.
 * 
 * @author Andre Gohr (bioinf (nospam:.) ag (nospam:@) googlemail (nospam:.)
 *         com)
 * 
 */
public class RepeatedSubSamplingAssessParameterSet extends ClassifierAssessmentAssessParameterSet {

	/**
	 * Constructs a new {@link RepeatedSubSamplingAssessParameterSet} with empty
	 * parameter values. This constructor should only be used to create
	 * &quot;filled&quot; {@link RepeatedSubSamplingAssessParameterSet}s, i.e.
	 * to create {@link RepeatedSubSamplingAssessParameterSet}s from a set of
	 * values and not to fill it from the platform user interface.
	 * 
	 * @throws UnsupportedOperationException
	 *             if the {@link RepeatedSubSamplingAssessParameterSet} could
	 *             not be constructed or the parameters could not be loaded
	 * @throws CloneNotSupportedException if the parameter for the subsampled sequences could not be created
	 * @throws IllegalValueException if the parameters could not be filled with the default values
	 * 
	 * @see ClassifierAssessmentAssessParameterSet#ClassifierAssessmentAssessParameterSet()
	 */
	public RepeatedSubSamplingAssessParameterSet() throws UnsupportedOperationException, CloneNotSupportedException, IllegalValueException {
		super();
		addParameters();
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link RepeatedSubSamplingAssessParameterSet} out of its XML
	 * representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * @throws NonParsableException
	 *             if the {@link RepeatedSubSamplingAssessParameterSet} could
	 *             not be reconstructed out of the XML representation (the
	 *             {@link StringBuffer} <code>representation</code> could not be
	 *             parsed)
	 * 
	 * @see ClassifierAssessmentAssessParameterSet#ClassifierAssessmentAssessParameterSet(StringBuffer)
	 */
	public RepeatedSubSamplingAssessParameterSet( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}

	/**
	 * Constructs a new {@link RepeatedSubSamplingAssessParameterSet} with given
	 * parameter values.
	 * 
	 * @param elementLength
	 *            defines the length of elements (sequences) the classifiers to
	 *            be assessed are able to classify
	 * @param exceptionIfMPNotComputable
	 *            a {@link RepeatedSubSamplingAssessParameterSet} is used in
	 *            combination with an
	 *            {@link de.jstacs.classifiers.performanceMeasures.AbstractPerformanceMeasure}-object to call
	 *            <code>assess( ... )</code>-methods of
	 *            {@link RepeatedSubSamplingExperiment}s. If
	 *            <code>exceptionIfMPNotComputable==true</code> an
	 *            {@link Exception} is thrown in case of a user selected measure
	 *            parameters that could not be computed.
	 * @param repeats
	 *            the number of repeats of each iteration (subsample test and
	 *            train data sets from user supplied data, train classifiers
	 *            using train data sets and test them using test data sets) of
	 *            that {@link RepeatedHoldOutExperiment} this
	 *            {@link RepeatedSubSamplingAssessParameterSet} is used with
	 * @param trainNumbers
	 *            an array containing for each class (the classifiers to be
	 *            assessed are capable to distinguish) the number of elements
	 *            the subsampled train data sets should contain
	 * @param testNumbers
	 *            an array containing for each class (the classifiers to be
	 *            assessed are capable to distinguish) the number of elements
	 *            the subsampled test data sets should contain
	 * 
	 * @throws IllegalValueException
	 *             in case of out-of-range or invalid given parameters
	 * @throws CloneNotSupportedException if the parameter for the subsampled sequences could not be created
	 * 
	 * @see ClassifierAssessmentAssessParameterSet#ClassifierAssessmentAssessParameterSet(int,
	 *      boolean)
	 */
	public RepeatedSubSamplingAssessParameterSet( int elementLength, boolean exceptionIfMPNotComputable, int repeats, double[] trainNumbers,
													double[] testNumbers ) throws IllegalValueException, CloneNotSupportedException {
		super( elementLength, exceptionIfMPNotComputable );
		addParameters();
		
		this.parameters.get( "repeats" ).setValue( repeats );

		ParameterSet[] tempPSA = new ParameterSet[trainNumbers.length];
		for( int i = 0; i < tempPSA.length; tempPSA[i] = getParameterSetContainingASingleDoubleValue( trainNumbers[i++], "training" ) );

		( (ExpandableParameterSet)( ( (ParameterSetContainer)( this.parameters.get( "trainDataNumbers" ) ) ).getValue() ) ).replaceContentWith( tempPSA );

		ParameterSet[] tempPSA2 = new ParameterSet[testNumbers.length];
		for( int i = 0; i < tempPSA2.length; tempPSA2[i] = getParameterSetContainingASingleDoubleValue( testNumbers[i++], "testing" ) );

		( (ExpandableParameterSet)( ( (ParameterSetContainer)( this.parameters.get( "testDataNumbers" ) ) ).getValue() ) ).replaceContentWith( tempPSA2 );

	}

	//	**********************
	//	member methods
	//	**********************

	/**
	 * Creates a new {@link ParameterSet} containing a single <code>double</code>-
	 * {@link SimpleParameter}. This {@link ParameterSet} is used as a part of
	 * the {@link ExpandableParameterSet} that contains the test data for a
	 * specific class.
	 * @throws DatatypeNotValidException 
	 */
	private ParameterSet getParameterSetContainingASingleDoubleValue( double num, final String train_test ) throws IllegalValueException {

		ParameterSet ret;
		try {
			ret = new SimpleParameterSet(new SimpleParameter( DataType.DOUBLE,
							"number",
							"Defines a number of elements of data used as " + train_test
									+ " items (class-specific) during a SubSamplingAssessment",
							true,
							new NumberValidator<Double>( 0d, Double.MAX_VALUE ) ) );
		} catch ( DatatypeNotValidException doesnothappen ) {
			throw new RuntimeException( doesnothappen );
		}

		ret.getParameterAt( 0 ).setValue( num );

		return ret;
	}

	private void addParameters() throws CloneNotSupportedException, IllegalValueException {
		//2-k
		try {
			this.parameters.add( new SimpleParameter( DataType.INT,
					"repeats",
					"Determines how often the procedure of " + "train/test classifers with random created "
							+ "train- and test-data should be repeated.",
					true,
					new NumberValidator<Integer>( 1, Integer.MAX_VALUE ) ) );
		} catch ( DatatypeNotValidException doesnothappen ) {
			throw new RuntimeException( doesnothappen );
		}

		//3-percents
		this.parameters.add( new ParameterSetContainer( "trainDataNumbers",
				"A RepeatedSubSamplingExperiment subsamples " + "the used train- and test-datasets from the given data (for each class) "
						+ "in each iteration. This ParameterSetcontainer "
						+ "contains an ExpandableParameterSet that contains for each class "
						+ "the number of the items (for each class), that should be subsampled "
						+ "and used as training-data.",
				new ExpandableParameterSet( getParameterSetContainingASingleDoubleValue( 1, "training" ),
						"number",
						"At pos i in this Expandable ParameterSet " + "defines the number of subsampled items, "
								+ "that should be used as train-data for class "
								+ "i in a RepeatedSubSamplingExperiment." )//new ExpandableParameterSet
		)//new ParameterSetContainer
		);//this.parameters.add(...)

		//4-percents
		this.parameters.add( new ParameterSetContainer( "testDataNumbers",
				"A RepeatedSubSamplingExperiment subsamples " + "the used train and test datasets from the given data (for each class) "
						+ "in each iteration. This ParameterSetcontainer "
						+ "contains an ExpandableParameterSet that contains for each class "
						+ "the number of the items (for each class), that should be subsampled "
						+ "and used as test-data.",
				new ExpandableParameterSet( getParameterSetContainingASingleDoubleValue( 1, "testing" ),
						"number",
						"At pos i in this Expandable ParameterSet " + "defines the number of subsampled items, "
								+ "that should be used as test-data for class "
								+ "i in a RepeatedSubSamplingExperiment." )//new ExpandableParameterSet
		)//new ParameterSetContainer
		);//this.parameters.add(...)
	}

	/**
	 * Returns the repeats defined by this
	 * {@link RepeatedSubSamplingAssessParameterSet} (repeats defines how many
	 * iterations (train and test classifiers) of that
	 * {@link RepeatedSubSamplingExperiment} this
	 * {@link RepeatedSubSamplingAssessParameterSet} is used with are
	 * performed).
	 * 
	 * @return the repeats defined by this
	 *         {@link RepeatedSubSamplingAssessParameterSet} (repeats defines
	 *         how many iterations (train and test classifiers) of that
	 *         {@link RepeatedSubSamplingExperiment} this
	 *         {@link RepeatedSubSamplingAssessParameterSet} is used with are
	 *         performed)
	 */
	public int getRepeats() {
		return ( (Integer)( this.getParameterForName( "repeats" ).getValue() ) ).intValue();
	}

	/**
	 * Returns an array containing the number of elements the subsampled (train
	 * | test) data sets should consist of.
	 * 
	 * @param train_case
	 *            if <code>true</code> then <code>(train | test)=train</code>
	 *            else <code>(train | test)=test</code>
	 * 
	 * @return an array class-wise containing the number of elements the
	 *         subsampled (train | test) data sets should consist of
	 */
	public double[] getTrain_TestNumbers( boolean train_case ) {

		String pos;
		if( train_case ) {
			pos = "trainDataNumbers";
		} else {
			pos = "testDataNumbers";
		}

		ExpandableParameterSet tempEPS = (ExpandableParameterSet)( this.getParameterForName( pos ).getValue() );

		double[] ret = new double[tempEPS.getNumberOfParameters()];

		for( int i = 0; i < ret.length; i++ ) {
			ret[i] = (Double)( ( (ParameterSet)( tempEPS.getParameterAt( i ).getValue() ) ).getParameterAt( 0 ).getValue() );
		};

		return ret;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.assessment.ClassifierAssessmentAssessParameterSet#getAnnotation()
	 */
	@Override
	public Collection<Result> getAnnotation() {
		ArrayList<Result> l = new ArrayList<Result>( 3 );
		l.add( new NumericalResult( "repeats", "The number of iterations", getRepeats() ) );
		l.add( new CategoricalResult( "number of elements for train",
				"The number of the items (for each class), that was subsampled and used as training-data.",
				Arrays.toString( getTrain_TestNumbers( true ) ) ) );
		l.add( new CategoricalResult( "number of elements for assessment",
				"The number of the items (for each class), that should be subsampled and used as assessment-data.",
				Arrays.toString( getTrain_TestNumbers( false ) ) ) );
		return l;
	}

}
