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

package de.jstacs.classifier.assessment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SimpleParameter;
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
	 * 
	 * @see ClassifierAssessmentAssessParameterSet#ClassifierAssessmentAssessParameterSet()
	 */
	public RepeatedSubSamplingAssessParameterSet() throws UnsupportedOperationException {
		super();
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
	protected RepeatedSubSamplingAssessParameterSet( StringBuffer representation ) throws NonParsableException {
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
	 *            combination with a
	 *            {@link de.jstacs.classifier.MeasureParameters}-object to call
	 *            <code>assess( ... )</code>-methods of
	 *            {@link RepeatedSubSamplingExperiment}s. If
	 *            <code>exceptionIfMPNotComputable==true</code> an
	 *            {@link Exception} is thrown in case of a user selected measure
	 *            parameters that could not be computed.
	 * @param repeats
	 *            the number of repeats of each iteration (subsample test and
	 *            train datasets from user supplied data, train classifiers
	 *            using train datasets and test them using test datasets) of
	 *            that {@link RepeatedHoldOutExperiment} this
	 *            {@link RepeatedSubSamplingAssessParameterSet} is used with
	 * @param trainNumbers
	 *            an array containing for each class (the classifiers to be
	 *            assessed are capable to distinguish) the number of elements
	 *            the subsampled train datasets should contain
	 * @param testNumbers
	 *            an array containing for each class (the classifiers to be
	 *            assessed are capable to distinguish) the number of elements
	 *            the subsampled test datasets should contain
	 * 
	 * @throws IllegalValueException
	 *             in case of out-of-range or invalid given parameters
	 * 
	 * @see ClassifierAssessmentAssessParameterSet#ClassifierAssessmentAssessParameterSet(int,
	 *      boolean)
	 */
	public RepeatedSubSamplingAssessParameterSet( int elementLength, boolean exceptionIfMPNotComputable, int repeats, int[] trainNumbers,
													int[] testNumbers ) throws IllegalValueException {
		super( elementLength, exceptionIfMPNotComputable );

		this.parameters.get( 2 ).setValue( new Integer( repeats ) );

		ParameterSet[] tempPSA = new ParameterSet[trainNumbers.length];
		for( int i = 0; i < tempPSA.length; tempPSA[i] = getParameterSetContainingASingleIntValue( trainNumbers[i++], "training" ) );

		( (ExpandableParameterSet)( ( (ParameterSetContainer)( this.parameters.get( 3 ) ) ).getValue() ) ).replaceContentWith( tempPSA );

		ParameterSet[] tempPSA2 = new ParameterSet[testNumbers.length];
		for( int i = 0; i < tempPSA2.length; tempPSA2[i] = getParameterSetContainingASingleIntValue( testNumbers[i++], "testing" ) );

		( (ExpandableParameterSet)( ( (ParameterSetContainer)( this.parameters.get( 4 ) ) ).getValue() ) ).replaceContentWith( tempPSA2 );

	}

	//	**********************
	//	member methods
	//	**********************

	/**
	 * Creates a new {@link ParameterSet} containing a single <code>int</code>-
	 * {@link SimpleParameter}. This {@link ParameterSet} is used as a part of
	 * the {@link ExpandableParameterSet} that contains the test data for a
	 * specific class.
	 */
	private ParameterSet getParameterSetContainingASingleIntValue( int num, final String train_test ) throws IllegalValueException {

		ParameterSet ret = new ParameterSet() {

			@Override
			protected void loadParameters() throws Exception {
				initParameterList( 1 );
				this.parameters.add( new SimpleParameter( DataType.INT,
						"number",
						"Defines a number of elements of data used as " + train_test
								+ " items (class-specific) during a SubSamplingAssessment",
						true,
						new NumberValidator<Integer>( 1, Integer.MAX_VALUE ) ) );
			}
		};

		ret.getParameterAt( 0 ).setValue( new Integer( num ) );

		return ret;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifier.assessment.ClassifierAssessmentAssessParameterSet#initializeMyParametersArrayList()
	 */
	@Override
	protected void initializeMyParametersArrayList() {
		initParameterList( 6 );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifier.assessment.ClassifierAssessmentAssessParameterSet#loadParameters()
	 */
	@Override
	protected void loadParameters() throws Exception {
		super.loadParameters();

		//2-k
		this.parameters.add( new SimpleParameter( DataType.INT,
				"repeats",
				"Determines how often the procedure of " + "train/test classifers with random created "
						+ "train- and test-data should be repeated.",
				true,
				new NumberValidator<Integer>( 1, Integer.MAX_VALUE ) ) );

		//3-percents
		this.parameters.add( new ParameterSetContainer( "trainDataNumbers",
				"A RepeatedSubSamplingExperiment subsamples " + "the used train- and test-datasets from the given data (for each class) "
						+ "in each iteration. This ParameterSetcontainer "
						+ "contains an ExpandableParameterSet that contains for each class "
						+ "the number of the items (for each class), that should be subsampled "
						+ "and used as training-data.",
				new ExpandableParameterSet( getParameterSetContainingASingleIntValue( 1, "training" ),
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
				new ExpandableParameterSet( getParameterSetContainingASingleIntValue( 1, "testing" ),
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
		return ( (Integer)( this.getParameterAt( 2 ).getValue() ) ).intValue();
	}

	/**
	 * Returns an array containing the number of elements the subsampled (train
	 * | test) datasets should consist of.
	 * 
	 * @param train_case
	 *            if <code>true</code> then <code>(train | test)=train</code>
	 *            else <code>(train | test)=test</code>
	 * 
	 * @return an array class-wise containing the number of elements the
	 *         subsampled (train | test) datasets should consist of
	 */
	public int[] getTrain_TestNumbers( boolean train_case ) {

		int pos;
		if( train_case ) {
			pos = 3;
		} else {
			pos = 4;
		}

		ExpandableParameterSet tempEPS = (ExpandableParameterSet)( this.getParameterAt( pos ).getValue() );

		int[] ret = new int[tempEPS.getNumberOfParameters()];

		for( int i = 0; i < ret.length; i++ ) {
			//holy shit, that's really unsexy
			ret[i] = ( (Integer)( ( (ParameterSet)( tempEPS.getParameterAt( i ).getValue() ) ).getParameterAt( 0 ).getValue() ) ).intValue();
		};

		return ret;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifier.assessment.ClassifierAssessmentAssessParameterSet#getAnnotation()
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
