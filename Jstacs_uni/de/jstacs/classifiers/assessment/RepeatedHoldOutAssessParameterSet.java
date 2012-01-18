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
import de.jstacs.data.DataSet.PartitionMethod;
import de.jstacs.io.NonParsableException;
import de.jstacs.parameters.EnumParameter;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.Result;

/**
 * This class implements a {@link ClassifierAssessmentAssessParameterSet} that
 * must be used to call method <code>assess( ... )</code> of a
 * {@link RepeatedHoldOutExperiment}. It contains user specific parameters
 * necessary for a run of a {@link RepeatedHoldOutExperiment}.
 * 
 * @author Andre Gohr (bioinf (nospam:.) ag (nospam:@) googlemail (nospam:.)
 *         com)
 * 
 */
public class RepeatedHoldOutAssessParameterSet extends ClassifierAssessmentAssessParameterSet {

	/**
	 * Constructs a new {@link RepeatedHoldOutAssessParameterSet} with empty
	 * parameter values. This constructor should only be used to create
	 * &quot;filled&quot; {@link RepeatedHoldOutAssessParameterSet}s, i.e. to
	 * create {@link RepeatedHoldOutAssessParameterSet}s from a set of values
	 * and not to fill it from the platform user interface.
	 * @throws ParameterException if the parameters could not be created
	 * @throws CloneNotSupportedException if the parameter for the percentages could not be created
	 * 
	 * @see ClassifierAssessmentAssessParameterSet#ClassifierAssessmentAssessParameterSet()
	 */
	public RepeatedHoldOutAssessParameterSet() throws ParameterException, CloneNotSupportedException {
		super();
		addParameters();
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link RepeatedHoldOutAssessParameterSet} out of its XML
	 * representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * @throws NonParsableException
	 *             if the {@link RepeatedHoldOutAssessParameterSet} could not be
	 *             reconstructed out of the XML representation (the
	 *             {@link StringBuffer} <code>representation</code> could not be
	 *             parsed)
	 * 
	 * @see ClassifierAssessmentAssessParameterSet#ClassifierAssessmentAssessParameterSet(StringBuffer)
	 * @see de.jstacs.Storable
	 */
	public RepeatedHoldOutAssessParameterSet( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}

	/**
	 * Constructs a new {@link RepeatedHoldOutAssessParameterSet} with given
	 * parameter values.
	 * 
	 * @param dataSplitMethod
	 *            defines the method used to split user supplied data into
	 *            <code>k</code> mutually exclusive random-splits (available
	 *            options are:
	 *            {@link PartitionMethod#PARTITION_BY_NUMBER_OF_ELEMENTS} and
	 *            {@link PartitionMethod#PARTITION_BY_NUMBER_OF_SYMBOLS})
	 * @param elementLength
	 *            defines the length of elements (sequences) the classifiers to
	 *            be assessed are able to classify
	 * @param exceptionIfMPNotComputable
	 *            a {@link RepeatedHoldOutAssessParameterSet} is used in
	 *            combination with an
	 *            {@link de.jstacs.classifiers.performanceMeasures.AbstractPerformanceMeasure}-object to call
	 *            <code>assess( ... )</code>-methods of
	 *            {@link RepeatedHoldOutExperiment}s. If
	 *            <code>exceptionIfMPNotComputable==true</code> an exception is
	 *            thrown in case of a user selected measure parameters that
	 *            could not be computed.
	 * 
	 * @param repeats
	 *            the number of repeats of each iteration (mutually exclusive,
	 *            randomly split data to obtain test and train datasets, train
	 *            classifiers using train datasets and test them using test
	 *            datasets) of that {@link RepeatedHoldOutExperiment} this
	 *            {@link RepeatedHoldOutAssessParameterSet} is used with
	 * 
	 * @param percents
	 *            this array contains class-wise the percentage of the user
	 *            supplied data that should be used as test data in each
	 *            iteration of that {@link RepeatedHoldOutExperiment} this
	 *            {@link RepeatedHoldOutAssessParameterSet} is used with
	 * 
	 * @throws ParameterException if the parameters could not be created
	 * @throws CloneNotSupportedException if the parameter for the percentages could not be created
	 * 
	 * @see ClassifierAssessmentAssessParameterSet#ClassifierAssessmentAssessParameterSet(int,
	 *      boolean)
	 * @see de.jstacs.data.DataSet.PartitionMethod
	 */
	public RepeatedHoldOutAssessParameterSet( PartitionMethod dataSplitMethod, int elementLength, boolean exceptionIfMPNotComputable,
												int repeats, double[] percents ) throws ParameterException, CloneNotSupportedException {
		super( elementLength, exceptionIfMPNotComputable );
		addParameters();
		this.parameters.get( 2 ).setValue( new Integer( repeats ) );

		ParameterSet[] tempPSA = new ParameterSet[percents.length];
		for( int i = 0; i < tempPSA.length; tempPSA[i] = getParameterSetContainingASingleDoubleValue( percents[i++] ) );

		( (ExpandableParameterSet)( ( (ParameterSetContainer)( this.parameters.get( 3 ) ) ).getValue() ) ).replaceContentWith( tempPSA );

		( this.parameters.get( 4 ) ).setValue( dataSplitMethod );

	}

	//	**********************
	//	member methods
	//	**********************

	/**
	 * Creates a new {@link ParameterSet} containing a single
	 * <code>double</code>-<code>SimpleParameter</code>. This
	 * {@link ParameterSet} is used as a part of the
	 * {@link ExpandableParameterSet} that contains the test data percent for a
	 * specific class. <br>
	 * 
	 * @param percent
	 *            the <code>double</code>-value to be contained in the returned
	 *            {@link ParameterSet}. If <code>percent==Double.NaN</code> no
	 *            values are contained in the returned {@link ParameterSet}.
	 *            (The {@link SimpleParameter} contained in the returned
	 *            {@link ParameterSet} contains no value).
	 * 
	 * @throws IllegalValueException
	 *             if something went wrong
	 * @throws DatatypeNotValidException 
	 */
	private ParameterSet getParameterSetContainingASingleDoubleValue( double percent ) throws IllegalValueException, DatatypeNotValidException {

		ParameterSet ret = new SimpleParameterSet( new SimpleParameter( DataType.DOUBLE,
						"percent",
						"Defines the percentage of the entire given data (for a " + "specific class) should be used as test-data in a "
								+ "RepeatedHoldOutExperiment.",
						true,
						new NumberValidator<Double>( 0d, 1d ) ) );
		if( !Double.isNaN( percent ) ) {
			ret.getParameterAt( 0 ).setValue( new Double( percent ) );
		}

		return ret;
	}
	
	private void addParameters() throws ParameterException, CloneNotSupportedException {
		//2-k
		this.parameters.add( new SimpleParameter( DataType.INT,
				"repeats",
				"Determines how often the procedure of " + "train/test classifers using random created "
						+ "train- and test-data should be repeated.",
				true,
				new NumberValidator<Integer>( 1, Integer.MAX_VALUE ) ) );

		//3-percents
		this.parameters.add( new ParameterSetContainer( "testDataPercentage",
				"A RepeatedHoldOutExperiment splits " + "the given data (for each class) in each iteration into "
						+ "a test-part and a train-part. This ParameterSetcontainer "
						+ "contains an ExpandableParameterSet that contains for each class "
						+ "the percent of the entire data (for each class), that should be used "
						+ "as test-data. (1-percent) defines the percent of train-data.",
				new ExpandableParameterSet( getParameterSetContainingASingleDoubleValue( Double.NaN ),
						"percent",
						"At pos i in this Expandable ParameterSet " + "defines the percent of all given data, "
								+ "that should be used as test-data for class "
								+ "i in a RepeatedHoldOutExperiment." )//new ExpandableParameterSet
		)//new ParameterSetContainer
		);//this.parameters.add(...)

		//4-dataSplitMethod
		this.parameters.add( new EnumParameter( PartitionMethod.class, "The method used to compute the percentages of the partitions", true ) );
	}

	/**
	 * Returns the repeats defined by this
	 * {@link RepeatedHoldOutAssessParameterSet} (repeats define how many
	 * iterations (train and test classifiers) of that
	 * {@link RepeatedHoldOutExperiment} this
	 * {@link RepeatedHoldOutAssessParameterSet} is used with are performed).
	 * 
	 * @return the repeats defined by this
	 *         {@link RepeatedHoldOutAssessParameterSet} (repeats define how
	 *         many iterations (train and test classifiers) of that
	 *         {@link RepeatedHoldOutExperiment} this
	 *         {@link RepeatedHoldOutAssessParameterSet} is used with are
	 *         performed)
	 */
	public int getRepeats() {
		return ( (Integer)( this.getParameterAt( 2 ).getValue() ) ).intValue();
	}

	/**
	 * Returns an array containing for each class the percentage of user
	 * supplied data that is used in each iteration as test dataset.
	 * 
	 * @return an array containing for each class the percentage of user
	 *         supplied data that is used in each iteration as test dataset
	 */
	public double[] getPercents() {

		ExpandableParameterSet tempEPS = (ExpandableParameterSet)( this.getParameterAt( 3 ).getValue() );

		double[] ret = new double[tempEPS.getNumberOfParameters()];

		for( int i = 0; i < ret.length; i++ ) {
			//holy shit, that's realy unsexy
			ret[i] = ( (Double)( ( (ParameterSet)( tempEPS.getParameterAt( i ).getValue() ) ).getParameterAt( 0 ).getValue() ) ).doubleValue();
		}

		return ret;
	}

	/**
	 * Returns the {@link PartitionMethod} defining how the mutually exclusive
	 * random-splits of user supplied data are generated.
	 * 
	 * @return the {@link PartitionMethod} defining how the mutually exclusive
	 *         random-splits of user supplied data are generated
	 * 
	 * @see de.jstacs.data.DataSet.PartitionMethod
	 */
	public PartitionMethod getDataSplitMethod() {
		return (PartitionMethod)( (EnumParameter)getParameterAt( 4 ) ).getValue();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.assessment.ClassifierAssessmentAssessParameterSet#getAnnotation()
	 */
	@Override
	public Collection<Result> getAnnotation() {
		ArrayList<Result> l = new ArrayList<Result>( 3 );
		l.add( new NumericalResult( "repeats", "The number of iterations", getRepeats() ) );
		l.add( new CategoricalResult( "percentage",
				"The percentage of the entire data (for each class), that was used in an assessment",
				Arrays.toString( getPercents() ) ) );
		l.add( new CategoricalResult( "dataSplitMethod",
				"Describes how data should be splitted in ClassiefierAssessment.evaluateClassifier())",
				getDataSplitMethod().name() ) );
		return l;
	}

}
