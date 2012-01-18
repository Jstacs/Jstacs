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
import java.util.Collection;

import de.jstacs.DataType;
import de.jstacs.data.DataSet.PartitionMethod;
import de.jstacs.io.NonParsableException;
import de.jstacs.parameters.EnumParameter;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.Result;

/**
 * This class implements a {@link ClassifierAssessmentAssessParameterSet} that
 * must be used to call the method <code>assess( ... )</code> of a
 * {@link Sampled_RepeatedHoldOutExperiment}. It contains user specific
 * parameters necessary for a run of a {@link Sampled_RepeatedHoldOutExperiment}
 * .
 * 
 * @author Jens Keilwagen
 */
public class Sampled_RepeatedHoldOutAssessParameterSet extends ClassifierAssessmentAssessParameterSet {

	/**
	 * Constructs a new {@link Sampled_RepeatedHoldOutAssessParameterSet} with
	 * empty parameter values. This constructor should only be used to create
	 * &quot;filled&quot; {@link Sampled_RepeatedHoldOutAssessParameterSet}s,
	 * i.e. to create {@link Sampled_RepeatedHoldOutAssessParameterSet}s from a
	 * set of values and not to fill it from the platform user-interface.
	 * 
	 * @throws UnsupportedOperationException
	 *             if the {@link Sampled_RepeatedHoldOutAssessParameterSet}
	 *             could not be constructed or the parameters could not be
	 *             loaded
	 * @throws ParameterException 
	 * 
	 * @see ClassifierAssessmentAssessParameterSet#ClassifierAssessmentAssessParameterSet()
	 */
	public Sampled_RepeatedHoldOutAssessParameterSet() throws UnsupportedOperationException, ParameterException {
		super();
		addParameters();
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link Sampled_RepeatedHoldOutAssessParameterSet} out of its
	 * XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * @throws NonParsableException
	 *             if the {@link Sampled_RepeatedHoldOutAssessParameterSet}
	 *             could not be reconstructed out of the XML representation (the
	 *             {@link StringBuffer} <code>representation</code> could not be
	 *             parsed)
	 * 
	 * @see ClassifierAssessmentAssessParameterSet#ClassifierAssessmentAssessParameterSet(StringBuffer)
	 */
	public Sampled_RepeatedHoldOutAssessParameterSet( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}

	/**
	 * Constructs a new {@link Sampled_RepeatedHoldOutAssessParameterSet} with
	 * given parameter values.
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
	 *            a {@link Sampled_RepeatedHoldOutAssessParameterSet} is used in
	 *            combination with an
	 *            {@link de.jstacs.classifier.performanceMeasures.AbstractPerformanceMeasure}-object to call
	 *            <code>assess( ... )</code>-methods of
	 *            {@link Sampled_RepeatedHoldOutExperiment}s. If
	 *            <code>exceptionIfMPNotComputable==true</code> an
	 *            {@link Exception} is thrown in case of a user selected measure
	 *            parameters that could not be computed.
	 * @param repeats
	 *            the number of repeats of each iteration (mutually exclusive,
	 *            randomly split data to obtain test and train datasets, train
	 *            classifiers using train datasets and test them using test
	 *            datasets) of that {@link RepeatedHoldOutExperiment} this
	 *            {@link Sampled_RepeatedHoldOutAssessParameterSet} is used with
	 * @param referenceClass
	 *            the index of the class for which the complete data set is
	 *            used, typically this should be the smallest data set (to meet
	 *            all constraints)
	 * @param percentage
	 *            the percentage of the <code>referenceClass</code> data that
	 *            should be used as test data in each iteration
	 * @param sameLength
	 *            if <code>true</code> for test and train dataset the sequences
	 *            of the non-reference classes have the same length as the
	 *            corresponding sequence of the reference class
	 * @throws ParameterException 
	 * 
	 * @see de.jstacs.data.DataSet.PartitionMethod
	 */
	public Sampled_RepeatedHoldOutAssessParameterSet( PartitionMethod dataSplitMethod, int elementLength,
														boolean exceptionIfMPNotComputable, int repeats, int referenceClass,
														double percentage, boolean sameLength ) throws ParameterException {
		super( elementLength, exceptionIfMPNotComputable );
		addParameters();
		this.parameters.get( 2 ).setValue( repeats );
		this.parameters.get( 3 ).setValue( referenceClass );
		this.parameters.get( 4 ).setValue( percentage );
		this.parameters.get( 5 ).setValue( dataSplitMethod );
		this.parameters.get( 6 ).setValue( sameLength );
	}

	private void addParameters() throws ParameterException {
		// 2-k
		this.parameters.add( new SimpleParameter( DataType.INT,
				"repeats",
				"Determines how often the procedure of train/test classifers using random created train- and test-data should be repeated.",
				true,
				new NumberValidator<Integer>( 1, Integer.MAX_VALUE ) ) );

		// 3-reference class
		this.parameters.add( new SimpleParameter( DataType.INT,
				"reference class",
				"the index of the class for which the complete data set is partitioned in an iteration; typically this should be the smallest data set (to meet all constraints)",
				true,
				new NumberValidator<Integer>( 0, Integer.MAX_VALUE ) ) );

		// 4-percentage
		this.parameters.add( new SimpleParameter( DataType.DOUBLE,
				"testDataPercentage",
				"the percentage of of the foreground data set  that is used in each iterations as test-data-set",
				true,
				new NumberValidator<Double>( 0d, 1d ) ) );

		// 5-dataSplitMethod
		this.parameters.add( new EnumParameter( PartitionMethod.class, "The method used to compute the percentages of the partitions", true ) );

		// 6-sameLength
		this.parameters.add( new SimpleParameter( DataType.BOOLEAN,
				"sameLength",
				"if true then for test respectively train data set the sequences of the non-reference" + "classes have the same length as the corresponding sequence of the reference class",
				true ) );
	}

	/**
	 * Returns the repeats defined by this
	 * {@link Sampled_RepeatedHoldOutAssessParameterSet} (repeats defines how
	 * many iterations (train and test classifiers) of that
	 * {@link Sampled_RepeatedHoldOutExperiment} this
	 * {@link Sampled_RepeatedHoldOutAssessParameterSet} is used with are
	 * performed).
	 * 
	 * @return the repeats defined by this
	 *         {@link Sampled_RepeatedHoldOutAssessParameterSet} (repeats
	 *         defines how many iterations (train and test classifiers) of that
	 *         {@link Sampled_RepeatedHoldOutExperiment} this
	 *         {@link Sampled_RepeatedHoldOutAssessParameterSet} is used with
	 *         are performed)
	 */
	public int getRepeats() {
		return (Integer)( this.getParameterAt( 2 ).getValue() );
	}

	/**
	 * Returns the index of the reference class.
	 * 
	 * @return the index of the reference class
	 */
	public int getReferenceClass() {
		return (Integer)( this.getParameterAt( 3 ).getValue() );
	}

	/**
	 * Returns the percentage of user supplied data that is used in each
	 * iteration as test dataset.
	 * 
	 * @return the percentage of user supplied data that is used in each
	 *         iteration as test dataset
	 */
	public double getPercent() {
		return (Double)( this.getParameterAt( 4 ).getValue() );
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
		return (PartitionMethod)( (EnumParameter)getParameterAt( 5 ) ).getValue();
	}

	/**
	 * Returns <code>true</code> if for test and train dataset the sequences of
	 * the non-reference classes have the same length as the corresponding
	 * sequence of the reference class.
	 * 
	 * @return returns <code>true</code> if for test and train data set the
	 *         sequences of the non-reference classes have the same length as
	 *         the corresponding sequence of the reference class
	 */
	public boolean sameLength() {
		return (Boolean)( getParameterAt( 6 ).getValue() );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifier.assessment.ClassifierAssessmentAssessParameterSet#getAnnotation()
	 */
	@Override
	public Collection<Result> getAnnotation() {
		ArrayList<Result> l = new ArrayList<Result>( 3 );
		l.add( new NumericalResult( "repeats", "The number of iterations", getRepeats() ) );
		l.add( new NumericalResult( "reference class",
				"The index of the class for which the complete data set is partitioned in an iteration",
				getReferenceClass() ) );
		l.add( new NumericalResult( "percentage",
				"The percentage of the entire data of the reference class, that was used in an assessment",
				getPercent() ) );
		l.add( new CategoricalResult( "dataSplitMethod",
				"Describes how data should be splitted in ClassifierAssessment.evaluateClassifier())",
				getDataSplitMethod().name() ) );
		l.add( new CategoricalResult( "sameLength",
				"if true then for test respectively train data set the sequences of the non-reference" + "classes have the same length as the corresponding sequence of the reference class",
				sameLength() ) );
		return l;
	}

}
