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

package de.jstacs.classifiers;

import java.util.LinkedList;

import de.jstacs.Storable;
import de.jstacs.classifiers.performanceMeasures.AbstractPerformanceMeasure;
import de.jstacs.classifiers.performanceMeasures.PerformanceMeasureParameterSet;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.DataSet.ElementEnumerator;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.StorableResult;

/**
 * The super class for any classifier.
 * 
 * <br>
 * <br>
 * 
 * <a name="order"> The order of the classes is never changed inside the
 * classifier. The samples you put in the methods like <code>train, test</code>
 * and <code>evaluate</code> should always have the same order that you have
 * used while instantiation of the object.</a>
 * 
 * <br>
 * <br>
 * 
 * <b>For two classes it is highly recommended to set the foreground as first
 * class and the second class as background.</b>
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public abstract class AbstractClassifier implements Storable, Cloneable {

	/**
	 * The underlying alphabets.
	 */
	private AlphabetContainer alphabets;

	/**
	 * The underlying length.
	 */
	private int length;

	/**
	 * The constructor for a homogeneous classifier. Such a classifier can
	 * handle sequences of arbitrary length.
	 * 
	 * @param abc
	 *            the alphabets that are used
	 * 
	 * @see AbstractClassifier#AbstractClassifier(AlphabetContainer, int)
	 */
	public AbstractClassifier( AlphabetContainer abc ) {
		this( abc, 0 );
	}

	/**
	 * The constructor for an inhomogeneous classifier. Such a classifier can
	 * handle sequences of fixed length.
	 * 
	 * @param abc
	 *            the alphabets that are used
	 * @param length
	 *            the length of the sequences that can be classified
	 * 
	 * @throws IllegalArgumentException
	 *             if the length and the possible length of the
	 *             {@link AlphabetContainer} does not match
	 */
	public AbstractClassifier( AlphabetContainer abc, int length ) throws IllegalArgumentException {
		int l = abc.getPossibleLength();
		if( l != 0 && l != length ) {
			throw new IllegalArgumentException( "The length and the possible length of the AlphabetContainer does not match." );
		}
		alphabets = abc;
		this.length = length;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link AbstractClassifier} out of its XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link AbstractClassifier} could not be reconstructed
	 *             out of the XML representation (the {@link StringBuffer} could
	 *             not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 */
	public AbstractClassifier( StringBuffer xml ) throws NonParsableException {
		alphabets = null;
		length = -1;
		fromXML( xml );
		if( length < 0 || alphabets == null ) {
			throw new NonParsableException( "The alphabets or the length were not set." );
		}
	}

	/**
	 * This method classifies a sequence and returns the index <code>i</code> of
	 * the class to which the sequence is assigned with
	 * <code>0 &lt; i &lt; getNumberOfClasses()</code>.
	 * 
	 * <br>
	 * <br>
	 * 
	 * This method should check that the sequence is defined over the underlying
	 * alphabet and length.
	 * 
	 * @param seq
	 *            the sequence to be classified
	 * 
	 * @return the index of the class to which the sequence is assigned
	 * 
	 * @throws Exception
	 *             if the classifier is not trained or something is wrong with
	 *             the sequence
	 */
	public abstract byte classify( Sequence seq ) throws Exception;

	/**
	 * This method classifies all sequences of a sample and returns an array of
	 * indices of the classes to which the respective sequences are assigned
	 * with for each index <code>i</code> in the array
	 * <code>0 &lt; i &lt; getNumberOfClasses()</code>.
	 * 
	 * @param s
	 *            the sample to be classified
	 * 
	 * @return an array of class assignments
	 * 
	 * @throws Exception
	 *             if something went wrong during the classification
	 */
	public byte[] classify( DataSet s ) throws Exception {
		byte[] clazz = new byte[s.getNumberOfElements()];
		ElementEnumerator ei = new ElementEnumerator( s );
		for( int i = 0; i < clazz.length; i++ ) {
			clazz[i] = classify( ei.nextElement() );
		}
		return clazz;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#clone()
	 */
	@Override
	public AbstractClassifier clone() throws CloneNotSupportedException {
		return (AbstractClassifier)super.clone();
	}

	/**
	 * This method evaluates the classifier and computes, for instance, the sensitivity for a given specificity, the
	 * area under the ROC curve and so on. This method should be used in any
	 * kind of {@link de.jstacs.classifiers.assessment.ClassifierAssessment} as, for instance, crossvalidation, hold out
	 * sampling, ... .
	 * 
	 * <br>
	 * <br>
	 * 
	 * For two classes it is highly recommended to set the foreground as first
	 * class and the second class as background, i.e. the first sample should be
	 * the foreground sample and the second should be background sample. See
	 * also <a href="#order">this comment</a>.
	 * 
	 * @param params
	 *            the current parameters defining the set of {@link AbstractPerformanceMeasure}s to be evaluated
	 * @param exceptionIfNotComputeable
	 *            indicates that the method throws an {@link Exception} if a measure
	 *            could not be computed
	 * @param s
	 *            the array of {@link DataSet}s
	 * 
	 * @return a set of results, if all results are scalars the return type is {@link NumericalResultSet}, otherwise {@link ResultSet}
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see #evaluate(PerformanceMeasureParameterSet, boolean, DataSet[], double[][])
	 */
	@SuppressWarnings( "unchecked" )
	public final ResultSet evaluate( PerformanceMeasureParameterSet params, boolean exceptionIfNotComputeable, DataSet... s ) throws Exception {
		return evaluate( params, exceptionIfNotComputeable, s, null );
	}
	
	/**
	 * This method evaluates the classifier and computes, for instance, the sensitivity for a given specificity, the
	 * area under the ROC curve and so on. This method should be used in any
	 * kind of {@link de.jstacs.classifiers.assessment.ClassifierAssessment} as, for instance, crossvalidation, hold out
	 * sampling, ... .
	 * 
	 * <br>
	 * <br>
	 * 
	 * For two classes it is highly recommended to set the foreground as first
	 * class and the second class as background, i.e. the first sample should be
	 * the foreground sample and the second should be background sample. See
	 * also <a href="#order">this comment</a>.
	 * 
	 * @param params
	 *            the current parameters defining the set of {@link AbstractPerformanceMeasure}s to be evaluated
	 * @param exceptionIfNotComputeable
	 *            indicates that the method throws an {@link Exception} if a measure
	 *            could not be computed
	 * @param s
	 *            the array of {@link DataSet}s
	 * @param weights
	 * 			  the weights of the sequences for each data set
	 * 
	 * @return a set of results, if all results are scalars the return type is {@link NumericalResultSet}, otherwise {@link ResultSet}
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see NumericalResultSet
	 * @see ResultSet
	 * @see #getResults(LinkedList, DataSet[], double[][], PerformanceMeasureParameterSet, boolean)
	 * 
	 * @see de.jstacs.classifiers.assessment.ClassifierAssessment
	 * @see de.jstacs.classifiers.assessment.ClassifierAssessment#assess(de.jstacs.classifiers.performanceMeasures.NumericalPerformanceMeasureParameterSet, de.jstacs.classifiers.assessment.ClassifierAssessmentAssessParameterSet, DataSet...)
	 */
	public final ResultSet evaluate( PerformanceMeasureParameterSet params, boolean exceptionIfNotComputeable, DataSet[] s, double[][] weights ) throws Exception {
		LinkedList list = new LinkedList();
		boolean isNumeric = getResults( list, s, weights, params, exceptionIfNotComputeable );
		if( isNumeric ) {
			return new NumericalResultSet( (LinkedList<NumericalResult>) list );
		} else {
			return new ResultSet( (LinkedList<Result>) list );
		}
	}

	/**
	 * This method computes the results for any evaluation of the classifier.
	 * 
	 * @param list 
	 *            a list adding the results
	 * @param s
	 *            the array of {@link DataSet}s
	 * @param weights
	 * 			  the weights of the sequences for each data set
	 * @param params
	 *            the current parameters
	 * @param exceptionIfNotComputeable
	 *            indicates the method throws an {@link Exception} if a measure
	 *            could not be computed
	 * 
	 * @return a boolean indicating if all results are numerical
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see #evaluate(PerformanceMeasureParameterSet, boolean, DataSet...)
	 * @see NumericalResult
	 * @see Result
	 */
	@SuppressWarnings( "unchecked" )
	protected boolean getResults( LinkedList list, DataSet[] s, double[][] weights, de.jstacs.classifiers.performanceMeasures.PerformanceMeasureParameterSet params, boolean exceptionIfNotComputeable ) throws Exception {
		if( s.length != getNumberOfClasses() ) {
			throw new ClassDimensionException();
		}
		double[][][] scores = getMultiClassScores( s );
		
		boolean isNumeric = true;
		AbstractPerformanceMeasure[] m = params.getAllMeasures();
		for( AbstractPerformanceMeasure current : m ) {
			ResultSet r = null;
			try {
				r = current.compute( scores, weights );
			} catch( Exception e ) {
				if( exceptionIfNotComputeable ) {
					throw e;
				}
			}
			if( r != null ) {
				isNumeric &= r instanceof NumericalResultSet;
				for( int j = 0; j < r.getNumberOfResults(); j++ ) {
					list.add( r.getResultAt(j) );
				}
			}
		}
		return isNumeric;
	}
	
	/**
	 * This method returns a multidimensional array with class specific scores. The first dimension is for the data set, the second for the sequences, and the third for the classes.
	 * The entry <code>result[d][n][c]</code> returns the score of class <code>c</code> for sequence <code>n</code> of the data set <code>s[d]</code>. The class with the maximum score
	 * for any sequence is the predicted class of the sequence. 
	 * 
	 * @param s the data sets
	 * @return a multidimensional array with class specific scores
	 * 
	 * @throws Exception if the scores can not be computed
	 * 
	 * @see #getResults(LinkedList, DataSet[], double[][], PerformanceMeasureParameterSet, boolean)
	 */
	protected double[][][] getMultiClassScores( DataSet[] s ) throws Exception {
		double[][][] scores = new double[getNumberOfClasses()][][];
		for( int d = 0; d < s.length; d++ ) {
			scores[d] = new double[s[d].getNumberOfElements()][scores.length];
			for( int n = 0; n < scores[d].length; n++ ) {
				scores[d][n][classify( s[d].getElementAt(n) )] = 1;
			}
		}
		return scores;
	}

	/**
	 * This method returns the container of alphabets that is used in the
	 * classifier.
	 * 
	 * @return the used container of alphabets
	 */
	public final AlphabetContainer getAlphabetContainer() {
		return alphabets;
	}

	/**
	 * Returns some information characterizing or describing the current
	 * instance of the classifier. This could be for instance the number of edges for
	 * a Bayesian network or an image showing some representation of the model
	 * of a class. The set of characteristics should always include the XML
	 * representation of the classifier. The corresponding result type is
	 * {@link StorableResult}.
	 * 
	 * @return the characteristics of the current instance of the classifier
	 * 
	 * @throws Exception
	 *             if some of the characteristics could not be defined
	 * 
	 * @see StorableResult
	 * @see AbstractClassifier#getNumericalCharacteristics()
	 * @see ResultSet#ResultSet(de.jstacs.results.Result[][])
	 */
	public ResultSet getCharacteristics() throws Exception {
		return new ResultSet( getNumericalCharacteristics().getResults(), new Result[]{ new StorableResult( "classifer",
				"the xml representation of the classifier",
				this ) } );
	}

	/**
	 * Returns a <b>short</b> description of the classifier.
	 * 
	 * @return a <b>short</b> description of the classifier
	 */
	public abstract String getInstanceName();

	/**
	 * Returns an array of {@link Result}s of dimension
	 * {@link #getNumberOfClasses()} that contains information about the
	 * classifier and for each class.<br>
	 * <br>
	 * 
	 * <code>
	 * res[0] = new CategoricalResult( "classifier", "the kind of classifier", getInstanceName() );<br>
	 * res[1] = new CategoricalResult( "class info 0", "some information about the class", "info0" );<br>
	 * res[2] = new CategoricalResult( "class info 1", "some information about the class", "info1" );<br>
	 * ...
	 * </code>
	 * 
	 * @return an array of {@link Result}s that contains information about the
	 *         classifier
	 */
	public abstract CategoricalResult[] getClassifierAnnotation();

	/**
	 * Returns the length of the sequences this classifier can handle or
	 * <code>0</code> for sequences of arbitrary length.
	 * 
	 * @return the length of the sequences the classifier can handle
	 */
	public final int getLength() {
		return length;
	}

	/**
	 * Returns the subset of numerical values that are also returned by
	 * {@link #getCharacteristics()}.
	 * 
	 * @return the numerical characteristics
	 * 
	 * @throws Exception
	 *             if some of the characteristics could not be defined
	 */
	public abstract NumericalResultSet getNumericalCharacteristics() throws Exception;

	/**
	 * Returns the number of classes that can be distinguished. For example if
	 * distinguishing between foreground and background this method should
	 * return 2, even if you use a mixture model for either foreground or
	 * background.
	 * 
	 * @return the number of classes that can be distinguished
	 */
	public abstract int getNumberOfClasses();

	/**
	 * This method gives information about the state of the classifier.
	 * 
	 * @return <code>true</code> if the classifier is initialized and therefore able
	 *         to classify sequences, otherwise <code>false</code>
	 */
	public abstract boolean isInitialized();

	/**
	 * Trains the {@link AbstractClassifier} object given the data as
	 * {@link DataSet}s.<br>
	 * This method should work non-incrementally. That means the result of the
	 * following series: <code>train(data1); train(data2);</code> should be a
	 * fully trained model over <code>data2</code> and not over
	 * <code>data1, data2</code>.
	 * 
	 * <br>
	 * <br>
	 * 
	 * This method should check that the {@link DataSet}s are defined over the
	 * underlying alphabet and length.
	 * 
	 * @param s
	 *            the data
	 *            <ul>
	 *            <li>either an array of {@link DataSet}s:
	 *            <code>train( new DataSet[]{s1,s2,s3})</code> or
	 *            <li>an enumeration of {@link DataSet}s:
	 *            <code>train(s1,s2,s3)</code>
	 *            </ul>
	 * 
	 * @throws Exception
	 *             if the training did not succeed
	 * 
	 * @see AbstractClassifier#train(DataSet[], double[][])
	 */
	public void train( DataSet... s ) throws Exception {
		train( s, new double[s.length][] );
	}

	/**
	 * This method trains a classifier over an array of weighted {@link DataSet}
	 * s. That is why the following has to be fulfilled:
	 * 
	 * <ul>
	 * <li> <code>s.length == weights.length</code>
	 * <li>and for all i:
	 * <code>weights[i] == null || s[i].getNumberOfElements() == weights[i].length</code>.
	 * </ul>
	 * 
	 * This method should work non-incrementally as the method
	 * {@link #train(DataSet...)}.
	 * 
	 * <br>
	 * <br>
	 * 
	 * This method should check that the {@link DataSet}s are defined over the
	 * underlying alphabet and length.
	 * 
	 * @param s
	 *            an array of {@link DataSet}s
	 * @param weights
	 *            the weights for the {@link DataSet}s
	 * 
	 * @throws Exception
	 *             if the weights are incorrect or the training did not succeed
	 * 
	 * @see AbstractClassifier#train(DataSet...)
	 */
	public abstract void train( DataSet[] s, double[][] weights ) throws Exception;

	// methods for Storable

	/**
	 * Returns the {@link String} that is used as tag for the XML representation
	 * of the classifier. This method is used by the methods
	 * {@link #fromXML(StringBuffer)} and {@link #toXML()}.
	 * 
	 * @return the {@link String} that is used as tag for the XML representation
	 *         of the classifier
	 */
	protected abstract String getXMLTag();

	private void fromXML( StringBuffer representation ) throws NonParsableException {
		StringBuffer xml = XMLParser.extractForTag( representation, getXMLTag() );
		alphabets = (AlphabetContainer) XMLParser.extractObjectForTags( xml, "alphabetcontainer" );
		length = XMLParser.extractObjectForTags( xml, "length", int.class );
		extractFurtherClassifierInfosFromXML( xml );
	}

	/**
	 * Extracts further information of a classifier from an XML representation.
	 * This method is used by the method {@link #fromXML(StringBuffer)} and
	 * should not be made public.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the information could not be parsed out of the XML
	 *             representation (the {@link StringBuffer} could not be parsed)
	 * 
	 * @see #fromXML(StringBuffer)
	 */
	protected abstract void extractFurtherClassifierInfosFromXML( StringBuffer xml ) throws NonParsableException;

	/* (non-Javadoc)
	 * @see de.jstacs.Storable#toXML()
	 */
	public final StringBuffer toXML() {
		StringBuffer xml = new StringBuffer( 100000 );
		XMLParser.appendObjectWithTags( xml, alphabets, "alphabetcontainer" );
		XMLParser.appendObjectWithTags( xml, length, "length" );
		xml.append( getFurtherClassifierInfos() );
		XMLParser.addTags( xml, getXMLTag() );
		return xml;
	}

	/**
	 * This method returns further information of a classifier as a
	 * {@link StringBuffer}. This method is used by the method {@link #toXML()}
	 * and should not be made public.
	 * 
	 * @return further information of a classifier as a {@link StringBuffer}
	 * 
	 * @see AbstractClassifier#toXML()
	 */
	protected abstract StringBuffer getFurtherClassifierInfos();
}
