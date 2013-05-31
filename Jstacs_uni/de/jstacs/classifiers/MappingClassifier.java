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

import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;

import de.jstacs.NotTrainedException;
import de.jstacs.classifiers.performanceMeasures.AbstractPerformanceMeasureParameterSet;
import de.jstacs.classifiers.performanceMeasures.PerformanceMeasure;
import de.jstacs.classifiers.performanceMeasures.PerformanceMeasureParameterSet;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.utils.IntList;
import de.jstacs.utils.Normalisation;

/**
 * This class allows the user to train the classifier on a given number of
 * classes and to evaluate the classifier on a smaller number of classes by
 * mapping classes together. For instance the user has a classifier for 3
 * classes, but likes to evaluate whether the classifier is able to discriminate
 * between class 1 and class 2 and 3. This is a good example where to use this
 * class. The user has to create its 3-class-classifier, create an instance of
 * this class using its classifier, map the test data sets together (
 * {@link MappingClassifier#mapDataSet(DataSet[])}) and invoke
 * {@link AbstractClassifier#evaluate(AbstractPerformanceMeasureParameterSet, boolean, DataSet...)}
 * with these mapped {@link DataSet}. Alternatively, the method
 * {@link AbstractClassifier#evaluate(AbstractPerformanceMeasureParameterSet, boolean, DataSet...)} can
 * be used directly and the {@link DataSet}s will be mapped internally
 * (i.e. inside the evaluate method).
 * 
 * @author Jens Keilwagen
 */
public class MappingClassifier extends AbstractScoreBasedClassifier {

	private AbstractScoreBasedClassifier classifier;

	private int[][] classMapping;

	private static int getNum( int[] mapping ) {
		HashSet<Integer> hash = new HashSet<Integer>();
		for( int i = 0; i < mapping.length; i++ ) {
			if( !hash.contains( mapping[i] ) ) {
				hash.add( mapping[i] );
			}
		}
		return hash.size();
	}

	/**
	 * Creates a new {@link MappingClassifier} from a given classifier and a
	 * class mapping.
	 * 
	 * @param classifier
	 *            the internal used classifier
	 * @param mapping
	 *            the mapping from the classes of the internal classifier to the
	 *            classes of this classifier
	 * 
	 * @throws CloneNotSupportedException
	 *             if something went wrong while cloning the classifier
	 */
	public MappingClassifier( AbstractScoreBasedClassifier classifier, int... mapping ) throws CloneNotSupportedException {
		super( classifier.getAlphabetContainer(), classifier.getLength(), getNum( mapping ) );

		if( mapping.length != classifier.getNumberOfClasses() ) {
			throw new IllegalArgumentException( "The length of the mapping is not correct." );
		}
		IntList[] il = new IntList[getNumberOfClasses()];
		for( int i = 0; i < il.length; i++ ) {
			il[i] = new IntList();
		}
		for( int i = 0; i < mapping.length; i++ ) {
			il[mapping[i]].add( i );
		}
		classMapping = new int[il.length][];
		for( int i = 0; i < il.length; i++ ) {
			if( il[i].length() == 0 ) {
				throw new IllegalArgumentException( "Mapping to class " + i + " is empty" );
			} else {
				classMapping[i] = il[i].toArray();
			}
		}
		this.classifier = classifier.clone();
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link MappingClassifier} out of its XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link MappingClassifier} could not be reconstructed
	 *             out of the XML representation (the {@link StringBuffer}
	 *             <code>representation</code> could not be parsed)
	 * 
	 * @see AbstractScoreBasedClassifier#AbstractScoreBasedClassifier(StringBuffer)
	 * @see de.jstacs.Storable
	 */
	public MappingClassifier( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractScoreBasedClassifier#extractFurtherClassifierInfosFromXML(java.lang.StringBuffer)
	 */
	@Override
	protected void extractFurtherClassifierInfosFromXML( StringBuffer xml ) throws NonParsableException {
		super.extractFurtherClassifierInfosFromXML( xml );
		classifier = XMLParser.extractObjectForTags( xml, "classifier", AbstractScoreBasedClassifier.class );
		classMapping = XMLParser.extractObjectForTags( xml, "mapping", int[][].class );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractScoreBasedClassifier#getFurtherClassifierInfos()
	 */
	@Override
	protected StringBuffer getFurtherClassifierInfos() {
		StringBuffer xml = super.getFurtherClassifierInfos();
		XMLParser.appendObjectWithTags( xml, classifier, "classifier" );
		XMLParser.appendObjectWithTags( xml, classMapping, "mapping" );
		return xml;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractScoreBasedClassifier#getScore(de.jstacs.data.Sequence, int, boolean)
	 */
	@Override
	protected double getScore( Sequence seq, int i, boolean check ) throws IllegalArgumentException, NotTrainedException, Exception {
		if( check ) {
			check( seq );
		}
		double res = classifier.getScore( seq, classMapping[i][0], true );
		for( int idx = 1; idx < classMapping[i].length; idx++ ) {
			res = Normalisation.getLogSum( res, classifier.getScore( seq, classMapping[i][1], false ) );
		}
		return res;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractClassifier#getClassifierAnnotation()
	 */
	@Override
	public CategoricalResult[] getClassifierAnnotation() {
		return classifier.getClassifierAnnotation();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractClassifier#getInstanceName()
	 */
	@Override
	public String getInstanceName() {
		return "MappingClassifier of " + classifier.getInstanceName();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractClassifier#getNumericalCharacteristics()
	 */
	@Override
	public NumericalResultSet getNumericalCharacteristics() throws Exception {
		return classifier.getNumericalCharacteristics();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractClassifier#getXMLTag()
	 */
	@Override
	protected String getXMLTag() {
		return getClass().getSimpleName();
	}

	/* 
	 * (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractClassifier#isInitialized()
	 */
	@Override
	public boolean isInitialized() {
		return classifier.isInitialized();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractClassifier#train(de.jstacs.data.DataSet[], double[][])
	 */
	@Override
	public void train( DataSet[] s, double[][] weights ) throws Exception {
		classifier.train( s, weights );
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractScoreBasedClassifier#getResults(java.util.LinkedList, de.jstacs.data.DataSet[], double[][], de.jstacs.classifiers.performanceMeasures.PerformanceMeasureParameterSet, boolean)
	 */
	@Override
	protected boolean getResults( LinkedList list, DataSet[] s, double[][] weights, AbstractPerformanceMeasureParameterSet<? extends PerformanceMeasure> params, boolean exceptionIfNotComputeable ) throws Exception {
		if( s.length == getNumberOfClasses() ) {
			return super.getResults( list, s, weights, params, exceptionIfNotComputeable );
		} else {
			
			
			return super.getResults( list, mapDataSet( s ), mapWeights(weights), params, exceptionIfNotComputeable );
		}
	}

	/**
	 * This method maps the given {@link DataSet}s to the internal classes.
	 * 
	 * @param s
	 *            the array of {@link DataSet}s
	 * 
	 * @return the array of data sets corresponding to the classes
	 */
	public DataSet[] mapDataSet( DataSet[] s ) {
		boolean[] in = new boolean[classifier.getNumberOfClasses()];
		DataSet[] mapped = new DataSet[classMapping.length];
		try {
			for( int j, i = 0; i < mapped.length; i++ ) {
				Arrays.fill( in, false );
				for( j = 0; j < classMapping[i].length; j++ ) {
					in[classMapping[i][j]] = true;
				}
				mapped[i] = DataSet.union( s, in );
			}
		} catch ( Exception e ) {
			// does not happen
			throw new RuntimeException();
		}
		return mapped;
	}
	
	/**
	 * This method maps the given {@link Sequence} weights to the internal classes.
	 * 
	 * @param w
	 *            the array of weights
	 * 
	 * @return the array of weights corresponding to the classes
	 */
	public double[][] mapWeights( double[][] w ) {
		double[][] mapped = new double[classMapping.length][];
		for( int j, i = 0; i < mapped.length; i++ ) {
			int anz = 0;
			for( j = 0; j < classMapping[i].length; j++ ) {
				anz += w[classMapping[i][j]].length;
			}
			mapped[i] = new double[anz];
			anz = 0;
			for( j = 0; j < classMapping[i].length; j++ ) {
				int current = w[classMapping[i][j]].length;
				System.arraycopy( w[classMapping[i][j]], 0, mapped[i], anz, current );
				anz += current;
			}
		}
		return mapped;
	}
}
