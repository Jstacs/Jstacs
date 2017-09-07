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

package de.jstacs.classifiers.trainSMBased;

import java.util.LinkedList;

import javax.naming.OperationNotSupportedException;

import de.jstacs.NotTrainedException;
import de.jstacs.classifiers.AbstractScoreBasedClassifier;
import de.jstacs.classifiers.ClassDimensionException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.StorableResult;
import de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModel;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.ToolBox;

/**
 * Classifier that works on {@link TrainableStatisticalModel}s for each of the different classes.
 * By calling the {@link TrainSMBasedClassifier#train(DataSet...)} method of this classifier,
 * all internal {@link TrainableStatisticalModel}s are learned on the {@link DataSet} for the corresponding class
 * using their own {@link TrainableStatisticalModel#train(DataSet)} method. In addition, the a-priori class
 * probabilities are estimated.
 * 
 * After training, the method {@link TrainSMBasedClassifier#getScore(Sequence, int)} returns the joint
 * probability (likelihood) of the provided {@link Sequence} and the specified class.
 * 
 * For two-class problems, the method {@link TrainSMBasedClassifier#getScores(DataSet)} returns the log-likelihood ratios for
 * all {@link Sequence}s in the provided {@link DataSet}.
 * 
 * The methods {@link TrainSMBasedClassifier#classify(Sequence)} and {@link TrainSMBasedClassifier#classify(DataSet)} use the likelihoods of {@link Sequence} 
 * and class and report the class yielding the maximum likelihood.
 * 
 * @author Jens Keilwagen
 * 
 * @see TrainableStatisticalModel
 */
public class TrainSMBasedClassifier extends AbstractScoreBasedClassifier {

	/**
	 * The internal {@link TrainableStatisticalModel}s. {@link TrainableStatisticalModel} 0 handles class 0;
	 * {@link TrainableStatisticalModel} 1 handles class 1 ... etc.
	 */
	protected TrainableStatisticalModel[] models;

	/**
	 * This method returns the possible length of a classifier that would use
	 * the given {@link TrainableStatisticalModel}s.
	 * 
	 * @param models
	 *            the {@link TrainableStatisticalModel}s that will be tested
	 * 
	 * @return the length of a classifier that would use the given models
	 * 
	 * @throws IllegalArgumentException
	 *             if no classifier could be created since the {@link TrainableStatisticalModel}s
	 *             have incompatible lengths
	 */
	public static int getPossibleLength( TrainableStatisticalModel... models ) throws IllegalArgumentException {
		int length = 0, l, i = 0;
		while( i < models.length ) {
			l = models[i++].getLength();
			if( l != 0 && l != length ) {
				if( length == 0 ) {
					length = l;
				} else {
					throw new IllegalArgumentException( "The models can't be used for one classifier. Since at least one model has length " + length
														+ ", while another has length "
														+ l
														+ "." );
				}
			}
		}
		return length;
	}

	/**
	 * This constructor creates a new instance with the given {@link TrainableStatisticalModel}s and
	 * clones these if necessary.
	 * 
	 * @param cloneModels
	 *            a switch to decide whether to clone the {@link TrainableStatisticalModel} or not
	 * @param models
	 *            the {@link TrainableStatisticalModel}s
	 * 
	 * @throws IllegalArgumentException
	 *             if the {@link TrainableStatisticalModel}s do not describe a common domain of
	 *             sequences
	 * @throws CloneNotSupportedException
	 *             if at least one {@link TrainableStatisticalModel} could not be cloned
	 * @throws ClassDimensionException
	 *             if the number of classes is below 2
	 * 
	 * @see AbstractScoreBasedClassifier#AbstractScoreBasedClassifier(AlphabetContainer,
	 *      int, int, double)
	 */
	protected TrainSMBasedClassifier( boolean cloneModels, TrainableStatisticalModel... models ) throws IllegalArgumentException, CloneNotSupportedException,
																			ClassDimensionException {
		super( models[0].getAlphabetContainer(), getPossibleLength( models ), models.length, -Math.log( (double)models.length ) );

		int i = checkAndSetModels( models, cloneModels );
		if( i <= 0 ) {
			throw new IllegalArgumentException( "Check length and AlphabetContainer of model " + ( -1 * i ) + "." );
		}
	}

	/**
	 * The default constructor that creates a new instance with the given
	 * {@link TrainableStatisticalModel}s.
	 * 
	 * @param models
	 *            the {@link TrainableStatisticalModel}s
	 * 
	 * @throws IllegalArgumentException
	 *             if the {@link TrainableStatisticalModel}s do not describe a common domain of
	 *             sequences
	 * @throws CloneNotSupportedException
	 *             if at least one {@link TrainableStatisticalModel} could not be cloned
	 * @throws ClassDimensionException
	 *             if the number of classes is below 2
	 * 
	 * @see TrainSMBasedClassifier#TrainSMBasedClassifier(boolean, TrainableStatisticalModel...)
	 */
	public TrainSMBasedClassifier( TrainableStatisticalModel... models ) throws IllegalArgumentException, CloneNotSupportedException, ClassDimensionException {
		this( true, models );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link TrainSMBasedClassifier} out of its XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link TrainSMBasedClassifier} could not be
	 *             reconstructed out of the XML representation (the
	 *             {@link StringBuffer} could not be parsed)
	 * 
	 * @see AbstractScoreBasedClassifier#AbstractScoreBasedClassifier(StringBuffer)
	 * @see de.jstacs.Storable
	 */
	public TrainSMBasedClassifier( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractScoreBasedClassifier#clone()
	 */
	@Override
	public TrainSMBasedClassifier clone() throws CloneNotSupportedException {
		TrainSMBasedClassifier clone = (TrainSMBasedClassifier)super.clone();
		clone.models = ArrayHandler.clone( models );
		return clone;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractClassifier#getCharacteristics()
	 */
	@Override
	public ResultSet getCharacteristics() throws Exception {
		ResultSet set;
		LinkedList<Result> list = new LinkedList<Result>();
		int i = 0, j;
		for( ; i < models.length; i++ ) {
			set = models[i].getCharacteristics();
			if( set != null && set.getNumberOfResults() > 0 ) {
				list.add( new NumericalResult( "class index", "the index of the class that produces the following results", i ) );
				for( j = 0; j < set.getNumberOfResults(); j++ ) {
					list.add( set.getResultAt( j ) );
				}
			}
		}
		list.add( new StorableResult( "classifer", "the xml representation of the classifier", this ) );
		return new ResultSet( list );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractClassifier#getInstanceName()
	 */
	@Override
	public String getInstanceName() {
		return "model-based classifier";
	}

	/**
	 * Returns a clone of the {@link TrainableStatisticalModel} for a specified class.
	 * 
	 * @param classIndex
	 *            the index of the specified class
	 * 
	 * @return a clone of the {@link TrainableStatisticalModel} of the specified class
	 * 
	 * @throws CloneNotSupportedException
	 *             if the {@link TrainableStatisticalModel} could not be cloned
	 * 
	 * @see TrainableStatisticalModel#clone()
	 */
	public TrainableStatisticalModel getModel( int classIndex ) throws CloneNotSupportedException {
		return models[classIndex].clone();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractClassifier#getNumericalCharacteristics()
	 */
	@Override
	public NumericalResultSet getNumericalCharacteristics() throws Exception {
		NumericalResultSet set;
		LinkedList<NumericalResult> list = new LinkedList<NumericalResult>();
		int i = 0, j;
		for( ; i < models.length; i++ ) {
			set = models[i].getNumericalCharacteristics();
			if( set != null && set.getNumberOfResults() > 0 ) {
				list.add( new NumericalResult( "class index", "the index of the class that produces the following results", i ) );
				for( j = 0; j < set.getNumberOfResults(); j++ ) {
					list.add( set.getResultAt( j ) );
				}
			}
		}
		return new NumericalResultSet( list );
	}

	/* 
	 * (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractClassifier#isInitialized()
	 */
	@Override
	public boolean isInitialized() {
		int i = 0;
		while( i < models.length && models[i].isInitialized() ) {
			i++;
		}
		return i == models.length;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractClassifier#train(de.jstacs.data.DataSet[], double[][])
	 */
	@Override
	public void train( DataSet[] s, double[][] weights ) throws Exception {
		if( weights != null && s.length != weights.length ) {
			throw new IllegalArgumentException( "data and weights do not match" );
		}
		if( models.length != s.length ) {
			throw new ClassDimensionException();
		}
		double[] c = new double[models.length];
		for( int i = 0; i < models.length; i++ ) {
			// estimate P(seq|class = i,\lambda)
			if( weights == null || weights[i] == null ) {
				models[i].train( s[i] );
			} else {
				models[i].train( s[i], weights[i] );
			}
			
			// estimate P(class = i|\lambda)
			if(getLength() > 0){
				c[i] = Math.log( s[i].getNumberOfElementsWithLength( getLength(), weights == null ? null : weights[i] ) );//XXX + models[i].getESS() );
			}else{
				c[i] = weights==null || weights[i] == null ? Math.log(s[i].getNumberOfElements()) : Math.log( ToolBox.sum( weights[i] ) );
			}
		}
		
		double norm = Normalisation.getLogSum( c );
		for(int i=0;i<c.length;i++){
			c[i] -= norm;
		}
		
		setClassWeights( false, c );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractScoreBasedClassifier#getFurtherClassifierInfos()
	 */
	@Override
	protected StringBuffer getFurtherClassifierInfos() {
		StringBuffer xml = super.getFurtherClassifierInfos();
		XMLParser.appendObjectWithTags( xml, models, "models" );
		return xml;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractScoreBasedClassifier#getScore(de.jstacs.data.Sequence, int, boolean)
	 */
	@Override
	protected double getScore( Sequence seq, int i, boolean check ) throws Exception {
		if( check ) {
			check( seq );
		}
		return models[i].getLogProbFor( seq ) + getClassWeight( i );
	}
	
	/**
	 * Returns the log likelihood ratios along the sequence <code>seq</code> for all
	 * sliding windows of length {@link #getLength()}. 
	 * @param seq the sequence
	 * @return the log likelihood ratios
	 * @throws Exception if the classifier has not been trained, it is no binary classifier, or the {@link AlphabetContainer}s of sequence and classifier do not match
	 */
	public double[] getLogLikelihoodRatio( Sequence seq ) throws Exception {
		if( !isInitialized() ) {
			throw new NotTrainedException( "The classifier is not trained yet." );
		}
		if( getNumberOfClasses()!= 2 ) {
			throw new IllegalArgumentException( "This method can only be used for binary classifiers." );
		}
		if( !getAlphabetContainer().checkConsistency( seq.getAlphabetContainer() ) ) {
			throw new IllegalArgumentException( "The sequence is not defined over the correct alphabets." );
		}
		
		double[] score = new double[seq.getLength()-getLength()+1];
		double constant = getClassWeight(0) - getClassWeight(1);
		for( int p = 0; p < score.length; p++ ) {
			score[p] = models[0].getLogProbFor(seq, p) - models[1].getLogProbFor(seq, p) + constant;
		}
		return score;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractScoreBasedClassifier#getScores(de.jstacs.data.DataSet)
	 */
	@Override
	public double[] getScores( DataSet s ) throws Exception {
		if( getNumberOfClasses() != 2 ) {
			throw new OperationNotSupportedException( "This method is only for 2-class-classifiers." );
		}
		if( s == null ) {
			return new double[0];
		}
		check( s );
		double[] score0 = models[0].getLogScoreFor( s );
		double[] score1 = models[1].getLogScoreFor( s );
		double c0 = getClassWeight( 0 ), c1 = getClassWeight( 1 );
		for( int i = 0; i < score0.length; i++ ) {
			score0[i] += c0 - ( score1[i] + c1 );
		}
		return score0;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractClassifier#classify(de.jstacs.data.DataSet)
	 */
	@Override
	public byte[] classify( DataSet s ) throws Exception {
		check( s );
		double[] best = models[0].getLogScoreFor( s ), current = new double[best.length];
		byte[] clazz = new byte[best.length];
		double cw = getClassWeight( 0 );
		int i = 0;
		for( ; i < best.length; i++ ) {
			best[i] += cw;
		}
		for( byte j = 1; j < getNumberOfClasses(); j++ ) {
			cw = getClassWeight( j );
			models[j].getLogScoreFor( s, current );
			for( i = 0; i < best.length; i++ ) {
				if( current[i] + cw > best[i] ) {
					best[i] = current[i] + cw;
					clazz[i] = j;
				}
			}
		}
		return clazz;
	}

	private static final String XML_TAG = "TrainSMBasedClassifier";

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractClassifier#getXMLTag()
	 */
	@Override
	protected String getXMLTag() {
		return XML_TAG;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractScoreBasedClassifier#extractFurtherClassifierInfosFromXML(java.lang.StringBuffer)
	 */
	@Override
	protected void extractFurtherClassifierInfosFromXML( StringBuffer xml ) throws NonParsableException {
		super.extractFurtherClassifierInfosFromXML( xml );
		int i;
		try {
			i = checkAndSetModels( XMLParser.extractObjectForTags( xml, "models", TrainableStatisticalModel[].class ), false );
		} catch ( CloneNotSupportedException e ) {
			NonParsableException n = new NonParsableException( "Clone not supported: " + e.getMessage() );
			n.setStackTrace( e.getStackTrace() );
			throw n;
		}
		if( i <= 0 ) {
			throw new NonParsableException( "Check length and AlphabetContainer of model " + ( -1 * i ) + "." );
		}
	}

	private int checkAndSetModels( TrainableStatisticalModel[] models, boolean clone ) throws CloneNotSupportedException {
		int i = 0, l, length = getLength();
		AlphabetContainer abc = getAlphabetContainer();
		this.models = new TrainableStatisticalModel[models.length];
		while( i < models.length ) {
			l = models[i].getLength();
			if( ( l != 0 && length != l ) || !models[i].getAlphabetContainer().checkConsistency( abc ) ) {
				return -i;
			}
			if( clone ) {
				this.models[i] = models[i++].clone();
			} else {
				this.models[i] = models[i++];
			}
		}
		return i;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.AbstractClassifier#getClassifierAnnotation()
	 */
	@Override
	public CategoricalResult[] getClassifierAnnotation() {
		CategoricalResult[] res = new CategoricalResult[models.length + 1];
		res[0] = new CategoricalResult( "classifier", "a <b>short</b> description of the classifier", getInstanceName() );
		int i = 0;
		while( i < models.length ) {
			res[i + 1] = new CategoricalResult( "class info " + i, "some information about the class", models[i++].getInstanceName() );
		}
		return res;
	}
	
	public String toString() {
		StringBuffer sb = new StringBuffer( models.length * 5000 );
		String heading = "model ";
		for( int i =0; i < models.length; i++ ) {
			sb.append( heading + i );
			sb.append( "\n" + models[i].toString() + "\n" );
		}
		sb.append( "class weights: " );
		for(int i=0;i<getNumberOfClasses();i++){
			sb.append( getClassWeight( i )+" " );
		}
		sb.append( "\n" );
		return sb.toString();
	}
}
