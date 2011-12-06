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

import java.util.Random;

import de.jstacs.WrongAlphabetException;
import de.jstacs.classifier.AbstractClassifier;
import de.jstacs.classifier.ClassDimensionException;
import de.jstacs.classifier.performanceMeasures.NumericalPerformanceMeasureParameters;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.Sequence;
import de.jstacs.data.DataSet.PartitionMethod;
import de.jstacs.models.Model;
import de.jstacs.utils.ProgressUpdater;

/**
 * This class is a special {@link ClassifierAssessment} that partitions the data
 * of a user-specified reference class (typically the smallest class) and
 * samples non-overlapping for all other classes, so that one gets the same
 * number of sequences (and the same lengths of the sequences) in each train and
 * test dataset.
 * 
 * @author Jens Keilwagen
 * 
 * @see Sampled_RepeatedHoldOutAssessParameterSet
 */
public class Sampled_RepeatedHoldOutExperiment extends ClassifierAssessment {

	/**
	 * Creates a new {@link Sampled_RepeatedHoldOutExperiment} from an array of
	 * {@link AbstractClassifier}s and a two-dimensional array of {@link Model}
	 * s, which are combined to additional classifiers. If
	 * <code>buildClassifiersByCrossProduct</code> is <code>true</code>, the
	 * cross-product of all {@link Model}s in <code>aMs</code> is built to
	 * obtain these classifiers.
	 * 
	 * @param aCs
	 *            the predefined classifiers
	 * @param aMs
	 *            the {@link Model}s that are used to build additional
	 *            classifiers
	 * @param buildClassifiersByCrossProduct
	 *            Determines how classifiers are constructed using the given
	 *            models. Suppose a k-class problem. In this case, each
	 *            classifier is supposed to consist of k models, one responsible
	 *            for each class. <br>
	 *            Let <code>S_i</code> be the set of all models in
	 *            <code>aMs[i]</code>. Let <code>S</code> be the set
	 *            <code>S_1 x S_2 x ... x S_k</code> (cross-product).<br>
	 * <br>
	 *            <code>true</code>: all possible classifiers consisting of a
	 *            subset (set of k models) of <code>S</code> are constructed <br>
	 *            <code>false</code>: one classifier consisting of the models
	 *            <code>aMs[0][i]</code>,<code>aMs[1][i]</code>,...,
	 *            <code>aMs[k][i]</code> for a fixed <code>i</code> is
	 *            constructed. In this case, all second dimensions of
	 *            <code>aMs</code> have to be equal, say <code>m</code>. In
	 *            total <code>m</code> classifiers are constructed.
	 * @param checkAlphabetConsistencyAndLength
	 *            indicates if alphabets and lengths shall be checked for
	 *            consistency
	 * 
	 * @throws IllegalArgumentException
	 *             if the classifiers have different lengths
	 * @throws WrongAlphabetException
	 *             if the classifiers use different alphabets
	 * @throws CloneNotSupportedException
	 *             if something went wrong while cloning
	 * @throws ClassDimensionException
	 *             if there is something wrong with the class dimension of the
	 *             classifier
	 * 
	 * @see ClassifierAssessment#ClassifierAssessment(AbstractClassifier[],
	 *      Model[][], boolean, boolean)
	 */
	protected Sampled_RepeatedHoldOutExperiment( AbstractClassifier[] aCs, Model[][] aMs, boolean buildClassifiersByCrossProduct,
													boolean checkAlphabetConsistencyAndLength ) throws IllegalArgumentException,
																								WrongAlphabetException,
																								CloneNotSupportedException,
																								ClassDimensionException {
		super( aCs, aMs, buildClassifiersByCrossProduct, checkAlphabetConsistencyAndLength );
	}

	/**
	 * Creates a new {@link Sampled_RepeatedHoldOutExperiment} from a set of
	 * {@link AbstractClassifier}s.
	 * 
	 * @param aCs
	 *            contains the classifiers to be assessed,<br>
	 *            If model based classifiers are trained, the order of models in
	 *            classifiers determines, which model will be trained using
	 *            which sample in method <code>assess( ... )</code>.<br>
	 *            For a two-class problem, it is recommended
	 *            <ul>
	 *            <li>to initiate the classifiers with models in order
	 *            (foreground model (positive class), background model (negative
	 *            class))
	 *            <li>to initiate an assessment object using models in order
	 *            (foreground model (positive class), background model (negative
	 *            class))
	 *            <li>to give data <code>s</code> in order (<code>s[0]</code>
	 *            contains foreground data, <code>s[1]</code> contains
	 *            background data)
	 *            </ul>
	 * @throws IllegalArgumentException
	 *             if the classifiers have different lengths
	 * @throws WrongAlphabetException
	 *             if not all given classifiers are defined on the same
	 *             {@link AlphabetContainer}
	 * @throws CloneNotSupportedException
	 *             if something went wrong while cloning
	 * @throws ClassDimensionException
	 *             if there is something wrong with the class dimension of the
	 *             classifier
	 * 
	 * @see ClassifierAssessment#ClassifierAssessment(AbstractClassifier...)
	 */
	public Sampled_RepeatedHoldOutExperiment( AbstractClassifier... aCs ) throws IllegalArgumentException, WrongAlphabetException,
																			CloneNotSupportedException, ClassDimensionException {
		super( aCs );
	}

	/**
	 * Creates a new {@link Sampled_RepeatedHoldOutExperiment} from a set of
	 * {@link Model}s. The argument <code>buildClassifiersByCrossProduct</code>
	 * determines how these {@link Model}s are combined to classifiers.
	 * 
	 * @param buildClassifiersByCrossProduct
	 * <br>
	 *            Determines how classifiers are constructed using the given
	 *            models. Suppose a k-class problem. In this case, each
	 *            classifier is supposed to consist of k models, one responsible
	 *            for each class. <br>
	 *            Let <code>S_i</code> be the set of all models in
	 *            <code>aMs[i]</code>. Let <code>S</code> be the set
	 *            <code>S_1 x S_2 x ... x S_k</code> (cross-product).<br>
	 * <br>
	 *            <code>true</code>: all possible classifiers consisting of a
	 *            subset (set of k models) of <code>S</code> are constructed <br>
	 *            <code>false</code>: one classifier consisting of the models
	 *            <code>aMs[0][i]</code>,<code>aMs[1][i]</code>,...,
	 *            <code>aMs[k][i]</code> for a fixed <code>i</code> is
	 *            constructed. In this case, all second dimensions of
	 *            <code>aMs</code> have to be equal, say <code>m</code>. In
	 *            total <code>m</code> classifiers are constructed.
	 * @param aMs
	 * <br>
	 *            Contains the models in the following way (suppose a k-class
	 *            problem): the first dimension encodes the class (here it is
	 *            k), the second dimension (<code>aMs[i]</code>) contains the
	 *            models according to class <code>i</code>.<br>
	 *            If models are trained directly (during assessment), the order
	 *            of given models during initiation of this assessment object
	 *            determines, which sample will be used for training which
	 *            model. In general the first model will be trained using the
	 *            first sample in <code>s</code>... . <br>
	 *            For a two-class problem, it is recommended
	 *            <ul>
	 *            <li>to initiate the classifiers with models in order
	 *            (foreground model (positive class), background model (negative
	 *            class))
	 *            <li>to initiate an assessment object using models in order
	 *            (foreground model (positive class), background model (negative
	 *            class))
	 *            <li>to give data <code>s</code> in order (<code>s[0]</code>
	 *            contains foreground data, <code>s[1]</code> contains
	 *            background data)
	 *            </ul>
	 * 
	 * @throws IllegalArgumentException
	 *             if the classifiers have different lengths
	 * @throws WrongAlphabetException
	 *             if not all given classifiers are defined on the same
	 *             {@link AlphabetContainer}
	 * @throws CloneNotSupportedException
	 *             if something went wrong while cloning
	 * @throws ClassDimensionException
	 *             if there is something wrong with the class dimension of the
	 *             classifier
	 * 
	 * @see ClassifierAssessment#ClassifierAssessment(boolean, Model[][])
	 */
	public Sampled_RepeatedHoldOutExperiment( boolean buildClassifiersByCrossProduct, Model[]... aMs ) throws IllegalArgumentException,
																										WrongAlphabetException,
																										CloneNotSupportedException,
																										ClassDimensionException {
		super( buildClassifiersByCrossProduct, aMs );
	}

	/**
	 * This constructor allows to assess a collection of given
	 * {@link AbstractClassifier}s and those constructed using the given
	 * {@link de.jstacs.models.Model}s by a
	 * {@link Sampled_RepeatedHoldOutExperiment}. <br>
	 * 
	 * @param aCs
	 *            contains some {@link AbstractClassifier} that should be
	 *            assessed in addition to the {@link AbstractClassifier}s
	 *            constructed using the given
	 *            {@link de.jstacs.models.Model}s
	 * @param buildClassifiersByCrossProduct
	 * <br>
	 *            Determines how classifiers are constructed using the given
	 *            models. Suppose a k-class problem. In this case, each
	 *            classifier is supposed to consist of k models, one responsible
	 *            for each class. <br>
	 *            Let <code>S_i</code> be the set of all models in
	 *            <code>aMs[i]</code>. Let <code>S</code> be the set
	 *            <code>S_1 x S_2 x ... x S_k</code> (cross-product).<br>
	 * <br>
	 *            <code>true</code>: all possible classifiers consisting of a
	 *            subset (set of k models) of <code>S</code> are constructed <br>
	 *            <code>false</code>: one classifier consisting of the models
	 *            <code>aMs[0][i]</code>,<code>aMs[1][i]</code>,...,
	 *            <code>aMs[k][i]</code> for a fixed <code>i</code> is
	 *            constructed. In this case, all second dimensions of
	 *            <code>aMs</code> have to be equal, say <code>m</code>. In
	 *            total <code>m</code> classifiers are constructed.
	 * @param aMs
	 * <br>
	 *            Contains the models in the following way (suppose a k-class
	 *            problem): the first dimension encodes the class (here it is
	 *            k), the second dimension (<code>aMs[i]</code>) contains the
	 *            models according to class <code>i</code>.<br>
	 *            If models are trained directly (during assessment), the order
	 *            of given models during initiation of this assessment object
	 *            determines, which sample will be used for training which
	 *            model. In general the first model will be trained using the
	 *            first sample in <code>s</code>... . <br>
	 *            For a two-class problem, it is recommended
	 *            <ul>
	 *            <li>to initiate the classifiers with models in order
	 *            (foreground model (positive class), background model (negative
	 *            class))
	 *            <li>to initiate a assessment object using models in order
	 *            (foreground model (positive class), background model (negative
	 *            class))
	 *            <li>to give data <code>s</code> in order (<code>s[0]</code>
	 *            contains foreground data, <code>s[1]</code> contains
	 *            background data)
	 *            </ul>
	 * 
	 * @throws IllegalArgumentException
	 *             if the classifiers have different lengths
	 * @throws WrongAlphabetException
	 *             if not all given classifiers are defined on the same
	 *             {@link AlphabetContainer}
	 * @throws CloneNotSupportedException
	 *             if something went wrong while cloning
	 * @throws ClassDimensionException
	 *             if there is something wrong with the class dimension of the
	 *             classifier
	 * 
	 * @see ClassifierAssessment#ClassifierAssessment(AbstractClassifier[],
	 *      boolean, Model[][])
	 */
	public Sampled_RepeatedHoldOutExperiment( AbstractClassifier[] aCs, boolean buildClassifiersByCrossProduct, Model[]... aMs )
																																throws IllegalArgumentException,
																																WrongAlphabetException,
																																CloneNotSupportedException,
																																ClassDimensionException {
		super( aCs, buildClassifiersByCrossProduct, aMs );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifier.assessment.ClassifierAssessment#evaluateClassifier(de.jstacs.classifier.MeasureParameters, de.jstacs.classifier.assessment.ClassifierAssessmentAssessParameterSet, de.jstacs.data.Sample[], de.jstacs.utils.ProgressUpdater)
	 */
	@Override
	protected void evaluateClassifier( NumericalPerformanceMeasureParameters mp, ClassifierAssessmentAssessParameterSet assessPS, DataSet[] s, ProgressUpdater pU ) throws IllegalArgumentException,
			Exception {
		if( s.length != 2 ) {
			throw new IllegalArgumentException( "This class can only handle two classes" );
		}
		Sampled_RepeatedHoldOutAssessParameterSet tempAssessPS = null;

		try {
			tempAssessPS = (Sampled_RepeatedHoldOutAssessParameterSet)assessPS;
		} catch ( ClassCastException e ) {
			throw new IllegalArgumentException( "Given AssessParameterSet assessPS is not of type " + "Sampled_RepeatedHoldOutAssessParameterSet." );
		}

		PartitionMethod splitMethod = tempAssessPS.getDataSplitMethod();
		int subSeqL = tempAssessPS.getElementLength();
		boolean exceptionIfMPNotComputable = tempAssessPS.getExceptionIfMPNotComputable(), sameLength = tempAssessPS.sameLength();
		int repeats = tempAssessPS.getRepeats(), referenceClass = tempAssessPS.getReferenceClass();
		double percents = tempAssessPS.getPercent();

		DataSet[][] sTrainTestClassWise = new DataSet[2][s.length];
		DataSet[] temp;

		pU.setMax( repeats );

		int j = 0, k = 0, clazz, l, iteration = 0, index, length;
		int[] number = new int[2];
		Random r = new Random();
		Sequence[][] all = new Sequence[s.length][];
		Sequence[][][] part = new Sequence[2][s.length][];
		AlphabetContainer[] abc = new AlphabetContainer[s.length];

		// cache all elements and alphabets
		for( ; j < all.length; j++ ) {
			if( j != referenceClass ) {
				all[j] = s[j].getAllElements();
				abc[j] = s[j].getAlphabetContainer();
			} else {
				all[j] = null;
				abc[j] = null;
			}
		}

		for( ; iteration < repeats; iteration++ ) {
			// split data of reference class
			temp = s[referenceClass].partition( percents, splitMethod, subSeqL );
			sTrainTestClassWise[0][referenceClass] = temp[0];
			sTrainTestClassWise[1][referenceClass] = temp[1];

			// create arrays for other classes
			if( iteration == 0 ) {
				number[0] = temp[0].getNumberOfElements();
				number[1] = temp[1].getNumberOfElements();
				for( ; k < s.length; k++ ) {
					if( k != referenceClass ) {
						part[0][k] = new Sequence[number[0]];
						part[1][k] = new Sequence[number[1]];
					}
				}
			}

			// create samples for other classes
			for( l = 0, j = 0; j < 2; j++ ) // train/test dataset
			{
				for( k = 0; k < number[j]; k++, l++ ) // number of sequences
				{
					length = temp[j].getElementAt( k ).getLength();
					for( clazz = 0; clazz < s.length; clazz++ ) // class of data
					{
						if( clazz != referenceClass ) {
							index = r.nextInt( all[clazz].length - l );
							part[j][clazz][k] = all[clazz][index];
							all[clazz][index] = all[clazz][all[clazz].length - l - 1];
							all[clazz][all[clazz].length - l - 1] = part[j][clazz][k];

							// truncate sequences if necessary
							if( sameLength ) {
								part[j][clazz][k] = part[j][clazz][k].getSubSequence( abc[clazz],
										r.nextInt( part[j][clazz][k].getLength() - length + 1 ),
										length );
							} else if( j == 1 && subSeqL != 0 ) {
								part[j][clazz][k] = part[j][clazz][k].getSubSequence( abc[clazz],
										r.nextInt( part[j][clazz][k].getLength() - subSeqL + 1 ),
										subSeqL );
							}
						}
					}
				}
				// create samples
				for( clazz = 0; clazz < s.length; clazz++ ) {
					if( clazz != referenceClass ) {
						sTrainTestClassWise[j][clazz] = new DataSet( "sampled data set", part[j][clazz] );
					}
				}
			}

			// train and test on the partitions
			train( sTrainTestClassWise[0] );
			test( mp, exceptionIfMPNotComputable, sTrainTestClassWise[1] );

			// progress updater
			pU.setValue( iteration + 1 );
			if( pU.isCancelled() ) {
				break;
			}
		}
	}
}
