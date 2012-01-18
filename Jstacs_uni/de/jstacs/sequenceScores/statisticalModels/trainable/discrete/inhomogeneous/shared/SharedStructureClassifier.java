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

package de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.shared;

import de.jstacs.algorithms.graphs.tensor.SymmetricTensor;
import de.jstacs.classifiers.ClassDimensionException;
import de.jstacs.classifiers.trainSMBased.TrainSMBasedClassifier;
import de.jstacs.data.DataSet;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.results.CategoricalResult;
import de.jstacs.sequenceScores.statisticalModels.trainable.AbstractTrainSM;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.FSDAGTrainSM;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.StructureLearner;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.StructureLearner.LearningType;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.StructureLearner.ModelType;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.parameters.BayesianNetworkTrainSMParameterSet;

/**
 * This class enables you to learn the structure on all classes of the
 * classifier together. A special case is, for instance, a <b>T</b>ree
 * <b>A</b>ugmented <b>N</b>aive Bayes (TAN).
 * 
 * @author Jens Keilwagen
 */
public class SharedStructureClassifier extends TrainSMBasedClassifier {

	private ModelType model;

	private byte order;

	private LearningType method;

	private StructureLearner sl;

	/**
	 * Creates a new {@link SharedStructureClassifier} from given
	 * {@link FSDAGTrainSM}s. This is the main constructor.
	 * 
	 * @param length
	 *            the sequence length
	 * @param model
	 *            the type of the model
	 * @param order
	 *            the order of the model
	 * @param method
	 *            the learning method
	 * @param models
	 *            the class models
	 * 
	 * @throws IllegalArgumentException
	 *             if <code>order</code> is below 0
	 * @throws CloneNotSupportedException
	 *             if at least one model could not be cloned
	 * @throws ClassDimensionException
	 *             if the class dimension is wrong (below 2)
	 * 
	 * @see ModelType
	 * @see LearningType
	 * @see TrainSMBasedClassifier#TrainSMBasedClassifier(boolean, de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModel...)
	 */
	public SharedStructureClassifier( int length, ModelType model, byte order, LearningType method, FSDAGTrainSM... models )
																															throws IllegalArgumentException,
																															CloneNotSupportedException,
																															ClassDimensionException {
		super( true, (AbstractTrainSM[])models );
		this.model = model;
		if( order < 0 ) {
			throw new IllegalArgumentException( "The value of order has to be non-negative." );
		}
		this.order = order;
		this.method = method;
		sl = new StructureLearner( getAlphabetContainer(), length );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link SharedStructureClassifier} out of its XML
	 * representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link SharedStructureClassifier} could not be
	 *             reconstructed out of the XML representation (the
	 *             {@link StringBuffer} could not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 * @see TrainSMBasedClassifier#TrainSMBasedClassifier(StringBuffer)
	 */
	public SharedStructureClassifier( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.trainSMBased.TrainSMBasedClassifier#clone()
	 */
	@Override
	public SharedStructureClassifier clone() throws CloneNotSupportedException {
		SharedStructureClassifier clone = (SharedStructureClassifier)super.clone();
		clone.sl = new StructureLearner( getAlphabetContainer(), getLength() );
		return clone;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.trainSMBased.TrainSMBasedClassifier#train(de.jstacs.data.DataSet[], double[][])
	 */
	@Override
	public void train( DataSet[] data, double[][] weights ) throws IllegalArgumentException, Exception {
		int dimension = models.length;
		SymmetricTensor[] parts = new SymmetricTensor[dimension];
		double[] w = new double[dimension];
		for( int i = 0; i < dimension; i++ ) {
			sl.setESS( ( (FSDAGTrainSM)models[i] ).getESS() );
			parts[i] = sl.getTensor( data[i], weights[i], order, method );
			w[i] = 1d;
		}
		FSDAGTrainSM.train( models, StructureLearner.getStructure( new SymmetricTensor( parts, w ), model, order ), weights, data );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.trainSMBased.TrainSMBasedClassifier#getInstanceName()
	 */
	@Override
	public String getInstanceName() {
		return "shared-structure classifier";
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.trainSMBased.TrainSMBasedClassifier#extractFurtherClassifierInfosFromXML(java.lang.StringBuffer)
	 */
	@Override
	protected void extractFurtherClassifierInfosFromXML( StringBuffer xml ) throws NonParsableException {
		super.extractFurtherClassifierInfosFromXML( xml );
		model = XMLParser.extractObjectForTags( xml, "model", ModelType.class );
		order = XMLParser.extractObjectForTags( xml, "order", byte.class );
		method = XMLParser.extractObjectForTags( xml, "method", LearningType.class );
		sl = new StructureLearner( getAlphabetContainer(), getLength() );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.trainSMBased.TrainSMBasedClassifier#getFurtherClassifierInfos()
	 */
	@Override
	protected StringBuffer getFurtherClassifierInfos() {
		StringBuffer xml = super.getFurtherClassifierInfos();
		XMLParser.appendObjectWithTags( xml, model, "model" );
		XMLParser.appendObjectWithTags( xml, order, "order" );
		XMLParser.appendObjectWithTags( xml, method, "method" );
		return xml;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.trainSMBased.TrainSMBasedClassifier#getClassifierAnnotation()
	 */
	@Override
	public CategoricalResult[] getClassifierAnnotation() {
		CategoricalResult[] res = new CategoricalResult[models.length + 1];
		res[0] = new CategoricalResult( "classifier", "a <b>short</b> description of the classifier", getInstanceName() );
		int i = 0;
		while( i < models.length ) {
			res[i + 1] = new CategoricalResult( "class info " + i,
					"some information about the class",
					BayesianNetworkTrainSMParameterSet.getModelInstanceName( model, order, method, ( (FSDAGTrainSM)models[i++] ).getESS() ) );
		}
		return res;
	}
}
