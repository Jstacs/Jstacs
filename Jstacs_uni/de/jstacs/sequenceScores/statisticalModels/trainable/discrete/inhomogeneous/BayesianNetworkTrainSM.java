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

package de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous;

import java.util.Arrays;

import de.jstacs.NonParsableException;
import de.jstacs.data.DataSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.DGTrainSMParameterSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.StructureLearner.LearningType;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.StructureLearner.ModelType;
import de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.parameters.BayesianNetworkTrainSMParameterSet;

/**
 * The class implements a Bayesian network (
 * {@link StructureLearner.ModelType#BN} ) of fixed order. It allows the user to
 * specify some kinds of specializations of BNs including inhomogeneous Markov
 * models ( {@link StructureLearner.ModelType#IMM} ) and permuted Markov models
 * ( {@link StructureLearner.ModelType#PMM} ).
 * 
 * @author Jens Keilwagen
 * 
 * @see ModelType
 */
public class BayesianNetworkTrainSM extends DAGTrainSM {

	private StructureLearner sl;

	/**
	 * Creates a new {@link BayesianNetworkTrainSM} from a given
	 * {@link BayesianNetworkTrainSMParameterSet}.
	 * 
	 * @param params
	 *            the given parameter set
	 * 
	 * @throws CloneNotSupportedException
	 *             if the parameter set could not be cloned
	 * @throws IllegalArgumentException
	 *             if the parameter set is not instantiated
	 * @throws NonParsableException
	 *             if the parameter set is not parsable
	 * 
	 * @see DAGTrainSM#DAGTrainSM(de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.parameters.IDGTrainSMParameterSet)
	 */
	public BayesianNetworkTrainSM( BayesianNetworkTrainSMParameterSet params ) throws CloneNotSupportedException, IllegalArgumentException,
																			NonParsableException {
		super( params );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link BayesianNetworkTrainSM} out of its XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link BayesianNetworkTrainSM} could not be
	 *             reconstructed out of the XML representation (the
	 *             {@link StringBuffer} could not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 * @see DAGTrainSM#DAGTrainSM(StringBuffer)
	 */
	public BayesianNetworkTrainSM( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.DAGTrainSM#clone()
	 */
	@Override
	public BayesianNetworkTrainSM clone() throws CloneNotSupportedException {
		BayesianNetworkTrainSM clone = (BayesianNetworkTrainSM)super.clone();
		clone.sl = new StructureLearner( alphabets, length, getESS() );
		return clone;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.TrainableStatisticalModel#getInstanceName()
	 */
	public String getInstanceName() {
		return ( (BayesianNetworkTrainSMParameterSet)params ).getModelInstanceName();
	}

	private static final String XML_TAG = "BayesianNetworkTrainSM";

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.DiscreteGraphicalTrainSM#getXMLTag()
	 */
	@Override
	protected String getXMLTag() {
		return XML_TAG;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.TrainableStatisticalModel#train(de.jstacs.data.DataSet, double[])
	 */
	public void train( DataSet data, double[] weights ) throws Exception {
		createConstraints( sl.getStructure( data, weights, getModel(), getMaximalMarkovOrder(), getMethod() ) );
		estimateParameters( data, weights );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.DAGTrainSM#getLogPriorTerm()
	 */
	@Override
	public double getLogPriorTerm() throws Exception {
		if( getMethod() == LearningType.BMA ) {
			throw new UnsupportedOperationException();
		}
		return super.getLogPriorTerm();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.trainableStatisticalModels.AbstractTrainSM#getMaximalMarkovOrder()
	 */
	@Override
	public byte getMaximalMarkovOrder() {
		return (byte)Math.min( length - 1, (Byte)params.getParameterAt( 3 ).getValue() );
	}

	private LearningType getMethod() {
		return (LearningType)params.getParameterAt( 4 ).getValue();
	}

	private ModelType getModel() {
		return (ModelType)params.getParameterAt( 2 ).getValue();
	}

	/**
	 * Counts the occurrence of the different indegrees and checks if the
	 * conventions are met.
	 * 
	 * @param structure
	 *            the structure
	 * @param maxOrder
	 *            the maximal order
	 * 
	 * @return an <code>int</code>-array containing the occurrence of indegrees
	 */
	protected int[] count( int[][] structure, byte maxOrder ) {
		int[] counts = new int[maxOrder + 1], help;
		int i = 0, j;
		for( ; i < length; i++ ) {
			if( i != structure[i][structure[i].length - 1] ) {
				throw new IllegalArgumentException( "The structure is not correct. Look at the parents of node " + i );
			}
			if( structure[i].length - 1 > maxOrder ) {
				throw new IllegalArgumentException( "The structure is not correct. There is at least one node (" + i
													+ ") that has to many parents." );
			}
			counts[structure[i].length - 1]++;
			help = structure[i].clone();
			Arrays.sort( help );
			j = 1;
			while( j < help.length && help[j - 1] < help[j] ) {
				j++;
			}
			if( j != help.length ) {
				throw new IllegalArgumentException( "The structure is not correct. Look at the parents of node " + i );
			}
		}
		return counts;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous.InhomogeneousDGTrainSM#set(de.jstacs.sequenceScores.statisticalModels.trainable.discrete.DGTrainSMParameterSet, boolean)
	 */
	@Override
	protected void set( DGTrainSMParameterSet parameter, boolean trained ) throws CloneNotSupportedException, NonParsableException {
		super.set( parameter, trained );
		sl = new StructureLearner( alphabets, length, getESS() );
	}
}
