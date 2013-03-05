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

package de.jstacs.motifDiscovery;

import de.jstacs.data.DataSet;

/**
 * This is the interface that any tool for de-novo motif discovery should implement that allows any modify-operations like shift, shrink and expand.
 * These operations are possible if the motif is mutable.  
 * 
 * @author Jan Grau, Jens Keilwagen
 * 
 * @see de.jstacs.motifDiscovery.Mutable
 */
public interface MutableMotifDiscoverer extends MotifDiscoverer {

	/**
	 * Manually modifies the motif model with index <code>motifIndex</code>. The two offsets <code>offsetLeft</code> and <code>offsetRight</code>
	 * define how many positions the left or right border positions shall be moved. Negative numbers indicate moves to the left while positive
	 * numbers correspond to moves to the right. The distribution for sequences to the left and right side of the motif shall be computed internally.
	 * 
	 * @param motifIndex the index of the motif in the motif discoverer
	 * @param offsetLeft the offset on the left side
	 * @param offsetRight the offset on the right side
	 * 
	 * @return <code>true</code> if the motif model was modified otherwise <code>false</code>
	 * 
	 * @throws Exception if some unexpected error occurred during the modification
	 * 
	 * @see MutableMotifDiscoverer#modifyMotif(int, int, int)
	 * @see Mutable#modify(int, int)
	 */
	public boolean modifyMotif( int motifIndex, int offsetLeft, int offsetRight ) throws Exception;
	
	/**
	 * This method allows to initialize the model of a motif manually using a weighted data set.
	 * 
	 * @param motifIndex the index of the motif in the motif discoverer
	 * @param data the data set of sequences
	 * @param weights either <code>null</code> or an array of length <code>data.getNumberofElements()</code> with non-negative weights.
	 * 
	 * @throws Exception if initialize was not possible
	 */
	public void initializeMotif( int motifIndex, DataSet data, double[] weights ) throws Exception;
	
	/**
	 * This method initializes the motif with index <code>motif</code> randomly using for instance {@link de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#initializeFunctionRandomly(boolean)}.
	 * Furthermore, if available, it also initializes the positional distribution.
	 *  
	 * @param motif the index of the motif
	 * 
	 * @throws Exception either if the index is wrong or if it is thrown by the method {@link de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel#initializeFunctionRandomly(boolean)}
	 */
	public void initializeMotifRandomly( int motif ) throws Exception;
	
	/**
	 * Adjusts all hidden parameters including duration and mixture parameters according to the current values of the remaining parameters.
	 * 
	 * @param index the index of the class of this {@link MutableMotifDiscoverer}
	 * @param data the array of data for all classes
	 * @param weights the weights for all sequences in data
	 * 
	 * @throws Exception thrown if the hidden parameters could not be adjusted
	 */
	public void adjustHiddenParameters( int index, DataSet[] data, double[][] weights ) throws Exception;
}
