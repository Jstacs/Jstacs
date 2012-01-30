/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Jstacs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.sequenceScores.statisticalModels.trainable;

import de.jstacs.NotTrainedException;
import de.jstacs.data.DataSet;
import de.jstacs.data.DataSet.WeightedDataSetFactory;
import de.jstacs.data.DataSet.WeightedDataSetFactory.SortOperation;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.results.NumericalResultSet;

/**
 * This class allows to train any {@link TrainableStatisticalModel} on {@link DataSet}s of {@link Sequence}s with
 * variable length if each individual length is at least {@link TrainableStatisticalModel#getLength()}. All other methods
 * are piped to the internally used {@link TrainableStatisticalModel}.
 * 
 * <br><br>
 * 
 * This class might be useful in any {@link de.jstacs.classifiers.assessment.ClassifierAssessment}.
 * 
 * @author Jens Keilwagen
 * 
 * @see de.jstacs.data.DataSet.WeightedDataSetFactory#DataSet.WeightedDataSetFactory(DataSet.WeightedDataSetFactory.SortOperation, DataSet, double[], int)
 * @see de.jstacs.classifiers.assessment.ClassifierAssessment
 */
public class VariableLengthWrapperTrainSM extends AbstractTrainableStatisticalModel {

	private TrainableStatisticalModel m;
	
	/**
	 * This is the main constructor that creates an instance from any {@link TrainableStatisticalModel}.
	 * 
	 * @param m the model
	 * 
	 * @throws CloneNotSupportedException if the mode <code>m</code> could not be cloned
	 */
	public VariableLengthWrapperTrainSM( TrainableStatisticalModel m ) throws CloneNotSupportedException {
		super( m.getAlphabetContainer(), m.getLength() );
		this.m = m.clone();
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link VariableLengthWrapperTrainSM} out of a {@link StringBuffer}.
	 * 
	 * @param stringBuff
	 *            the {@link StringBuffer} to be parsed
	 * 
	 * @throws NonParsableException
	 *             is thrown if the {@link StringBuffer} could not be parsed
	 */
	public VariableLengthWrapperTrainSM( StringBuffer stringBuff ) throws NonParsableException {
		super( stringBuff );
	}
	
	public VariableLengthWrapperTrainSM clone() throws CloneNotSupportedException {
		VariableLengthWrapperTrainSM clone = (VariableLengthWrapperTrainSM) super.clone();
		clone.m = m.clone();
		return clone;
	}

	@Override
	protected void fromXML( StringBuffer xml ) throws NonParsableException {
		StringBuffer content = XMLParser.extractForTag( xml, getInstanceName() );
		m = XMLParser.extractObjectForTags( content, "model", TrainableStatisticalModel.class );
		alphabets = m.getAlphabetContainer();
		length = m.getLength();
	}

	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, m, "model" );
		XMLParser.addTags( xml, getInstanceName() );
		return xml;
	}

	public String getInstanceName() {
		return m.getInstanceName();
	}

	public double getLogPriorTerm() throws Exception {
		return m.getLogPriorTerm();
	}

	public NumericalResultSet getNumericalCharacteristics() throws Exception {
		return m.getNumericalCharacteristics();
	}

	public double getLogProbFor( Sequence sequence, int startpos, int endpos ) throws NotTrainedException, Exception {
		return m.getLogProbFor( sequence, startpos, endpos );
	}

	public boolean isInitialized() {
		return m.isInitialized();
	}

	public void train( DataSet data, double[] weights ) throws Exception {
		WeightedDataSetFactory wsf = new WeightedDataSetFactory(SortOperation.NO_SORT,data,weights,length);
		m.train( wsf.getDataSet(), wsf.getWeights() );
	}
}
