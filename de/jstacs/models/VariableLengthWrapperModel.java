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

package de.jstacs.models;

import de.jstacs.NonParsableException;
import de.jstacs.NotTrainedException;
import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;
import de.jstacs.data.Sample.WeightedSampleFactory;
import de.jstacs.data.Sample.WeightedSampleFactory.SortOperation;
import de.jstacs.io.XMLParser;
import de.jstacs.results.NumericalResultSet;

/**
 * This class allows to train any {@link Model} on {@link Sample}s of {@link Sequence}s with
 * variable length if each individual length is at least {@link Model#getLength()}. All other methods
 * are piped to the internally used {@link Model}.
 * 
 * <br><br>
 * 
 * This class might be useful in any {@link de.jstacs.classifier.assessment.ClassifierAssessment}.
 * 
 * @author Jens Keilwagen
 * 
 * @see de.jstacs.data.Sample.WeightedSampleFactory#Sample.WeightedSampleFactory(SortOperation, Sample, double[], int)
 * @see de.jstacs.classifier.assessment.ClassifierAssessment
 */
public class VariableLengthWrapperModel extends AbstractModel {

	private Model m;
	
	/**
	 * This is the main constructor that creates an instance from any {@link Model}.
	 * 
	 * @param m the model
	 * 
	 * @throws CloneNotSupportedException if the mode <code>m</code> could not be cloned
	 */
	public VariableLengthWrapperModel( Model m ) throws CloneNotSupportedException {
		super( m.getAlphabetContainer(), m.getLength() );
		this.m = m.clone();
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link VariableLengthWrapperModel} out of a {@link StringBuffer}.
	 * 
	 * @param stringBuff
	 *            the {@link StringBuffer} to be parsed
	 * 
	 * @throws NonParsableException
	 *             is thrown if the {@link StringBuffer} could not be parsed
	 */
	public VariableLengthWrapperModel( StringBuffer stringBuff ) throws NonParsableException {
		super( stringBuff );
	}
	
	public VariableLengthWrapperModel clone() throws CloneNotSupportedException {
		VariableLengthWrapperModel clone = (VariableLengthWrapperModel) super.clone();
		clone.m = m.clone();
		return clone;
	}

	@Override
	protected void fromXML( StringBuffer xml ) throws NonParsableException {
		StringBuffer content = XMLParser.extractForTag( xml, getInstanceName() );
		m = XMLParser.extractObjectForTags( content, "model", Model.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
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

	public void train( Sample data, double[] weights ) throws Exception {
		WeightedSampleFactory wsf = new WeightedSampleFactory(SortOperation.NO_SORT,data,weights,length);
		m.train( wsf.getSample(), wsf.getWeights() );
	}
}
