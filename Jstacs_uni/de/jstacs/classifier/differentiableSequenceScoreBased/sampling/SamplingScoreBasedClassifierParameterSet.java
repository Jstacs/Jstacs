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

package de.jstacs.classifier.differentiableSequenceScoreBased.sampling;

import de.jstacs.DataType;
import de.jstacs.classifier.differentiableSequenceScoreBased.sampling.SamplingScoreBasedClassifier.SamplingScheme;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.parameters.EnumParameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SequenceScoringParameterSet;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;

/**
 * {@link ParameterSet} to instantiate a {@link SamplingScoreBasedClassifier}.
 * @author Jan Grau
 *
 */
public abstract class SamplingScoreBasedClassifierParameterSet extends SequenceScoringParameterSet {
	
	/**
	 * Create a new {@link SamplingScoreBasedClassifierParameterSet}.
	 * @param instanceClass the class, which must be a subclass of {@link SamplingScoreBasedClassifier}
	 * @param alphabet the alphabet of the classifier
	 * @param length the length of the sequences that can be classified
	 * @param numStarts the number of independent samplings
	 * @param scheme the sampling scheme
	 * @param testSamplings the number of samplings between checks for the lenght of the burn-in phase
	 * @param stationarySamplings the number of samplings in the stationary phase
	 * @param freeParameters if only free parameters shall be used
	 * @param adaptVariance if the variance shall be adapted to the size of the event space of a random variable
	 * @param outfilePrefix the prefix of the temporary files for storing the parameters
	 * @throws Exception if the parameters could not be created
	 * @see SamplingScoreBasedClassifier
	 */
	public SamplingScoreBasedClassifierParameterSet(Class<? extends SamplingScoreBasedClassifier> instanceClass, AlphabetContainer alphabet, int length, 
			int numStarts, SamplingScheme scheme, int testSamplings, int stationarySamplings, boolean freeParameters, boolean adaptVariance, String outfilePrefix) throws Exception{
		this(instanceClass,alphabet,length,numStarts,testSamplings,stationarySamplings,outfilePrefix);
		this.getParameterAt( 1 ).setValue( scheme );
		this.getParameterAt( 4 ).setValue( freeParameters );
		this.getParameterAt( 5 ).setValue( adaptVariance );
	}
	
	/**
	 * Create a new {@link SamplingScoreBasedClassifierParameterSet} with a grouped sampling scheme, sampling all parameters
	 * (and not only the free ones), and adaption of the variance.
	 * @param instanceClass the class, which must be a subclass of {@link SamplingScoreBasedClassifier}
	 * @param alphabet the alphabet of the classifier
	 * @param length the length of the sequences that can be classified
	 * @param numStarts the number of independent samplings
	 * @param testSamplings the number of samplings between checks for the lenght of the burn-in phase
	 * @param stationarySamplings the number of samplings in the stationary phase
	 * @param outfilePrefix the prefix of the temporary files for storing the parameters
	 * @throws Exception if the parameters could not be created
	 * @see SamplingScoreBasedClassifier
	 */
	public SamplingScoreBasedClassifierParameterSet(Class<? extends SamplingScoreBasedClassifier> instanceClass, AlphabetContainer alphabet, int length, 
			int numStarts, int testSamplings, int stationarySamplings, String outfilePrefix) throws Exception{
		super( instanceClass, alphabet, length, length == 0 );
		
		parameters.add( new SimpleParameter( DataType.INT, "Number of starts", "Number of parallel sampling starts", true ) );
		parameters.add( new EnumParameter( SamplingScheme.class, "The sampling scheme", true, SamplingScheme.GROUPED.name() ));
		parameters.add( new SimpleParameter( DataType.INT, "Test samplings", "The number of sampling steps between to burn-in tests", true, 100 ) );
		parameters.add( new SimpleParameter( DataType.INT, "Stationary samplings", "The number of sampling steps in the stationary phase", true, 1000 ) );
		parameters.add( new SimpleParameter( DataType.BOOLEAN, "Free parameters", "Use only free parameters", true, false ) );
		parameters.add( new SimpleParameter( DataType.BOOLEAN, "Adapt variance", "Adapt the variance to the size of event spaces for each random variable", true, true ) );
		parameters.add( new SimpleParameter( DataType.STRING, "Outfile prefix", "The prefix of the outfiles where the parameters are stored", true ) );
		
		this.getParameterAt( 0 ).setValue( numStarts );
		this.getParameterAt( 2 ).setValue( testSamplings );
		this.getParameterAt( 3 ).setValue( stationarySamplings );
		this.getParameterAt( 6 ).setValue( outfilePrefix );
	}
	
	public SamplingScoreBasedClassifierParameterSet clone() throws CloneNotSupportedException{
		return (SamplingScoreBasedClassifierParameterSet) super.clone();
	}

	/**
	 * Returns the number of independent sampling starts
	 * @return the number of starts
	 */
	public int getNumberOfStarts(){
		return (Integer) this.getParameterAt( 0 ).getValue();
	}
	
	/**
	 * Returns the sampling scheme
	 * @return the scheme
	 */
	public SamplingScheme getSamplingScheme(){
		return (SamplingScheme) this.getParameterAt( 1 ).getValue();
	}
	
	/**
	 * Returns the number of samplings between checks for the stationary phase
	 * @return te number of samplings
	 */
	public int getNumberOfTestSamplings(){
		return (Integer) this.getParameterAt( 2 ).getValue();
	}
	
	/**
	 * Returns the number of samplings steps in the stationary phase
	 * @return the number of steps
	 */
	public int getNumberOfStationarySamplings(){
		return (Integer) this.getParameterAt( 3 ).getValue();
	}
	
	/**
	 * Returns <code>true</code> if only free parameters shall be used
	 * @return if only free parameters shall be used
	 */
	public boolean getFreeParameters(){
		return (Boolean) this.getParameterAt( 4 ).getValue();
	}
	
	/**
	 * Returns true if the sampling variance shall be adapted to the size 
	 * of the event space of a random variable
	 * @return if the sampling variance shall be adapted
	 * @see SamplingScoreBasedClassifier
	 */
	public boolean getAdaptVariance(){
		return (Boolean) this.getParameterAt( 5 ).getValue();
	}
	
	/**
	 * Returns the prefix of the temporary files for storing sampled
	 * parameter values
	 * @return the prefix
	 */
	public String getOutfilePrefix(){
		return (String) this.getParameterAt( 6 ).getValue();
	}

	/**
	 * Sets the number of starts to <code>i</code>
	 * @param i the new number of starts
	 * @throws IllegalValueException if this value is not allowed
	 */
	public void setNumberOfStarts( int i ) throws IllegalValueException {
		this.getParameterAt( 0 ).setValue( i );
	}
	
}
