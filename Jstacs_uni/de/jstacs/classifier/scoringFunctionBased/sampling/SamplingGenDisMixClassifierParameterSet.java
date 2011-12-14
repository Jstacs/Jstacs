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

package de.jstacs.classifier.scoringFunctionBased.sampling;

import de.jstacs.DataType;
import de.jstacs.classifier.scoringFunctionBased.gendismix.LogGenDisMixFunction;
import de.jstacs.classifier.scoringFunctionBased.sampling.SamplingScoreBasedClassifier.SamplingScheme;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.NumberValidator;

/**
 * {@link ParameterSet} to instantiate a {@link SamplingGenDisMixClassifier}.
 * @author Jan Grau
 *
 */
public class SamplingGenDisMixClassifierParameterSet extends SamplingScoreBasedClassifierParameterSet {

	/**
	 * Create a new {@link SamplingGenDisMixClassifierParameterSet}.
	 * @param alphabet the alphabet of the classifier
	 * @param length the length of the sequences that can be classified
	 * @param numStarts the number of independent samplings
	 * @param scheme the sampling scheme
	 * @param testSamplings the number of samplings between checks for the lenght of the burn-in phase
	 * @param stationarySamplings the number of samplings in the stationary phase
	 * @param freeParameters if only free parameters shall be used
	 * @param adaptVariance if the variance shall be adapted to the size of the event space of a random variable
	 * @param outfilePrefix the prefix of the temporary files for storing the parameters
	 * @param threads the number of threads for evaluating the {@link LogGenDisMixFunction}
	 * @throws Exception if the parameters could not be created
	 */
	public SamplingGenDisMixClassifierParameterSet( AlphabetContainer alphabet, int length, int numStarts, SamplingScheme scheme,
			int testSamplings, int stationarySamplings, boolean freeParameters,
			boolean adaptVariance, String outfilePrefix, int threads ) throws Exception {
		this( SamplingGenDisMixClassifier.class,
				alphabet,
				length,
				numStarts,
				scheme,
				testSamplings,
				stationarySamplings,
				freeParameters,
				adaptVariance,
				outfilePrefix, threads );
	}
	
	
	/**
	 * Create a new {@link SamplingGenDisMixClassifierParameterSet}.
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
	 * @param threads the number of threads for evaluating the {@link LogGenDisMixFunction}
	 * @throws Exception if the parameters could not be created
	 */
	protected SamplingGenDisMixClassifierParameterSet( Class<? extends SamplingScoreBasedClassifier> instanceClass,
													AlphabetContainer alphabet, int length, int numStarts, SamplingScheme scheme,
													int testSamplings, int stationarySamplings, boolean freeParameters,
													boolean adaptVariance, String outfilePrefix, int threads ) throws Exception {
		super( instanceClass,
				alphabet,
				length,
				numStarts,
				scheme,
				testSamplings,
				stationarySamplings,
				freeParameters,
				adaptVariance,
				outfilePrefix );
		parameters.add( new SimpleParameter( DataType.INT, "Threads", "The number of threads used for computation", true, new NumberValidator<Integer>( 1, 128 ),1 ) );
		this.getParameterForName( "Threads" ).setValue( threads );
	}

	/**
	 * Create a new {@link SamplingGenDisMixClassifierParameterSet} with a grouped sampling scheme, sampling all parameters
	 * (and not only the free ones), and adaption of the variance.
	 * @param alphabet the alphabet of the classifier
	 * @param length the length of the sequences that can be classified
	 * @param numStarts the number of independent samplings
	 * @param testSamplings the number of samplings between checks for the lenght of the burn-in phase
	 * @param stationarySamplings the number of samplings in the stationary phase
	 * @param outfilePrefix the prefix of the temporary files for storing the parameters
	 * @param threads the number of threads for evaluating the {@link LogGenDisMixFunction}
	 * @throws Exception if the parameters could not be created
	 */
	public SamplingGenDisMixClassifierParameterSet( AlphabetContainer alphabet, int length, int numStarts, int testSamplings,
													int stationarySamplings, String outfilePrefix, int threads ) throws Exception {
		super( SamplingGenDisMixClassifier.class, alphabet, length, numStarts, testSamplings, stationarySamplings, outfilePrefix );
		this.getParameterAt( parameters.size()-1 ).setValue( threads );
	}
	
	
	
	@Override
	public String getInstanceComment() {
		return null;
	}

	@Override
	public String getInstanceName() {
		return null;
	}
	
	/**
	 * Returns the number of threads for evaluating the {@link LogGenDisMixFunction}
	 * @return the number of threads
	 */
	public int getNumberOfThreads() {
		return (Integer) this.getParameterAt( parameters.size()-1 ).getValue();
	}

}
