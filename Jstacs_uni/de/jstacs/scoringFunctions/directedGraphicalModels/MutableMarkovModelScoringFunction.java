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

package de.jstacs.scoringFunctions.directedGraphicalModels;

import java.util.ArrayList;

import de.jstacs.NonParsableException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sample;
import de.jstacs.io.XMLParser;
import de.jstacs.motifDiscovery.Mutable;
import de.jstacs.scoringFunctions.SamplingScoringFunction;
import de.jstacs.scoringFunctions.directedGraphicalModels.structureLearning.measures.InhomogeneousMarkov;
import de.jstacs.scoringFunctions.mix.motifSearch.DurationScoringFunction;

/**
 * This class implements a {@link de.jstacs.scoringFunctions.AbstractNormalizableScoringFunction} for an inhomogeneous Markov model.
 * The modeled length can be modified which might be very important for de-novo motif discovery.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class MutableMarkovModelScoringFunction extends BayesianNetworkScoringFunction implements Mutable, SamplingScoringFunction{
	
	private DurationScoringFunction lengthPenalty;

	/**
	 * This constructor creates an instance with an prior for the modeled length.
	 * 
	 * @param alphabet the {@link AlphabetContainer} of the {@link MutableMarkovModelScoringFunction}
	 * @param length the initial length of the modeled sequences
	 * @param ess the equivalent sample size
	 * @param plugInParameters a switch whether to use plug-in parameters of not
	 * @param order the order of the Markov model
	 * @param lengthPenalty the prior on the modeled sequence length
	 * 
	 * @throws Exception if super class constructor throws an {@link Exception} or if the <code>lengthPenalty</code> does not allow the initial length
	 */
	public MutableMarkovModelScoringFunction(AlphabetContainer alphabet,
			int length, double ess, boolean plugInParameters,
			int order, DurationScoringFunction lengthPenalty ) throws Exception {
		this( alphabet, length, ess, plugInParameters, new InhomogeneousMarkov( order ), lengthPenalty );
	}
	
	/**
	 * This constructor creates an instance without any prior for the modeled length.
	 * 
	 * @param alphabet the {@link AlphabetContainer} of the {@link MutableMarkovModelScoringFunction}
	 * @param length the initial length of the modeled sequences
	 * @param ess the equivalent sample size
	 * @param plugInParameters a switch whether to use plug-in parameters of not
	 * @param structureMeasure an {@link InhomogeneousMarkov} {@link de.jstacs.scoringFunctions.directedGraphicalModels.structureLearning.measures.Measure} for the structure
	 * 
	 * @throws Exception if super class constructor throws an {@link Exception}
	 */
	public MutableMarkovModelScoringFunction( AlphabetContainer alphabet, int length, double ess, boolean plugInParameters,
			InhomogeneousMarkov structureMeasure ) throws Exception {
		super( alphabet, length, ess, plugInParameters, structureMeasure );
	}

	/**
	 * This constructor creates an instance with an prior for the modeled length.
	 * 
	 * @param alphabet the {@link AlphabetContainer} of the {@link MutableMarkovModelScoringFunction}
	 * @param length the initial length of the modeled sequences
	 * @param ess the equivalent sample size
	 * @param plugInParameters a switch whether to use plug-in parameters of not
	 * @param structureMeasure a {@link InhomogeneousMarkov} {@link de.jstacs.scoringFunctions.directedGraphicalModels.structureLearning.measures.Measure} for the structure
	 * @param lengthPenalty the prior on the modeled sequence length
	 * 
	 * @throws Exception if super class constructor throws an {@link Exception} or if the <code>lengthPenalty</code> does not allow the initial length
	 */
	public MutableMarkovModelScoringFunction( AlphabetContainer alphabet, int length, double ess, boolean plugInParameters,
			InhomogeneousMarkov structureMeasure, DurationScoringFunction lengthPenalty ) throws Exception {
		this( alphabet, length, ess, plugInParameters, structureMeasure );
		this.lengthPenalty = lengthPenalty;
		if( lengthPenalty!= null && !lengthPenalty.isPossible( length ) ) {
			throw new IllegalArgumentException( "This motif length is not possible: " + length );
		}
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Recreates a {@link MutableMarkovModelScoringFunction} from its XML
	 * representation as saved by the method {@link #toXML()}.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML code could not be parsed
	 */
	public MutableMarkovModelScoringFunction( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}	
	
	private static final String XML_TAG = "MutableMarkovModelScoringFunction";
	
	protected void fromXML(StringBuffer source) throws NonParsableException {
		StringBuffer sb = XMLParser.extractForTag( source, XML_TAG );
		lengthPenalty = XMLParser.extractObjectForTags( sb, "lengthPenalty", DurationScoringFunction.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		super.fromXML( sb );
	}
	
	public StringBuffer toXML() {
		StringBuffer sb = super.toXML();
		XMLParser.appendObjectWithTags( sb, lengthPenalty, "lengthPenalty" );
		XMLParser.addTags( sb, XML_TAG );
		return sb;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.scoringFunctions.directedGraphicalModels.BayesianNetworkScoringFunction#getLogPriorTerm()
	 */
	@Override
	public double getLogPriorTerm() {
		if(lengthPenalty != null){
			return super.getLogPriorTerm() + lengthPenalty.getLogScore( length );
		}else{
			return super.getLogPriorTerm();
		}
	}
	
	/**
	 * Returns the order of the inhomogeneous Markov model.
	 * @return the order
	 */
	public int getOrder(){
		return ((InhomogeneousMarkov)structureMeasure).getOrder();
	}

	public boolean modify( int offsetLeft, int offsetRight ) {

		if(! getAlphabetContainer().isSimple() ){
			return false;
		}
		if(offsetLeft == 0 && offsetRight == 0){
			return true;
		}else{
			this.precomputeNormalization();
			this.normalizeParameters();
			this.precomputeNormalization();
			
			ParameterTree[] backTrees = trees;
			int backLength = this.length;
			try{
				this.length = this.length - offsetLeft + offsetRight;
				if( lengthPenalty != null && !lengthPenalty.isPossible( length ) ) {
					throw new IllegalArgumentException( "This motif length is not possible: " + length );
				}
				this.createTrees( new Sample[]{null,null}, new double[][]{null,null} );
				int indexNew = 0, indexOld = 0;

				//left side
				if( offsetLeft >= 0 )
				{
					indexOld = offsetLeft;
				}
				else
				{
					indexNew = -offsetLeft;
				}

				//copy
				for( ;indexOld<backTrees.length && indexNew<trees.length; indexNew++, indexOld++ ){
					trees[indexNew].copy( backTrees[indexOld] );
				}

				//right side
				if( indexNew < trees.length ) {
					indexNew = trees.length;
				}
				logNormalizationConstant = null;
				return true;
			}catch(Exception e){
				this.length = backLength;
				trees = backTrees;
				logNormalizationConstant = null;
				return false;
			}
		}
	}

	
	private void normalizeParameters() {
		for(int i=0;i<trees.length;i++){
			trees[i].normalizeParameters();
		}
	}

	@Override
	public int[][] getSamplingGroups( int parameterOffset ) {
		ArrayList<int[]> list = new ArrayList<int[]>();
		for(int i=0;i<trees.length;i++){
			for(int j=0;j<trees[i].getNumberOfSamplingSteps();j++){
				list.add( trees[i].getParameterIndexesForSamplingStep( j, parameterOffset ) );
			}
			parameterOffset += trees[i].getNumberOfParameters();
		}
		return list.toArray( new int[0][0] );
	}
	
	
}
