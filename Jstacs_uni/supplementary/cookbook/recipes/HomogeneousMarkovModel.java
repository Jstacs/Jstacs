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

package supplementary.cookbook.recipes;

import java.text.NumberFormat;
import java.util.Arrays;

import de.jstacs.NotTrainedException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.sequenceScores.statisticalModels.trainable.AbstractTrainableStatisticalModel;
import de.jstacs.utils.ToolBox;
 
 
 
public class HomogeneousMarkovModel extends AbstractTrainableStatisticalModel {
 
	private double[] logProbs;//array for the parameters, i.e. the probabilities for each symbol
 
	public HomogeneousMarkovModel( AlphabetContainer alphabets ) throws Exception {
		super( alphabets, 0 ); //we have a homogeneous TrainableStatisticalModel, hence the length is set to 0
		//a homogeneous TrainableStatisticalModel can only handle simple alphabets
		if(! (alphabets.isSimple() && alphabets.isDiscrete()) ){
			throw new Exception("Only simple and discrete alphabets allowed");
		}
		//initialize parameter array
		this.logProbs = new double[(int) alphabets.getAlphabetLengthAt( 0 )];
		Arrays.fill( logProbs, -Math.log(logProbs.length) );
	}
 
	public HomogeneousMarkovModel( StringBuffer stringBuff ) throws NonParsableException { 
        super( stringBuff ); 
    }
 
	protected void fromXML( StringBuffer xml ) throws NonParsableException {
		//extract our XML-code
		xml = XMLParser.extractForTag( xml, "homogeneousMarkovModel" );
		//extract all the variables using XMLParser
		alphabets = (AlphabetContainer) XMLParser.extractObjectForTags( xml, "alphabets" );
		length = XMLParser.extractObjectForTags( xml, "length", int.class );
		logProbs = XMLParser.extractObjectForTags( xml, "logProbs", double[].class );
	}
 
	public StringBuffer toXML() {
		StringBuffer buf = new StringBuffer();
		//pack all the variables using XMLParser
		XMLParser.appendObjectWithTags( buf, alphabets, "alphabets" );
		XMLParser.appendObjectWithTags( buf, length, "length" );
		XMLParser.appendObjectWithTags( buf, logProbs, "logProbs" );
		//add our own tag
		XMLParser.addTags( buf, "homogeneousMarkovModel" );
		return buf;
	}
 
	public String getInstanceName() { 
            return "Homogeneous Markov model of order 0"; 
        }
 
	public double getLogPriorTerm() throws Exception { 
            //we use ML-estimation, hence no prior term
            return 0; 
        } 
 
	public NumericalResultSet getNumericalCharacteristics() throws Exception {
		//we do not have much to tell here
		return new NumericalResultSet(new NumericalResult("Number of parameters","The number of parameters this model uses",logProbs.length));
	}
 
	public double getLogProbFor( Sequence sequence, int startpos, int endpos ) throws NotTrainedException, Exception {
		double seqLogProb = 0.0;
		//compute the log-probability of the sequence between startpos and endpos (inclusive)
		//as sum of the single symbol log-probabilities
		for(int i=startpos;i<=endpos;i++){
			//directly access the array by the numerical representation of the symbols
			seqLogProb += logProbs[sequence.discreteVal( i )];
		}
		return seqLogProb;
	}
 
	public boolean isInitialized() {
        return true; 
    }
 
	public void train( DataSet data, double[] weights ) throws Exception {
		//reset the parameter array
		Arrays.fill( logProbs, 0.0 );
		//default sequence weight
		double w = 1;
		//for each sequence in the data set
		for(int i=0;i<data.getNumberOfElements();i++){
			//retrieve sequence
			Sequence seq = data.getElementAt( i );
			//if we do have any weights, use them
			if(weights != null){
				w = weights[i];
			}
			//for each position in the sequence
			for(int j=0;j<seq.getLength();j++){
				//count symbols, weighted by weights
				logProbs[ seq.discreteVal( j ) ] += w;
			}
		}
		//compute normalization
		double norm = 0.0;
		for(int i=0;i<logProbs.length;i++){ norm += logProbs[i]; }
		//normalize probs to obtain proper probabilities
		for(int i=0;i<logProbs.length;i++){ logProbs[i] = Math.log( logProbs[i]/norm ); }
	}

	@Override
	public String toString(NumberFormat nf) {
		return ToolBox.toString(logProbs, nf);
	} 
}