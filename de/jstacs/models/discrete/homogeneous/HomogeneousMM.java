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

package de.jstacs.models.discrete.homogeneous;

import java.util.Map;
import java.util.Random;
import java.util.TreeMap;

import de.jstacs.NonParsableException;
import de.jstacs.NotTrainedException;
import de.jstacs.WrongAlphabetException;
import de.jstacs.data.DataSet;
import de.jstacs.data.Sequence;
import de.jstacs.data.DataSet.ElementEnumerator;
import de.jstacs.data.sequences.IntSequence;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.io.XMLParser;
import de.jstacs.models.discrete.DGMParameterSet;
import de.jstacs.models.discrete.homogeneous.parameters.HomMMParameterSet;
import de.jtem.numericalMethods.calculus.specialFunctions.Gamma;

/**
 * This class implements homogeneous Markov models (hMM) of arbitrary order.
 * 
 * @author Jens Keilwagen
 * 
 * @see HomMMParameterSet
 */
public class HomogeneousMM extends HomogeneousModel {

	private HomCondProb[] condProb;

	/**
	 * Creates a new homogeneous Markov model from a parameter set.
	 * 
	 * @param params
	 *            the parameter set
	 * 
	 * @throws CloneNotSupportedException
	 *             if the parameter set could not be cloned
	 * @throws IllegalArgumentException
	 *             if the parameter set is not instantiated
	 * @throws NonParsableException
	 *             if the parameter set is not parsable
	 * 
	 * @see HomMMParameterSet
	 * @see HomogeneousModel#HomogeneousModel(de.jstacs.models.discrete.homogeneous.parameters.HomogeneousModelParameterSet)
	 */
	public HomogeneousMM( HomMMParameterSet params ) throws CloneNotSupportedException, IllegalArgumentException, NonParsableException {
		super( params );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link HomogeneousMM} out of its XML representation.
	 * 
	 * @param stringBuff
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link HomogeneousMM} could not be reconstructed out
	 *             of the XML representation (the {@link StringBuffer} could not
	 *             be parsed)
	 * 
	 * @see de.jstacs.Storable
	 */
	public HomogeneousMM( StringBuffer stringBuff ) throws NonParsableException {
		super( stringBuff );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.discrete.DiscreteGraphicalModel#clone()
	 */
	@Override
	public HomogeneousMM clone() throws CloneNotSupportedException {
		HomogeneousMM clone = (HomogeneousMM)super.clone();
		clone.condProb = clone.cloneHomProb( condProb );
		return clone;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.discrete.homogeneous.HomogeneousModel#getRandomSequence(java.util.Random, int)
	 */
	@Override
	protected Sequence getRandomSequence( Random r, int length ) throws WrongAlphabetException, WrongSequenceTypeException {
		int[] seq = new int[length];
		int j = 0, val = 0;
		for( j = 0; j < order && j < length; j++ ) {
			seq[j] = chooseFromDistr( condProb[j], val, val + powers[1] - 1, r.nextDouble() );
			val = ( val + seq[j] ) * powers[1];
		}
		while( j < length ) {
			seq[j] = chooseFromDistr( condProb[order], val, val + powers[1] - 1, r.nextDouble() );
			val = ( ( val + seq[j++] ) % powers[order] ) * powers[1];
		}
		return new IntSequence( alphabets, seq );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.Model#getInstanceName()
	 */
	public String getInstanceName() {
		return "hMM(" + getMaximalMarkovOrder() + ") " + ( getESS() == 0 ? "ML" : "MAP" );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.Model#getLogPriorTerm()
	 */
	public double getLogPriorTerm() throws Exception {
		if( !trained ) {
			throw new NotTrainedException();
		}
		double p = 0, pot, ess = getESS();
		if( ess != 0 ) {
			int counter1 = 0, counter2, anz1;
			for( ; counter1 < condProb.length; counter1++ ) {
				anz1 = condProb[counter1].getNumberOfSpecificConstraints();
				pot = ess / (double)anz1;
				p += anz1 * ( Gamma.logOfGamma( powers[1] * pot ) / powers[1] - Gamma.logOfGamma( pot ) );
				for( counter2 = 0; counter2 < anz1; counter2++ ) {
					p += pot * condProb[counter1].getLnFreq( counter2 );
				}
			}
		}
		return p;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.discrete.homogeneous.HomogeneousModel#logProbFor(de.jstacs.data.Sequence, int, int)
	 */
	@Override
	protected double logProbFor( Sequence sequence, int startpos, int endpos ) {
		if( endpos < startpos ) {
			return 0;
		}
		int idx = sequence.discreteVal( startpos++ );
		double erg = condProb[0].getLnFreq( idx );
		for( int i = 1; i < order && startpos <= endpos; i++ ) {
			idx = idx * powers[1] + sequence.discreteVal( startpos++ );
			erg += condProb[i].getLnFreq( idx );
		}
		while( startpos <= endpos ) {
			idx = ( idx % powers[order] ) * powers[1] + sequence.discreteVal( startpos++ );
			erg += condProb[order].getLnFreq( idx );
		}
		return erg;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.discrete.DiscreteGraphicalModel#toString()
	 */
	@Override
	public String toString() {
		String erg = "description: " + getDescription();
		if( trained ) {
			erg += "\n\nprobabilities:\n";
			for( int i = 0; i <= order; i++ ) {
				erg += condProb[i].getFreqInfo( alphabets ) + "\n";
			}
		}
		return erg;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.Model#train(de.jstacs.data.Sample, double[])
	 */
	public void train( DataSet data, double[] weights ) throws Exception {
		train( new DataSet[]{ data }, new double[][]{ weights } );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.discrete.homogeneous.HomogeneousModel#train(de.jstacs.data.Sample[], double[][])
	 */
	@Override
	public void train( DataSet[] data, double[][] weights ) throws Exception {
		// check
		if( data.length != weights.length ) {
			throw new IllegalArgumentException( "The constraint data.length == weights.length is not fulfilled." );
		}

		int i = 0;

		// reset container of counter
		while( i < condProb.length ) {
			condProb[i++].reset();
		}

		// count
		for( i = 0; i < data.length; i++ ) {
			if( data[i] != null ) {
				countHomogeneous( data[i], weights[i] );
			}
		}
		// estimate
		estimate();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.discrete.DiscreteGraphicalModel#getFurtherModelInfos()
	 */
	@Override
	protected StringBuffer getFurtherModelInfos() {
		if( condProb != null ) {
			StringBuffer xml = new StringBuffer( 1000 );
			XMLParser.appendObjectWithTags(xml, condProb, "condProb");
			/*TODO remove
			int l = condProb.length;
			StringBuffer source = new StringBuffer( 25 + l * 500 );
			XMLParser.appendObjectWithTags( source, l, "length" );
			for( int j = 0; j < l; j++ ) {
				XMLParser.appendObjectWithTagsAndAttributes( source, condProb[j].toXML().toString(), "pos", "val=\"" + j + "\"" );
			}
			XMLParser.addTags( source, "condProb" );
			xml.append( source );
			*/
			return xml;
		} else {
			return null;
		}
	}

	private static final String XML_TAG = "HomogeneousMarkovModel";

	/* (non-Javadoc)
	 * @see de.jstacs.models.discrete.DiscreteGraphicalModel#getXMLTag()
	 */
	@Override
	protected String getXMLTag() {
		return XML_TAG;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.discrete.homogeneous.HomogeneousModel#set(de.jstacs.models.discrete.DGMParameterSet, boolean)
	 */
	@Override
	protected void set( DGMParameterSet params, boolean trained ) throws CloneNotSupportedException, NonParsableException {
		super.set( params, trained );
		if( !trained ) {
			byte j, i = 0;
			int k;
			condProb = new HomCondProb[order + 1];
			int[] pos;
			for( ; i < condProb.length; i++ ) {
				pos = new int[i + 1];
				for( j = 0; j <= i; j++ ) {
					pos[j] = j;
				}
				k = powers[i] * powers[1];
				condProb[i] = new HomCondProb( pos, k );
			}
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.discrete.DiscreteGraphicalModel#setFurtherModelInfos(java.lang.StringBuffer)
	 */
	@Override
	protected void setFurtherModelInfos( StringBuffer xml ) throws NonParsableException {
		if( trained ) {
			condProb = XMLParser.extractObjectAndAttributesForTags(xml, "condProb", null, null, HomCondProb[].class, HomogeneousModel.class, this);
			/*TODO remove
			StringBuffer help = new StringBuffer(XMLParser.extractObjectForTags( xml, "condProb", String.class ) );
			condProb = new HomCondProb[XMLParser.extractObjectForTags( help, "length", int.class )];			Map<String, String> filter = new TreeMap<String, String>();
			for( int j = 0; j < condProb.length; j++ ) {
				filter.clear();
				filter.put( "val", ""+j );
				condProb[j] = new HomCondProb(  new StringBuffer(XMLParser.extractObjectAndAttributesForTags( help, "pos", null, filter, String.class)) );
			}
			*/
		}
	}

	/**
	 * Counts homogeneously the symbols in the {@link DataSet} for this
	 * {@link HomogeneousMM}.
	 * 
	 * @param data
	 *            the {@link DataSet}
	 * @param weights
	 *            the weights of the {@link Sequence} in the {@link DataSet} or
	 *            <code>null</code>
	 * 
	 * @throws IllegalArgumentException
	 *             if the weights do not have the right dimension
	 * @throws WrongAlphabetException
	 *             if the alphabets do not match
	 */
	private void countHomogeneous( DataSet data, double[] weights ) throws WrongAlphabetException {
		int d = data.getNumberOfElements(), counter1, lengthCounter, l;
		Sequence seq;

		// check some constraints
		if( weights != null && d != weights.length ) {
			throw new IllegalArgumentException( "The weights are not suitable for the data (wrong dimension)." );
		}
		if( !alphabets.checkConsistency( data.getAlphabetContainer() ) ) {
			throw new WrongAlphabetException( "The alphabets of the model and the Sample are not suitable." );
		}

		// fill the constraints with the absolute frequency in the data
		int idx;
		ElementEnumerator ei = new ElementEnumerator( data );
		double w = 1;
		for( counter1 = 0; counter1 < d; counter1++ ) {
			seq = ei.nextElement();
			l = seq.getLength();
			idx = 0;
			if( weights != null ) {
				w = weights[counter1];
			}
			for( lengthCounter = 0; lengthCounter < order && lengthCounter < l; lengthCounter++ ) {
				idx = idx * powers[1] + seq.discreteVal( lengthCounter );
				condProb[lengthCounter].add( idx, w );
			}
			condProb[order].addAll( seq, w, lengthCounter, idx );
			/*while( lengthCounter < l ) {
				idx = ( idx % powers[order] ) * powers[1] + seq.discreteVal( lengthCounter++ );
				condProb[order].add( idx, w );
			}*/
		}
	}

	private void estimate() {
		double ess = getESS();
		for( int i = 0; i < condProb.length; i++ ) {
			condProb[i].estimate( ess );
		}
		trained = true;
	}
}
