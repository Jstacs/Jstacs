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
package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.emissions.discrete;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.ReferenceSequenceAnnotation;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;

/**
 * This class implements a discrete emission that depends on some {@link ReferenceSequenceAnnotation}
 * at a certain reference position. This class can be used in so-called conditional profile HMM for
 * the match states.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class ReferenceSequenceDiscreteEmission extends AbstractConditionalDiscreteEmission {

	private AlphabetContainer refCon;
	/**
	 * The index in the reference sequence
	 */
	protected int refIdx;
	
	/**
	 * Returns the reference sequence annotated to <code>seq</code>.
	 * @param seq the sequence
	 * @return the reference sequence
	 * @see ReferenceSequenceAnnotation
	 */
	protected static final Sequence getReferenceSequence(Sequence seq){
		return ((ReferenceSequenceAnnotation) seq.getSequenceAnnotationByType( "reference", 0 )).getReferenceSequence();
	}
	
	/**
	 * This is a simple constructor for a {@link ReferenceSequenceDiscreteEmission} based on the equivalent sample size.
	 * 
	 * @param con the {@link AlphabetContainer} of this emission
	 * @param refCon the {@link AlphabetContainer} of the reference
	 * @param refIdx the index in the reference sequence
	 * @param ess the equivalent sample size (ess) of this emission that is equally distributed over all parameters
	 * 
	 * @throws WrongAlphabetException if the {@link AlphabetContainer} is not discrete or simple
	 * 
	 * @see #ReferenceSequenceDiscreteEmission(AlphabetContainer, AlphabetContainer, int, double[][])
	 */
	public ReferenceSequenceDiscreteEmission( AlphabetContainer con, AlphabetContainer refCon, int refIdx, double ess ) throws WrongAlphabetException {
		this( con, refCon, refIdx, getHyperParams(ess, (int) refCon.getAlphabetLengthAt(refIdx), (int) con.getAlphabetLengthAt(0)) );
	}
	
	/**
	 * This is a simple constructor for a {@link ReferenceSequenceDiscreteEmission} based on the equivalent sample size.
	 * 
	 * @param con the {@link AlphabetContainer} of this emission
	 * @param refCon the {@link AlphabetContainer} of the reference
	 * @param refIdx the index in the reference sequence
	 * @param ess the equivalent sample size (ess) of this emission that is equally distributed over all parameters
	 * @param initHyperParams the individual hyper parameters for each parameter used in {@link #initializeFunctionRandomly()}
	 * 
	 * @throws WrongAlphabetException if the {@link AlphabetContainer} is not discrete or simple
	 * 
	 * @see #ReferenceSequenceDiscreteEmission(AlphabetContainer, AlphabetContainer, int, double[][])
	 */
	public ReferenceSequenceDiscreteEmission( AlphabetContainer con, AlphabetContainer refCon, int refIdx, double ess, double[][] initHyperParams ) throws WrongAlphabetException {
		this( con, refCon, refIdx, getHyperParams(ess, (int) refCon.getAlphabetLengthAt(refIdx), (int) con.getAlphabetLengthAt(0)), initHyperParams );
	}

	/**
	 * This constructor creates a {@link ReferenceSequenceDiscreteEmission} defining the individual hyper parameters. 
	 * 
	 * @param con the {@link AlphabetContainer} of this emission
	 * @param refCon the {@link AlphabetContainer} of the reference
	 * @param refIdx the index in the reference sequence
	 * @param hyperParams the individual hyper parameters for each parameter
	 * 
	 * @throws IllegalArgumentException if the dimension of the hyper parameters and the size of the alphabet defined by reference {@link AlphabetContainer} do not match
	 * @throws WrongAlphabetException if the {@link AlphabetContainer} is not discrete or simple
	 */
	public ReferenceSequenceDiscreteEmission( AlphabetContainer con, AlphabetContainer refCon, int refIdx, double[][] hyperParams ) throws IllegalArgumentException, WrongAlphabetException {
		super( con, hyperParams );
		if(refCon.getAlphabetLengthAt( refIdx ) != hyperParams.length){
			throw new IllegalArgumentException("Hyper-parameters do not match length of alphabet");
		}
		this.refCon = refCon;
		this.refIdx = refIdx;
	}
	
	/**
	 * This constructor creates a {@link ReferenceSequenceDiscreteEmission} defining the individual hyper parameters. 
	 * 
	 * @param con the {@link AlphabetContainer} of this emission
	 * @param refCon the {@link AlphabetContainer} of the reference
	 * @param refIdx the index in the reference sequence
	 * @param hyperParams the individual hyper parameters for each parameter
	 * @param initHyperParams the individual hyper parameters for each parameter used in {@link #initializeFunctionRandomly()}
	 * 
	 * @throws IllegalArgumentException if the dimension of the hyper parameters and the size of the alphabet defined by reference {@link AlphabetContainer} do not match
	 * @throws WrongAlphabetException if the {@link AlphabetContainer} is not discrete or simple
	 */
	public ReferenceSequenceDiscreteEmission( AlphabetContainer con, AlphabetContainer refCon, int refIdx, double[][] hyperParams, double[][] initHyperParams ) throws IllegalArgumentException, WrongAlphabetException {
		super( con, hyperParams,initHyperParams );
		if(refCon.getAlphabetLengthAt( refIdx ) != hyperParams.length){
			throw new IllegalArgumentException("Hyper-parameters do not match length of alphabet");
		}
		this.refCon = refCon;
		this.refIdx = refIdx;
	}

	/**
	 * Creates a {@link ReferenceSequenceDiscreteEmission} from its XML representation.
	 * @param xml the XML representation.
	 * @throws NonParsableException if the XML representation could not be parsed
	 */
	public ReferenceSequenceDiscreteEmission( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	@Override
	protected int getConditionIndex( boolean forward, int seqPos, Sequence seq ) {
		Sequence ref = getReferenceSequence( seq );
		if(refIdx < ref.getLength()){
			return ref.discreteVal( refIdx );
		}else{
			return -1;
		}
	}

	@Override
	protected String getCondition(int i) {
		return " | R_" + refIdx + "=" + refCon.getSymbol(0, i);
	}

	@Override
	protected void appendFurtherInformation( StringBuffer xml ) {
		XMLParser.appendObjectWithTags( xml, refCon, "refCon" );
		XMLParser.appendObjectWithTags( xml, refIdx, "refIdx" );
	}

	@Override
	protected void extractFurtherInformation( StringBuffer xml ) throws NonParsableException {
		refCon = (AlphabetContainer) XMLParser.extractObjectForTags( xml, "refCon" );
		refIdx = XMLParser.extractObjectForTags( xml, "refIdx", int.class );
	}
}