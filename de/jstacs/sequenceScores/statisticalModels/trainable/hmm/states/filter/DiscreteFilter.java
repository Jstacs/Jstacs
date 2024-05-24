package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.filter;

import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;

/**
 * A filter based on the discrete values of {@link Sequence} and a position.
 * 
 * @author Jens Keilwagen
 */
public class DiscreteFilter implements Filter {

	/**
	 * The result for each possible discrete value
	 */
	private boolean[] result;
	/**
	 * The relative offset
	 */
	private int offset;
	
	/**
	 * The main constructor.
	 * 
	 * @param offset the relative offset from the anchor
	 * @param result the result for each possible discrete value
	 */
	public DiscreteFilter( int offset, boolean[] result ) {
		this.offset=offset;
		this.result = result.clone();
	}
	
	/**
	 * The constructor for the Storable interface.
	 * 
	 * @param xml the xml representation
	 * 
	 * @throws NonParsableException if an error occurs while parsing
	 * 
	 * @see XMLParser
	 */
	public DiscreteFilter( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag(xml, "DiscreteFilter");
		offset = (Integer) XMLParser.extractObjectForTags(xml, "offset");
		result = (boolean[]) XMLParser.extractObjectForTags(xml, "result");
		
	}
	
	@Override
	public boolean isAccepted( int anchor, Sequence seq ) {
		int pos = anchor+offset;
		if( pos >= 0 && pos < seq.getLength() ) {
			return result[seq.discreteVal(anchor+offset)];
		} else {
			return false;
		}
	}

	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags(xml, offset, "offset");
		XMLParser.appendObjectWithTags(xml, result, "result");
		XMLParser.addTags(xml, "DiscreteFilter");
		return xml;
	}
}