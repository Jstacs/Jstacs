package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.filter;

import de.jstacs.Storable;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;

/**
 * A filter that combines other filters.
 * 
 * @author Jens Keilwagen
 */
public class AndFilter implements Filter {

	private Filter[] filter;
	
	/**
	 * Main constructor.
	 * 
	 * @param filter individual filters to be combined
	 */
	public AndFilter( Filter... filter ) {
		this.filter=filter.clone();
	}
	
	/**
	 * The constructor for the {@link Storable} interface.
	 * 
	 * @param xml the xml representation of an object
	 * 
	 * @throws NonParsableException if an error occurs while parsing
	 */
	public AndFilter( StringBuffer xml ) throws NonParsableException {
		xml=XMLParser.extractForTag(xml, "AndFilter");
		filter = (Filter[]) XMLParser.extractObjectForTags(xml, "filter");
	}

	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags(xml, filter, "filter");
		XMLParser.addTags(xml, "AndFilter");
		return xml;
	}

	@Override
	public boolean isAccepted(int anchor, Sequence seq) {
		int i = 0;
		while( i < filter.length && filter[i].isAccepted(anchor, seq) ) {
			i++;
		}
		return i == filter.length;
	}
}