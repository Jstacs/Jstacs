package de.jstacs.sequenceScores.statisticalModels.trainable.hmm.states.filter;

import de.jstacs.Storable;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;

/**
 * A filter that negates the result of another filter.
 * 
 * @author Jens Keilwagen
 */
public class NotFilter implements Filter {

	private Filter filter;
	
	/**
	 * Main constructor.
	 * 
	 * @param filter individual filters to be combined
	 */
	public NotFilter( Filter filter ) {
		this.filter=filter;
	}
	
	/**
	 * The constructor for the {@link Storable} interface.
	 * 
	 * @param xml the xml representation of an object
	 * 
	 * @throws NonParsableException if an error occurs while parsing
	 */
	public NotFilter( StringBuffer xml ) throws NonParsableException {
		xml=XMLParser.extractForTag(xml, "NotFilter");
		filter = (Filter) XMLParser.extractObjectForTags(xml, "filter");
	}

	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags(xml, filter, "filter");
		XMLParser.addTags(xml, "NotFilter");
		return xml;
	}
	
	@Override
	public boolean isAccepted(int anchor, Sequence seq) {
		return !filter.isAccepted(anchor, seq);
	}
}